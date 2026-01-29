"""
CompGene Adapter Runner - Unified Execution Framework.

This module provides the AdapterRunner class for executing external tools
through the BaseAdapter interface with unified timeout handling and retry logic.

Source: ADR-002 Tool Adaptation Layer, ADR-003 Error Handling
"""

import os
import random
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.logging import DualLogger
from workflow.lib.audit import collect_run_metadata, write_run_json
from workflow.adapters.base import BaseAdapter, AdapterContext, RunResult


# Maximum characters of stderr to include in error details
MAX_STDERR_DETAILS = 500


# =============================================================================
# RetryConfig - Retry Configuration
# =============================================================================

@dataclass
class RetryConfig:
    """
    Configuration for retry behavior with exponential backoff.

    Attributes:
        max_retries: Maximum number of retry attempts (default: 3).
        base_delay: Initial delay in seconds before first retry (default: 1.0).
        max_delay: Maximum delay in seconds between retries (default: 60.0).
        jitter: Fraction of delay to add as random jitter (default: 0.1).

    Example:
        >>> config = RetryConfig(max_retries=3, base_delay=1.0)
        >>> # First retry waits ~1s, second ~2s, third ~4s (with jitter)
    """

    max_retries: int = 3
    base_delay: float = 1.0
    max_delay: float = 60.0
    jitter: float = 0.1


# =============================================================================
# AdapterRunner - Main Execution Class
# =============================================================================

class AdapterRunner:
    """
    Unified runner for executing external tools through BaseAdapter interface.

    Provides a consistent execution flow:
    1. Version check
    2. Input validation
    3. Command execution with timeout
    4. Output verification
    5. Output parsing
    6. Audit record generation

    Also handles:
    - Graceful timeout with SIGTERM → SIGKILL escalation
    - Automatic retry for retryable errors with exponential backoff

    Attributes:
        retry_config: Configuration for retry behavior.
        logger: Optional DualLogger for logging execution progress.
        grace_period: Seconds to wait after SIGTERM before SIGKILL (default: 10).
        meta_dir: Directory for audit metadata files.

    Example:
        >>> runner = AdapterRunner()
        >>> result = runner.run(adapter, ctx)
    """

    def __init__(
        self,
        retry_config: Optional[RetryConfig] = None,
        logger: Optional[DualLogger] = None,
        grace_period: float = 10.0,
        meta_dir: Path = Path("results/meta")
    ):
        """
        Initialize the AdapterRunner.

        Args:
            retry_config: Retry configuration (default: RetryConfig()).
            logger: Optional logger for execution progress.
            grace_period: Seconds to wait after SIGTERM (default: 10).
            meta_dir: Directory for audit metadata files.
        """
        self.retry_config = retry_config or RetryConfig()
        self.logger = logger
        self.grace_period = grace_period
        self.meta_dir = meta_dir

    def run(self, adapter: BaseAdapter, ctx: AdapterContext) -> RunResult:
        """
        Execute an adapter with full execution flow.

        Executes the following steps in order:
        1. Version check - verify tool version compatibility
        2. Input validation - validate all inputs
        3. Command execution - run with timeout and retry
        4. Output verification - check expected outputs exist
        5. Output parsing - extract results from outputs
        6. Audit recording - write .run.json metadata

        Args:
            adapter: The adapter instance to execute.
            ctx: The execution context with inputs, outputs, config.

        Returns:
            RunResult from adapter.parse_outputs().

        Raises:
            CompGeneError: On version mismatch, validation failure,
                execution error, or missing outputs.
        """
        start_time = time.time()
        rule = adapter.spec.name
        wildcards = ctx.wildcards

        # Use context logger if available, otherwise use runner logger
        logger = ctx.logger or self.logger

        try:
            # Step 1: Version check
            version = self._check_version(adapter, logger)

            # Step 2: Input validation
            self._validate_inputs(adapter, ctx, logger)

            # Step 3: Build command
            cmd = adapter.build_command(ctx)
            if logger:
                logger.info(f"Built command: {' '.join(cmd)}")

            # Step 4: Execute with timeout and retry
            timeout = adapter.timeout_seconds(ctx)
            returncode, stdout, stderr = self._execute_with_retry(
                cmd, timeout, adapter, ctx, logger
            )

            # Step 5: Verify outputs
            self._verify_outputs(adapter, ctx, logger)

            # Step 6: Parse outputs
            result = adapter.parse_outputs(ctx)
            if logger:
                logger.info(f"Parsed outputs: {len(result.outputs)} files")

            # Step 7: Write audit record (success)
            runtime = time.time() - start_time
            self._write_audit(
                rule=rule,
                wildcards=wildcards,
                cmd=cmd,
                version=version,
                ctx=ctx,
                runtime=runtime,
                exit_code=returncode,
                logger=logger
            )

            return result

        except CompGeneError as e:
            # Write audit record (failure)
            runtime = time.time() - start_time
            self._write_audit(
                rule=rule,
                wildcards=wildcards,
                cmd=adapter.build_command(ctx) if hasattr(adapter, 'build_command') else [],
                version=getattr(self, '_cached_version', 'unknown'),
                ctx=ctx,
                runtime=runtime,
                exit_code=e.to_exit_code(),
                error_code=e.error_code.value,
                error_message=str(e),
                logger=logger
            )
            raise

    def _check_version(
        self,
        adapter: BaseAdapter,
        logger: Optional[DualLogger]
    ) -> str:
        """
        Check tool version and verify compatibility.

        Args:
            adapter: The adapter to check.
            logger: Optional logger.

        Returns:
            Version string.

        Raises:
            CompGeneError: If version is incompatible (E_TOOL_VERSION).
        """
        if logger:
            logger.info(f"Checking {adapter.spec.name} version...")

        version = adapter.check_version()
        self._cached_version = version

        if not adapter.spec.check_version_compatible(version):
            msg = (
                f"Version {version} not compatible with "
                f"{adapter.spec.name} requirements "
                f"(min: {adapter.spec.min_version}, max: {adapter.spec.max_version})"
            )
            if logger:
                logger.error(msg)
            raise CompGeneError(ErrorCode.E_TOOL_VERSION, msg)

        if logger:
            logger.info(f"Version {version} is compatible")

        return version

    def _validate_inputs(
        self,
        adapter: BaseAdapter,
        ctx: AdapterContext,
        logger: Optional[DualLogger]
    ) -> None:
        """
        Validate inputs using the adapter.

        Args:
            adapter: The adapter to use for validation.
            ctx: The execution context.
            logger: Optional logger.

        Raises:
            CompGeneError: If validation fails.
        """
        if logger:
            logger.info("Validating inputs...")

        adapter.validate_inputs(ctx)

        if logger:
            logger.info("Input validation passed")

    def _verify_outputs(
        self,
        adapter: BaseAdapter,
        ctx: AdapterContext,
        logger: Optional[DualLogger]
    ) -> None:
        """
        Verify all expected outputs exist.

        Args:
            adapter: The adapter.
            ctx: The execution context.
            logger: Optional logger.

        Raises:
            CompGeneError: If any output is missing (E_OUTPUT_MISSING).
        """
        expected = adapter.expected_outputs(ctx)
        missing = [p for p in expected if not p.exists()]

        if missing:
            paths_str = ", ".join(str(p) for p in missing)
            msg = f"Expected outputs not found: {paths_str}"
            if logger:
                logger.error(msg)
            raise CompGeneError(
                ErrorCode.E_OUTPUT_MISSING,
                msg,
                details=f"Missing {len(missing)} of {len(expected)} expected outputs"
            )

        if logger:
            logger.info(f"All {len(expected)} expected outputs verified")

    def _execute_with_retry(
        self,
        cmd: list[str],
        timeout: int,
        adapter: BaseAdapter,
        ctx: AdapterContext,
        logger: Optional[DualLogger]
    ) -> tuple[int, str, str]:
        """
        Execute command with timeout and retry logic.

        Args:
            cmd: Command to execute.
            timeout: Timeout in seconds.
            adapter: The adapter for error classification.
            ctx: The execution context.
            logger: Optional logger.

        Returns:
            Tuple of (returncode, stdout, stderr).

        Raises:
            CompGeneError: On non-retryable error or max retries exceeded.
        """
        last_error: Optional[CompGeneError] = None

        for attempt in range(self.retry_config.max_retries + 1):
            if attempt > 0 and logger:
                logger.warning(f"Retry attempt {attempt}/{self.retry_config.max_retries}")

            try:
                returncode, stdout, stderr = self._execute_once(cmd, timeout, logger)

                if returncode == 0:
                    return returncode, stdout, stderr

                # Classify the error
                error_code, is_retryable = adapter.classify_error(ctx, returncode, stderr)

                if not is_retryable:
                    if logger:
                        logger.error(f"Non-retryable error: {error_code.value}")
                    raise CompGeneError(
                        error_code,
                        f"Command failed with exit code {returncode}",
                        details=stderr[:MAX_STDERR_DETAILS] if stderr else None
                    )

                # Retryable error - check if we have retries left
                if attempt >= self.retry_config.max_retries:
                    if logger:
                        logger.error(f"Max retries ({self.retry_config.max_retries}) exceeded")
                    raise CompGeneError(
                        error_code,
                        f"Command failed after {self.retry_config.max_retries + 1} attempts",
                        details=stderr[:MAX_STDERR_DETAILS] if stderr else None
                    )

                # Calculate backoff and wait
                delay = self._calculate_backoff(attempt)
                if logger:
                    logger.warning(f"Retryable error ({error_code.value}), waiting {delay:.1f}s...")
                time.sleep(delay)

            except CompGeneError as e:
                last_error = e
                if not e.is_retryable or attempt >= self.retry_config.max_retries:
                    raise

                delay = self._calculate_backoff(attempt)
                if logger:
                    logger.warning(f"Retryable error, waiting {delay:.1f}s before retry...")
                time.sleep(delay)

        # Should not reach here, but just in case
        if last_error:
            raise last_error
        raise CompGeneError(
            ErrorCode.E_NONZERO_EXIT,
            "Unexpected execution flow in retry logic"
        )

    def _execute_once(
        self,
        cmd: list[str],
        timeout: int,
        logger: Optional[DualLogger]
    ) -> tuple[int, str, str]:
        """
        Execute command once with timeout handling.

        Args:
            cmd: Command to execute.
            timeout: Timeout in seconds.
            logger: Optional logger.

        Returns:
            Tuple of (returncode, stdout, stderr).

        Raises:
            CompGeneError: On timeout (E_TIMEOUT).
        """
        if logger:
            logger.info(f"Executing: {cmd[0]} (timeout: {timeout}s)")

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True  # Allow killing entire process group
        )

        try:
            stdout, stderr = process.communicate(timeout=timeout)
            return process.returncode, stdout, stderr

        except subprocess.TimeoutExpired:
            # Graceful termination sequence
            return self._handle_timeout(process, logger)

    def _handle_timeout(
        self,
        process: subprocess.Popen,
        logger: Optional[DualLogger]
    ) -> tuple[int, str, str]:
        """
        Handle process timeout with graceful termination.

        Sends SIGTERM, waits for grace period, then SIGKILL if needed.

        Args:
            process: The timed-out process.
            logger: Optional logger.

        Raises:
            CompGeneError: Always raises E_TIMEOUT.
        """
        if logger:
            logger.warning("Process timed out, sending SIGTERM...")

        # Send SIGTERM to process group
        # Use try-except as process may have already exited
        try:
            pgid = os.getpgid(process.pid)
            os.killpg(pgid, signal.SIGTERM)
        except (ProcessLookupError, OSError) as e:
            if logger:
                logger.warning(f"Could not send SIGTERM to process group: {e}")
            # Process may have already exited, try to get output
            stdout, stderr = process.communicate()
            raise CompGeneError(
                ErrorCode.E_TIMEOUT,
                "Process timed out and was terminated",
                details=stderr[:MAX_STDERR_DETAILS] if stderr else None
            )

        # Wait for grace period
        try:
            stdout, stderr = process.communicate(timeout=self.grace_period)
            if logger:
                logger.info("Process terminated gracefully after SIGTERM")

        except subprocess.TimeoutExpired:
            # Force kill
            if logger:
                logger.warning(f"Grace period ({self.grace_period}s) expired, sending SIGKILL...")
            try:
                os.killpg(pgid, signal.SIGKILL)
            except (ProcessLookupError, OSError) as e:
                if logger:
                    logger.warning(f"Could not send SIGKILL to process group: {e}")
            stdout, stderr = process.communicate()
            if logger:
                logger.info("Process killed with SIGKILL")

        raise CompGeneError(
            ErrorCode.E_TIMEOUT,
            "Process timed out and was terminated",
            details=stderr[:MAX_STDERR_DETAILS] if stderr else None
        )

    def _calculate_backoff(self, attempt: int) -> float:
        """
        Calculate exponential backoff delay with jitter.

        Args:
            attempt: Current attempt number (0-indexed).

        Returns:
            Delay in seconds.
        """
        base = self.retry_config.base_delay * (2 ** attempt)
        delay = min(base, self.retry_config.max_delay)

        # Add jitter
        jitter_range = delay * self.retry_config.jitter
        delay += random.uniform(-jitter_range, jitter_range)

        return max(0, delay)  # Ensure non-negative

    def _write_audit(
        self,
        rule: str,
        wildcards: dict[str, str],
        cmd: list[str],
        version: str,
        ctx: AdapterContext,
        runtime: float,
        exit_code: int,
        error_code: Optional[str] = None,
        error_message: Optional[str] = None,
        logger: Optional[DualLogger] = None
    ) -> Path:
        """
        Write audit metadata to .run.json file.

        Args:
            rule: Rule name.
            wildcards: Wildcard values.
            cmd: Executed command.
            version: Tool version.
            ctx: Execution context.
            runtime: Runtime in seconds.
            exit_code: Process exit code.
            error_code: Optional error code if failed.
            error_message: Optional error message if failed.
            logger: Optional logger.

        Returns:
            Path to written audit file.
        """
        # Prepare output paths for checksums (only on success)
        output_paths = None
        if exit_code == 0:
            output_paths = ctx.outputs

        metadata = collect_run_metadata(
            rule=rule,
            wildcards=wildcards,
            cmd=cmd,
            tool_version=version,
            input_paths=ctx.inputs,
            threads=ctx.threads,
            runtime_seconds=runtime,
            exit_code=exit_code,
            output_paths=output_paths,
            error_code=error_code,
            error_message=error_message
        )

        audit_path = write_run_json(metadata, self.meta_dir)

        if logger:
            logger.info(f"Audit written to: {audit_path}")

        return audit_path
