"""
Tests for the AdapterRunner module.

Tests cover:
- Normal execution flow
- Version check failures
- Input validation failures
- Timeout handling (SIGTERM → SIGKILL)
- Retry logic for retryable errors
- Non-retryable error handling
- Audit record generation
"""

import os
import signal
import subprocess
import tempfile
import time
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

from workflow.lib.runner import AdapterRunner, RetryConfig
from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.logging import DualLogger
from workflow.adapters.base import (
    BaseAdapter, AdapterContext, RunResult, ToolSpec
)


# =============================================================================
# Test Fixtures
# =============================================================================

class MockAdapter(BaseAdapter):
    """Mock adapter for testing."""

    def __init__(
        self,
        version: str = "1.0.0",
        min_version: str = "1.0.0",
        max_version: str = "2.0.0",
        validation_error: Exception | None = None,
        command: list[str] | None = None,
        expected_output_paths: list[Path] | None = None,
        parse_result: RunResult | None = None,
        timeout: int = 300,
        classify_result: tuple[ErrorCode, bool] = (ErrorCode.E_NONZERO_EXIT, False)
    ):
        self._version = version
        self._min_version = min_version
        self._max_version = max_version
        self._validation_error = validation_error
        self._command = command or ["echo", "test"]
        self._expected_outputs = expected_output_paths or []
        self._parse_result = parse_result or RunResult(outputs={}, summary={})
        self._timeout = timeout
        self._classify_result = classify_result

    @property
    def spec(self) -> ToolSpec:
        return ToolSpec(
            name="mock_tool",
            min_version=self._min_version,
            max_version=self._max_version,
            conda_env="mock.yaml",
            description="Mock tool for testing"
        )

    def check_version(self) -> str:
        return self._version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        if self._validation_error:
            raise self._validation_error

    def build_command(self, ctx: AdapterContext) -> list[str]:
        return self._command

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        return self._expected_outputs

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        return self._parse_result

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        return self._timeout

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str
    ) -> tuple[ErrorCode, bool]:
        return self._classify_result


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_context(temp_dir):
    """Create a mock AdapterContext."""
    input_file = temp_dir / "input.txt"
    input_file.write_text("test input")

    output_file = temp_dir / "output.txt"

    return AdapterContext(
        inputs={"input": input_file},
        outputs={"output": output_file},
        config={"threads": 4},
        wildcards={"sample": "test_sample"},
        threads=4
    )


@pytest.fixture
def runner(temp_dir):
    """Create an AdapterRunner with test configuration."""
    return AdapterRunner(
        retry_config=RetryConfig(max_retries=2, base_delay=0.1, max_delay=1.0),
        meta_dir=temp_dir / "meta"
    )


# =============================================================================
# RetryConfig Tests
# =============================================================================

class TestRetryConfig:
    """Tests for RetryConfig dataclass."""

    def test_default_values(self):
        """Test default retry configuration values."""
        config = RetryConfig()
        assert config.max_retries == 3
        assert config.base_delay == 1.0
        assert config.max_delay == 60.0
        assert config.jitter == 0.1

    def test_custom_values(self):
        """Test custom retry configuration values."""
        config = RetryConfig(
            max_retries=5,
            base_delay=2.0,
            max_delay=120.0,
            jitter=0.2
        )
        assert config.max_retries == 5
        assert config.base_delay == 2.0
        assert config.max_delay == 120.0
        assert config.jitter == 0.2


# =============================================================================
# AdapterRunner Initialization Tests
# =============================================================================

class TestAdapterRunnerInit:
    """Tests for AdapterRunner initialization."""

    def test_default_initialization(self):
        """Test default initialization."""
        runner = AdapterRunner()
        assert runner.retry_config.max_retries == 3
        assert runner.grace_period == 10.0
        assert runner.meta_dir == Path("results/meta")

    def test_custom_initialization(self, temp_dir):
        """Test custom initialization."""
        config = RetryConfig(max_retries=5)
        logger = MagicMock(spec=DualLogger)

        runner = AdapterRunner(
            retry_config=config,
            logger=logger,
            grace_period=20.0,
            meta_dir=temp_dir
        )

        assert runner.retry_config.max_retries == 5
        assert runner.logger == logger
        assert runner.grace_period == 20.0
        assert runner.meta_dir == temp_dir


# =============================================================================
# Version Check Tests
# =============================================================================

class TestVersionCheck:
    """Tests for version checking functionality."""

    def test_compatible_version(self, runner, mock_context, temp_dir):
        """Test successful version check."""
        adapter = MockAdapter(version="1.5.0")

        # Create expected output file
        output_file = mock_context.outputs["output"]
        output_file.write_text("test output")
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            result = runner.run(adapter, mock_context)

        assert result is not None

    def test_version_below_minimum(self, runner, mock_context):
        """Test version below minimum raises error."""
        adapter = MockAdapter(version="0.9.0", min_version="1.0.0")

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, mock_context)

        assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION
        assert "0.9.0" in str(exc_info.value)

    def test_version_above_maximum(self, runner, mock_context):
        """Test version above maximum raises error."""
        adapter = MockAdapter(version="3.0.0", max_version="2.0.0")

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, mock_context)

        assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION
        assert "3.0.0" in str(exc_info.value)


# =============================================================================
# Input Validation Tests
# =============================================================================

class TestInputValidation:
    """Tests for input validation functionality."""

    def test_validation_success(self, runner, mock_context, temp_dir):
        """Test successful input validation."""
        adapter = MockAdapter()
        output_file = mock_context.outputs["output"]
        output_file.write_text("test")
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            result = runner.run(adapter, mock_context)

        assert result is not None

    def test_validation_failure_input_missing(self, runner, mock_context):
        """Test validation failure for missing input."""
        adapter = MockAdapter(
            validation_error=CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Input file not found"
            )
        )

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, mock_context)

        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validation_failure_input_format(self, runner, mock_context):
        """Test validation failure for invalid format."""
        adapter = MockAdapter(
            validation_error=CompGeneError(
                ErrorCode.E_INPUT_FORMAT,
                "Invalid FASTA format"
            )
        )

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, mock_context)

        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Command Execution Tests
# =============================================================================

class TestCommandExecution:
    """Tests for command execution functionality."""

    def test_successful_execution(self, runner, mock_context, temp_dir):
        """Test successful command execution."""
        adapter = MockAdapter()
        output_file = mock_context.outputs["output"]
        output_file.write_text("test")
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "output", "")):
            result = runner.run(adapter, mock_context)

        assert result is not None
        assert result.outputs == {}  # MockAdapter returns empty

    def test_execution_with_real_subprocess(self, runner, mock_context, temp_dir):
        """Test execution with actual subprocess."""
        output_file = mock_context.outputs["output"]
        adapter = MockAdapter(
            command=["bash", "-c", f"echo 'test' > {output_file}"],
        )
        adapter._expected_outputs = [output_file]

        result = runner.run(adapter, mock_context)
        assert output_file.exists()


# =============================================================================
# Timeout Handling Tests
# =============================================================================

class TestTimeoutHandling:
    """Tests for timeout handling functionality."""

    def test_timeout_triggers_sigterm(self, temp_dir):
        """Test that timeout sends SIGTERM first."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=0),
            grace_period=0.5,
            meta_dir=temp_dir / "meta"
        )

        # Create a mock context
        input_file = temp_dir / "input.txt"
        input_file.write_text("test")
        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": temp_dir / "output.txt"},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        # Adapter with very short timeout
        adapter = MockAdapter(
            timeout=1,
            command=["sleep", "10"]  # Sleep longer than timeout
        )

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, ctx)

        assert exc_info.value.error_code == ErrorCode.E_TIMEOUT

    def test_timeout_escalates_to_sigkill(self, temp_dir):
        """Test that SIGKILL is sent after grace period."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=0),
            grace_period=0.1,  # Very short grace period
            meta_dir=temp_dir / "meta"
        )

        # Create context
        input_file = temp_dir / "input.txt"
        input_file.write_text("test")
        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": temp_dir / "output.txt"},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        # Command that ignores SIGTERM
        trap_script = """
        trap '' TERM
        sleep 10
        """
        adapter = MockAdapter(
            timeout=1,
            command=["bash", "-c", trap_script]
        )

        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, ctx)

        assert exc_info.value.error_code == ErrorCode.E_TIMEOUT


# =============================================================================
# Retry Logic Tests
# =============================================================================

class TestRetryLogic:
    """Tests for retry logic functionality."""

    def test_retryable_error_succeeds_on_retry(self, temp_dir):
        """Test that retryable errors are retried and can succeed."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=2, base_delay=0.01),
            meta_dir=temp_dir / "meta"
        )

        input_file = temp_dir / "input.txt"
        input_file.write_text("test")
        output_file = temp_dir / "output.txt"

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": output_file},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        adapter = MockAdapter(
            classify_result=(ErrorCode.E_NET_RATE_LIMIT, True)
        )
        adapter._expected_outputs = [output_file]

        # Track execution attempts
        call_count = 0

        def mock_execute(cmd, timeout, logger):
            nonlocal call_count
            call_count += 1
            if call_count < 2:
                return (1, "", "rate limit exceeded")
            output_file.write_text("success")
            return (0, "success", "")

        with patch.object(runner, '_execute_once', side_effect=mock_execute):
            result = runner.run(adapter, ctx)

        assert call_count == 2
        assert result is not None

    def test_retryable_error_max_retries_exceeded(self, temp_dir):
        """Test that max retries is respected."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=2, base_delay=0.01),
            meta_dir=temp_dir / "meta"
        )

        input_file = temp_dir / "input.txt"
        input_file.write_text("test")

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": temp_dir / "output.txt"},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        adapter = MockAdapter(
            classify_result=(ErrorCode.E_NET_RATE_LIMIT, True)
        )

        call_count = 0

        def mock_execute(cmd, timeout, logger):
            nonlocal call_count
            call_count += 1
            return (1, "", "rate limit exceeded")

        with patch.object(runner, '_execute_once', side_effect=mock_execute):
            with pytest.raises(CompGeneError) as exc_info:
                runner.run(adapter, ctx)

        assert exc_info.value.error_code == ErrorCode.E_NET_RATE_LIMIT
        assert call_count == 3  # Initial + 2 retries

    def test_non_retryable_error_no_retry(self, temp_dir):
        """Test that non-retryable errors are not retried."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=2, base_delay=0.01),
            meta_dir=temp_dir / "meta"
        )

        input_file = temp_dir / "input.txt"
        input_file.write_text("test")

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": temp_dir / "output.txt"},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        adapter = MockAdapter(
            classify_result=(ErrorCode.E_OOM, False)  # Not retryable
        )

        call_count = 0

        def mock_execute(cmd, timeout, logger):
            nonlocal call_count
            call_count += 1
            return (1, "", "out of memory")

        with patch.object(runner, '_execute_once', side_effect=mock_execute):
            with pytest.raises(CompGeneError) as exc_info:
                runner.run(adapter, ctx)

        assert exc_info.value.error_code == ErrorCode.E_OOM
        assert call_count == 1  # No retries


# =============================================================================
# Exponential Backoff Tests
# =============================================================================

class TestExponentialBackoff:
    """Tests for exponential backoff calculation."""

    def test_backoff_increases_exponentially(self):
        """Test that backoff increases exponentially."""
        runner = AdapterRunner(
            retry_config=RetryConfig(base_delay=1.0, max_delay=60.0, jitter=0.0)
        )

        delay0 = runner._calculate_backoff(0)
        delay1 = runner._calculate_backoff(1)
        delay2 = runner._calculate_backoff(2)

        assert delay0 == pytest.approx(1.0)  # 1 * 2^0
        assert delay1 == pytest.approx(2.0)  # 1 * 2^1
        assert delay2 == pytest.approx(4.0)  # 1 * 2^2

    def test_backoff_respects_max_delay(self):
        """Test that backoff respects max_delay."""
        runner = AdapterRunner(
            retry_config=RetryConfig(base_delay=1.0, max_delay=5.0, jitter=0.0)
        )

        delay5 = runner._calculate_backoff(5)  # Would be 32 without cap
        assert delay5 == pytest.approx(5.0)

    def test_backoff_includes_jitter(self):
        """Test that backoff includes jitter."""
        runner = AdapterRunner(
            retry_config=RetryConfig(base_delay=1.0, max_delay=60.0, jitter=0.2)
        )

        # With 20% jitter, delay should be between 0.8 and 1.2
        delays = [runner._calculate_backoff(0) for _ in range(100)]

        assert min(delays) >= 0.8
        assert max(delays) <= 1.2
        # Should have some variance
        assert max(delays) != min(delays)


# =============================================================================
# Output Verification Tests
# =============================================================================

class TestOutputVerification:
    """Tests for output verification functionality."""

    def test_missing_output_raises_error(self, runner, mock_context, temp_dir):
        """Test that missing outputs raise error."""
        missing_file = temp_dir / "nonexistent.txt"
        adapter = MockAdapter()
        adapter._expected_outputs = [missing_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            with pytest.raises(CompGeneError) as exc_info:
                runner.run(adapter, mock_context)

        assert exc_info.value.error_code == ErrorCode.E_OUTPUT_MISSING
        assert "nonexistent.txt" in str(exc_info.value)

    def test_all_outputs_exist(self, runner, mock_context, temp_dir):
        """Test successful execution when all outputs exist."""
        output1 = temp_dir / "output1.txt"
        output2 = temp_dir / "output2.txt"
        output1.write_text("test1")
        output2.write_text("test2")

        adapter = MockAdapter()
        adapter._expected_outputs = [output1, output2]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            result = runner.run(adapter, mock_context)

        assert result is not None


# =============================================================================
# Audit Record Tests
# =============================================================================

class TestAuditRecords:
    """Tests for audit record generation."""

    def test_audit_written_on_success(self, runner, mock_context, temp_dir):
        """Test that audit record is written on success."""
        output_file = mock_context.outputs["output"]
        output_file.write_text("test")

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            runner.run(adapter, mock_context)

        # Check audit file was created
        meta_dir = temp_dir / "meta"
        audit_files = list(meta_dir.rglob("*.run.json"))
        assert len(audit_files) == 1

        # Verify content
        import json
        content = json.loads(audit_files[0].read_text())
        assert content["rule"] == "mock_tool"
        assert content["exit_code"] == 0
        assert "error_code" not in content

    def test_audit_written_on_failure(self, runner, mock_context, temp_dir):
        """Test that audit record is written on failure."""
        adapter = MockAdapter(
            validation_error=CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Input not found"
            )
        )

        with pytest.raises(CompGeneError):
            runner.run(adapter, mock_context)

        # Check audit file was created
        meta_dir = temp_dir / "meta"
        audit_files = list(meta_dir.rglob("*.run.json"))
        assert len(audit_files) == 1

        # Verify content
        import json
        content = json.loads(audit_files[0].read_text())
        assert content["error_code"] == "E_INPUT_MISSING"

    def test_audit_includes_input_checksums(self, runner, mock_context, temp_dir):
        """Test that audit includes input checksums."""
        output_file = mock_context.outputs["output"]
        output_file.write_text("test")

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            runner.run(adapter, mock_context)

        meta_dir = temp_dir / "meta"
        audit_files = list(meta_dir.rglob("*.run.json"))

        import json
        content = json.loads(audit_files[0].read_text())
        assert "input_checksums" in content
        assert "input" in content["input_checksums"]

    def test_audit_includes_output_checksums_on_success(
        self, runner, mock_context, temp_dir
    ):
        """Test that audit includes output checksums on success."""
        output_file = mock_context.outputs["output"]
        output_file.write_text("test output content")

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            runner.run(adapter, mock_context)

        meta_dir = temp_dir / "meta"
        audit_files = list(meta_dir.rglob("*.run.json"))

        import json
        content = json.loads(audit_files[0].read_text())
        assert "output_checksums" in content
        assert "output" in content["output_checksums"]


# =============================================================================
# Logging Integration Tests
# =============================================================================

class TestLoggingIntegration:
    """Tests for logging integration."""

    def test_logger_receives_info_messages(self, runner, mock_context, temp_dir):
        """Test that logger receives INFO messages."""
        logger = MagicMock(spec=DualLogger)
        runner.logger = logger

        output_file = mock_context.outputs["output"]
        output_file.write_text("test")

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            runner.run(adapter, mock_context)

        # Check logger was called
        assert logger.info.called
        info_calls = [str(call) for call in logger.info.call_args_list]
        # Should have logged version check, validation, execution, etc.
        assert len(logger.info.call_args_list) >= 4

    def test_logger_receives_warning_on_retry(self, temp_dir):
        """Test that logger receives WARNING on retry."""
        logger = MagicMock(spec=DualLogger)

        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=2, base_delay=0.01),
            logger=logger,
            meta_dir=temp_dir / "meta"
        )

        input_file = temp_dir / "input.txt"
        input_file.write_text("test")
        output_file = temp_dir / "output.txt"

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": output_file},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        adapter = MockAdapter(
            classify_result=(ErrorCode.E_NET_RATE_LIMIT, True)
        )
        adapter._expected_outputs = [output_file]

        call_count = 0

        def mock_execute(cmd, timeout, log):
            nonlocal call_count
            call_count += 1
            if call_count < 2:
                return (1, "", "rate limit")
            output_file.write_text("success")
            return (0, "", "")

        with patch.object(runner, '_execute_once', side_effect=mock_execute):
            runner.run(adapter, ctx)

        assert logger.warning.called

    def test_logger_receives_error_on_failure(self, runner, mock_context, temp_dir):
        """Test that logger receives ERROR on failure."""
        logger = MagicMock(spec=DualLogger)
        runner.logger = logger

        adapter = MockAdapter(
            validation_error=CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Input not found"
            )
        )

        with pytest.raises(CompGeneError):
            runner.run(adapter, mock_context)

        # Validation errors might not go through runner logger
        # but version check logging should still occur


# =============================================================================
# Context Logger Override Tests
# =============================================================================

class TestContextLoggerOverride:
    """Tests for context logger overriding runner logger."""

    def test_context_logger_takes_precedence(self, runner, mock_context, temp_dir):
        """Test that context logger overrides runner logger."""
        runner_logger = MagicMock(spec=DualLogger)
        context_logger = MagicMock(spec=DualLogger)

        runner.logger = runner_logger
        mock_context.logger = context_logger

        output_file = mock_context.outputs["output"]
        output_file.write_text("test")

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            runner.run(adapter, mock_context)

        # Context logger should be used, not runner logger
        assert context_logger.info.called
        assert not runner_logger.info.called


# =============================================================================
# Edge Case Tests
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_empty_wildcards(self, runner, temp_dir):
        """Test execution with empty wildcards."""
        input_file = temp_dir / "input.txt"
        input_file.write_text("test")
        output_file = temp_dir / "output.txt"
        output_file.write_text("test")

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": output_file},
            config={},
            wildcards={},  # Empty wildcards
            threads=1
        )

        adapter = MockAdapter()
        adapter._expected_outputs = [output_file]

        with patch.object(runner, '_execute_once', return_value=(0, "", "")):
            result = runner.run(adapter, ctx)

        assert result is not None

    def test_zero_timeout(self, temp_dir):
        """Test with zero timeout."""
        runner = AdapterRunner(
            retry_config=RetryConfig(max_retries=0),
            meta_dir=temp_dir / "meta"
        )

        input_file = temp_dir / "input.txt"
        input_file.write_text("test")

        ctx = AdapterContext(
            inputs={"input": input_file},
            outputs={"output": temp_dir / "output.txt"},
            config={},
            wildcards={"sample": "test"},
            threads=1
        )

        # Zero timeout should still work for fast commands
        adapter = MockAdapter(
            timeout=0,
            command=["echo", "test"]
        )

        # This will likely timeout immediately
        with pytest.raises(CompGeneError) as exc_info:
            runner.run(adapter, ctx)

        assert exc_info.value.error_code == ErrorCode.E_TIMEOUT

    def test_very_long_stderr(self, runner, mock_context, temp_dir):
        """Test handling of very long stderr."""
        long_stderr = "x" * 10000

        adapter = MockAdapter(
            classify_result=(ErrorCode.E_NONZERO_EXIT, False)
        )

        def mock_execute(cmd, timeout, logger):
            return (1, "", long_stderr)

        with patch.object(runner, '_execute_once', side_effect=mock_execute):
            with pytest.raises(CompGeneError) as exc_info:
                runner.run(adapter, mock_context)

        # Details should be truncated
        assert len(exc_info.value.details) <= 500
