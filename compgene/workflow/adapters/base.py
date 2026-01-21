"""
CompGene Adapter Framework - Base Classes and Data Structures.

This module provides the foundation for external tool integration:
- ToolSpec: Tool specification with version constraints
- RunResult: Standardized execution results
- AdapterContext: Runtime context for adapter methods
- BaseAdapter: Abstract base class for all tool adapters

Source: ADR-002 Tool Adaptation Layer
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
import json
import re

from workflow.lib.errors import ErrorCode
from workflow.lib.logging import DualLogger


# =============================================================================
# ToolSpec - Tool Specification
# =============================================================================

@dataclass
class ToolSpec:
    """
    Specification for an external tool adapter.

    Defines the tool's identity, version constraints, and conda environment
    for consistent tool management across the pipeline.

    Attributes:
        name: Tool name (e.g., "orthofinder", "busco").
        min_version: Minimum supported version (inclusive).
        max_version: Maximum supported version (inclusive), None for no upper limit.
        conda_env: Name of conda environment file (e.g., "orthofinder.yaml").
        description: Optional human-readable description of the tool.

    Example:
        >>> spec = ToolSpec(
        ...     name="orthofinder",
        ...     min_version="2.5.0",
        ...     max_version="2.5.99",
        ...     conda_env="orthofinder.yaml",
        ...     description="OrthoFinder for orthology inference"
        ... )
    """

    name: str
    min_version: str
    max_version: Optional[str] = None
    conda_env: str = ""
    description: str = ""

    def check_version_compatible(self, version: str) -> bool:
        """
        Check if a version string is compatible with this spec.

        Args:
            version: Version string to check (e.g., "2.5.5").

        Returns:
            True if version is within [min_version, max_version] range.
        """
        parsed = self._parse_version(version)
        min_parsed = self._parse_version(self.min_version)

        if parsed < min_parsed:
            return False

        if self.max_version is not None:
            max_parsed = self._parse_version(self.max_version)
            if parsed > max_parsed:
                return False

        return True

    def _parse_version(self, version: str) -> tuple[int, ...]:
        """
        Parse a version string into a comparable tuple.

        Handles versions like "2.5.5", "1.6.3", "5.4.7".
        Non-numeric parts are stripped.

        Args:
            version: Version string to parse.

        Returns:
            Tuple of integers for comparison.
        """
        # Extract numeric parts only
        parts = re.findall(r'\d+', version)
        return tuple(int(p) for p in parts) if parts else (0,)

    def to_dict(self) -> dict[str, Any]:
        """
        Convert to dictionary for serialization.

        Returns:
            Dictionary with all spec fields. Empty optional fields are omitted.
        """
        result: dict[str, Any] = {
            "name": self.name,
            "min_version": self.min_version,
        }
        if self.max_version is not None:
            result["max_version"] = self.max_version
        if self.conda_env:
            result["conda_env"] = self.conda_env
        if self.description:
            result["description"] = self.description
        return result


# =============================================================================
# RunResult - Execution Result
# =============================================================================

@dataclass
class RunResult:
    """
    Standardized result from an adapter's parse_outputs method.

    Captures execution outputs, summary statistics, and optional metrics
    for consistent result handling across all adapters.

    Attributes:
        outputs: Dictionary mapping output names to file paths.
        summary: Dictionary with summary statistics (tool-specific).
        metrics: Optional dictionary with performance metrics.
        warnings: Optional list of warning messages.

    Example:
        >>> result = RunResult(
        ...     outputs={"orthogroups": Path("results/Orthogroups.tsv")},
        ...     summary={"total_orthogroups": 12345, "species_count": 5},
        ...     metrics={"runtime_seconds": 3600.5},
        ...     warnings=["Low memory detected during execution"]
        ... )
    """

    outputs: dict[str, Path]
    summary: dict[str, Any]
    metrics: Optional[dict[str, Any]] = None
    warnings: Optional[list[str]] = None

    def to_dict(self) -> dict[str, Any]:
        """
        Convert to dictionary for serialization.

        Converts Path objects to strings for JSON compatibility.

        Returns:
            Dictionary with all result fields.
        """
        result: dict[str, Any] = {
            "outputs": {k: str(v) for k, v in self.outputs.items()},
            "summary": self.summary,
        }
        if self.metrics is not None:
            result["metrics"] = self.metrics
        if self.warnings is not None:
            result["warnings"] = self.warnings
        return result

    def to_json(self, indent: int = 2) -> str:
        """
        Convert to JSON string.

        Args:
            indent: JSON indentation level.

        Returns:
            JSON string representation.
        """
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)


# =============================================================================
# AdapterContext - Runtime Context
# =============================================================================

@dataclass
class AdapterContext:
    """
    Runtime context passed to adapter methods.

    Provides all necessary information for adapter methods to execute,
    including input/output paths, configuration, and logging.

    Attributes:
        inputs: Dictionary mapping input names to paths.
        outputs: Dictionary mapping output names to expected paths.
        config: Configuration dictionary from Snakemake.
        wildcards: Snakemake wildcards for this rule execution.
        threads: Number of threads allocated for execution.
        logger: Optional DualLogger for structured logging.
        temp_dir: Optional temporary directory for intermediate files.
            (Extension to ADR-002 for adapter-specific temp file management)

    Example:
        >>> ctx = AdapterContext(
        ...     inputs={"proteins": Path("data/proteins.fa")},
        ...     outputs={"orthogroups": Path("results/orthogroups.tsv")},
        ...     config={"threads": 8, "memory_mb": 16000},
        ...     wildcards={"species_set": "lemur_macaque"},
        ...     threads=8
        ... )
    """

    inputs: dict[str, Path]
    outputs: dict[str, Path]
    config: dict[str, Any]
    wildcards: dict[str, str]
    threads: int
    logger: Optional[DualLogger] = None
    temp_dir: Optional[Path] = None

    def get_input(self, name: str) -> Path:
        """
        Get an input path by name.

        Args:
            name: Input name.

        Returns:
            Path to the input.

        Raises:
            KeyError: If input name not found.
        """
        if name not in self.inputs:
            raise KeyError(f"Input '{name}' not found. Available: {list(self.inputs.keys())}")
        return self.inputs[name]

    def get_output(self, name: str) -> Path:
        """
        Get an output path by name.

        Args:
            name: Output name.

        Returns:
            Path to the output.

        Raises:
            KeyError: If output name not found.
        """
        if name not in self.outputs:
            raise KeyError(f"Output '{name}' not found. Available: {list(self.outputs.keys())}")
        return self.outputs[name]

    def get_config(self, key: str, default: Any = None) -> Any:
        """
        Get a configuration value.

        Args:
            key: Configuration key.
            default: Default value if key not found.

        Returns:
            Configuration value or default.
        """
        return self.config.get(key, default)

    def get_wildcard(self, name: str) -> str:
        """
        Get a wildcard value.

        Args:
            name: Wildcard name.

        Returns:
            Wildcard value.

        Raises:
            KeyError: If wildcard name not found.
        """
        if name not in self.wildcards:
            raise KeyError(f"Wildcard '{name}' not found. Available: {list(self.wildcards.keys())}")
        return self.wildcards[name]


# =============================================================================
# BaseAdapter - Abstract Base Class
# =============================================================================

class BaseAdapter(ABC):
    """
    Abstract base class for all external tool adapters.

    Defines the 8-method interface that all adapters must implement
    for consistent tool integration across the CompGene pipeline.

    Required Methods:
        - spec: Tool specification property
        - check_version: Get tool version string
        - validate_inputs: Validate inputs before execution
        - build_command: Build the command line arguments
        - expected_outputs: List expected output files
        - parse_outputs: Parse outputs after execution
        - timeout_seconds: Get execution timeout
        - classify_error: Classify errors for recovery

    Example:
        >>> class OrthoFinderAdapter(BaseAdapter):
        ...     @property
        ...     def spec(self) -> ToolSpec:
        ...         return ToolSpec(name="orthofinder", min_version="2.5.0")
        ...     # ... implement other methods
    """

    @property
    @abstractmethod
    def spec(self) -> ToolSpec:
        """
        Get the tool specification.

        Returns:
            ToolSpec instance describing the tool.
        """
        ...

    @abstractmethod
    def check_version(self) -> str:
        """
        Get the installed tool version.

        Should execute the tool's version command and parse the output.

        Returns:
            Version string (e.g., "2.5.5").

        Raises:
            CompGeneError: If tool not found or version cannot be determined.
        """
        ...

    @abstractmethod
    def validate_inputs(self, ctx: AdapterContext) -> None:
        """
        Validate inputs before execution.

        Check that all required inputs exist and are in the correct format.

        Args:
            ctx: Adapter context with inputs and configuration.

        Raises:
            CompGeneError: If validation fails (E_INPUT_MISSING, E_INPUT_FORMAT).
        """
        ...

    @abstractmethod
    def build_command(self, ctx: AdapterContext) -> list[str]:
        """
        Build the command line arguments for execution.

        Args:
            ctx: Adapter context with inputs, outputs, and configuration.

        Returns:
            List of command line arguments (first element is the executable).
        """
        ...

    @abstractmethod
    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        """
        List the expected output files after successful execution.

        Used to verify that the tool completed successfully.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            List of expected output file paths.
        """
        ...

    @abstractmethod
    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        """
        Parse outputs after successful execution.

        Extract summary statistics and metrics from tool outputs.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            RunResult with outputs, summary, and optional metrics.
        """
        ...

    @abstractmethod
    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """
        Get the execution timeout in seconds.

        May vary based on input size or configuration.

        Args:
            ctx: Adapter context for calculating timeout.

        Returns:
            Timeout in seconds.
        """
        ...

    @abstractmethod
    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str
    ) -> tuple[ErrorCode, bool]:
        """
        Classify an execution error for recovery handling.

        Analyze the return code and stderr to determine the error type
        and whether automatic retry is appropriate.

        Args:
            ctx: Adapter context.
            returncode: Process return code.
            stderr: Standard error output.

        Returns:
            Tuple of (ErrorCode, is_retryable).
        """
        ...


# =============================================================================
# Error Classification Helpers
# =============================================================================

def classify_common_errors(returncode: int, stderr: str) -> tuple[ErrorCode, bool]:
    """
    Classify common error patterns shared across tools.

    Helper function for adapter implementations to handle common cases.

    Args:
        returncode: Process return code.
        stderr: Standard error output.

    Returns:
        Tuple of (ErrorCode, is_retryable).
    """
    stderr_lower = stderr.lower()

    # Tool not found (command not found)
    if returncode == 127:
        return (ErrorCode.E_TOOL_NOT_FOUND, False)

    # Permission denied
    if returncode == 126:
        return (ErrorCode.E_TOOL_NOT_FOUND, False)

    # Timeout indicators
    if "timeout" in stderr_lower or "timed out" in stderr_lower:
        return (ErrorCode.E_TIMEOUT, True)

    # Out of memory
    if "out of memory" in stderr_lower or ("memory" in stderr_lower and "error" in stderr_lower):
        return (ErrorCode.E_OOM, False)
    if "cannot allocate" in stderr_lower or "bad_alloc" in stderr_lower:
        return (ErrorCode.E_OOM, False)

    # Disk space
    if "no space left" in stderr_lower or "disk full" in stderr_lower:
        return (ErrorCode.E_DISK_FULL, False)

    # Network/rate limit
    if "rate limit" in stderr_lower or "too many requests" in stderr_lower:
        return (ErrorCode.E_NET_RATE_LIMIT, True)
    if "connection" in stderr_lower and ("refused" in stderr_lower or "timeout" in stderr_lower):
        return (ErrorCode.E_NET_RATE_LIMIT, True)

    # Input format errors
    if "format" in stderr_lower and "error" in stderr_lower:
        return (ErrorCode.E_INPUT_FORMAT, False)
    if "invalid" in stderr_lower and ("input" in stderr_lower or "file" in stderr_lower):
        return (ErrorCode.E_INPUT_FORMAT, False)

    # Input missing
    if "not found" in stderr_lower or "no such file" in stderr_lower:
        return (ErrorCode.E_INPUT_MISSING, False)

    # Default: non-zero exit
    return (ErrorCode.E_NONZERO_EXIT, False)
