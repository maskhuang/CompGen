"""
CompGene Error Code System.

This module provides standardized error codes, recovery suggestions,
and a custom exception class for the CompGene pipeline.

Source: ADR-003 Error Handling and Recovery
"""

from enum import Enum
import sys
from typing import NoReturn, Optional


class ErrorCode(str, Enum):
    """
    Standardized error codes for the CompGene pipeline.

    Inherits from str and Enum for JSON serialization compatibility.

    Error Code Categories:
    - E_INPUT_*: Input-related errors (missing files, format issues)
    - E_TOOL_*: External tool errors (not found, version mismatch)
    - E_*: Resource and runtime errors (timeout, memory, disk, network)
    """

    E_INPUT_MISSING = "E_INPUT_MISSING"
    """Input file or directory not found. Check the file path."""

    E_INPUT_FORMAT = "E_INPUT_FORMAT"
    """Input file format error. Validate GFF/FASTA format."""

    E_TOOL_NOT_FOUND = "E_TOOL_NOT_FOUND"
    """Required external tool not found. Activate conda environment."""

    E_TOOL_VERSION = "E_TOOL_VERSION"
    """Tool version mismatch. Update to required version."""

    E_TIMEOUT = "E_TIMEOUT"
    """Operation timed out. Increase timeout or reduce input size."""

    E_OOM = "E_OOM"
    """Out of memory error. Reduce threads or increase available memory."""

    E_NET_RATE_LIMIT = "E_NET_RATE_LIMIT"
    """Network rate limit exceeded. Wait and retry automatically."""

    E_DISK_FULL = "E_DISK_FULL"
    """Disk space exhausted. Free up disk space before retrying."""

    E_NONZERO_EXIT = "E_NONZERO_EXIT"
    """Tool returned non-zero exit code. Check logs for details."""


# =============================================================================
# Error Recovery Mapping
# =============================================================================

ERROR_RECOVERY: dict[ErrorCode, tuple[bool, str]] = {
    ErrorCode.E_INPUT_MISSING: (False, "检查输入文件路径是否正确"),
    ErrorCode.E_INPUT_FORMAT: (False, "验证 GFF/FASTA 格式是否符合规范"),
    ErrorCode.E_TOOL_NOT_FOUND: (False, "激活对应的 conda 环境"),
    ErrorCode.E_TOOL_VERSION: (False, "更新工具版本至要求范围"),
    ErrorCode.E_TIMEOUT: (True, "增加超时时间或减少输入数据规模"),
    ErrorCode.E_OOM: (False, "减少 threads 或增加可用内存"),
    ErrorCode.E_NET_RATE_LIMIT: (True, "等待后自动重试"),
    ErrorCode.E_DISK_FULL: (False, "清理磁盘空间后重试"),
    ErrorCode.E_NONZERO_EXIT: (False, "查看日志了解详细错误"),
}


def get_recovery(code: ErrorCode) -> tuple[bool, str]:
    """
    Get recovery information for an error code.

    Args:
        code: The ErrorCode to look up.

    Returns:
        Tuple of (is_retryable, recovery_suggestion).
        - is_retryable: True if the operation can be automatically retried.
        - recovery_suggestion: Human-readable suggestion in Chinese.
    """
    return ERROR_RECOVERY.get(code, (False, "未知错误，请查看日志"))


# =============================================================================
# Exit Code Mapping
# =============================================================================

# Exit codes follow Unix conventions and FR49 requirements
EXIT_CODES: dict[ErrorCode, int] = {
    ErrorCode.E_INPUT_MISSING: 3,    # Input error
    ErrorCode.E_INPUT_FORMAT: 3,     # Input error
    ErrorCode.E_TOOL_NOT_FOUND: 4,   # Tool error
    ErrorCode.E_TOOL_VERSION: 4,     # Tool error
    ErrorCode.E_TIMEOUT: 5,          # Resource error
    ErrorCode.E_OOM: 5,              # Resource error
    ErrorCode.E_DISK_FULL: 5,        # Resource error
    ErrorCode.E_NET_RATE_LIMIT: 6,   # Network error
    ErrorCode.E_NONZERO_EXIT: 1,     # General error
}

# Special exit codes not tied to ErrorCode
EXIT_SUCCESS = 0
EXIT_GENERAL_ERROR = 1
EXIT_CONFIG_ERROR = 2


# =============================================================================
# CompGeneError Exception Class
# =============================================================================

class CompGeneError(Exception):
    """
    Base exception class for CompGene pipeline errors.

    Provides standardized error handling with error codes, recovery suggestions,
    and process exit code mapping.

    Attributes:
        error_code: The ErrorCode enum value identifying the error type.
        message: Human-readable error description.
        details: Optional additional error details (e.g., file path, tool output).
        is_retryable: Whether the operation can be automatically retried.
        recovery_suggestion: Human-readable suggestion for resolving the error.

    Example:
        >>> raise CompGeneError(
        ...     ErrorCode.E_INPUT_MISSING,
        ...     "Input file not found",
        ...     details="/path/to/missing/file.fa"
        ... )
    """

    def __init__(
        self,
        error_code: ErrorCode,
        message: str,
        details: Optional[str] = None
    ):
        """
        Initialize a CompGeneError.

        Args:
            error_code: The ErrorCode identifying the error type.
            message: Human-readable error description.
            details: Optional additional context (file paths, tool output, etc.).
        """
        self.error_code = error_code
        self.message = message
        self.details = details

        # Look up recovery information
        self.is_retryable, self.recovery_suggestion = get_recovery(error_code)

        super().__init__(self._format_message())

    def _format_message(self) -> str:
        """Format the exception message."""
        return f"[{self.error_code.value}] {self.message}"

    def __str__(self) -> str:
        """Return formatted error string."""
        return self._format_message()

    def to_exit_code(self) -> int:
        """
        Get the process exit code for this error.

        Returns:
            Integer exit code following Unix conventions:
            - 0: Success (not applicable for errors)
            - 1: General error
            - 2: Configuration error
            - 3: Input error
            - 4: Tool error
            - 5: Resource error
            - 6: Network error
        """
        return EXIT_CODES.get(self.error_code, EXIT_GENERAL_ERROR)

    def __reduce__(self):
        """Support pickle serialization for multiprocessing compatibility."""
        return (
            self.__class__,
            (self.error_code, self.message, self.details)
        )


# =============================================================================
# Error Formatting and Output
# =============================================================================

def format_error_message(error: CompGeneError, use_color: bool = True) -> str:
    """
    Format a CompGeneError for terminal output.

    Args:
        error: The CompGeneError to format.
        use_color: Whether to use ANSI color codes (default: True).

    Returns:
        Formatted error message string with error code, message,
        recovery suggestion, and optional details.
    """
    # ANSI color codes
    RED = "\033[91m" if use_color else ""
    YELLOW = "\033[93m" if use_color else ""
    CYAN = "\033[96m" if use_color else ""
    RESET = "\033[0m" if use_color else ""
    BOLD = "\033[1m" if use_color else ""

    lines = []

    # Error header
    lines.append(f"{RED}{BOLD}ERROR [{error.error_code.value}]{RESET}")

    # Error message
    lines.append(f"  {error.message}")

    # Details if present
    if error.details:
        lines.append(f"  {CYAN}Details:{RESET} {error.details}")

    # Recovery suggestion
    retry_indicator = "[可重试]" if error.is_retryable else "[不可重试]"
    lines.append(f"  {YELLOW}Recovery:{RESET} {error.recovery_suggestion} {retry_indicator}")

    return "\n".join(lines)


def exit_with_error(error: CompGeneError, use_color: bool = True) -> NoReturn:
    """
    Print formatted error message and exit with appropriate exit code.

    This function prints the error to stderr and calls sys.exit()
    with the error's exit code. This function never returns.

    Args:
        error: The CompGeneError to report.
        use_color: Whether to use ANSI color codes (default: True).

    Raises:
        SystemExit: Always raised with the error's exit code.
    """
    print(format_error_message(error, use_color=use_color), file=sys.stderr)
    sys.exit(error.to_exit_code())


# =============================================================================
# Configuration Error (Special Case for Exit Code 2)
# =============================================================================

class ConfigurationError(CompGeneError):
    """
    Special exception for configuration validation errors.

    Always returns exit code 2 regardless of the underlying error code.
    Used when config validation fails (from Story 1.2).

    Note:
        This class uses E_INPUT_FORMAT as the internal error_code for
        recovery suggestion lookup, but to_exit_code() is overridden
        to always return EXIT_CONFIG_ERROR (2) as required by FR49.
    """

    def __init__(self, message: str, details: Optional[str] = None):
        """
        Initialize a ConfigurationError.

        Args:
            message: Description of the configuration error.
            details: Optional details about which config field is invalid.
        """
        # Use E_INPUT_FORMAT as the base error code, but override exit code
        super().__init__(ErrorCode.E_INPUT_FORMAT, message, details)

    def to_exit_code(self) -> int:
        """Configuration errors always return exit code 2."""
        return EXIT_CONFIG_ERROR

    def __reduce__(self):
        """Support pickle serialization for multiprocessing compatibility."""
        return (
            self.__class__,
            (self.message, self.details)
        )
