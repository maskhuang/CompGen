"""
CompGene Shared Library

This package provides shared utilities for the CompGene pipeline.
"""

from workflow.lib.errors import (
    ErrorCode,
    ERROR_RECOVERY,
    EXIT_CODES,
    EXIT_SUCCESS,
    EXIT_GENERAL_ERROR,
    EXIT_CONFIG_ERROR,
    CompGeneError,
    ConfigurationError,
    get_recovery,
    format_error_message,
    exit_with_error,
)

__all__ = [
    "ErrorCode",
    "ERROR_RECOVERY",
    "EXIT_CODES",
    "EXIT_SUCCESS",
    "EXIT_GENERAL_ERROR",
    "EXIT_CONFIG_ERROR",
    "CompGeneError",
    "ConfigurationError",
    "get_recovery",
    "format_error_message",
    "exit_with_error",
]
