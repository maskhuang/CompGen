"""
CompGene Shared Library

This package provides shared utilities for the CompGene pipeline.

Modules:
    errors: Error codes and exception classes
    io: Atomic file writing utilities
    checksum: File and directory checksum computation
    logging: Dual-format logging (text + JSON Lines)
    audit: Rule execution audit metadata collection
    runner: Unified adapter execution with retry and timeout handling
"""

# Error handling (Story 1.3)
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

# I/O utilities (Story 1.4)
from workflow.lib.io import (
    atomic_write,
    atomic_write_json,
    atomic_append,
)

# Checksum utilities (Story 1.4)
from workflow.lib.checksum import (
    compute_file_checksum,
    compute_dir_checksum,
    compute_checksums,
)

# Logging (Story 1.4)
from workflow.lib.logging import (
    DualLogger,
    get_log_paths,
    create_logger,
)

# Audit (Story 1.4)
from workflow.lib.audit import (
    RunMetadata,
    get_audit_path,
    collect_run_metadata,
    write_run_json,
    create_and_write_audit,
)

# Runner (Story 1.6) - Import directly from workflow.lib.runner to avoid circular import
# Usage: from workflow.lib.runner import AdapterRunner, RetryConfig

__all__ = [
    # errors
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
    # io
    "atomic_write",
    "atomic_write_json",
    "atomic_append",
    # checksum
    "compute_file_checksum",
    "compute_dir_checksum",
    "compute_checksums",
    # logging
    "DualLogger",
    "get_log_paths",
    "create_logger",
    # audit
    "RunMetadata",
    "get_audit_path",
    "collect_run_metadata",
    "write_run_json",
    "create_and_write_audit",
    # runner - import directly from workflow.lib.runner
]
