"""
CompGene Audit Metadata Collection.

This module provides functionality for collecting and writing
audit metadata for rule executions.

Source: ADR-004 Audit Granularity (Rule-level)
"""

import json
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from workflow.lib.io import atomic_write_json
from workflow.lib.checksum import compute_checksums


# =============================================================================
# RunMetadata Data Class
# =============================================================================

@dataclass
class RunMetadata:
    """
    Audit metadata for a rule execution.

    Contains all information needed to reproduce and trace
    the execution of a Snakemake rule.

    Attributes:
        rule: The rule name.
        wildcards: Dictionary of wildcard values.
        cmd: Command executed as list of strings.
        tool_version: Version string of the external tool.
        input_checksums: Checksums of input files/directories.
        threads: Number of threads used.
        runtime_seconds: Execution time in seconds.
        exit_code: Process exit code (0 = success).
        timestamp: ISO 8601 timestamp of completion.
        output_checksums: Optional checksums of output files.
        error_code: Optional ErrorCode value if error occurred.
        error_message: Optional error message if error occurred.
    """
    rule: str
    wildcards: dict[str, str]
    cmd: list[str]
    tool_version: str
    input_checksums: dict[str, str]
    threads: int
    runtime_seconds: float
    exit_code: int
    timestamp: str
    output_checksums: Optional[dict[str, str]] = None
    error_code: Optional[str] = None
    error_message: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        """
        Convert to a JSON-serializable dictionary.

        Excludes None values for cleaner output.

        Returns:
            Dictionary with non-None fields only.
        """
        return {k: v for k, v in asdict(self).items() if v is not None}

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
# Metadata Collection Functions
# =============================================================================

def get_audit_path(
    rule: str,
    wildcards: dict[str, str],
    meta_dir: Path = Path("results/meta")
) -> Path:
    """
    Generate the audit file path for a rule execution.

    Args:
        rule: The rule name.
        wildcards: Dictionary of wildcard values.
        meta_dir: Base directory for metadata (default: "results/meta").

    Returns:
        Path to the .run.json file.

    Example:
        >>> path = get_audit_path("orthofinder", {"species_set": "lemur"})
        >>> print(path)
        results/meta/orthofinder/species_set=lemur.run.json
    """
    if wildcards:
        wildcard_str = "_".join(f"{k}={v}" for k, v in sorted(wildcards.items()))
    else:
        wildcard_str = "default"

    return meta_dir / rule / f"{wildcard_str}.run.json"


def collect_run_metadata(
    rule: str,
    wildcards: dict[str, str],
    cmd: list[str],
    tool_version: str,
    input_paths: dict[str, Path],
    threads: int,
    runtime_seconds: float,
    exit_code: int,
    output_paths: Optional[dict[str, Path]] = None,
    error_code: Optional[str] = None,
    error_message: Optional[str] = None
) -> RunMetadata:
    """
    Collect all metadata for a rule execution.

    Computes input checksums and optionally output checksums.

    Args:
        rule: The rule name.
        wildcards: Dictionary of wildcard values.
        cmd: Command executed as list of strings.
        tool_version: Version string of the external tool.
        input_paths: Dictionary mapping names to input file/dir paths.
        threads: Number of threads used.
        runtime_seconds: Execution time in seconds.
        exit_code: Process exit code.
        output_paths: Optional dictionary mapping names to output paths.
        error_code: Optional ErrorCode value if error occurred.
        error_message: Optional error message if error occurred.

    Returns:
        Populated RunMetadata instance.

    Example:
        >>> metadata = collect_run_metadata(
        ...     rule="orthofinder",
        ...     wildcards={"species_set": "lemur"},
        ...     cmd=["orthofinder", "-f", "proteins/"],
        ...     tool_version="2.5.5",
        ...     input_paths={"proteins": Path("proteins/")},
        ...     threads=8,
        ...     runtime_seconds=3600.5,
        ...     exit_code=0
        ... )
    """
    # Generate timestamp
    timestamp = datetime.now(timezone.utc).isoformat()

    # Compute input checksums
    input_checksums = compute_checksums(input_paths)

    # Compute output checksums if provided
    output_checksums = None
    if output_paths:
        output_checksums = compute_checksums(output_paths)

    return RunMetadata(
        rule=rule,
        wildcards=wildcards,
        cmd=cmd,
        tool_version=tool_version,
        input_checksums=input_checksums,
        threads=threads,
        runtime_seconds=runtime_seconds,
        exit_code=exit_code,
        timestamp=timestamp,
        output_checksums=output_checksums,
        error_code=error_code,
        error_message=error_message
    )


def write_run_json(
    metadata: RunMetadata,
    meta_dir: Path = Path("results/meta")
) -> Path:
    """
    Write audit metadata to a .run.json file.

    Args:
        metadata: The RunMetadata to write.
        meta_dir: Base directory for metadata (default: "results/meta").

    Returns:
        Path to the written file.

    Example:
        >>> path = write_run_json(metadata)
        >>> print(f"Audit written to: {path}")
    """
    audit_path = get_audit_path(metadata.rule, metadata.wildcards, meta_dir)
    atomic_write_json(audit_path, metadata.to_dict())
    return audit_path


# =============================================================================
# Convenience Functions
# =============================================================================

def create_and_write_audit(
    rule: str,
    wildcards: dict[str, str],
    cmd: list[str],
    tool_version: str,
    input_paths: dict[str, Path],
    threads: int,
    runtime_seconds: float,
    exit_code: int,
    meta_dir: Path = Path("results/meta"),
    output_paths: Optional[dict[str, Path]] = None,
    error_code: Optional[str] = None,
    error_message: Optional[str] = None
) -> tuple[RunMetadata, Path]:
    """
    Collect metadata and write to file in one step.

    Convenience function that combines collect_run_metadata and write_run_json.

    Args:
        rule: The rule name.
        wildcards: Dictionary of wildcard values.
        cmd: Command executed as list of strings.
        tool_version: Version string of the external tool.
        input_paths: Dictionary mapping names to input file/dir paths.
        threads: Number of threads used.
        runtime_seconds: Execution time in seconds.
        exit_code: Process exit code.
        meta_dir: Base directory for metadata.
        output_paths: Optional dictionary mapping names to output paths.
        error_code: Optional ErrorCode value if error occurred.
        error_message: Optional error message if error occurred.

    Returns:
        Tuple of (RunMetadata, Path to written file).
    """
    metadata = collect_run_metadata(
        rule=rule,
        wildcards=wildcards,
        cmd=cmd,
        tool_version=tool_version,
        input_paths=input_paths,
        threads=threads,
        runtime_seconds=runtime_seconds,
        exit_code=exit_code,
        output_paths=output_paths,
        error_code=error_code,
        error_message=error_message
    )

    audit_path = write_run_json(metadata, meta_dir)

    return metadata, audit_path
