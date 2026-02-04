"""
CompGene BUSCO Adapter - Quality Control Tool Integration.

This module provides the adapter for BUSCO (Benchmarking Universal Single-Copy
Orthologs) tool integration, enabling annotation completeness assessment.

Source: Story 2.4 BUSCO QC, ADR-002 Tool Adaptation Layer
"""

from pathlib import Path
import re
import subprocess

from .base import BaseAdapter, ToolSpec, RunResult, AdapterContext, classify_common_errors
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# BUSCO Summary Parsing Patterns
# =============================================================================

BUSCO_SUMMARY_PATTERN = re.compile(
    r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
)

LINEAGE_PATTERN = re.compile(r'The lineage dataset is: (\S+)')
VERSION_PATTERN = re.compile(r'BUSCO version is: (\S+)')


# =============================================================================
# Helper Functions
# =============================================================================

def parse_short_summary(content: str) -> dict:
    """
    Parse BUSCO short_summary.txt content.

    Extracts completeness statistics, lineage, and version information
    from BUSCO summary output.

    Args:
        content: Raw text content of short_summary.txt.

    Returns:
        Dictionary containing:
        - complete_pct: Complete BUSCOs percentage
        - single_copy_pct: Single-copy BUSCOs percentage
        - duplicated_pct: Duplicated BUSCOs percentage
        - fragmented_pct: Fragmented BUSCOs percentage
        - missing_pct: Missing BUSCOs percentage
        - total: Total BUSCO groups searched
        - lineage: Lineage dataset name (if found)
        - busco_version: BUSCO version (if found)

    Raises:
        ValueError: If BUSCO summary notation cannot be parsed.
    """
    match = BUSCO_SUMMARY_PATTERN.search(content)
    if not match:
        raise ValueError("Cannot parse BUSCO summary: notation not found")

    result = {
        "complete_pct": float(match.group(1)),
        "single_copy_pct": float(match.group(2)),
        "duplicated_pct": float(match.group(3)),
        "fragmented_pct": float(match.group(4)),
        "missing_pct": float(match.group(5)),
        "total": int(match.group(6)),
    }

    # Extract lineage
    lineage_match = LINEAGE_PATTERN.search(content)
    result["lineage"] = lineage_match.group(1) if lineage_match else "unknown"

    # Extract version
    version_match = VERSION_PATTERN.search(content)
    result["busco_version"] = version_match.group(1) if version_match else "unknown"

    return result


# =============================================================================
# BuscoAdapter Class
# =============================================================================

class BuscoAdapter(BaseAdapter):
    """
    Adapter for BUSCO (Benchmarking Universal Single-Copy Orthologs).

    BUSCO provides quantitative measures for the assessment of genome assembly,
    gene set, and transcriptome completeness based on evolutionarily-informed
    expectations of gene content.

    Supported Versions: 5.x, 6.x

    Example:
        >>> adapter = BuscoAdapter()
        >>> version = adapter.check_version()
        >>> print(f"BUSCO {version}")
    """

    SUPPORTED_VERSIONS = (5, 6)  # Support BUSCO 5.x and 6.x
    DEFAULT_TIMEOUT = 1800  # 30 minutes

    @property
    def spec(self) -> ToolSpec:
        """
        Get the BUSCO tool specification.

        Returns:
            ToolSpec with BUSCO version constraints and conda environment.
        """
        return ToolSpec(
            name="busco",
            min_version="5.0.0",
            max_version="6.99.99",
            conda_env="busco.yaml",
            description="BUSCO annotation completeness assessment",
        )

    def check_version(self) -> str:
        """
        Get the installed BUSCO version.

        Executes `busco --version` and parses the output to extract
        the version string. Validates that the version is supported (5.x or 6.x).

        Returns:
            Version string (e.g., "5.4.7", "6.0.0").

        Raises:
            CompGeneError: E_TOOL_NOT_FOUND if BUSCO is not installed or
                          version cannot be detected.
            CompGeneError: E_TOOL_VERSION if version is not 5.x or 6.x.
        """
        try:
            result = subprocess.run(
                ["busco", "--version"],
                capture_output=True,
                text=True,
                timeout=30,
            )
        except FileNotFoundError:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "BUSCO not found. Ensure BUSCO is installed and in PATH.",
            )
        except subprocess.TimeoutExpired:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "BUSCO version check timed out.",
            )

        # Parse version from stdout or stderr
        # BUSCO 5.x outputs: "BUSCO 5.4.7"
        # BUSCO 6.x outputs: "BUSCO 6.0.0"
        combined_output = result.stdout + result.stderr
        match = re.search(r'BUSCO (\d+\.\d+\.\d+)', combined_output)

        if not match:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Cannot parse BUSCO version from output.",
                details=combined_output[:200] if combined_output else "No output",
            )

        version = match.group(1)
        major = int(version.split('.')[0])

        if major not in self.SUPPORTED_VERSIONS:
            raise CompGeneError(
                ErrorCode.E_TOOL_VERSION,
                f"BUSCO {version} not supported. Requires version 5.x or 6.x.",
                details=f"Detected major version: {major}",
            )

        return version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        """
        Validate inputs before BUSCO execution.

        Checks:
        1. Protein sequence file exists
        2. BUSCO lineage is configured

        Args:
            ctx: Adapter context with inputs and configuration.

        Raises:
            CompGeneError: E_INPUT_MISSING if proteins file not found
                          or lineage not configured.
        """
        # Check protein file exists
        proteins_path = ctx.inputs.get("proteins")
        if proteins_path is None:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Protein input not specified in context.",
            )

        if not proteins_path.exists():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Protein file not found: {proteins_path}",
            )

        # Check lineage configuration
        busco_config = ctx.config.get("busco", {})
        lineage = busco_config.get("lineage", "")

        if not lineage:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "BUSCO lineage not specified in config. "
                "Add 'busco.lineage' to your configuration (e.g., 'primates_odb10').",
            )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        """
        Build the BUSCO command line arguments.

        Args:
            ctx: Adapter context with inputs, outputs, and configuration.

        Returns:
            List of command line arguments for BUSCO execution.
        """
        busco_config = ctx.config.get("busco", {})

        cmd = [
            "busco",
            "-i", str(ctx.inputs["proteins"]),
            "-l", busco_config["lineage"],
            "-o", ctx.wildcards["species"],
            "-m", busco_config.get("mode", "proteins"),
            "-c", str(ctx.threads),
            "-f",  # Force overwrite existing results
        ]

        # Optional: specify download path for lineage databases
        if download_path := busco_config.get("download_path"):
            cmd.extend(["--download_path", str(download_path)])

        # Optional: offline mode (don't download lineage)
        if busco_config.get("offline", False):
            cmd.append("--offline")

        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        """
        List expected output files after successful BUSCO execution.

        Per AC1, BUSCO outputs include:
        - short_summary.txt
        - full_table.tsv
        - missing_busco_list.tsv

        Args:
            ctx: Adapter context with output paths.

        Returns:
            List of expected output file paths.
        """
        outputs = []

        # Primary output: short_summary.txt
        if "summary" in ctx.outputs:
            outputs.append(ctx.outputs["summary"])

        # Full table with per-gene results
        if "full_table" in ctx.outputs:
            outputs.append(ctx.outputs["full_table"])

        # Missing BUSCO list
        if "missing_list" in ctx.outputs:
            outputs.append(ctx.outputs["missing_list"])

        return outputs

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        """
        Parse BUSCO outputs after successful execution.

        Extracts completeness statistics from short_summary.txt.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            RunResult with outputs, summary statistics, and optional metrics.
        """
        # Locate summary file
        if "summary" in ctx.outputs:
            summary_path = ctx.outputs["summary"]
        elif "output_dir" in ctx.outputs:
            summary_path = ctx.outputs["output_dir"] / "short_summary.txt"
        else:
            raise CompGeneError(
                ErrorCode.E_OUTPUT_MISSING,
                "Cannot determine BUSCO summary file path.",
            )

        # Parse summary content
        content = summary_path.read_text()
        summary = parse_short_summary(content)

        return RunResult(
            outputs={"summary": summary_path},
            summary=summary,
            metrics=None,
            warnings=None,
        )

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """
        Get the execution timeout in seconds.

        Returns the configured timeout or default (30 minutes).

        Args:
            ctx: Adapter context for calculating timeout.

        Returns:
            Timeout in seconds.
        """
        busco_config = ctx.config.get("busco", {})
        return busco_config.get("timeout", self.DEFAULT_TIMEOUT)

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str,
    ) -> tuple[ErrorCode, bool]:
        """
        Classify a BUSCO execution error.

        Analyzes return code and stderr to determine error type
        and whether retry is appropriate.

        Args:
            ctx: Adapter context.
            returncode: Process return code.
            stderr: Standard error output.

        Returns:
            Tuple of (ErrorCode, is_retryable).
        """
        stderr_lower = stderr.lower()

        # BUSCO-specific: lineage database not found
        if "lineage" in stderr_lower and "not found" in stderr_lower:
            return (ErrorCode.E_INPUT_MISSING, False)

        # Memory errors
        if "memory" in stderr_lower or "oom" in stderr_lower:
            return (ErrorCode.E_OOM, False)

        # SIGKILL typically indicates timeout
        if returncode == -9:
            return (ErrorCode.E_TIMEOUT, True)

        # Fall back to common error classification
        return classify_common_errors(returncode, stderr)
