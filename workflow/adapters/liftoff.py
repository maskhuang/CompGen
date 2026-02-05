"""
CompGene Liftoff Adapter - Annotation Lifting Tool Integration.

This module provides the adapter for Liftoff tool integration, enabling
annotation transfer from a reference genome to a target genome assembly.

Source: Story 5.1 Liftoff Adapter, ADR-002 Tool Adaptation Layer
"""

from pathlib import Path
import re
import subprocess

from .base import BaseAdapter, ToolSpec, RunResult, AdapterContext, classify_common_errors
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Liftoff Output Parsing
# =============================================================================

def parse_liftoff_stats(lifted_gff: Path, unmapped_file: Path) -> dict:
    """
    Parse Liftoff output files to extract statistics.

    Args:
        lifted_gff: Path to the lifted GFF3 output file.
        unmapped_file: Path to the unmapped features file.

    Returns:
        Dictionary containing:
        - lifted_genes: Number of successfully lifted genes
        - lifted_features: Number of all lifted features (genes, mRNAs, exons, etc.)
        - unmapped_genes: Number of unmapped genes
        - total_genes: Total genes processed
        - lift_rate: Gene-level success rate (0.0 to 1.0)
    """
    lifted_genes = 0
    lifted_features = 0
    unmapped_count = 0

    # Count lifted features from GFF - distinguish genes from other features
    if lifted_gff.exists():
        with open(lifted_gff, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    lifted_features += 1
                    # Count gene-level features for accurate lift_rate
                    parts = line.split('\t')
                    if len(parts) >= 3 and parts[2].lower() == 'gene':
                        lifted_genes += 1

    # Count unmapped features (unmapped file contains gene IDs)
    if unmapped_file.exists():
        with open(unmapped_file, 'r') as f:
            for line in f:
                if line.strip():
                    unmapped_count += 1

    total_genes = lifted_genes + unmapped_count
    lift_rate = lifted_genes / total_genes if total_genes > 0 else 0.0

    return {
        "lifted_genes": lifted_genes,
        "lifted_features": lifted_features,
        "unmapped_genes": unmapped_count,
        "total_genes": total_genes,
        "lift_rate": round(lift_rate, 4),
    }


def parse_liftoff_gff_coverage(lifted_gff: Path) -> dict:
    """
    Parse coverage and identity statistics from lifted GFF attributes.

    Liftoff adds coverage and sequence_ID attributes to lifted features.

    Args:
        lifted_gff: Path to the lifted GFF3 output file.

    Returns:
        Dictionary with coverage statistics:
        - mean_coverage: Average coverage of lifted features
        - min_coverage: Minimum coverage
        - max_coverage: Maximum coverage
        - feature_count: Number of features with coverage info
    """
    coverages = []

    if lifted_gff.exists():
        coverage_pattern = re.compile(r'coverage=([0-9.]+)')
        with open(lifted_gff, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                match = coverage_pattern.search(line)
                if match:
                    try:
                        coverages.append(float(match.group(1)))
                    except ValueError:
                        pass

    if not coverages:
        return {
            "mean_coverage": 0.0,
            "min_coverage": 0.0,
            "max_coverage": 0.0,
            "feature_count": 0,
        }

    return {
        "mean_coverage": round(sum(coverages) / len(coverages), 4),
        "min_coverage": round(min(coverages), 4),
        "max_coverage": round(max(coverages), 4),
        "feature_count": len(coverages),
    }


# =============================================================================
# LiftoffAdapter Class
# =============================================================================

class LiftoffAdapter(BaseAdapter):
    """
    Adapter for Liftoff annotation lifting tool.

    Liftoff is a tool that accurately maps annotations from one genome
    assembly to another. It is particularly useful for transferring
    annotations from a well-annotated reference genome to a newly
    assembled target genome.

    Supported Versions: 1.6.x

    Example:
        >>> adapter = LiftoffAdapter()
        >>> version = adapter.check_version()
        >>> print(f"Liftoff {version}")
    """

    SUPPORTED_MAJOR_VERSION = 1
    SUPPORTED_MINOR_MIN = 6
    DEFAULT_TIMEOUT = 3600  # 60 minutes

    @property
    def spec(self) -> ToolSpec:
        """
        Get the Liftoff tool specification.

        Returns:
            ToolSpec with Liftoff version constraints and conda environment.
        """
        return ToolSpec(
            name="liftoff",
            min_version="1.6.0",
            max_version="1.99.99",
            conda_env="liftoff.yaml",
            description="Liftoff annotation lifting tool",
        )

    def check_version(self) -> str:
        """
        Get the installed Liftoff version.

        Executes `liftoff --version` and parses the output to extract
        the version string. Validates that the version is supported (1.6.x).

        Returns:
            Version string (e.g., "1.6.3").

        Raises:
            CompGeneError: E_TOOL_NOT_FOUND if Liftoff is not installed or
                          version cannot be detected.
            CompGeneError: E_TOOL_VERSION if version is not 1.6.x.
        """
        try:
            result = subprocess.run(
                ["liftoff", "--version"],
                capture_output=True,
                text=True,
                timeout=30,
            )
        except FileNotFoundError:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Liftoff not found. Ensure Liftoff is installed and in PATH.",
            )
        except subprocess.TimeoutExpired:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Liftoff version check timed out.",
            )

        # Parse version from stdout or stderr
        # Liftoff outputs: "liftoff 1.6.3"
        combined_output = result.stdout + result.stderr
        match = re.search(r'liftoff\s+(\d+\.\d+\.\d+)', combined_output, re.IGNORECASE)

        if not match:
            # Try alternate format without "liftoff" prefix
            match = re.search(r'^(\d+\.\d+\.\d+)', combined_output.strip())

        if not match:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Cannot parse Liftoff version from output.",
                details=combined_output[:200] if combined_output else "No output",
            )

        version = match.group(1)
        parts = version.split('.')
        major, minor = int(parts[0]), int(parts[1])

        if major != self.SUPPORTED_MAJOR_VERSION or minor < self.SUPPORTED_MINOR_MIN:
            raise CompGeneError(
                ErrorCode.E_TOOL_VERSION,
                f"Liftoff {version} not supported. Requires version 1.6.x or higher.",
                details=f"Detected version: {version}",
            )

        return version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        """
        Validate inputs before Liftoff execution.

        Checks:
        1. Reference GFF file exists
        2. Reference FASTA file exists
        3. Target FASTA file exists

        Args:
            ctx: Adapter context with inputs and configuration.

        Raises:
            CompGeneError: E_INPUT_MISSING if any required input is missing.
        """
        # Check reference GFF
        ref_gff = ctx.inputs.get("reference_gff")
        if ref_gff is None:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Reference GFF not specified in context.",
            )
        if not ref_gff.exists():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Reference GFF file not found: {ref_gff}",
            )

        # Check reference FASTA
        ref_fa = ctx.inputs.get("reference_fa")
        if ref_fa is None:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Reference FASTA not specified in context.",
            )
        if not ref_fa.exists():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Reference FASTA file not found: {ref_fa}",
            )

        # Check target FASTA
        target_fa = ctx.inputs.get("target_fa")
        if target_fa is None:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Target FASTA not specified in context.",
            )
        if not target_fa.exists():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Target FASTA file not found: {target_fa}",
            )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        """
        Build the Liftoff command line arguments.

        Args:
            ctx: Adapter context with inputs, outputs, and configuration.

        Returns:
            List of command line arguments for Liftoff execution.
        """
        liftoff_config = ctx.config.get("liftoff", {})

        cmd = [
            "liftoff",
            "-g", str(ctx.inputs["reference_gff"]),
            "-o", str(ctx.outputs["lifted_gff"]),
            "-u", str(ctx.outputs["unmapped"]),
            "-dir", str(ctx.outputs["intermediate_dir"]),
            "-p", str(ctx.threads),
        ]

        # Optional: coverage threshold
        min_coverage = liftoff_config.get("min_coverage")
        if min_coverage is not None:
            cmd.extend(["-a", str(min_coverage)])

        # Optional: sequence identity threshold
        min_identity = liftoff_config.get("min_identity")
        if min_identity is not None:
            cmd.extend(["-s", str(min_identity)])

        # Optional: allow multiple copies
        if liftoff_config.get("copies", False):
            cmd.append("-copies")

        # Optional: flank sequence
        flank = liftoff_config.get("flank")
        if flank is not None and flank > 0:
            cmd.extend(["-flank", str(flank)])

        # Positional arguments: target genome, then reference genome
        cmd.append(str(ctx.inputs["target_fa"]))
        cmd.append(str(ctx.inputs["reference_fa"]))

        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        """
        List expected output files after successful Liftoff execution.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            List of expected output file paths.
        """
        outputs = []

        # Primary output: lifted GFF3
        if "lifted_gff" in ctx.outputs:
            outputs.append(ctx.outputs["lifted_gff"])

        # Unmapped features list (may be empty but should exist)
        if "unmapped" in ctx.outputs:
            outputs.append(ctx.outputs["unmapped"])

        # Statistics file (generated by run_liftoff.py script)
        if "stats" in ctx.outputs:
            outputs.append(ctx.outputs["stats"])

        return outputs

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        """
        Parse Liftoff outputs after successful execution.

        Extracts statistics from the lifted GFF and unmapped features file.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            RunResult with outputs, summary statistics, and optional metrics.
        """
        lifted_gff = ctx.outputs.get("lifted_gff")
        unmapped = ctx.outputs.get("unmapped")

        if lifted_gff is None or unmapped is None:
            raise CompGeneError(
                ErrorCode.E_OUTPUT_MISSING,
                "Cannot determine Liftoff output paths.",
            )

        # Parse basic statistics
        stats = parse_liftoff_stats(lifted_gff, unmapped)

        # Parse coverage statistics from GFF
        coverage_stats = parse_liftoff_gff_coverage(lifted_gff)

        # Combine statistics
        summary = {
            **stats,
            "coverage_stats": coverage_stats,
        }

        return RunResult(
            outputs={
                "lifted_gff": lifted_gff,
                "unmapped": unmapped,
            },
            summary=summary,
            metrics=None,
            warnings=None,
        )

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """
        Get the execution timeout in seconds.

        Returns the configured timeout or default (60 minutes).

        Args:
            ctx: Adapter context for calculating timeout.

        Returns:
            Timeout in seconds.
        """
        liftoff_config = ctx.config.get("liftoff", {})
        return liftoff_config.get("timeout", self.DEFAULT_TIMEOUT)

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str,
    ) -> tuple[ErrorCode, bool]:
        """
        Classify a Liftoff execution error.

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

        # Liftoff-specific: input file not found
        if "no such file" in stderr_lower or "not found" in stderr_lower:
            return (ErrorCode.E_INPUT_MISSING, False)

        # GFF format errors
        if "gff" in stderr_lower and ("error" in stderr_lower or "invalid" in stderr_lower):
            return (ErrorCode.E_INPUT_FORMAT, False)

        # FASTA format errors
        if "fasta" in stderr_lower and ("error" in stderr_lower or "invalid" in stderr_lower):
            return (ErrorCode.E_INPUT_FORMAT, False)

        # Memory errors
        if "memory" in stderr_lower or "oom" in stderr_lower:
            return (ErrorCode.E_OOM, False)
        if "cannot allocate" in stderr_lower or "bad_alloc" in stderr_lower:
            return (ErrorCode.E_OOM, False)

        # SIGKILL typically indicates timeout or OOM killer
        if returncode == -9:
            return (ErrorCode.E_TIMEOUT, True)

        # Fall back to common error classification
        return classify_common_errors(returncode, stderr)
