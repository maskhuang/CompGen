"""
CompGene OrthoFinder Adapter - Orthology Inference Tool Integration.

This module provides the adapter for OrthoFinder tool integration,
enabling orthogroup inference across multiple species.

Source: Story 3.1 OrthoFinder Adapter, ADR-002 Tool Adaptation Layer
"""

from pathlib import Path
import re
import subprocess

from .base import BaseAdapter, ToolSpec, RunResult, AdapterContext, classify_common_errors
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# OrthoFinder Output Parsing
# =============================================================================

def parse_gene_count_tsv(gene_count_path: Path) -> dict:
    """
    Parse Orthogroups.GeneCount.tsv to extract statistics.

    Args:
        gene_count_path: Path to Orthogroups.GeneCount.tsv.

    Returns:
        Dictionary with:
        - total_orthogroups: Number of orthogroups
        - species_names: List of species names from header
        - single_copy_count: Orthogroups with exactly 1 gene per species

    Raises:
        ValueError: If file cannot be parsed.
    """
    if not gene_count_path.exists():
        return {
            "total_orthogroups": 0,
            "species_names": [],
            "single_copy_count": 0,
        }

    total_ogs = 0
    single_copy = 0
    species_names = []

    with open(gene_count_path) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if i == 0:
                # Header: Orthogroup, species1, species2, ..., Total
                species_names = parts[1:-1]  # Exclude Orthogroup and Total columns
                continue

            total_ogs += 1

            # Check if single-copy: all species have exactly 1 gene
            if species_names:
                counts = parts[1:-1]  # Exclude Orthogroup and Total columns
                if all(c.strip() == '1' for c in counts if c.strip()):
                    single_copy += 1

    return {
        "total_orthogroups": total_ogs,
        "species_names": species_names,
        "single_copy_count": single_copy,
    }


def count_unassigned_genes(unassigned_path: Path) -> int:
    """
    Count unassigned genes from Orthogroups_UnassignedGenes.tsv.

    Args:
        unassigned_path: Path to Orthogroups_UnassignedGenes.tsv.

    Returns:
        Number of unassigned gene entries (rows minus header).
    """
    if not unassigned_path.exists():
        return 0

    count = 0
    with open(unassigned_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # Skip header
            if line.strip():
                count += 1
    return count


# =============================================================================
# OrthoFinderAdapter Class
# =============================================================================

class OrthoFinderAdapter(BaseAdapter):
    """
    Adapter for OrthoFinder orthology inference tool.

    OrthoFinder infers orthogroups (gene families) from protein sequences
    of multiple species using an all-vs-all DIAMOND/BLAST search followed
    by MCL clustering.

    Supported Versions: 2.5.x

    Example:
        >>> adapter = OrthoFinderAdapter()
        >>> version = adapter.check_version()
        >>> print(f"OrthoFinder {version}")
    """

    SUPPORTED_MAJOR_VERSION = 2
    SUPPORTED_MINOR_MIN = 5
    DEFAULT_TIMEOUT = 14400  # 4 hours

    @property
    def spec(self) -> ToolSpec:
        """
        Get the OrthoFinder tool specification.

        Returns:
            ToolSpec with OrthoFinder version constraints and conda environment.
        """
        return ToolSpec(
            name="orthofinder",
            min_version="2.5.0",
            max_version="2.5.99",
            conda_env="orthofinder.yaml",
            description="OrthoFinder orthology inference tool",
        )

    def check_version(self) -> str:
        """
        Get the installed OrthoFinder version.

        Executes `orthofinder --version` and parses the output to extract
        the version string. Validates that the version is 2.5.x.

        Returns:
            Version string (e.g., "2.5.5").

        Raises:
            CompGeneError: E_TOOL_NOT_FOUND if OrthoFinder is not installed.
            CompGeneError: E_TOOL_VERSION if version is not 2.5.x.
        """
        try:
            result = subprocess.run(
                ["orthofinder", "--version"],
                capture_output=True,
                text=True,
                timeout=30,
            )
        except FileNotFoundError:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "OrthoFinder not found. Ensure OrthoFinder is installed and in PATH.",
            )
        except subprocess.TimeoutExpired:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "OrthoFinder version check timed out.",
            )

        # OrthoFinder version 2.5.5
        combined_output = result.stdout + result.stderr
        match = re.search(r'OrthoFinder version (\d+\.\d+\.\d+)', combined_output)

        if not match:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Cannot parse OrthoFinder version from output.",
                details=combined_output[:200] if combined_output else "No output",
            )

        version = match.group(1)
        parts = version.split('.')
        major, minor = int(parts[0]), int(parts[1])

        if major != self.SUPPORTED_MAJOR_VERSION or minor < self.SUPPORTED_MINOR_MIN:
            raise CompGeneError(
                ErrorCode.E_TOOL_VERSION,
                f"OrthoFinder {version} not supported. Requires version 2.5.x.",
                details=f"Detected: {major}.{minor}",
            )

        return version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        """
        Validate inputs before OrthoFinder execution.

        Checks:
        1. Proteins directory exists
        2. Directory contains at least 2 .fa files (minimum for comparison)

        Args:
            ctx: Adapter context with inputs and configuration.

        Raises:
            CompGeneError: E_INPUT_MISSING if directory not found or
                          insufficient species files.
        """
        proteins_dir = ctx.inputs.get("proteins_dir")
        if proteins_dir is None:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "Proteins directory not specified in context.",
            )

        if not proteins_dir.is_dir():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Proteins directory not found: {proteins_dir}",
            )

        fa_files = list(proteins_dir.glob("*.fa"))
        if len(fa_files) < 2:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Need at least 2 species .fa files in {proteins_dir}, found {len(fa_files)}. "
                "OrthoFinder requires protein sequences from multiple species.",
            )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        """
        Build the OrthoFinder command line arguments.

        Args:
            ctx: Adapter context with inputs, outputs, and configuration.

        Returns:
            List of command line arguments for OrthoFinder execution.
        """
        of_config = ctx.config.get("orthofinder", {})

        cmd = [
            "orthofinder",
            "-f", str(ctx.inputs["proteins_dir"].resolve()),
            "-t", str(ctx.threads),
            "-a", str(min(ctx.threads, 4)),
        ]

        # Output directory
        if "output_dir" in ctx.outputs:
            cmd.extend(["-o", str(ctx.outputs["output_dir"].resolve())])

        # Search method (diamond, diamond_ultra_sens, blast)
        search_method = of_config.get("search_method", "diamond")
        cmd.extend(["-S", search_method])

        # Fixed result suffix for deterministic paths
        cmd.extend(["-n", "compgene"])

        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        """
        List expected output files after successful OrthoFinder execution.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            List of expected output file paths.
        """
        results_dir = ctx.outputs["results_dir"]
        return [
            results_dir / "Orthogroups" / "Orthogroups.tsv",
            results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv",
            results_dir / "Species_Tree" / "SpeciesTree_rooted.txt",
        ]

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        """
        Parse OrthoFinder outputs after successful execution.

        Extracts orthogroup statistics from output files.

        Args:
            ctx: Adapter context with output paths.

        Returns:
            RunResult with outputs, summary statistics, and metrics.
        """
        results_dir = ctx.outputs["results_dir"]

        orthogroups_tsv = results_dir / "Orthogroups" / "Orthogroups.tsv"
        gene_count_tsv = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
        unassigned_tsv = results_dir / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
        species_tree = results_dir / "Species_Tree" / "SpeciesTree_rooted.txt"

        # Parse gene count statistics
        gc_stats = parse_gene_count_tsv(gene_count_tsv)

        # Count unassigned genes
        unassigned_count = count_unassigned_genes(unassigned_tsv)

        summary = {
            "total_orthogroups": gc_stats["total_orthogroups"],
            "species_count": len(gc_stats["species_names"]),
            "species_names": gc_stats["species_names"],
            "single_copy_orthogroups": gc_stats["single_copy_count"],
            "unassigned_genes": unassigned_count,
        }

        outputs = {
            "orthogroups_tsv": orthogroups_tsv,
            "gene_count_tsv": gene_count_tsv,
            "species_tree": species_tree,
        }

        return RunResult(
            outputs=outputs,
            summary=summary,
        )

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """
        Get the execution timeout in seconds.

        OrthoFinder can be slow for large species sets:
        - 3 species: ~1-2 hours
        - 7+ species: several hours

        Default: 4 hours (14400 seconds).

        Args:
            ctx: Adapter context for calculating timeout.

        Returns:
            Timeout in seconds.
        """
        of_config = ctx.config.get("orthofinder", {})
        return of_config.get("timeout", self.DEFAULT_TIMEOUT)

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str,
    ) -> tuple[ErrorCode, bool]:
        """
        Classify an OrthoFinder execution error.

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

        # DIAMOND/BLAST search errors
        if "diamond" in stderr_lower and "error" in stderr_lower:
            return (ErrorCode.E_NONZERO_EXIT, False)

        # No input sequences
        if "no fasta files" in stderr_lower or "no sequences" in stderr_lower:
            return (ErrorCode.E_INPUT_MISSING, False)

        # MCL clustering errors
        if "mcl" in stderr_lower and "error" in stderr_lower:
            return (ErrorCode.E_NONZERO_EXIT, False)

        # Memory errors
        if "memory" in stderr_lower or "oom" in stderr_lower:
            return (ErrorCode.E_OOM, False)

        # SIGKILL typically indicates timeout or OOM kill
        if returncode == -9:
            return (ErrorCode.E_TIMEOUT, True)

        # Fall back to common error classification
        return classify_common_errors(returncode, stderr)
