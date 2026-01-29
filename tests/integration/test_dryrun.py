"""
Integration tests for dry-run mode functionality (Story 1.7).

Tests that dry-run mode correctly previews execution without modifying files.
Note: These tests verify the underlying mechanisms that support Snakemake's
--dry-run flag rather than testing Snakemake directly.
"""

import os
from pathlib import Path

import pytest

from workflow.lib.io import atomic_write, atomic_write_json, cleanup_temp_files
from workflow.lib.audit import (
    mark_rule_complete,
    is_rule_complete,
    get_checkpoint_path,
)
from workflow.lib.checksum import ChecksumCache, compute_file_checksum


# =============================================================================
# Test Dry-Run Simulation Support
# =============================================================================

class TestDryRunSimulationSupport:
    """Tests for dry-run mode simulation support."""

    def test_can_check_completion_without_modification(self, tmp_path: Path) -> None:
        """[P1] Test that checking rule completion doesn't modify state."""
        # Given
        meta_dir = tmp_path / "meta"
        rule = "test_rule"
        wildcards = {"key": "value"}

        # Take snapshot of directory state
        files_before = set()
        if meta_dir.exists():
            files_before = set(meta_dir.rglob("*"))

        # When - check completion (should be read-only)
        result = is_rule_complete(rule, wildcards, meta_dir)

        # Then - no files created
        files_after = set()
        if meta_dir.exists():
            files_after = set(meta_dir.rglob("*"))

        assert result is False
        assert files_before == files_after

    def test_checksum_computation_is_read_only(self, tmp_path: Path) -> None:
        """[P1] Test that checksum computation doesn't modify files."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("Test content")
        original_mtime = file_path.stat().st_mtime

        # When
        checksum = compute_file_checksum(file_path)

        # Then
        new_mtime = file_path.stat().st_mtime
        assert original_mtime == new_mtime
        assert checksum.startswith("sha256:")

    def test_cache_get_is_read_only(self, tmp_path: Path) -> None:
        """[P1] Test that cache get operation doesn't modify cache file."""
        # Given
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        cache_file = cache_dir / ".checksum_cache.json"

        # Create empty cache
        cache = ChecksumCache(cache_dir)
        cache.save()

        # Record cache file state
        original_size = cache_file.stat().st_size if cache_file.exists() else 0

        # When - only read operations
        cache2 = ChecksumCache(cache_dir)
        result = cache2.get(tmp_path / "nonexistent.txt")

        # Then - no modification
        assert result is None
        new_size = cache_file.stat().st_size if cache_file.exists() else 0
        assert original_size == new_size


# =============================================================================
# Test Output State Inspection
# =============================================================================

class TestOutputStateInspection:
    """Tests for inspecting output state without modification."""

    def test_determine_required_rules(self, tmp_path: Path) -> None:
        """[P1] Test determining which rules need to run (dry-run logic)."""
        # Given
        meta_dir = tmp_path / "meta"
        output_dir = tmp_path / "results"
        output_dir.mkdir()

        # Complete some rules
        mark_rule_complete("step1", {"species": "mmur"}, meta_dir)
        atomic_write(output_dir / "step1_mmur.txt", "done")

        mark_rule_complete("step2", {"species": "mmur"}, meta_dir)
        atomic_write(output_dir / "step2_mmur.txt", "done")

        # Step 3 not complete

        # When - check which rules need to run
        rules = ["step1", "step2", "step3"]
        wildcards = {"species": "mmur"}
        needs_run = [r for r in rules if not is_rule_complete(r, wildcards, meta_dir)]

        # Then
        assert needs_run == ["step3"]

    def test_inspect_outputs_without_execution(self, tmp_path: Path) -> None:
        """[P1] Test inspecting expected outputs without running rules."""
        # Given - simulate expected outputs for rules
        expected_outputs = {
            "standardize": [
                tmp_path / "results" / "standardized" / "mmur" / "genome.fa.gz",
                tmp_path / "results" / "standardized" / "mmur" / "annotation.gff3.gz",
            ],
            "qc_busco": [
                tmp_path / "results" / "qc" / "mmur" / "busco" / "summary.txt",
            ],
        }

        # Create only standardize outputs
        for output in expected_outputs["standardize"]:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text("mock content")

        # When - check what exists
        existing = {}
        missing = {}
        for rule, outputs in expected_outputs.items():
            existing[rule] = [o for o in outputs if o.exists()]
            missing[rule] = [o for o in outputs if not o.exists()]

        # Then
        assert len(existing["standardize"]) == 2
        assert len(missing["standardize"]) == 0
        assert len(existing["qc_busco"]) == 0
        assert len(missing["qc_busco"]) == 1

    def test_dry_run_leaves_no_artifacts(self, tmp_path: Path) -> None:
        """[P1] Test that dry-run analysis leaves no artifacts."""
        # Given
        meta_dir = tmp_path / "meta"
        output_dir = tmp_path / "results"

        # Record initial state
        def count_files(d: Path) -> int:
            if not d.exists():
                return 0
            return len(list(d.rglob("*")))

        initial_meta = count_files(meta_dir)
        initial_output = count_files(output_dir)

        # When - simulate dry-run operations (read-only checks)
        rules = ["step1", "step2", "step3"]
        for rule in rules:
            is_rule_complete(rule, {}, meta_dir)

        # Then - no new files created
        assert count_files(meta_dir) == initial_meta
        assert count_files(output_dir) == initial_output


# =============================================================================
# Test Execution Plan Generation
# =============================================================================

class TestExecutionPlanGeneration:
    """Tests for generating execution plans (dry-run output)."""

    def test_generate_execution_order(self, tmp_path: Path) -> None:
        """[P2] Test generating correct execution order for rules."""
        # Given - define rule dependencies
        rule_deps = {
            "standardize": [],
            "qc_busco": ["standardize"],
            "orthofinder": ["standardize"],
            "annotation": ["orthofinder"],
            "report": ["qc_busco", "annotation"],
        }

        meta_dir = tmp_path / "meta"

        # Complete some rules
        mark_rule_complete("standardize", {}, meta_dir)
        mark_rule_complete("qc_busco", {}, meta_dir)

        # When - find rules that need to run
        def get_pending_rules(rules: dict, completed: set) -> list:
            """Get rules to run in order based on dependencies."""
            pending = []
            for rule, deps in rules.items():
                if rule not in completed:
                    if all(d in completed for d in deps):
                        pending.append(rule)
            return pending

        completed = {"standardize", "qc_busco"}
        pending = get_pending_rules(rule_deps, completed)

        # Then - orthofinder should be next (depends only on standardize)
        assert "orthofinder" in pending
        assert "annotation" not in pending  # Depends on orthofinder
        assert "report" not in pending  # Depends on annotation

    def test_identify_all_remaining_rules(self, tmp_path: Path) -> None:
        """[P2] Test identifying all rules that would run."""
        # Given
        meta_dir = tmp_path / "meta"
        all_rules = ["step1", "step2", "step3", "step4", "step5"]

        # Complete first two
        mark_rule_complete("step1", {}, meta_dir)
        mark_rule_complete("step2", {}, meta_dir)

        # When
        would_run = [r for r in all_rules if not is_rule_complete(r, {}, meta_dir)]

        # Then
        assert would_run == ["step3", "step4", "step5"]
        assert len(would_run) == 3


# =============================================================================
# Test DAG Inspection Support
# =============================================================================

class TestDAGInspectionSupport:
    """Tests for DAG inspection support (--dag flag)."""

    def test_can_determine_rule_dependencies(self, tmp_path: Path) -> None:
        """[P2] Test that we can track rule dependencies for DAG visualization."""
        # Given - simulate tracking of rule inputs/outputs
        rule_io = {
            "standardize": {
                "inputs": ["raw_data/genome.fa", "raw_data/annotation.gff3"],
                "outputs": ["results/standardized/genome.fa.gz"],
            },
            "orthofinder": {
                "inputs": ["results/standardized/genome.fa.gz"],
                "outputs": ["results/orthology/orthogroups.tsv"],
            },
        }

        # When - find dependencies
        def get_dependent_rules(rule: str, rules: dict) -> list:
            """Find rules that depend on this rule's outputs."""
            outputs = set(rules.get(rule, {}).get("outputs", []))
            dependents = []
            for r, io in rules.items():
                if r != rule:
                    inputs = set(io.get("inputs", []))
                    if inputs & outputs:
                        dependents.append(r)
            return dependents

        deps = get_dependent_rules("standardize", rule_io)

        # Then
        assert "orthofinder" in deps
