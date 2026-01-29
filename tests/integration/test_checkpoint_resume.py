"""
Integration tests for checkpoint and resume functionality (Story 1.7).

Tests the --rerun-incomplete behavior and checkpoint marker system.
"""

import os
import time
from pathlib import Path

import pytest

from workflow.lib.io import atomic_write, cleanup_temp_files, get_temp_path
from workflow.lib.audit import (
    mark_rule_complete,
    is_rule_complete,
    clear_checkpoint,
    get_checkpoint_path,
)
from workflow.lib.checksum import ChecksumCache


# =============================================================================
# Test Checkpoint Marker Integration
# =============================================================================

class TestCheckpointMarkerIntegration:
    """Integration tests for checkpoint markers."""

    def test_full_checkpoint_lifecycle(self, tmp_path: Path) -> None:
        """[P1] Test complete checkpoint lifecycle: create, check, clear."""
        # Given
        rule = "orthofinder"
        wildcards = {"species_set": "lemur_macaque"}
        meta_dir = tmp_path / "meta"

        # Initially not complete
        assert not is_rule_complete(rule, wildcards, meta_dir)

        # When - mark complete
        path = mark_rule_complete(rule, wildcards, meta_dir)

        # Then - verify complete
        assert path.exists()
        assert is_rule_complete(rule, wildcards, meta_dir)

        # When - clear checkpoint
        result = clear_checkpoint(rule, wildcards, meta_dir)

        # Then - back to incomplete
        assert result is True
        assert not is_rule_complete(rule, wildcards, meta_dir)

    def test_multiple_rules_independent_checkpoints(self, tmp_path: Path) -> None:
        """[P1] Test that checkpoints for different rules are independent."""
        # Given
        meta_dir = tmp_path / "meta"
        rules = ["standardize", "qc_busco", "orthofinder"]
        wildcards = {"species": "mmur"}

        # Mark only first two as complete
        mark_rule_complete(rules[0], wildcards, meta_dir)
        mark_rule_complete(rules[1], wildcards, meta_dir)

        # Then
        assert is_rule_complete(rules[0], wildcards, meta_dir)
        assert is_rule_complete(rules[1], wildcards, meta_dir)
        assert not is_rule_complete(rules[2], wildcards, meta_dir)

    def test_checkpoint_survives_process_restart(self, tmp_path: Path) -> None:
        """[P1] Test that checkpoints persist across process simulated restarts."""
        # Given
        meta_dir = tmp_path / "meta"
        rule = "test_rule"
        wildcards = {"key": "value"}

        # Simulate first process
        mark_rule_complete(rule, wildcards, meta_dir)
        checkpoint_path = get_checkpoint_path(rule, wildcards, meta_dir)

        # Then - verify file exists on disk
        assert checkpoint_path.exists()

        # Simulate process restart by re-checking
        assert is_rule_complete(rule, wildcards, meta_dir)


# =============================================================================
# Test Atomic Write for Checkpoint Safety
# =============================================================================

class TestAtomicWriteCheckpointSafety:
    """Integration tests for atomic write ensuring checkpoint safety."""

    def test_atomic_write_no_partial_file_on_success(self, tmp_path: Path) -> None:
        """[P1] Test that successful atomic write leaves no temp files."""
        # Given
        output_path = tmp_path / "output.txt"
        content = "Complete content"

        # When
        atomic_write(output_path, content)

        # Then
        assert output_path.exists()
        assert not get_temp_path(output_path).exists()
        assert output_path.read_text() == content

    def test_cleanup_removes_orphaned_temp_files(self, tmp_path: Path) -> None:
        """[P1] Test that cleanup removes orphaned .tmp files."""
        # Given - simulate interrupted writes
        temp1 = tmp_path / "file1.txt.tmp"
        temp2 = tmp_path / "subdir" / "file2.json.tmp"
        temp2.parent.mkdir(parents=True)
        temp1.write_text("interrupted")
        temp2.write_text("also interrupted")

        # When
        removed = cleanup_temp_files(tmp_path)

        # Then
        assert len(removed) == 2
        assert not temp1.exists()
        assert not temp2.exists()

    def test_completed_files_not_affected_by_cleanup(self, tmp_path: Path) -> None:
        """[P1] Test that cleanup doesn't affect completed files."""
        # Given - mix of completed and temp files
        completed = tmp_path / "completed.txt"
        orphaned = tmp_path / "orphaned.txt.tmp"
        completed.write_text("done")
        orphaned.write_text("interrupted")

        # When
        cleanup_temp_files(tmp_path)

        # Then
        assert completed.exists()
        assert completed.read_text() == "done"
        assert not orphaned.exists()


# =============================================================================
# Test Cache Performance for Resume
# =============================================================================

class TestCachePerformanceResume:
    """Integration tests for cache hit performance (NFR4: < 1 second)."""

    def test_checksum_cache_hit_under_1_second(self, tmp_path: Path) -> None:
        """[P1] Test that checksum cache hits are fast enough for NFR4."""
        # Given - create multiple files and cache their checksums
        files = []
        for i in range(100):
            f = tmp_path / f"file_{i:03d}.txt"
            f.write_text(f"Content for file {i}")
            files.append(f)

        cache = ChecksumCache(tmp_path / "cache")

        # First pass - compute checksums
        for f in files:
            cache.get_or_compute(f)
        cache.save()

        # When - measure cache hit time for all files
        start = time.perf_counter()
        for f in files:
            cache.get_or_compute(f)
        elapsed = time.perf_counter() - start

        # Then - should be well under 1 second
        assert elapsed < 1.0, f"100 cache hits took {elapsed:.3f}s (NFR4: < 1s)"

    def test_checkpoint_check_under_1_second(self, tmp_path: Path) -> None:
        """[P1] Test that checkpoint checks are fast enough for NFR4."""
        # Given - create many checkpoints
        meta_dir = tmp_path / "meta"
        checkpoints = []
        for i in range(100):
            rule = f"rule_{i:03d}"
            wildcards = {"species": f"sp_{i:03d}"}
            mark_rule_complete(rule, wildcards, meta_dir)
            checkpoints.append((rule, wildcards))

        # When - check all checkpoints
        start = time.perf_counter()
        for rule, wildcards in checkpoints:
            is_rule_complete(rule, wildcards, meta_dir)
        elapsed = time.perf_counter() - start

        # Then - should be well under 1 second
        assert elapsed < 1.0, f"100 checkpoint checks took {elapsed:.3f}s (NFR4: < 1s)"


# =============================================================================
# Test Simulated Interrupt and Resume
# =============================================================================

class TestSimulatedInterruptResume:
    """Integration tests simulating process interruption and resume."""

    def test_resume_after_simulated_interrupt(self, tmp_path: Path) -> None:
        """[P1] Test resume scenario after simulated interruption."""
        # Given - simulate partial completion
        meta_dir = tmp_path / "meta"
        output_dir = tmp_path / "results"
        output_dir.mkdir()

        # Complete first two steps
        step1_output = output_dir / "step1.txt"
        step2_output = output_dir / "step2.txt"
        step3_output = output_dir / "step3.txt"

        atomic_write(step1_output, "Step 1 complete")
        mark_rule_complete("step1", {}, meta_dir)

        atomic_write(step2_output, "Step 2 complete")
        mark_rule_complete("step2", {}, meta_dir)

        # Step 3 interrupted - leave orphaned temp file
        step3_temp = get_temp_path(step3_output)
        step3_temp.write_text("Step 3 interrupted")

        # When - simulate resume: cleanup temp and check status
        cleanup_temp_files(output_dir)

        # Then - verify state
        assert step1_output.exists()  # Completed
        assert step2_output.exists()  # Completed
        assert not step3_output.exists()  # Needs re-run
        assert not step3_temp.exists()  # Cleaned up

        assert is_rule_complete("step1", {}, meta_dir)
        assert is_rule_complete("step2", {}, meta_dir)
        assert not is_rule_complete("step3", {}, meta_dir)

    def test_file_change_invalidates_downstream(self, tmp_path: Path) -> None:
        """[P2] Test that file changes are detected for cache invalidation."""
        # Given
        cache_dir = tmp_path / "cache"
        input_file = tmp_path / "input.txt"
        input_file.write_text("Original content")

        cache = ChecksumCache(cache_dir)
        original_checksum = cache.get_or_compute(input_file)
        cache.save()

        # When - modify file
        time.sleep(0.01)  # Ensure mtime changes
        input_file.write_text("Modified content")

        # Reload cache and compute
        cache2 = ChecksumCache(cache_dir)
        new_checksum = cache2.get_or_compute(input_file)

        # Then - checksums should differ (change detected)
        assert original_checksum != new_checksum
