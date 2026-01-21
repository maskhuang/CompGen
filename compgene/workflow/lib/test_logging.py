"""
Unit tests for workflow/lib/logging.py

Tests dual-format logging functionality.
"""

import json
import logging
from pathlib import Path

import pytest

from workflow.lib.logging import (
    DualLogger,
    get_log_paths,
    format_wildcards_for_path,
    create_logger,
    ANSIColors,
)


# =============================================================================
# Test get_log_paths
# =============================================================================

class TestGetLogPaths:
    """Tests for get_log_paths function."""

    def test_generates_both_paths(self, tmp_path: Path) -> None:
        """[P1] Given rule and wildcards, then both .log and .jsonl paths returned."""
        # Given
        rule = "orthofinder"
        wildcards = {"species_set": "lemur"}

        # When
        log_path, jsonl_path = get_log_paths(rule, wildcards, tmp_path)

        # Then
        assert log_path.suffix == ".log"
        assert jsonl_path.suffix == ".jsonl"

    def test_includes_rule_in_path(self, tmp_path: Path) -> None:
        """[P1] Given rule name, then rule included in path."""
        # Given
        rule = "qc_busco"
        wildcards = {}

        # When
        log_path, _ = get_log_paths(rule, wildcards, tmp_path)

        # Then
        assert rule in str(log_path)

    def test_includes_wildcards_in_filename(self, tmp_path: Path) -> None:
        """[P1] Given wildcards, then wildcards included in filename."""
        # Given
        rule = "test_rule"
        wildcards = {"species": "mmur", "sample": "s1"}

        # When
        log_path, _ = get_log_paths(rule, wildcards, tmp_path)

        # Then
        filename = log_path.stem
        assert "species=mmur" in filename
        assert "sample=s1" in filename

    def test_uses_default_for_empty_wildcards(self, tmp_path: Path) -> None:
        """[P2] Given empty wildcards, then 'default' used in filename."""
        # Given
        rule = "test_rule"
        wildcards: dict[str, str] = {}

        # When
        log_path, _ = get_log_paths(rule, wildcards, tmp_path)

        # Then
        assert "default" in log_path.stem


# =============================================================================
# Test format_wildcards_for_path
# =============================================================================

class TestFormatWildcardsForPath:
    """Tests for format_wildcards_for_path function."""

    def test_formats_single_wildcard(self) -> None:
        """[P2] Given single wildcard, then formatted correctly."""
        # Given
        wildcards = {"species": "mmur"}

        # When
        result = format_wildcards_for_path(wildcards)

        # Then
        assert result == "species=mmur"

    def test_formats_multiple_wildcards_sorted(self) -> None:
        """[P2] Given multiple wildcards, then sorted and joined."""
        # Given
        wildcards = {"z_last": "value1", "a_first": "value2"}

        # When
        result = format_wildcards_for_path(wildcards)

        # Then
        assert result == "a_first=value2_z_last=value1"

    def test_returns_default_for_empty(self) -> None:
        """[P2] Given empty wildcards, then 'default' returned."""
        # Given
        wildcards: dict[str, str] = {}

        # When
        result = format_wildcards_for_path(wildcards)

        # Then
        assert result == "default"


# =============================================================================
# Test DualLogger
# =============================================================================

class TestDualLogger:
    """Tests for DualLogger class."""

    def test_creates_log_files(self, tmp_path: Path) -> None:
        """[P1] Given logger initialized, when logging, then both files created."""
        # Given
        logger = DualLogger("test_rule", {"species": "mmur"}, tmp_path)

        # When
        logger.info("Test message")

        # Then
        assert logger.log_path.exists()
        assert logger.jsonl_path.exists()

    def test_writes_to_text_log(self, tmp_path: Path) -> None:
        """[P1] Given logger, when info called, then message in text log."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, use_color=False)

        # When
        logger.info("Test message")

        # Then
        content = logger.log_path.read_text()
        assert "Test message" in content
        assert "[INFO]" in content
        assert "test_rule" in content

    def test_writes_to_jsonl_log(self, tmp_path: Path) -> None:
        """[P1] Given logger, when info called, then valid JSON in jsonl log."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path)

        # When
        logger.info("Test message")

        # Then
        content = logger.jsonl_path.read_text().strip()
        record = json.loads(content)
        assert record["msg"] == "Test message"
        assert record["level"] == "INFO"
        assert record["rule"] == "test_rule"

    def test_includes_extra_fields_in_jsonl(self, tmp_path: Path) -> None:
        """[P1] Given extra kwargs, when logging, then included in JSON."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path)

        # When
        logger.info("Starting", species_count=3, memory_mb=4096)

        # Then
        content = logger.jsonl_path.read_text().strip()
        record = json.loads(content)
        assert record["species_count"] == 3
        assert record["memory_mb"] == 4096

    def test_respects_log_level(self, tmp_path: Path) -> None:
        """[P1] Given WARNING level, when debug called, then not logged."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, level=logging.WARNING)

        # When
        logger.debug("Debug message")
        logger.warning("Warning message")

        # Then
        content = logger.log_path.read_text()
        assert "Debug message" not in content
        assert "Warning message" in content

    def test_all_log_levels_work(self, tmp_path: Path) -> None:
        """[P2] Given logger, when all levels called, then all logged appropriately."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, level=logging.DEBUG, use_color=False)

        # When
        logger.debug("Debug")
        logger.info("Info")
        logger.warning("Warning")
        logger.error("Error")

        # Then
        content = logger.log_path.read_text()
        assert "[DEBUG]" in content
        assert "[INFO]" in content
        assert "[WARNING]" in content
        assert "[ERROR]" in content

    def test_includes_wildcards_in_jsonl(self, tmp_path: Path) -> None:
        """[P2] Given wildcards, when logging, then included in JSON."""
        # Given
        wildcards = {"species": "mmur", "sample": "s1"}
        logger = DualLogger("test_rule", wildcards, tmp_path)

        # When
        logger.info("Test")

        # Then
        content = logger.jsonl_path.read_text().strip()
        record = json.loads(content)
        assert record["wildcards"] == wildcards

    def test_timestamp_format_in_text_log(self, tmp_path: Path) -> None:
        """[P2] Given logger, when logging, then timestamp in correct format."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, use_color=False)

        # When
        logger.info("Test")

        # Then
        content = logger.log_path.read_text()
        # Should match YYYY-MM-DD HH:MM:SS format
        import re
        assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", content)

    def test_timestamp_is_iso8601_in_jsonl(self, tmp_path: Path) -> None:
        """[P2] Given logger, when logging, then ISO8601 timestamp in JSON."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path)

        # When
        logger.info("Test")

        # Then
        content = logger.jsonl_path.read_text().strip()
        record = json.loads(content)
        # ISO8601 format check
        assert "T" in record["ts"]
        assert ":" in record["ts"]

    def test_set_level_changes_threshold(self, tmp_path: Path) -> None:
        """[P2] Given logger, when set_level called, then threshold changes."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, level=logging.WARNING)

        # When
        logger.debug("Before")
        logger.set_level(logging.DEBUG)
        logger.debug("After")

        # Then
        content = logger.log_path.read_text()
        assert "Before" not in content
        assert "After" in content

    def test_get_level_from_string(self, tmp_path: Path) -> None:
        """[P2] Given level string, when get_level_from_string, then correct constant."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path)

        # When/Then
        assert logger.get_level_from_string("DEBUG") == logging.DEBUG
        assert logger.get_level_from_string("INFO") == logging.INFO
        assert logger.get_level_from_string("WARNING") == logging.WARNING
        assert logger.get_level_from_string("ERROR") == logging.ERROR
        assert logger.get_level_from_string("invalid") == logging.INFO  # Default

    def test_color_output_disabled(self, tmp_path: Path) -> None:
        """[P2] Given use_color=False, then no ANSI codes in output."""
        # Given
        logger = DualLogger("test_rule", {}, tmp_path, use_color=False)

        # When
        logger.info("Test")

        # Then
        content = logger.log_path.read_text()
        assert "\033[" not in content  # No ANSI escape codes


# =============================================================================
# Test create_logger
# =============================================================================

class TestCreateLogger:
    """Tests for create_logger convenience function."""

    def test_creates_dual_logger(self, tmp_path: Path) -> None:
        """[P1] Given parameters, when create_logger called, then DualLogger returned."""
        # Given/When
        logger = create_logger("test_rule", {"species": "mmur"}, tmp_path)

        # Then
        assert isinstance(logger, DualLogger)
        assert logger.rule == "test_rule"

    def test_accepts_level_as_string(self, tmp_path: Path) -> None:
        """[P2] Given level as string, when create_logger called, then level set."""
        # Given/When
        logger = create_logger("test_rule", level="DEBUG", log_dir=tmp_path)

        # Then
        assert logger.level == logging.DEBUG

    def test_uses_defaults(self, tmp_path: Path) -> None:
        """[P2] Given minimal params, when create_logger called, then defaults used."""
        # Given/When
        logger = create_logger("test_rule", log_dir=tmp_path)

        # Then
        assert logger.wildcards == {}
        assert logger.use_color is True
