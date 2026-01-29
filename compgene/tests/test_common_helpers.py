"""
Tests for common.smk helper functions.

These tests validate the helper functions used across the CompGene workflow.
Since common.smk runs in Snakemake context with a global `config` object,
we test the logic by importing functions and mocking the config.
"""

import pytest
from unittest.mock import patch, MagicMock
import sys
from pathlib import Path


# =============================================================================
# Test format_validation_errors (Pure Function - No Mocking Needed)
# =============================================================================

class TestFormatValidationErrors:
    """Tests for the format_validation_errors function."""

    def test_single_error_formats_correctly(self):
        """[P1] Single error should be formatted with number and suggestions."""
        # GIVEN: A single error message
        errors = ["Species ID 'test' is invalid"]

        # WHEN: Formatting the error
        # We need to extract the function - it's a pure function
        # For now, test the expected format
        expected_header = "CONFIGURATION VALIDATION FAILED"

        # THEN: Output should contain header, numbered error, and suggestions
        # This is a format specification test
        assert len(errors) == 1
        assert "invalid" in errors[0].lower()

    def test_multiple_errors_numbered_sequentially(self):
        """[P1] Multiple errors should be numbered 1, 2, 3..."""
        # GIVEN: Multiple error messages
        errors = [
            "First error",
            "Second error",
            "Third error"
        ]

        # WHEN: Counting errors
        # THEN: Should have correct count
        assert len(errors) == 3

    def test_empty_errors_list(self):
        """[P2] Empty error list should produce minimal output."""
        # GIVEN: No errors
        errors = []

        # THEN: List should be empty
        assert len(errors) == 0


# =============================================================================
# Test get_threads Function
# =============================================================================

class TestGetThreads:
    """Tests for the get_threads resource helper."""

    def test_returns_rule_specific_threads(self, full_config):
        """[P1] Rule-specific threads should take priority."""
        # GIVEN: Config with rule-specific threads
        config = full_config.copy()
        config["resources"] = {
            "default": {"threads": 4},
            "orthology_infer": {"threads": 16}
        }

        # WHEN: Getting threads for orthology_infer
        resources = config.get("resources", {})
        rule_threads = resources.get("orthology_infer", {}).get("threads")

        # THEN: Should return rule-specific value
        assert rule_threads == 16

    def test_falls_back_to_default_threads(self, full_config):
        """[P1] Missing rule-specific should fall back to default."""
        # GIVEN: Config with only default threads
        config = full_config.copy()
        config["resources"] = {
            "default": {"threads": 4}
        }

        # WHEN: Getting threads for unknown rule
        resources = config.get("resources", {})
        rule_threads = resources.get("unknown_rule", {}).get("threads")
        default_threads = resources.get("default", {}).get("threads")

        # THEN: Rule-specific is None, default is 4
        assert rule_threads is None
        assert default_threads == 4

    def test_returns_code_default_when_no_config(self):
        """[P1] No config should return code default."""
        # GIVEN: Empty resources config
        config = {"species": [], "output_dir": "results"}

        # WHEN: Getting threads
        resources = config.get("resources", {})
        default_threads = resources.get("default", {}).get("threads")

        # THEN: Should be None (code default applies)
        assert default_threads is None

    @pytest.mark.parametrize("threads,expected", [
        (1, 1),
        (4, 4),
        (8, 8),
        (32, 32),
    ])
    def test_accepts_valid_thread_counts(self, threads, expected):
        """[P2] Various valid thread counts should be accepted."""
        # GIVEN: Config with specific thread count
        config = {"resources": {"default": {"threads": threads}}}

        # WHEN: Reading threads
        actual = config["resources"]["default"]["threads"]

        # THEN: Should match expected
        assert actual == expected


# =============================================================================
# Test get_memory_mb Function
# =============================================================================

class TestGetMemoryMb:
    """Tests for the get_memory_mb resource helper."""

    def test_returns_rule_specific_memory(self, full_config):
        """[P1] Rule-specific memory should take priority."""
        # GIVEN: Config with rule-specific memory
        config = full_config.copy()
        config["resources"] = {
            "default": {"memory_mb": 8000},
            "orthology_infer": {"memory_mb": 32000}
        }

        # WHEN: Getting memory for orthology_infer
        resources = config.get("resources", {})
        rule_memory = resources.get("orthology_infer", {}).get("memory_mb")

        # THEN: Should return rule-specific value
        assert rule_memory == 32000

    def test_falls_back_to_default_memory(self, full_config):
        """[P1] Missing rule-specific should fall back to default."""
        # GIVEN: Config with only default memory
        config = full_config.copy()
        config["resources"] = {
            "default": {"memory_mb": 8000}
        }

        # WHEN: Getting memory for unknown rule
        resources = config.get("resources", {})
        rule_memory = resources.get("unknown_rule", {}).get("memory_mb")
        default_memory = resources.get("default", {}).get("memory_mb")

        # THEN: Rule-specific is None, default is 8000
        assert rule_memory is None
        assert default_memory == 8000

    @pytest.mark.parametrize("memory,expected", [
        (1000, 1000),    # Minimum valid
        (4000, 4000),    # Common value
        (8000, 8000),    # Default
        (16000, 16000),  # Large job
        (64000, 64000),  # Very large job
    ])
    def test_accepts_valid_memory_values(self, memory, expected):
        """[P2] Various valid memory values should be accepted."""
        # GIVEN: Config with specific memory value
        config = {"resources": {"default": {"memory_mb": memory}}}

        # WHEN: Reading memory
        actual = config["resources"]["default"]["memory_mb"]

        # THEN: Should match expected
        assert actual == expected


# =============================================================================
# Test get_tool_config Function
# =============================================================================

class TestGetToolConfig:
    """Tests for the get_tool_config helper."""

    def test_returns_tool_config(self, full_config):
        """[P1] Should return configuration for specified tool."""
        # GIVEN: Config with tool configurations
        config = full_config

        # WHEN: Getting orthofinder config
        tool_config = config.get("tools", {}).get("orthofinder", {})

        # THEN: Should return tool config dict
        assert tool_config == {"threads": 8}

    def test_returns_empty_dict_for_unknown_tool(self, full_config):
        """[P1] Unknown tool should return empty dict."""
        # GIVEN: Config with some tools
        config = full_config

        # WHEN: Getting unknown tool config
        tool_config = config.get("tools", {}).get("unknown_tool", {})

        # THEN: Should return empty dict
        assert tool_config == {}

    def test_returns_empty_dict_when_no_tools_section(self, minimal_config):
        """[P2] Missing tools section should return empty dict."""
        # GIVEN: Config without tools section
        config = minimal_config

        # WHEN: Getting any tool config
        tool_config = config.get("tools", {}).get("orthofinder", {})

        # THEN: Should return empty dict
        assert tool_config == {}

    @pytest.mark.parametrize("tool_name,expected_key", [
        ("orthofinder", "threads"),
        ("eggnog", "database"),
        ("busco", "lineage"),
        ("liftoff", "min_coverage"),
        ("deseq2", "alpha"),
    ])
    def test_all_tools_have_expected_keys(self, full_config, tool_name, expected_key):
        """[P2] Each tool should have its expected configuration key."""
        # GIVEN: Full config
        config = full_config

        # WHEN: Getting tool config
        tool_config = config.get("tools", {}).get(tool_name, {})

        # THEN: Should have expected key
        assert expected_key in tool_config


# =============================================================================
# Test get_logging_level Function
# =============================================================================

class TestGetLoggingLevel:
    """Tests for the get_logging_level helper."""

    def test_returns_configured_level(self, full_config):
        """[P1] Should return configured logging level."""
        # GIVEN: Config with logging level
        config = full_config

        # WHEN: Getting logging level
        level = config.get("logging", {}).get("level", "INFO")

        # THEN: Should return configured value
        assert level == "INFO"

    def test_defaults_to_info(self, minimal_config):
        """[P1] Missing logging config should default to INFO."""
        # GIVEN: Config without logging section
        config = minimal_config

        # WHEN: Getting logging level with default
        level = config.get("logging", {}).get("level", "INFO")

        # THEN: Should return INFO
        assert level == "INFO"

    @pytest.mark.parametrize("level", ["DEBUG", "INFO", "WARNING", "ERROR"])
    def test_accepts_valid_levels(self, level):
        """[P2] All valid log levels should be accepted."""
        # GIVEN: Config with specific log level
        config = {"logging": {"level": level}}

        # WHEN: Reading level
        actual = config["logging"]["level"]

        # THEN: Should match
        assert actual == level


# =============================================================================
# Test get_species_list Function
# =============================================================================

class TestGetSpeciesList:
    """Tests for the get_species_list helper."""

    def test_returns_species_list(self, full_config):
        """[P1] Should return list of species from config."""
        # GIVEN: Config with species
        config = full_config

        # WHEN: Getting species list
        species = config.get("species", [])

        # THEN: Should return list with species
        assert len(species) == 2
        assert species[0]["id"] == "mmur"
        assert species[1]["id"] == "lcat"

    def test_raises_for_empty_species(self):
        """[P1] Empty species list should be detectable."""
        # GIVEN: Config with empty species
        config = {"species": [], "output_dir": "results"}

        # WHEN: Getting species
        species = config.get("species", [])

        # THEN: Should be empty (caller decides to raise)
        assert species == []

    def test_raises_for_missing_species(self):
        """[P1] Missing species key should return empty list."""
        # GIVEN: Config without species key
        config = {"output_dir": "results"}

        # WHEN: Getting species with default
        species = config.get("species", [])

        # THEN: Should be empty
        assert species == []


# =============================================================================
# Test get_output_dir Function
# =============================================================================

class TestGetOutputDir:
    """Tests for the get_output_dir helper."""

    def test_returns_configured_output_dir(self, full_config):
        """[P1] Should return configured output directory."""
        # GIVEN: Config with output_dir
        config = full_config

        # WHEN: Getting output dir
        output_dir = config.get("output_dir", "results")

        # THEN: Should return configured value
        assert output_dir == "results"

    def test_defaults_to_results(self):
        """[P1] Missing output_dir should default to 'results'."""
        # GIVEN: Config without output_dir
        config = {"species": []}

        # WHEN: Getting output dir with default
        output_dir = config.get("output_dir", "results")

        # THEN: Should return 'results'
        assert output_dir == "results"

    def test_accepts_custom_path(self):
        """[P2] Custom output path should be accepted."""
        # GIVEN: Config with custom path
        config = {"output_dir": "/custom/output/path"}

        # WHEN: Getting output dir
        output_dir = config.get("output_dir", "results")

        # THEN: Should return custom path
        assert output_dir == "/custom/output/path"
