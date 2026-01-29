"""
Tests for the validate_config.py CLI script.

These tests validate the command-line interface for configuration validation.
"""

import pytest
import subprocess
import sys
from pathlib import Path


# Paths
PROJECT_ROOT = Path(__file__).parent.parent
SCRIPT_PATH = PROJECT_ROOT / "workflow" / "scripts" / "validate_config.py"
FIXTURES_DIR = Path(__file__).parent / "fixtures"
SCHEMA_PATH = PROJECT_ROOT / "schemas" / "config.schema.yaml"


class TestValidateConfigCLI:
    """Tests for the validate_config.py command-line interface."""

    def test_valid_config_exits_zero(self):
        """[P1] Valid config should exit with code 0."""
        # GIVEN: Valid config file
        config_path = FIXTURES_DIR / "valid_config.yaml"

        # WHEN: Running validation script
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit with code 0
        assert result.returncode == 0
        assert "valid" in result.stdout.lower()

    def test_invalid_config_exits_nonzero(self):
        """[P1] Invalid config should exit with non-zero code."""
        # GIVEN: Invalid config file
        config_path = FIXTURES_DIR / "invalid_config.yaml"

        # WHEN: Running validation script
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit with code 1
        assert result.returncode == 1
        assert "FAILED" in result.stdout or "error" in result.stdout.lower()

    def test_missing_config_file_exits_nonzero(self, tmp_path):
        """[P1] Missing config file should exit with non-zero code."""
        # GIVEN: Non-existent config file
        config_path = tmp_path / "nonexistent.yaml"

        # WHEN: Running validation script
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit with code 1
        assert result.returncode == 1
        assert "not found" in result.stdout.lower() or "error" in result.stdout.lower()

    def test_quiet_mode_suppresses_success_message(self):
        """[P2] Quiet mode should not print success message."""
        # GIVEN: Valid config file
        config_path = FIXTURES_DIR / "valid_config.yaml"

        # WHEN: Running with --quiet flag
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), "-q", str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit 0 with no output
        assert result.returncode == 0
        assert result.stdout.strip() == ""

    def test_custom_schema_path(self):
        """[P2] Should accept custom schema path."""
        # GIVEN: Valid config and explicit schema path
        config_path = FIXTURES_DIR / "valid_config.yaml"

        # WHEN: Running with --schema flag
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH),
             "--schema", str(SCHEMA_PATH),
             str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit 0
        assert result.returncode == 0

    def test_invalid_yaml_syntax_exits_nonzero(self, tmp_path):
        """[P1] Invalid YAML syntax should exit with non-zero code."""
        # GIVEN: Config file with invalid YAML
        config_path = tmp_path / "invalid_yaml.yaml"
        config_path.write_text("""
species:
  - id: test
    name: "unclosed quote
        """)

        # WHEN: Running validation
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), str(config_path)],
            capture_output=True,
            text=True
        )

        # THEN: Should exit with code 1
        assert result.returncode == 1
        assert "yaml" in result.stdout.lower() or "parsing" in result.stdout.lower()


class TestValidateConfigIntegration:
    """Integration tests for validate_config.py functions."""

    def test_validate_config_function_with_valid_config(self):
        """[P1] validate_config function should return empty list for valid config."""
        # GIVEN: Import the validation function
        sys.path.insert(0, str(PROJECT_ROOT / "workflow" / "scripts"))
        from validate_config import validate_config

        # WHEN: Validating valid config
        errors = validate_config(
            FIXTURES_DIR / "valid_config.yaml",
            SCHEMA_PATH
        )

        # THEN: Should return empty list
        assert errors == []

    def test_validate_config_function_with_invalid_config(self):
        """[P1] validate_config function should return errors for invalid config."""
        # GIVEN: Import the validation function
        sys.path.insert(0, str(PROJECT_ROOT / "workflow" / "scripts"))
        from validate_config import validate_config

        # WHEN: Validating invalid config
        errors = validate_config(
            FIXTURES_DIR / "invalid_config.yaml",
            SCHEMA_PATH
        )

        # THEN: Should return non-empty list
        assert len(errors) > 0

    def test_validate_business_rules_detects_duplicates(self):
        """[P1] Business rules should detect duplicate species IDs."""
        # GIVEN: Import the business rule function
        sys.path.insert(0, str(PROJECT_ROOT / "workflow" / "scripts"))
        from validate_config import validate_business_rules

        # AND: Config with duplicate species IDs
        config = {
            "species": [
                {"id": "mmur", "name": "Species 1", "annotation": "a.gff3", "genome": "g.fa"},
                {"id": "mmur", "name": "Species 2", "annotation": "b.gff3", "genome": "h.fa"}
            ],
            "output_dir": "results"
        }

        # WHEN: Validating business rules
        errors = validate_business_rules(config)

        # THEN: Should detect duplicate
        assert len(errors) > 0
        assert any("duplicate" in e.lower() or "mmur" in e.lower() for e in errors)

    def test_validate_business_rules_passes_unique_ids(self):
        """[P1] Business rules should pass for unique species IDs."""
        # GIVEN: Import the business rule function
        sys.path.insert(0, str(PROJECT_ROOT / "workflow" / "scripts"))
        from validate_config import validate_business_rules

        # AND: Config with unique species IDs
        config = {
            "species": [
                {"id": "mmur", "name": "Species 1", "annotation": "a.gff3", "genome": "g.fa"},
                {"id": "lcat", "name": "Species 2", "annotation": "b.gff3", "genome": "h.fa"}
            ],
            "output_dir": "results"
        }

        # WHEN: Validating business rules
        errors = validate_business_rules(config)

        # THEN: Should pass
        assert errors == []
