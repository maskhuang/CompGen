"""
Tests for configuration validation.

Tests the JSON Schema validation and business rule validation
for CompGene configuration files.
"""

import pytest
import yaml
from pathlib import Path

# Try to import jsonschema, skip tests if not available
pytest.importorskip("jsonschema")
from jsonschema import Draft7Validator, ValidationError


# Paths
FIXTURES_DIR = Path(__file__).parent / "fixtures"
SCHEMA_PATH = Path(__file__).parent.parent / "schemas" / "config.schema.yaml"


def load_yaml(path: Path) -> dict:
    """Load a YAML file."""
    with open(path) as f:
        return yaml.safe_load(f)


@pytest.fixture
def schema():
    """Load the configuration schema."""
    return load_yaml(SCHEMA_PATH)


@pytest.fixture
def valid_config():
    """Load the valid test configuration."""
    return load_yaml(FIXTURES_DIR / "valid_config.yaml")


@pytest.fixture
def invalid_config():
    """Load the invalid test configuration."""
    return load_yaml(FIXTURES_DIR / "invalid_config.yaml")


class TestSchemaValidation:
    """Tests for JSON Schema validation."""

    def test_valid_config_passes(self, schema, valid_config):
        """Test that a valid configuration passes validation."""
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(valid_config))
        assert len(errors) == 0, f"Valid config should pass: {errors}"

    def test_invalid_config_fails(self, schema, invalid_config):
        """Test that an invalid configuration fails validation."""
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(invalid_config))
        assert len(errors) > 0, "Invalid config should have errors"

    def test_missing_species_fails(self, schema):
        """Test that missing species field fails validation."""
        config = {
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        # Should fail because 'species' is required
        assert any("species" in str(e.message) for e in errors)

    def test_missing_output_dir_fails(self, schema):
        """Test that missing output_dir field fails validation."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ]
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        # Should fail because 'output_dir' is required
        assert any("output_dir" in str(e.message) for e in errors)

    def test_empty_species_list_fails(self, schema):
        """Test that empty species list fails validation."""
        config = {
            "species": [],
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        # Should fail because minItems: 1
        assert len(errors) > 0

    def test_invalid_species_id_pattern_fails(self, schema):
        """Test that invalid species ID pattern fails validation."""
        config = {
            "species": [
                {
                    "id": "Invalid_ID",  # Should be lowercase
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        # Should fail because of pattern mismatch
        assert any("pattern" in str(e.message).lower() or "Invalid_ID" in str(e.message) for e in errors)

    def test_invalid_logging_level_fails(self, schema):
        """Test that invalid logging level fails validation."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results",
            "logging": {
                "level": "TRACE"  # Invalid - not in enum
            }
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        assert any("TRACE" in str(e.message) for e in errors)

    def test_invalid_thread_count_fails(self, schema):
        """Test that thread count < 1 fails validation."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results",
            "resources": {
                "default": {
                    "threads": 0  # Invalid - must be >= 1
                }
            }
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        assert len(errors) > 0

    def test_invalid_memory_mb_fails(self, schema):
        """Test that memory_mb < 1000 fails validation."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results",
            "resources": {
                "default": {
                    "memory_mb": 500  # Invalid - must be >= 1000
                }
            }
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        assert len(errors) > 0

    def test_deseq2_alpha_range(self, schema):
        """Test that deseq2 alpha must be between 0 and 1."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results",
            "tools": {
                "deseq2": {
                    "alpha": 1.5  # Invalid - must be <= 1
                }
            }
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        assert len(errors) > 0


class TestMinimalValidConfig:
    """Tests with minimal valid configurations."""

    def test_minimal_config_passes(self, schema):
        """Test that minimal required config passes."""
        config = {
            "species": [
                {
                    "id": "test",
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))
        assert len(errors) == 0

    def test_multiple_species_passes(self, schema):
        """Test that multiple species configuration passes."""
        config = {
            "species": [
                {
                    "id": "species1",
                    "name": "Species One",
                    "annotation": "sp1.gff3",
                    "genome": "sp1.fa"
                },
                {
                    "id": "species2",
                    "name": "Species Two",
                    "annotation": "sp2.gff3",
                    "genome": "sp2.fa"
                },
                {
                    "id": "species3",
                    "name": "Species Three",
                    "annotation": "sp3.gff3",
                    "genome": "sp3.fa"
                }
            ],
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))
        assert len(errors) == 0


class TestSpeciesIdPattern:
    """Tests for species ID pattern validation."""

    @pytest.mark.parametrize("species_id,should_pass", [
        ("mmur", True),           # Valid: lowercase
        ("species1", True),       # Valid: lowercase with number
        ("sp_test", True),        # Valid: lowercase with underscore
        ("a", True),              # Valid: single letter
        ("test123", True),        # Valid: letters then numbers
        ("Mmur", False),          # Invalid: starts with uppercase
        ("MMUR", False),          # Invalid: all uppercase
        ("1species", False),      # Invalid: starts with number
        ("sp-test", False),       # Invalid: contains hyphen
        ("sp.test", False),       # Invalid: contains dot
        ("sp test", False),       # Invalid: contains space
    ])
    def test_species_id_pattern(self, schema, species_id, should_pass):
        """Test various species ID patterns."""
        config = {
            "species": [
                {
                    "id": species_id,
                    "name": "Test Species",
                    "annotation": "test.gff3",
                    "genome": "test.fa"
                }
            ],
            "output_dir": "results"
        }
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(config))

        if should_pass:
            # Filter for pattern errors only
            pattern_errors = [e for e in errors if "pattern" in str(e.schema_path)]
            assert len(pattern_errors) == 0, f"'{species_id}' should be valid"
        else:
            assert len(errors) > 0, f"'{species_id}' should be invalid"


class TestBusinessRuleValidation:
    """Tests for business rule validation beyond JSON Schema."""

    def test_duplicate_species_ids_detected(self):
        """Test that duplicate species IDs are detected."""
        # Import the validation function
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
        from validate_config import validate_business_rules

        config = {
            "species": [
                {
                    "id": "mmur",
                    "name": "Species One",
                    "annotation": "sp1.gff3",
                    "genome": "sp1.fa"
                },
                {
                    "id": "mmur",  # Duplicate!
                    "name": "Species Two",
                    "annotation": "sp2.gff3",
                    "genome": "sp2.fa"
                }
            ],
            "output_dir": "results"
        }
        errors = validate_business_rules(config)
        assert len(errors) > 0, "Duplicate species IDs should be detected"
        assert any("duplicate" in e.lower() or "mmur" in e.lower() for e in errors)

    def test_unique_species_ids_pass(self):
        """Test that unique species IDs pass validation."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
        from validate_config import validate_business_rules

        config = {
            "species": [
                {
                    "id": "mmur",
                    "name": "Species One",
                    "annotation": "sp1.gff3",
                    "genome": "sp1.fa"
                },
                {
                    "id": "lcat",
                    "name": "Species Two",
                    "annotation": "sp2.gff3",
                    "genome": "sp2.fa"
                }
            ],
            "output_dir": "results"
        }
        errors = validate_business_rules(config)
        assert len(errors) == 0, "Unique species IDs should pass"
