"""Tests for Story 1.2: Configuration schema validation."""

import sys
from pathlib import Path

import pytest
import yaml

# Add workflow/lib to path for imports
workflow_lib = Path(__file__).resolve().parent.parent / "workflow" / "lib"
sys.path.insert(0, str(workflow_lib))

from config import (
    ConfigValidationError,
    get_config_value,
    set_config_value,
    merge_cli_config,
    format_validation_error,
    validate_config_value,
    load_config,
    _deep_copy_dict,
)


@pytest.fixture
def project_root():
    """Return the project root directory."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def valid_config():
    """Return a valid configuration dictionary."""
    return {
        "project": {
            "name": "test_analysis",
            "output_dir": "./results",
        },
        "species": [
            {
                "name": "Microcebus_murinus",
                "assembly": "/path/to/assembly.fa",
                "annotation": "/path/to/annotation.gff",
            }
        ],
        "analysis": {
            "orthology": {"method": "orthofinder"},
            "expression": {"enabled": False},
        },
        "resources": {
            "threads": 8,
            "memory": "16G",
        },
    }


@pytest.fixture
def schema_path(project_root):
    """Return path to the config schema."""
    return project_root / "schemas" / "config.schema.yaml"


class TestGetConfigValue:
    """Tests for get_config_value function."""

    def test_simple_key(self):
        config = {"key": "value"}
        assert get_config_value(config, "key") == "value"

    def test_nested_key(self):
        config = {"level1": {"level2": {"level3": "deep_value"}}}
        assert get_config_value(config, "level1.level2.level3") == "deep_value"

    def test_missing_key_returns_default(self):
        config = {"key": "value"}
        assert get_config_value(config, "missing", "default") == "default"

    def test_partial_path_returns_default(self):
        config = {"level1": {"level2": "value"}}
        assert get_config_value(config, "level1.missing.key", "default") == "default"


class TestSetConfigValue:
    """Tests for set_config_value function."""

    def test_simple_key(self):
        config = {"key": "old"}
        set_config_value(config, "key", "new")
        assert config["key"] == "new"

    def test_nested_key(self):
        config = {"level1": {"level2": "old"}}
        set_config_value(config, "level1.level2", "new")
        assert config["level1"]["level2"] == "new"

    def test_creates_missing_intermediate_keys(self):
        config = {}
        set_config_value(config, "a.b.c", "value")
        assert config["a"]["b"]["c"] == "value"


class TestMergeCliConfig:
    """Tests for merge_cli_config function."""

    def test_simple_override(self):
        base = {"key": "base_value"}
        cli = {"key": "cli_value"}
        result = merge_cli_config(base, cli)
        assert result["key"] == "cli_value"

    def test_nested_override(self):
        base = {"level1": {"level2": "base_value", "other": "keep"}}
        cli = {"level1": {"level2": "cli_value"}}
        result = merge_cli_config(base, cli)
        assert result["level1"]["level2"] == "cli_value"
        assert result["level1"]["other"] == "keep"

    def test_dot_notation_override(self):
        base = {"resources": {"threads": 8, "memory": "16G"}}
        cli = {"resources.threads": 16}
        result = merge_cli_config(base, cli)
        assert result["resources"]["threads"] == 16
        assert result["resources"]["memory"] == "16G"

    def test_does_not_modify_original(self):
        base = {"key": "original"}
        cli = {"key": "override"}
        merge_cli_config(base, cli)
        assert base["key"] == "original"


class TestValidateConfigValue:
    """Tests for validate_config_value function."""

    def test_valid_string(self):
        validate_config_value("field", "value", str)

    def test_valid_int(self):
        validate_config_value("field", 42, int)

    def test_valid_bool(self):
        validate_config_value("field", True, bool)

    def test_invalid_type_raises(self):
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config_value("threads", "not_an_int", int)
        assert "threads" in str(exc_info.value)
        assert "Expected int" in str(exc_info.value)

    def test_pattern_validation_passes(self):
        validate_config_value("name", "valid_name", str, pattern=r"^[a-zA-Z][a-zA-Z0-9_]*$")

    def test_pattern_validation_fails(self):
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config_value("name", "123invalid", str, pattern=r"^[a-zA-Z][a-zA-Z0-9_]*$")
        assert "pattern" in str(exc_info.value)

    def test_minimum_validation_passes(self):
        validate_config_value("threads", 8, int, minimum=1)

    def test_minimum_validation_fails(self):
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config_value("threads", 0, int, minimum=1)
        assert ">= 1" in str(exc_info.value)

    def test_maximum_validation_passes(self):
        validate_config_value("threads", 64, int, maximum=128)

    def test_maximum_validation_fails(self):
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config_value("threads", 256, int, maximum=128)
        assert "<= 128" in str(exc_info.value)

    def test_enum_validation_passes(self):
        validate_config_value("method", "orthofinder", str, enum=["orthofinder"])

    def test_enum_validation_fails(self):
        with pytest.raises(ConfigValidationError) as exc_info:
            validate_config_value("method", "invalid", str, enum=["orthofinder"])
        assert "orthofinder" in str(exc_info.value)


class TestFormatValidationError:
    """Tests for format_validation_error function."""

    def test_includes_error_message(self):
        result = format_validation_error("Some error")
        assert "Some error" in result

    def test_includes_schema_path(self):
        result = format_validation_error("Error", "project.name")
        assert "project.name" in result

    def test_provides_suggestion_for_pattern(self):
        result = format_validation_error("pattern mismatch")
        assert "Suggestion" in result

    def test_provides_suggestion_for_required(self):
        result = format_validation_error("required property")
        assert "required" in result.lower()


class TestSchemaValidation:
    """Tests for JSON Schema validation."""

    def test_schema_file_exists(self, schema_path):
        assert schema_path.is_file(), f"Schema file not found: {schema_path}"

    def test_schema_is_valid_yaml(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        assert schema is not None
        assert "$schema" in schema
        assert schema["type"] == "object"

    def test_schema_has_required_fields(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        assert "project" in schema.get("required", [])
        assert "species" in schema.get("required", [])

    def test_schema_has_project_name_pattern(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        project_props = schema["properties"]["project"]["properties"]
        assert "pattern" in project_props["name"]
        assert project_props["name"]["pattern"] == "^[a-zA-Z][a-zA-Z0-9_]*$"

    def test_schema_has_threads_constraints(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        threads = schema["properties"]["resources"]["properties"]["threads"]
        assert threads["type"] == "integer"
        assert threads["minimum"] == 1
        assert threads["maximum"] == 128

    def test_schema_has_memory_pattern(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        memory = schema["properties"]["resources"]["properties"]["memory"]
        assert "pattern" in memory

    def test_schema_has_orthology_method_enum(self, schema_path):
        with open(schema_path) as f:
            schema = yaml.safe_load(f)
        method = schema["properties"]["analysis"]["properties"]["orthology"]["properties"]["method"]
        assert "enum" in method
        assert "orthofinder" in method["enum"]


class TestConfigFileValidation:
    """Tests for actual config file validation."""

    def test_default_config_has_required_structure(self, project_root):
        config_path = project_root / "config" / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "project" in config
        assert "name" in config["project"]
        assert "output_dir" in config["project"]
        assert "species" in config
        assert "analysis" in config
        assert "resources" in config

    def test_default_config_has_valid_project_name(self, project_root):
        config_path = project_root / "config" / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        import re
        name = config["project"]["name"]
        assert re.match(r"^[a-zA-Z][a-zA-Z0-9_]*$", name), f"Invalid project name: {name}"

    def test_default_config_has_valid_resources(self, project_root):
        config_path = project_root / "config" / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        threads = config["resources"]["threads"]
        assert isinstance(threads, int)
        assert 1 <= threads <= 128

        import re
        memory = config["resources"]["memory"]
        assert re.match(r"^[1-9][0-9]*[GMK]$", memory), f"Invalid memory format: {memory}"


class TestLoadConfig:
    """Tests for load_config function."""

    def test_load_valid_config(self, project_root):
        config_path = project_root / "config" / "config.yaml"
        config = load_config(config_path)
        assert isinstance(config, dict)
        assert "project" in config

    def test_load_config_file_not_found(self, tmp_path):
        nonexistent = tmp_path / "nonexistent.yaml"
        with pytest.raises(FileNotFoundError) as exc_info:
            load_config(nonexistent)
        assert "not found" in str(exc_info.value)

    def test_load_config_empty_file(self, tmp_path):
        empty_file = tmp_path / "empty.yaml"
        empty_file.write_text("")
        config = load_config(empty_file)
        assert config == {}

    def test_load_config_with_valid_config_fixture(self, tmp_path, valid_config):
        """Test loading a config file with known valid structure."""
        config_file = tmp_path / "test_config.yaml"
        config_file.write_text(yaml.dump(valid_config))
        loaded = load_config(config_file)
        assert loaded["project"]["name"] == "test_analysis"
        assert loaded["resources"]["threads"] == 8


class TestDeepCopyDict:
    """Tests for _deep_copy_dict function."""

    def test_shallow_dict(self):
        original = {"a": 1, "b": 2}
        copy = _deep_copy_dict(original)
        assert copy == original
        copy["a"] = 999
        assert original["a"] == 1

    def test_nested_dict(self):
        original = {"level1": {"level2": {"level3": "value"}}}
        copy = _deep_copy_dict(original)
        copy["level1"]["level2"]["level3"] = "modified"
        assert original["level1"]["level2"]["level3"] == "value"

    def test_dict_with_list(self):
        original = {"items": [{"name": "item1"}, {"name": "item2"}]}
        copy = _deep_copy_dict(original)
        copy["items"][0]["name"] = "modified"
        assert original["items"][0]["name"] == "item1"

    def test_dict_with_mixed_types(self):
        original = {
            "string": "value",
            "number": 42,
            "boolean": True,
            "none": None,
            "list": [1, 2, 3],
            "nested": {"key": "value"},
        }
        copy = _deep_copy_dict(original)
        assert copy == original
        copy["nested"]["key"] = "modified"
        assert original["nested"]["key"] == "value"


class TestConfigHelperLogic:
    """Tests for helper function logic that mirrors common.smk functions.

    These tests verify the logic used by Snakemake helper functions
    without requiring the Snakemake runtime.
    """

    def test_get_threads_logic(self, valid_config):
        """Test thread retrieval logic."""
        threads = get_config_value(valid_config, "resources.threads", 8)
        assert threads == 8
        assert isinstance(threads, int)

    def test_get_memory_logic(self, valid_config):
        """Test memory retrieval logic."""
        memory = get_config_value(valid_config, "resources.memory", "16G")
        assert memory == "16G"

    def test_get_output_dir_logic(self, valid_config):
        """Test output directory retrieval logic."""
        output_dir = get_config_value(valid_config, "project.output_dir", "./results")
        assert output_dir == "./results"

    def test_get_species_list_logic(self, valid_config):
        """Test species list retrieval logic."""
        species = get_config_value(valid_config, "species", [])
        assert len(species) == 1
        assert species[0]["name"] == "Microcebus_murinus"

    def test_get_species_names_logic(self, valid_config):
        """Test species names extraction logic."""
        species_list = get_config_value(valid_config, "species", [])
        names = [s["name"] for s in species_list]
        assert names == ["Microcebus_murinus"]

    def test_is_expression_enabled_logic(self, valid_config):
        """Test expression enabled flag logic."""
        enabled = get_config_value(valid_config, "analysis.expression.enabled", False)
        assert enabled is False

    def test_get_orthology_method_logic(self, valid_config):
        """Test orthology method retrieval logic."""
        method = get_config_value(valid_config, "analysis.orthology.method", "orthofinder")
        assert method == "orthofinder"
