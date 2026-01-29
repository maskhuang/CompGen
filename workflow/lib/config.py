"""Configuration validation and utilities for CompGene pipeline.

This module provides functions for loading, validating, and accessing
configuration values with user-friendly error messages.
"""

import re
from pathlib import Path
from typing import Any

import yaml


class ConfigValidationError(Exception):
    """Raised when configuration validation fails."""

    def __init__(self, field: str, value: Any, message: str, suggestion: str = ""):
        self.field = field
        self.value = value
        self.message = message
        self.suggestion = suggestion
        super().__init__(self.format_error())

    def format_error(self) -> str:
        """Format the error as a user-friendly message."""
        lines = [
            "Configuration validation failed:",
            f"  Field: {self.field}",
            f"  Value: {self.value!r}",
            f"  Error: {self.message}",
        ]
        if self.suggestion:
            lines.append(f"  Suggestion: {self.suggestion}")
        return "\n".join(lines)


def get_config_value(config: dict, key: str, default: Any = None) -> Any:
    """Get a nested configuration value using dot notation.

    Args:
        config: Configuration dictionary
        key: Dot-separated key path (e.g., "analysis.orthology.method")
        default: Default value if key not found

    Returns:
        The configuration value or default

    Examples:
        >>> config = {"analysis": {"orthology": {"method": "orthofinder"}}}
        >>> get_config_value(config, "analysis.orthology.method")
        'orthofinder'
        >>> get_config_value(config, "missing.key", "default")
        'default'
    """
    keys = key.split(".")
    value = config
    for k in keys:
        if isinstance(value, dict) and k in value:
            value = value[k]
        else:
            return default
    return value


def set_config_value(config: dict, key: str, value: Any) -> None:
    """Set a nested configuration value using dot notation.

    Args:
        config: Configuration dictionary (modified in place)
        key: Dot-separated key path (e.g., "resources.threads")
        value: Value to set

    Examples:
        >>> config = {"resources": {"threads": 8}}
        >>> set_config_value(config, "resources.threads", 16)
        >>> config["resources"]["threads"]
        16
    """
    keys = key.split(".")
    current = config
    for k in keys[:-1]:
        if k not in current:
            current[k] = {}
        current = current[k]
    current[keys[-1]] = value


def merge_cli_config(base_config: dict, cli_overrides: dict) -> dict:
    """Merge CLI configuration overrides into base configuration.

    CLI overrides can use dot notation for nested keys.
    CLI values take precedence over base configuration.

    Args:
        base_config: Base configuration dictionary
        cli_overrides: Dictionary of CLI overrides (may use dot notation keys)

    Returns:
        Merged configuration dictionary

    Examples:
        >>> base = {"resources": {"threads": 8, "memory": "16G"}}
        >>> cli = {"resources": {"threads": 16}}
        >>> result = merge_cli_config(base, cli)
        >>> result["resources"]["threads"]
        16
    """
    result = _deep_copy_dict(base_config)

    for key, value in cli_overrides.items():
        if "." in key:
            set_config_value(result, key, value)
        elif isinstance(value, dict) and key in result and isinstance(result[key], dict):
            result[key] = merge_cli_config(result[key], value)
        else:
            result[key] = value

    return result


def _deep_copy_dict(d: dict) -> dict:
    """Create a deep copy of a dictionary."""
    result = {}
    for key, value in d.items():
        if isinstance(value, dict):
            result[key] = _deep_copy_dict(value)
        elif isinstance(value, list):
            result[key] = [_deep_copy_dict(v) if isinstance(v, dict) else v for v in value]
        else:
            result[key] = value
    return result


def format_validation_error(error_message: str, schema_path: str = "") -> str:
    """Format a JSON Schema validation error as a user-friendly message.

    Args:
        error_message: Raw error message from jsonschema
        schema_path: Path in the schema where validation failed

    Returns:
        Formatted error message with suggestions
    """
    suggestions = {
        "pattern": "Check the format matches the required pattern",
        "required": "Add the missing required field to your configuration",
        "type": "Ensure the value has the correct type",
        "enum": "Use one of the allowed values",
        "minimum": "Use a value greater than or equal to the minimum",
        "maximum": "Use a value less than or equal to the maximum",
        "minLength": "Provide a non-empty value",
        "additionalProperties": "Remove unknown fields from configuration",
    }

    suggestion = ""
    for key, msg in suggestions.items():
        if key in error_message.lower():
            suggestion = msg
            break

    lines = ["Configuration validation failed:"]
    if schema_path:
        lines.append(f"  Location: {schema_path}")
    lines.append(f"  Error: {error_message}")
    if suggestion:
        lines.append(f"  Suggestion: {suggestion}")

    return "\n".join(lines)


def load_config(config_path: str | Path) -> dict:
    """Load a YAML configuration file.

    Args:
        config_path: Path to the configuration file

    Returns:
        Parsed configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If YAML parsing fails
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path) as f:
        return yaml.safe_load(f) or {}


def validate_config_value(field: str, value: Any, expected_type: type, **constraints) -> None:
    """Validate a single configuration value.

    Args:
        field: Field name for error messages
        value: Value to validate
        expected_type: Expected Python type
        **constraints: Additional constraints (pattern, minimum, maximum, enum)

    Raises:
        ConfigValidationError: If validation fails
    """
    if not isinstance(value, expected_type):
        raise ConfigValidationError(
            field=field,
            value=value,
            message=f"Expected {expected_type.__name__}, got {type(value).__name__}",
            suggestion=f"Change the value to a {expected_type.__name__}",
        )

    if "pattern" in constraints and isinstance(value, str):
        if not re.match(constraints["pattern"], value):
            raise ConfigValidationError(
                field=field,
                value=value,
                message=f"Value does not match pattern: {constraints['pattern']}",
                suggestion="Check the format matches the required pattern",
            )

    if "minimum" in constraints and isinstance(value, (int, float)):
        if value < constraints["minimum"]:
            raise ConfigValidationError(
                field=field,
                value=value,
                message=f"Value must be >= {constraints['minimum']}",
                suggestion=f"Use a value of at least {constraints['minimum']}",
            )

    if "maximum" in constraints and isinstance(value, (int, float)):
        if value > constraints["maximum"]:
            raise ConfigValidationError(
                field=field,
                value=value,
                message=f"Value must be <= {constraints['maximum']}",
                suggestion=f"Use a value of at most {constraints['maximum']}",
            )

    if "enum" in constraints:
        if value not in constraints["enum"]:
            raise ConfigValidationError(
                field=field,
                value=value,
                message=f"Value must be one of: {constraints['enum']}",
                suggestion=f"Choose from: {', '.join(str(v) for v in constraints['enum'])}",
            )
