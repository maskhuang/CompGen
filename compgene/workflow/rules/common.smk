"""
Common utilities and shared functions for CompGene workflow.

This module provides helper functions used across multiple rules.
"""

import sys
import logging
from pathlib import Path

# Get configured logging level (defaults to INFO if not set)
_log_level_str = config.get("logging", {}).get("level", "INFO")
_log_level = getattr(logging, _log_level_str.upper(), logging.INFO)

# Configure logging with the configured level
logging.basicConfig(
    level=_log_level,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("compgene")


def validate_config_business_rules():
    """
    Perform custom business rule validation beyond JSON Schema.

    This validates constraints that cannot be expressed in JSON Schema,
    such as file existence checks and cross-field validations.

    Raises:
        ValueError: If any business rule validation fails
    """
    errors = []

    # Validate species IDs are unique
    species_ids = [s.get("id") for s in config.get("species", [])]
    if len(species_ids) != len(set(species_ids)):
        duplicates = [sid for sid in species_ids if species_ids.count(sid) > 1]
        errors.append(
            f"Duplicate species IDs found: {set(duplicates)}. "
            "Each species must have a unique ID."
        )

    # Note: Resource constraints (threads >= 1, memory_mb >= 1000) are validated
    # by JSON Schema. Business rules here focus on cross-field validations
    # that cannot be expressed in JSON Schema.

    # Report all errors
    if errors:
        error_msg = format_validation_errors(errors)
        raise ValueError(error_msg)

    logger.info("Configuration validation passed")


def format_validation_errors(errors: list) -> str:
    """
    Format validation errors into a clear, user-friendly message.

    Args:
        errors: List of error messages

    Returns:
        Formatted error string with suggestions
    """
    header = "\n" + "=" * 60 + "\n"
    header += "CONFIGURATION VALIDATION FAILED\n"
    header += "=" * 60 + "\n\n"

    body = "The following issues were found:\n\n"
    for i, error in enumerate(errors, 1):
        body += f"  {i}. {error}\n\n"

    footer = "-" * 60 + "\n"
    footer += "Suggestions:\n"
    footer += "  - Check config/config.yaml for typos\n"
    footer += "  - Ensure all required fields are present\n"
    footer += "  - Verify file paths exist and are accessible\n"
    footer += "  - Run 'snakemake --dry-run' to preview execution\n"
    footer += "=" * 60 + "\n"

    return header + body + footer


def log_config_override(key: str, file_value, cli_value):
    """
    Log when a configuration value is overridden via command line.

    Args:
        key: Configuration key that was overridden
        file_value: Original value from config file
        cli_value: New value from command line
    """
    logger.info(
        f"Config override: {key} = {cli_value} (was: {file_value})"
    )


def get_species_list():
    """
    Get list of species from configuration.

    Returns:
        list: Species identifiers from config

    Raises:
        ValueError: If species list is empty or missing
    """
    species = config.get("species", [])
    if not species:
        raise ValueError(
            "Configuration error: 'species' list is empty or missing. "
            "Please define at least one species in config/config.yaml"
        )
    return species


def get_output_dir():
    """
    Get the output directory from configuration.

    Returns:
        str: Path to output directory
    """
    return config.get("output_dir", "results")


def get_threads(rule_name: str, default: int = 1) -> int:
    """
    Get thread count for a specific rule from configuration.

    Priority order (highest to lowest):
    1. Rule-specific setting: resources.<rule_name>.threads
    2. Default setting: resources.default.threads
    3. Function default parameter

    Note: Command-line overrides (--config) are automatically merged
    by Snakemake before this function is called.

    Args:
        rule_name: Name of the rule
        default: Default thread count if not specified

    Returns:
        int: Number of threads to use
    """
    resources = config.get("resources", {})

    # Try rule-specific first
    rule_threads = resources.get(rule_name, {}).get("threads")
    if rule_threads is not None:
        return rule_threads

    # Fall back to default config
    default_threads = resources.get("default", {}).get("threads")
    if default_threads is not None:
        return default_threads

    # Use function default
    return default


def get_memory_mb(rule_name: str, default: int = 4000) -> int:
    """
    Get memory allocation for a specific rule from configuration.

    Priority order (highest to lowest):
    1. Rule-specific setting: resources.<rule_name>.memory_mb
    2. Default setting: resources.default.memory_mb
    3. Function default parameter

    Note: Command-line overrides (--config) are automatically merged
    by Snakemake before this function is called.

    Args:
        rule_name: Name of the rule
        default: Default memory in MB if not specified

    Returns:
        int: Memory allocation in MB
    """
    resources = config.get("resources", {})

    # Try rule-specific first
    rule_memory = resources.get(rule_name, {}).get("memory_mb")
    if rule_memory is not None:
        return rule_memory

    # Fall back to default config
    default_memory = resources.get("default", {}).get("memory_mb")
    if default_memory is not None:
        return default_memory

    # Use function default
    return default


def get_tool_config(tool_name: str) -> dict:
    """
    Get configuration for a specific tool.

    Args:
        tool_name: Name of the tool (e.g., 'orthofinder', 'busco')

    Returns:
        dict: Tool configuration, empty dict if not configured
    """
    return config.get("tools", {}).get(tool_name, {})


def get_logging_level() -> str:
    """
    Get the configured logging level.

    Returns:
        str: Logging level (DEBUG, INFO, WARNING, ERROR)
    """
    return config.get("logging", {}).get("level", "INFO")
