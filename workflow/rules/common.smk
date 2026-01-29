# Module: common
# Shared functions and configuration loading
# Implementation: Story 1.1, 1.2

import sys
from pathlib import Path

# Add workflow/lib to Python path for imports
workflow_lib = Path(workflow.basedir).parent / "lib"
if str(workflow_lib) not in sys.path:
    sys.path.insert(0, str(workflow_lib))


def get_config_value(key: str, default=None):
    """Get a nested configuration value using dot notation.

    Args:
        key: Dot-separated key path (e.g., "analysis.orthology.method")
        default: Default value if key not found

    Returns:
        The configuration value or default
    """
    keys = key.split(".")
    value = config
    for k in keys:
        if isinstance(value, dict) and k in value:
            value = value[k]
        else:
            return default
    return value


def get_threads(rule_name: str = None) -> int:
    """Get the number of threads to use.

    Returns the configured threads value, capped by workflow.cores if set.

    Args:
        rule_name: Optional rule name for future per-rule configuration

    Returns:
        Number of threads to use
    """
    configured = get_config_value("resources.threads", 8)
    if workflow.cores:
        return min(configured, workflow.cores)
    return configured


def get_memory(rule_name: str = None) -> str:
    """Get the memory allocation string.

    Args:
        rule_name: Optional rule name for future per-rule configuration

    Returns:
        Memory string (e.g., "16G")
    """
    return get_config_value("resources.memory", "16G")


def get_output_dir() -> Path:
    """Get the configured output directory.

    Returns:
        Path to the output directory
    """
    return Path(get_config_value("project.output_dir", "./results"))


def get_species_list() -> list:
    """Get the list of species to analyze.

    Returns:
        List of species configuration dictionaries
    """
    return get_config_value("species", [])


def get_species_names() -> list:
    """Get the list of species names.

    Returns:
        List of species name strings
    """
    return [s["name"] for s in get_species_list()]


def is_expression_enabled() -> bool:
    """Check if expression analysis is enabled.

    Returns:
        True if expression analysis should be run
    """
    return get_config_value("analysis.expression.enabled", False)


def get_orthology_method() -> str:
    """Get the orthology inference method.

    Returns:
        Orthology method name (currently only "orthofinder")
    """
    return get_config_value("analysis.orthology.method", "orthofinder")
