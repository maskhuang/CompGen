"""
Common utilities and shared functions for CompGene workflow.

This module provides helper functions used across multiple rules.
"""

import sys
from pathlib import Path


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

    Args:
        rule_name: Name of the rule
        default: Default thread count if not specified

    Returns:
        int: Number of threads to use
    """
    resources = config.get("resources", {})
    return resources.get(rule_name, {}).get("threads", default)


def get_memory_mb(rule_name: str, default: int = 4000) -> int:
    """
    Get memory allocation for a specific rule from configuration.

    Args:
        rule_name: Name of the rule
        default: Default memory in MB if not specified

    Returns:
        int: Memory allocation in MB
    """
    resources = config.get("resources", {})
    return resources.get(rule_name, {}).get("memory_mb", default)
