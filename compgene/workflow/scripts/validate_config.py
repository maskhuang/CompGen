#!/usr/bin/env python3
"""
Configuration validation script for CompGene.

This script validates a configuration file against the JSON Schema
and reports any errors in a user-friendly format.

Usage:
    python validate_config.py config/config.yaml
    python validate_config.py --schema schemas/config.schema.yaml config.yaml
"""

import argparse
import sys
from pathlib import Path

import yaml

try:
    import jsonschema
    from jsonschema import Draft7Validator, ValidationError
except ImportError:
    print("Error: jsonschema package required. Install with: pip install jsonschema")
    sys.exit(1)


def load_yaml(path: Path) -> dict:
    """Load a YAML file and return its contents."""
    with open(path) as f:
        return yaml.safe_load(f)


def format_validation_error(error: ValidationError) -> str:
    """Format a validation error into a readable message."""
    path = " -> ".join(str(p) for p in error.absolute_path) or "(root)"
    return f"  - Path: {path}\n    Error: {error.message}"


def validate_business_rules(config: dict) -> list[str]:
    """
    Validate business rules that cannot be expressed in JSON Schema.

    Args:
        config: Configuration dictionary

    Returns:
        List of error messages (empty if valid)
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

    return errors


def validate_config(config_path: Path, schema_path: Path) -> list[str]:
    """
    Validate a configuration file against a JSON Schema and business rules.

    Args:
        config_path: Path to the configuration YAML file
        schema_path: Path to the JSON Schema YAML file

    Returns:
        List of error messages (empty if valid)
    """
    try:
        config = load_yaml(config_path)
    except yaml.YAMLError as e:
        return [f"YAML parsing error: {e}"]
    except FileNotFoundError:
        return [f"Configuration file not found: {config_path}"]

    try:
        schema = load_yaml(schema_path)
    except yaml.YAMLError as e:
        return [f"Schema parsing error: {e}"]
    except FileNotFoundError:
        return [f"Schema file not found: {schema_path}"]

    # Validate against schema
    validator = Draft7Validator(schema)
    errors = []

    for error in sorted(validator.iter_errors(config), key=lambda e: str(e.path)):
        errors.append(format_validation_error(error))

    # Validate business rules (only if schema validation passed)
    if not errors:
        business_errors = validate_business_rules(config)
        for err in business_errors:
            errors.append(f"  - Business Rule Error: {err}")

    return errors


def main():
    """Main entry point for the validation script."""
    parser = argparse.ArgumentParser(
        description="Validate CompGene configuration file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_config.py config/config.yaml
  python validate_config.py --schema custom_schema.yaml config.yaml
        """,
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to configuration YAML file",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path(__file__).parent.parent.parent / "schemas" / "config.schema.yaml",
        help="Path to JSON Schema file (default: schemas/config.schema.yaml)",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Only output errors, no success message",
    )

    args = parser.parse_args()

    # Validate
    errors = validate_config(args.config, args.schema)

    if errors:
        print("=" * 60)
        print("CONFIGURATION VALIDATION FAILED")
        print("=" * 60)
        print(f"\nConfig file: {args.config}")
        print(f"Schema file: {args.schema}\n")
        print("Errors found:\n")
        for error in errors:
            print(error)
            print()
        print("-" * 60)
        print("Suggestions:")
        print("  - Check config file for typos")
        print("  - Ensure all required fields are present")
        print("  - Verify data types match schema requirements")
        print("=" * 60)
        sys.exit(1)
    else:
        if not args.quiet:
            print(f"Configuration valid: {args.config}")
        sys.exit(0)


if __name__ == "__main__":
    main()
