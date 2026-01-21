#!/usr/bin/env python3
"""
CompGene Adapter Runner CLI - Entry Point for Snakemake Rules.

This script provides a command-line interface for executing external tools
through the BaseAdapter framework. It is designed to be called from
Snakemake rules.

Usage:
    python run_adapter.py <adapter_name> --inputs '{"key": "path"}' --outputs '{"key": "path"}'
        --config '{"key": "value"}' --wildcards '{"key": "value"}' --threads 4
        [--meta-dir results/meta] [--log-dir logs]

Example:
    python run_adapter.py orthofinder \
        --inputs '{"proteins": "data/proteins/"}' \
        --outputs '{"orthogroups": "results/orthogroups.tsv"}' \
        --config '{}' \
        --wildcards '{"species_set": "mammals"}' \
        --threads 8

Source: ADR-002 Tool Adaptation Layer
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

from workflow.lib.errors import CompGeneError, exit_with_error
from workflow.lib.logging import create_logger
from workflow.lib.runner import AdapterRunner, RetryConfig
from workflow.adapters.base import AdapterContext


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Execute an adapter through the CompGene runner framework.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        "adapter_name",
        help="Name of the adapter to run (e.g., 'orthofinder', 'busco')"
    )

    parser.add_argument(
        "--inputs",
        required=True,
        help="JSON dict mapping input names to paths"
    )

    parser.add_argument(
        "--outputs",
        required=True,
        help="JSON dict mapping output names to paths"
    )

    parser.add_argument(
        "--config",
        default="{}",
        help="JSON dict of configuration values (default: {})"
    )

    parser.add_argument(
        "--wildcards",
        default="{}",
        help="JSON dict of wildcard values (default: {})"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use (default: 1)"
    )

    parser.add_argument(
        "--meta-dir",
        default="results/meta",
        help="Directory for audit metadata (default: results/meta)"
    )

    parser.add_argument(
        "--log-dir",
        default="logs",
        help="Directory for log files (default: logs)"
    )

    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Maximum retry attempts for retryable errors (default: 3)"
    )

    parser.add_argument(
        "--grace-period",
        type=float,
        default=10.0,
        help="Seconds to wait after SIGTERM before SIGKILL (default: 10.0)"
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)"
    )

    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI color output"
    )

    return parser.parse_args()


def load_adapter(adapter_name: str):
    """
    Dynamically load an adapter class by name.

    Args:
        adapter_name: Name of the adapter (e.g., 'orthofinder').

    Returns:
        Instantiated adapter.

    Raises:
        CompGeneError: If adapter module or class not found.
    """
    # Convert adapter_name to module and class name
    # e.g., "orthofinder" -> "workflow.adapters.orthofinder.OrthoFinderAdapter"
    module_name = f"workflow.adapters.{adapter_name}"
    class_name = f"{adapter_name.title().replace('_', '')}Adapter"

    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as e:
        from workflow.lib.errors import ErrorCode
        raise CompGeneError(
            ErrorCode.E_TOOL_NOT_FOUND,
            f"Adapter module not found: {module_name}",
            details=str(e)
        )

    try:
        adapter_class = getattr(module, class_name)
    except AttributeError:
        from workflow.lib.errors import ErrorCode
        raise CompGeneError(
            ErrorCode.E_TOOL_NOT_FOUND,
            f"Adapter class not found: {class_name} in {module_name}",
            details=f"Available: {dir(module)}"
        )

    return adapter_class()


def parse_json_arg(value: str, arg_name: str) -> dict:
    """
    Parse a JSON string argument into a dict.

    Args:
        value: JSON string.
        arg_name: Argument name for error messages.

    Returns:
        Parsed dictionary.

    Raises:
        CompGeneError: If JSON is invalid.
    """
    try:
        result = json.loads(value)
        if not isinstance(result, dict):
            from workflow.lib.errors import ErrorCode
            raise CompGeneError(
                ErrorCode.E_INPUT_FORMAT,
                f"Argument --{arg_name} must be a JSON object, not {type(result).__name__}"
            )
        return result
    except json.JSONDecodeError as e:
        from workflow.lib.errors import ErrorCode
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Invalid JSON in --{arg_name}",
            details=str(e)
        )


def paths_from_dict(d: dict) -> dict[str, Path]:
    """Convert string values to Path objects."""
    return {k: Path(v) for k, v in d.items()}


def main() -> int:
    """
    Main entry point for the run_adapter CLI.

    Returns:
        Exit code (0 for success, non-zero for errors).
    """
    args = parse_args()

    try:
        # Parse JSON arguments
        inputs = paths_from_dict(parse_json_arg(args.inputs, "inputs"))
        outputs = paths_from_dict(parse_json_arg(args.outputs, "outputs"))
        config = parse_json_arg(args.config, "config")
        wildcards = parse_json_arg(args.wildcards, "wildcards")

        # Load adapter
        adapter = load_adapter(args.adapter_name)

        # Create logger
        logger = create_logger(
            rule=adapter.spec.name,
            wildcards=wildcards,
            log_dir=Path(args.log_dir),
            level=args.log_level,
            use_color=not args.no_color
        )

        # Create context
        ctx = AdapterContext(
            inputs=inputs,
            outputs=outputs,
            config=config,
            wildcards=wildcards,
            threads=args.threads,
            logger=logger
        )

        # Create runner
        retry_config = RetryConfig(max_retries=args.max_retries)
        runner = AdapterRunner(
            retry_config=retry_config,
            logger=logger,
            grace_period=args.grace_period,
            meta_dir=Path(args.meta_dir)
        )

        # Execute
        result = runner.run(adapter, ctx)

        # Log success summary
        logger.info(
            f"Execution complete: {len(result.outputs)} outputs, "
            f"summary: {result.summary}"
        )

        return 0

    except CompGeneError as e:
        exit_with_error(e, use_color=not args.no_color)
        return e.to_exit_code()  # Never reached, but for type checker

    except Exception as e:
        # Unexpected error
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
