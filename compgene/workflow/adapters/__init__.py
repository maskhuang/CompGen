"""
CompGene Tool Adapters

This package provides adapter classes for external bioinformatics tools.
Each adapter implements the BaseAdapter interface for consistent tool integration.

Modules:
    base: Core adapter framework (BaseAdapter, ToolSpec, RunResult, AdapterContext)
"""

# Base adapter framework (Story 1.5)
from workflow.adapters.base import (
    ToolSpec,
    RunResult,
    AdapterContext,
    BaseAdapter,
    classify_common_errors,
)

__all__ = [
    # Data structures
    "ToolSpec",
    "RunResult",
    "AdapterContext",
    # Base class
    "BaseAdapter",
    # Helpers
    "classify_common_errors",
]
