"""Shared test fixtures for CompGene."""

from pathlib import Path

import pytest


@pytest.fixture
def project_root():
    """Return the project root directory."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def sample_config(project_root):
    """Return path to the sample config fixture."""
    return project_root / "tests" / "fixtures" / "sample_config.yaml"
