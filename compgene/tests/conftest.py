"""
Pytest configuration and shared fixtures for CompGene tests.

This module provides reusable test fixtures and configuration
for the entire test suite.
"""

import yaml
from pathlib import Path
from unittest.mock import MagicMock

# Import pytest only when available (allows factory functions to work standalone)
try:
    import pytest
    PYTEST_AVAILABLE = True
except ImportError:
    PYTEST_AVAILABLE = False


# =============================================================================
# Factory Functions (for generating test data) - No pytest required
# =============================================================================

def create_species(
    id: str = "test",
    name: str = "Test Species",
    annotation: str = "test.gff3",
    genome: str = "test.fa",
    **overrides
) -> dict:
    """
    Factory function to create a species configuration.

    Args:
        id: Species short code
        name: Species full name
        annotation: Path to annotation file
        genome: Path to genome file
        **overrides: Additional fields to override

    Returns:
        dict: Species configuration
    """
    species = {
        "id": id,
        "name": name,
        "annotation": annotation,
        "genome": genome
    }
    species.update(overrides)
    return species


def create_config(
    species_count: int = 1,
    output_dir: str = "results",
    include_tools: bool = False,
    include_resources: bool = False,
    include_logging: bool = False,
    **overrides
) -> dict:
    """
    Factory function to create a configuration dict.

    Args:
        species_count: Number of species to generate
        output_dir: Output directory path
        include_tools: Include tool configuration
        include_resources: Include resource configuration
        include_logging: Include logging configuration
        **overrides: Additional fields to override

    Returns:
        dict: Configuration dictionary
    """
    config = {
        "species": [
            create_species(id=f"sp{i}", name=f"Species {i}")
            for i in range(1, species_count + 1)
        ],
        "output_dir": output_dir
    }

    if include_tools:
        config["tools"] = {
            "orthofinder": {"threads": 8},
            "eggnog": {"database": "eukaryota"},
            "busco": {"lineage": "eukaryota_odb10"}
        }

    if include_resources:
        config["resources"] = {
            "default": {"threads": 4, "memory_mb": 8000}
        }

    if include_logging:
        config["logging"] = {"level": "INFO"}

    config.update(overrides)
    return config


# =============================================================================
# Pytest Fixtures (only defined when pytest is available)
# =============================================================================

if PYTEST_AVAILABLE:
    @pytest.fixture
    def project_root() -> Path:
        """Return the project root directory."""
        return Path(__file__).parent.parent

    @pytest.fixture
    def fixtures_dir() -> Path:
        """Return the test fixtures directory."""
        return Path(__file__).parent / "fixtures"

    @pytest.fixture
    def schema_path(project_root) -> Path:
        """Return the path to the config schema."""
        return project_root / "schemas" / "config.schema.yaml"

    @pytest.fixture
    def valid_config(fixtures_dir) -> dict:
        """Load the valid test configuration."""
        with open(fixtures_dir / "valid_config.yaml") as f:
            return yaml.safe_load(f)

    @pytest.fixture
    def invalid_config(fixtures_dir) -> dict:
        """Load the invalid test configuration."""
        with open(fixtures_dir / "invalid_config.yaml") as f:
            return yaml.safe_load(f)

    @pytest.fixture
    def minimal_config() -> dict:
        """Return a minimal valid configuration."""
        return {
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

    @pytest.fixture
    def full_config() -> dict:
        """Return a fully-populated configuration for testing."""
        return {
            "species": [
                {
                    "id": "mmur",
                    "name": "Microcebus murinus",
                    "annotation": "mmur.gff3",
                    "genome": "mmur.fa"
                },
                {
                    "id": "lcat",
                    "name": "Lemur catta",
                    "annotation": "lcat.gff3",
                    "genome": "lcat.fa"
                }
            ],
            "output_dir": "results",
            "tools": {
                "orthofinder": {"threads": 8},
                "eggnog": {"database": "eukaryota"},
                "busco": {"lineage": "eukaryota_odb10"},
                "liftoff": {"min_coverage": 0.5, "min_identity": 0.5},
                "deseq2": {"alpha": 0.05}
            },
            "resources": {
                "default": {"threads": 4, "memory_mb": 8000},
                "orthology_infer": {"threads": 8, "memory_mb": 16000}
            },
            "logging": {"level": "INFO"}
        }

    @pytest.fixture
    def mock_snakemake_config(full_config):
        """Create a mock for Snakemake's global config object."""
        return full_config

    @pytest.fixture
    def mock_logger():
        """Create a mock logger for testing logging calls."""
        return MagicMock()

    @pytest.fixture
    def schema(schema_path) -> dict:
        """Load the configuration schema."""
        with open(schema_path) as f:
            return yaml.safe_load(f)

    @pytest.fixture
    def species_factory():
        """Provide the species factory function as a fixture."""
        return create_species

    @pytest.fixture
    def config_factory():
        """Provide the config factory function as a fixture."""
        return create_config
