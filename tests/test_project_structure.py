"""Tests for Story 1.1: Project skeleton initialization."""

import pytest


class TestDirectoryStructure:
    """AC #1: Complete directory structure generated."""

    EXPECTED_DIRS = [
        "workflow",
        "workflow/rules",
        "workflow/scripts",
        "workflow/adapters",
        "workflow/lib",
        "workflow/envs",
        "config",
        "schemas",
        "profiles/local",
        "resources",
        "results/meta",
        "logs",
        "tests/integration/test_data",
        "tests/fixtures",
    ]

    def test_all_directories_exist(self, project_root):
        for d in self.EXPECTED_DIRS:
            path = project_root / d
            assert path.is_dir(), f"Missing directory: {d}"


class TestSnakefile:
    """AC #2: Snakefile contains min_version and validate."""

    def test_snakefile_exists(self, project_root):
        snakefile = project_root / "workflow" / "Snakefile"
        assert snakefile.is_file(), "workflow/Snakefile not found"

    def test_snakefile_has_min_version(self, project_root):
        content = (project_root / "workflow" / "Snakefile").read_text()
        assert 'min_version("9.0")' in content

    def test_snakefile_has_validate(self, project_root):
        content = (project_root / "workflow" / "Snakefile").read_text()
        assert "validate(" in content

    def test_snakefile_has_configfile(self, project_root):
        content = (project_root / "workflow" / "Snakefile").read_text()
        assert "configfile:" in content

    def test_snakefile_includes_all_rule_modules(self, project_root):
        content = (project_root / "workflow" / "Snakefile").read_text()
        modules = [
            "common",
            "standardize",
            "qc",
            "orthology",
            "annotation",
            "validation",
            "matrices",
            "expression",
            "reporting",
            "ingest",
        ]
        for module in modules:
            assert f'include: "rules/{module}.smk"' in content, (
                f"Missing include for {module}.smk"
            )

    def test_snakefile_has_rule_all(self, project_root):
        content = (project_root / "workflow" / "Snakefile").read_text()
        assert "rule all:" in content


class TestRuleFiles:
    """AC #2: Rule files exist for all modules."""

    EXPECTED_RULES = [
        "common.smk",
        "standardize.smk",
        "qc.smk",
        "orthology.smk",
        "annotation.smk",
        "validation.smk",
        "matrices.smk",
        "expression.smk",
        "reporting.smk",
        "ingest.smk",
    ]

    def test_all_rule_files_exist(self, project_root):
        rules_dir = project_root / "workflow" / "rules"
        for rule_file in self.EXPECTED_RULES:
            path = rules_dir / rule_file
            assert path.is_file(), f"Missing rule file: {rule_file}"


class TestPyprojectToml:
    """AC #3: pyproject.toml defines project metadata and dependencies."""

    def test_pyproject_exists(self, project_root):
        assert (project_root / "pyproject.toml").is_file()

    def test_pyproject_has_required_content(self, project_root):
        content = (project_root / "pyproject.toml").read_text()
        assert 'name = "compgene"' in content
        assert 'requires-python = ">=3.11"' in content
        assert "snakemake" in content
        assert "pytest" in content
        assert "ruff" in content


class TestPythonPackages:
    """AC #1: Python packages have __init__.py."""

    def test_adapters_package(self, project_root):
        assert (project_root / "workflow" / "adapters" / "__init__.py").is_file()

    def test_lib_package(self, project_root):
        assert (project_root / "workflow" / "lib" / "__init__.py").is_file()


class TestConfigFiles:
    """AC #1, #2: Configuration files exist and have required content."""

    def test_config_yaml_exists(self, project_root):
        assert (project_root / "config" / "config.yaml").is_file()

    def test_config_yaml_has_required_sections(self, project_root):
        content = (project_root / "config" / "config.yaml").read_text()
        assert "project:" in content
        assert "name:" in content
        assert "output_dir:" in content
        assert "species:" in content
        assert "analysis:" in content
        assert "resources:" in content

    def test_schema_yaml_exists(self, project_root):
        assert (project_root / "schemas" / "config.schema.yaml").is_file()

    def test_schema_yaml_has_required_structure(self, project_root):
        content = (project_root / "schemas" / "config.schema.yaml").read_text()
        assert "$schema:" in content
        assert "type: object" in content
        assert "required:" in content
        assert "project" in content
        assert "species" in content

    def test_gitignore_exists(self, project_root):
        assert (project_root / ".gitignore").is_file()

    def test_gitignore_excludes_runtime(self, project_root):
        content = (project_root / ".gitignore").read_text()
        assert "results/" in content
        assert "logs/" in content
        assert ".snakemake/" in content
        assert "*.pyc" in content

    def test_local_profile_exists(self, project_root):
        assert (project_root / "profiles" / "local" / "config.yaml").is_file()


class TestCondaEnvironments:
    """AC #1: Conda environment files exist."""

    EXPECTED_ENVS = [
        "base.yaml",
        "orthofinder.yaml",
        "eggnog.yaml",
        "liftoff.yaml",
        "busco.yaml",
        "deseq2.yaml",
    ]

    def test_all_env_files_exist(self, project_root):
        envs_dir = project_root / "workflow" / "envs"
        for env_file in self.EXPECTED_ENVS:
            path = envs_dir / env_file
            assert path.is_file(), f"Missing env file: {env_file}"

    def test_base_env_has_channels(self, project_root):
        content = (project_root / "workflow" / "envs" / "base.yaml").read_text()
        assert "conda-forge" in content
        assert "bioconda" in content

    def test_all_env_files_have_channels(self, project_root):
        envs_dir = project_root / "workflow" / "envs"
        for env_file in self.EXPECTED_ENVS:
            content = (envs_dir / env_file).read_text()
            assert "conda-forge" in content, f"{env_file} missing conda-forge channel"
            assert "bioconda" in content, f"{env_file} missing bioconda channel"

    def test_tool_env_files_have_dependencies(self, project_root):
        tool_envs = {
            "orthofinder.yaml": "orthofinder",
            "eggnog.yaml": "eggnog-mapper",
            "liftoff.yaml": "liftoff",
            "busco.yaml": "busco",
            "deseq2.yaml": "deseq2",
        }
        envs_dir = project_root / "workflow" / "envs"
        for env_file, tool in tool_envs.items():
            content = (envs_dir / env_file).read_text()
            assert tool in content, f"{env_file} missing {tool} dependency"
