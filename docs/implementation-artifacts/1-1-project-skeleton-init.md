# Story 1-1: Project Skeleton Initialization

**Status:** done
**Completed:** 2026-01-28

## Description

Create the initial Snakemake 9.x project structure with all required directories, configuration files, and testing infrastructure.

## Acceptance Criteria

- [x] Complete directory structure per architecture specification
- [x] Snakefile with min_version, configfile, validate directives
- [x] All rule modules exist (empty stubs)
- [x] Configuration schema and default config
- [x] Conda environment files for all tools
- [x] Local execution profile
- [x] pyproject.toml with project metadata
- [x] Test infrastructure validates project structure

## Tasks

All tasks completed - see implementation notes below.

## Dev Agent Record

### Implementation Notes
- Created full Snakemake 9.x project structure
- All 23 tests passing for project structure validation
- Code review fixes applied:
  - Added snakemake>=9.14 to pyproject.toml dependencies
  - Changed /tmp to $TMPDIR in profiles/local/config.yaml
  - Added content validation tests for config.yaml and conda env files
  - Removed duplicate fixture

### Files Created
- `workflow/Snakefile` - Main workflow entry point
- `workflow/rules/common.smk` - Shared functions
- `workflow/rules/downloads.smk` - Download rules stub
- `workflow/rules/orthology.smk` - Orthology rules stub
- `workflow/rules/expression.smk` - Expression rules stub
- `workflow/rules/reporting.smk` - Reporting rules stub
- `workflow/lib/__init__.py` - Library package
- `workflow/adapters/__init__.py` - Adapters package
- `config/config.yaml` - Default configuration
- `schemas/config.schema.yaml` - Configuration schema
- `profiles/local/config.yaml` - Local execution profile
- `envs/base.yaml` - Base conda environment
- `envs/orthofinder.yaml` - OrthoFinder environment
- `envs/star.yaml` - STAR aligner environment
- `envs/rsem.yaml` - RSEM environment
- `envs/deseq2.yaml` - DESeq2 environment
- `pyproject.toml` - Project metadata
- `.gitignore` - Git ignore patterns
- `tests/test_project_structure.py` - Structure validation tests

### Test Results
```
23 passed
```
