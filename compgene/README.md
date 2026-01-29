# CompGene

Comparative genomics pipeline for ortholog analysis.

## Overview

CompGene is a Snakemake-based workflow for comparative genomics analysis, including:

- Data standardization and quality control (BUSCO)
- Orthology inference (OrthoFinder)
- Functional annotation (eggNOG-mapper)
- Absence verification (Liftoff)
- Presence/absence matrix generation
- Expression analysis (DESeq2)

## Requirements

- Python 3.11+
- Snakemake 9.0+
- Conda/Mamba (for environment management)

## Installation

```bash
# Create conda environment
conda create -n compgene python=3.11 snakemake=9.14 -c conda-forge -c bioconda

# Activate environment
conda activate compgene

# Install package (development mode)
pip install -e ".[dev]"
```

## Quick Start

1. Copy and edit the configuration file:

```bash
cp config/config.yaml config/my_config.yaml
# Edit my_config.yaml with your species and paths
```

2. Run the pipeline (from the compgene root directory):

```bash
snakemake --snakefile workflow/Snakefile --cores 8 --use-conda
```

Or use the local profile:

```bash
snakemake --snakefile workflow/Snakefile --profile profiles/local
```

3. Preview execution plan (dry run):

```bash
snakemake --snakefile workflow/Snakefile --dry-run
```

## Project Structure

```
compgene/
├── workflow/
│   ├── Snakefile         # Main entry point
│   ├── rules/            # Snakemake rules
│   ├── scripts/          # Python scripts
│   ├── adapters/         # Tool adapters
│   ├── lib/              # Shared utilities
│   └── envs/             # Conda environments
├── config/               # Configuration files
├── schemas/              # Config validation schemas
├── profiles/             # Execution profiles
├── tests/                # Test suite
├── resources/            # Static resources
├── results/              # Output directory
└── logs/                 # Log files
```

## Configuration

See `config/config.yaml` for available options.

### Configuration Validation

The pipeline validates your configuration at startup using JSON Schema. If validation fails, you'll see a clear error message explaining what needs to be fixed.

Required configuration:
- `species` - At least one species with id, name, annotation, and genome paths
- `output_dir` - Output directory path

### Command-Line Overrides

You can override configuration values from the command line:

```bash
# Override output directory
snakemake --snakefile workflow/Snakefile --config output_dir=/custom/path

# Override logging level
snakemake --snakefile workflow/Snakefile --config "logging={level: DEBUG}"

# Override resource defaults
snakemake --snakefile workflow/Snakefile --config "resources={default: {threads: 16}}"

# Multiple overrides
snakemake --snakefile workflow/Snakefile \
  --config output_dir=/custom/path \
  --config "logging={level: DEBUG}"
```

Command-line values take precedence over values in config.yaml.

## Checkpoint Resume and Dry-Run

### Dry-Run Mode

Preview the execution plan without running any rules:

```bash
# Standard dry-run
snakemake --snakefile workflow/Snakefile --dry-run

# Quiet mode (less output)
snakemake --snakefile workflow/Snakefile --dry-run --quiet

# Generate DAG visualization
snakemake --snakefile workflow/Snakefile --dag | dot -Tsvg > dag.svg
```

### Checkpoint Resume

If a pipeline run is interrupted, you can resume from where it left off:

```bash
# Resume from incomplete rules
snakemake --snakefile workflow/Snakefile --rerun-incomplete

# Continue with other rules if some fail
snakemake --snakefile workflow/Snakefile --keep-going
```

### Caching and Rerun Triggers

By default, Snakemake uses file modification time (mtime) to determine if rules need to be re-run. For stricter checking based on file content:

```bash
# Use checksums instead of mtime
snakemake --snakefile workflow/Snakefile --rerun-triggers checksum
```

Performance note: Cache hit detection is optimized to complete in < 1 second even for large file sets.

### Force Re-run

To force re-run of specific rules:

```bash
# Force re-run a specific rule
snakemake --snakefile workflow/Snakefile --forcerun rule_name

# Force re-run all rules
snakemake --snakefile workflow/Snakefile --forceall
```

## FAQ

### How do I know if a rule needs to be re-run?

Use dry-run mode to preview what would run:
```bash
snakemake --snakefile workflow/Snakefile --dry-run
```

If nothing is listed, all rules are up-to-date.

### Why is my rule re-running even though inputs haven't changed?

This can happen if:
1. Output files were manually deleted
2. File timestamps were modified (e.g., by copying files)
3. You're using checksum mode and file content changed

Use `--rerun-triggers mtime` (default) for timestamp-based checking or `--rerun-triggers checksum` for content-based checking.

### How do I clean up after an interrupted run?

Interrupted runs may leave behind `.tmp` files. The pipeline includes utilities to clean these up.

From the `compgene/` directory (after running `pip install -e .`):
```python
from workflow.lib.io import cleanup_temp_files
from pathlib import Path

removed = cleanup_temp_files(Path("results/"))
print(f"Cleaned up {len(removed)} temp files")
```

## License

MIT License
