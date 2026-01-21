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

## License

MIT License
