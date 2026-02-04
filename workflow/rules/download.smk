# Module: download
# Rules for downloading data from external sources (NCBI)
# Implementation: Story 2.3
#
# This module handles downloading genomic data from NCBI:
#   results/downloads/ncbi/{accession}/
#   ├── genome.fna.gz           # Compressed genome sequence
#   ├── annotation.gff.gz       # Compressed GFF annotation
#   └── download_manifest.json  # Download metadata and checksums

import os
from pathlib import Path

from workflow.lib.ncbi import (
    download_ncbi_dataset,
    is_ncbi_accession,
    is_cached,
    parse_ncbi_accession,
)


# =============================================================================
# Helper Functions
# =============================================================================

def get_ncbi_species():
    """
    Get list of species that use NCBI accessions.

    Returns species configs where annotation starts with NCBI: or matches
    GCF_/GCA_ accession pattern.
    """
    ncbi_species = []
    for species in get_species_list():
        annotation = str(species.get("annotation", ""))
        if is_ncbi_accession(annotation):
            ncbi_species.append(species)
    return ncbi_species


def get_ncbi_accession_for_species(species_name: str) -> str:
    """
    Get the NCBI accession for a species.

    Args:
        species_name: The species identifier.

    Returns:
        The NCBI accession (without NCBI: prefix).
    """
    for species in get_species_list():
        if species["name"] == species_name:
            annotation = str(species.get("annotation", ""))
            if is_ncbi_accession(annotation):
                # Parse to get normalized accession
                parsed = parse_ncbi_accession(annotation)
                return parsed.full_accession
    return ""


def get_ncbi_download_dir() -> Path:
    """
    Get the NCBI download cache directory.

    Returns:
        Path to the downloads/ncbi directory.
    """
    return get_output_dir() / "downloads" / "ncbi"


def get_species_ncbi_genome(wildcards):
    """
    Get the downloaded genome path for an NCBI species.

    Args:
        wildcards: Snakemake wildcards with 'species' attribute.

    Returns:
        Path to downloaded genome file.
    """
    accession = get_ncbi_accession_for_species(wildcards.species)
    if accession:
        return str(get_ncbi_download_dir() / accession / "genome.fna.gz")
    return ""


def get_species_ncbi_annotation(wildcards):
    """
    Get the downloaded annotation path for an NCBI species.

    Args:
        wildcards: Snakemake wildcards with 'species' attribute.

    Returns:
        Path to downloaded annotation file.
    """
    accession = get_ncbi_accession_for_species(wildcards.species)
    if accession:
        return str(get_ncbi_download_dir() / accession / "annotation.gff.gz")
    return ""


def get_ncbi_download_targets():
    """
    Generate all NCBI download manifest targets.

    Returns:
        List of manifest file paths for all NCBI species.
    """
    targets = []
    for species in get_ncbi_species():
        annotation = str(species.get("annotation", ""))
        parsed = parse_ncbi_accession(annotation)
        accession = parsed.full_accession
        manifest_path = get_ncbi_download_dir() / accession / "download_manifest.json"
        targets.append(str(manifest_path))
    return targets


def get_ncbi_species_names():
    """
    Get names of species using NCBI data.

    Returns:
        List of species name strings.
    """
    return [s["name"] for s in get_ncbi_species()]


# =============================================================================
# Download Rules
# =============================================================================

rule ncbi_download:
    """
    Download genomic data for a single NCBI accession.

    Downloads genome, annotation, and protein files from NCBI FTP.
    Uses caching to avoid redundant downloads.

    AC2: Automatic download from NCBI
    AC3: Local cache mechanism
    AC4: Rate limit compliance
    AC5: Error handling and retry
    """
    output:
        manifest="{output_dir}/downloads/ncbi/{accession}/download_manifest.json",
        genome="{output_dir}/downloads/ncbi/{accession}/genome.fna.gz",
        annotation="{output_dir}/downloads/ncbi/{accession}/annotation.gff.gz",
    log:
        "{output_dir}/logs/ncbi_download/{accession}.log"
    params:
        force=get_config_value("ncbi.force_download", False),
        api_key=os.environ.get("NCBI_API_KEY")
    run:
        import logging
        from pathlib import Path

        # Set up file handler for this rule
        logger = logging.getLogger(f"ncbi_download.{wildcards.accession}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            output_dir = Path(output.manifest).parent.parent
            accession = wildcards.accession

            logger.info(f"Downloading NCBI dataset: {accession}")
            logger.info(f"  Output directory: {output_dir}")
            logger.info(f"  Force download: {params.force}")
            logger.info(f"  API key: {'set' if params.api_key else 'not set'}")

            result = download_ncbi_dataset(
                accession=accession,
                output_dir=output_dir,
                force=params.force,
                api_key=params.api_key,
                download_genome=True,
                download_annotation=True,
                download_protein=False
            )

            if result.from_cache:
                logger.info("Using cached data")
            else:
                logger.info("Download complete")

            logger.info(f"  Genome: {result.genome_path}")
            logger.info(f"  Annotation: {result.annotation_path}")
            if result.protein_path:
                logger.info(f"  Protein: {result.protein_path}")
            logger.info(f"  Manifest: {result.manifest_path}")

        except Exception as e:
            logger.error(f"Download failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


rule ncbi_download_all:
    """
    Download all NCBI datasets for configured species.

    This is an aggregate rule that triggers downloads for all species
    that use NCBI accessions.
    """
    input:
        manifests=get_ncbi_download_targets()
    output:
        touch("{output_dir}/downloads/ncbi/.downloads_complete")


# =============================================================================
# Standardization Integration Rules
# =============================================================================

rule standardize_ncbi_genome:
    """
    Standardize a downloaded NCBI genome to the standard output path.

    This rule is triggered for species with NCBI annotations and
    copies the downloaded genome to the standardized location.

    AC6: Integration with standardization
    """
    input:
        genome=get_species_ncbi_genome
    output:
        genome="{output_dir}/standardized/{species}/genome.fa.gz"
    log:
        "{output_dir}/logs/standardize_ncbi_genome/{species}.log"
    wildcard_constraints:
        species="|".join(get_ncbi_species_names()) if get_ncbi_species_names() else "NEVER_MATCH"
    run:
        import logging
        import shutil
        from pathlib import Path

        # Set up file handler
        logger = logging.getLogger(f"standardize_ncbi_genome.{wildcards.species}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            input_path = Path(input.genome)
            output_path = Path(output.genome)

            logger.info(f"Standardizing NCBI genome: {input_path} -> {output_path}")

            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Copy file (already gzipped from NCBI)
            shutil.copy2(input_path, output_path)

            logger.info("NCBI genome standardization complete")
        except Exception as e:
            logger.error(f"NCBI genome standardization failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


rule standardize_ncbi_annotation:
    """
    Standardize a downloaded NCBI annotation to the standard output path.

    This rule is triggered for species with NCBI annotations and
    copies the downloaded annotation to the standardized location.

    AC6: Integration with standardization
    """
    input:
        annotation=get_species_ncbi_annotation
    output:
        annotation="{output_dir}/standardized/{species}/annotation.gff3.gz"
    log:
        "{output_dir}/logs/standardize_ncbi_annotation/{species}.log"
    wildcard_constraints:
        species="|".join(get_ncbi_species_names()) if get_ncbi_species_names() else "NEVER_MATCH"
    run:
        import logging
        import shutil
        from pathlib import Path

        # Set up file handler
        logger = logging.getLogger(f"standardize_ncbi_annotation.{wildcards.species}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            input_path = Path(input.annotation)
            output_path = Path(output.annotation)

            logger.info(f"Standardizing NCBI annotation: {input_path} -> {output_path}")

            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Copy file (already gzipped from NCBI)
            shutil.copy2(input_path, output_path)

            logger.info("NCBI annotation standardization complete")
        except Exception as e:
            logger.error(f"NCBI annotation standardization failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


rule standardize_all_ncbi:
    """
    Standardize all species with NCBI data.

    This is a checkpoint rule that triggers standardization for all
    species that use NCBI accessions.
    """
    input:
        proteins=expand(
            "{output_dir}/standardized/{species}/proteins.longest.fa.gz",
            output_dir=get_output_dir(),
            species=get_ncbi_species_names()
        ) if get_ncbi_species_names() else [],
        mapping=expand(
            "{output_dir}/standardized/{species}/representative_transcripts.tsv",
            output_dir=get_output_dir(),
            species=get_ncbi_species_names()
        ) if get_ncbi_species_names() else []
    output:
        touch("{output_dir}/standardized/.standardize_ncbi_complete")
