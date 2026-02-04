# Module: standardize
# Rules for data standardization
# Implementation: Story 2.2
#
# This module standardizes local annotation files to a unified directory structure:
#   results/standardized/{species}/
#   ├── genome.fa.gz              # Compressed genome sequence
#   ├── annotation.gff3.gz        # Compressed GFF3 annotation
#   ├── proteins.longest.fa.gz    # Representative transcript proteins
#   └── representative_transcripts.tsv  # Representative transcript mapping

from pathlib import Path

from workflow.lib.standardize import (
    standardize_genome,
    standardize_annotation,
    extract_proteins,
)


# =============================================================================
# Helper Functions
# =============================================================================

def get_species_genome(wildcards):
    """Get the genome path for a species from config."""
    for species in get_species_list():
        if species["name"] == wildcards.species:
            return species.get("assembly", "")
    return ""


def get_species_annotation(wildcards):
    """Get the annotation path for a species from config."""
    for species in get_species_list():
        if species["name"] == wildcards.species:
            return species.get("annotation", "")
    return ""


def is_local_file(path_str: str) -> bool:
    """Check if a path is a local file (not an NCBI accession)."""
    if not path_str:
        return False
    # NCBI accessions start with GCF_ or GCA_
    if path_str.startswith("GCF_") or path_str.startswith("GCA_"):
        return False
    return True


def get_local_species():
    """Get list of species with local files (not NCBI accessions)."""
    local_species = []
    for species in get_species_list():
        assembly = species.get("assembly", "")
        annotation = species.get("annotation", "")
        # Species must have both local assembly and annotation
        if is_local_file(assembly) and is_local_file(annotation):
            local_species.append(species["name"])
    return local_species


# =============================================================================
# Standardization Rules
# =============================================================================

rule standardize_genome:
    """
    Standardize a genome FASTA file to the standard output path.

    Copies and compresses (if not already compressed) the genome file.
    Uses atomic write pattern for checkpoint safety.

    AC1: Standard directory structure
    AC4: File compression and atomic write
    """
    input:
        genome=lambda wildcards: get_species_genome(wildcards)
    output:
        genome="{output_dir}/standardized/{species}/genome.fa.gz"
    log:
        "{output_dir}/logs/standardize_genome/{species}.log"
    params:
        output_dir=get_output_dir()
    run:
        import logging
        from pathlib import Path

        # Set up file handler for this rule (avoids basicConfig conflicts)
        logger = logging.getLogger(f"standardize_genome.{wildcards.species}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            input_path = Path(input.genome)
            output_path = Path(output.genome)

            logger.info(f"Standardizing genome: {input_path} -> {output_path}")

            count = standardize_genome(input_path, output_path)

            logger.info(f"Genome standardization complete: {count} sequences")
        except Exception as e:
            logger.error(f"Genome standardization failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


rule standardize_annotation:
    """
    Standardize an annotation GFF3/GTF file to the standard output path.

    Copies and compresses the annotation file.
    Uses atomic write pattern for checkpoint safety.

    AC1: Standard directory structure
    AC4: File compression and atomic write
    """
    input:
        annotation=lambda wildcards: get_species_annotation(wildcards)
    output:
        annotation="{output_dir}/standardized/{species}/annotation.gff3.gz"
    log:
        "{output_dir}/logs/standardize_annotation/{species}.log"
    params:
        output_dir=get_output_dir()
    run:
        import logging
        from pathlib import Path

        # Set up file handler for this rule (avoids basicConfig conflicts)
        logger = logging.getLogger(f"standardize_annotation.{wildcards.species}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            input_path = Path(input.annotation)
            output_path = Path(output.annotation)

            logger.info(f"Standardizing annotation: {input_path} -> {output_path}")

            count = standardize_annotation(input_path, output_path)

            logger.info(f"Annotation standardization complete: {count} features")
        except Exception as e:
            logger.error(f"Annotation standardization failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


rule standardize_proteins:
    """
    Extract representative proteins from standardized annotation and genome.

    This rule:
    1. Parses the annotation and builds gene hierarchy
    2. Selects representative transcripts (longest CDS)
    3. Extracts CDS sequences from genome
    4. Translates to protein sequences
    5. Writes output files in OrthoFinder-compatible format

    AC2: Representative transcript selection
    AC3: Protein sequence extraction
    AC4: File compression and atomic write
    """
    input:
        genome="{output_dir}/standardized/{species}/genome.fa.gz",
        annotation="{output_dir}/standardized/{species}/annotation.gff3.gz"
    output:
        proteins="{output_dir}/standardized/{species}/proteins.longest.fa.gz",
        mapping="{output_dir}/standardized/{species}/representative_transcripts.tsv"
    log:
        "{output_dir}/logs/standardize_proteins/{species}.log"
    params:
        output_dir=get_output_dir(),
        min_protein_length=10
    run:
        import logging
        from pathlib import Path

        # Set up file handler for this rule (avoids basicConfig conflicts)
        logger = logging.getLogger(f"standardize_proteins.{wildcards.species}")
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log[0])
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)

        try:
            annotation_path = Path(input.annotation)
            genome_path = Path(input.genome)
            proteins_path = Path(output.proteins)
            mapping_path = Path(output.mapping)

            logger.info(f"Extracting proteins for species: {wildcards.species}")
            logger.info(f"  Genome: {genome_path}")
            logger.info(f"  Annotation: {annotation_path}")

            genes_count, proteins_count = extract_proteins(
                species_id=wildcards.species,
                annotation_path=annotation_path,
                genome_path=genome_path,
                proteins_output=proteins_path,
                mapping_output=mapping_path,
                min_protein_length=params.min_protein_length
            )

            logger.info(f"Protein extraction complete:")
            logger.info(f"  Genes: {genes_count}")
            logger.info(f"  Proteins: {proteins_count}")
        except Exception as e:
            logger.error(f"Protein extraction failed: {e}")
            raise
        finally:
            handler.close()
            logger.removeHandler(handler)


# =============================================================================
# Aggregate Rules
# =============================================================================

rule standardize_all_local:
    """
    Standardize all species with local files.

    This is a checkpoint rule that triggers standardization for all
    species that have local assembly and annotation files.
    """
    input:
        proteins=expand(
            "{output_dir}/standardized/{species}/proteins.longest.fa.gz",
            output_dir=get_output_dir(),
            species=get_local_species()
        ),
        mapping=expand(
            "{output_dir}/standardized/{species}/representative_transcripts.tsv",
            output_dir=get_output_dir(),
            species=get_local_species()
        )
    output:
        touch("{output_dir}/standardized/.standardize_local_complete")
    params:
        output_dir=get_output_dir()
