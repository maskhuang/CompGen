"""
CompGene Data Standardization Module.

This module provides utilities for standardizing local annotation files
into a unified directory structure for downstream analysis.

Key features:
- Standardize genome and annotation files to consistent paths
- Extract representative transcript proteins
- CDS to protein translation with codon tables
- Support for gzip compression and atomic writes

Source: Story 2.2 - 本地数据标准化
"""

import gzip
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.fasta import (
    FastaRecord,
    parse_fasta,
    write_fasta_gzip,
)
from workflow.lib.gff import (
    GeneModel,
    TranscriptModel,
    build_gene_hierarchy,
    parse_annotation,
    validate_gff,
)
from workflow.lib.io import atomic_write


# =============================================================================
# Constants
# =============================================================================

# Standard genetic code (NCBI Table 1)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Complement bases for reverse complement
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class SpeciesData:
    """
    Configuration for a single species' input data.

    Attributes:
        species_id: Short identifier for the species (e.g., "mmur").
        genome_path: Path to the genome FASTA file.
        annotation_path: Path to the annotation GFF3/GTF file.
        output_dir: Directory for standardized output files.
        name: Optional full species name.
    """

    species_id: str
    genome_path: Path
    annotation_path: Path
    output_dir: Path
    name: Optional[str] = None

    def __post_init__(self):
        """Convert string paths to Path objects."""
        if isinstance(self.genome_path, str):
            self.genome_path = Path(self.genome_path)
        if isinstance(self.annotation_path, str):
            self.annotation_path = Path(self.annotation_path)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)


@dataclass
class RepresentativeTranscript:
    """
    Information about a representative transcript.

    Attributes:
        gene_id: Gene identifier.
        transcript_id: Selected representative transcript ID.
        cds_length: Total CDS length in base pairs.
        exon_length: Total exon length in base pairs.
        is_representative: Always True for this record.
    """

    gene_id: str
    transcript_id: str
    cds_length: int
    exon_length: int
    is_representative: bool = True


# =============================================================================
# Input Validation
# =============================================================================

def validate_species_inputs(species: SpeciesData) -> None:
    """
    Validate that species input files exist and have valid format.

    Args:
        species: SpeciesData configuration.

    Raises:
        CompGeneError: If files are missing or have invalid format.
    """
    # Check genome file exists
    if not species.genome_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Genome file not found: {species.genome_path}",
            details=f"Check config.yaml species.{species.species_id}.genome"
        )

    # Check annotation file exists
    if not species.annotation_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Annotation file not found: {species.annotation_path}",
            details=f"Check config.yaml species.{species.species_id}.annotation"
        )

    # Validate annotation format
    is_valid, error = validate_gff(species.annotation_path)
    if not is_valid:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Invalid annotation format: {species.annotation_path}",
            details=error
        )


def validate_has_cds_features(annotation_path: Path) -> None:
    """
    Validate that annotation contains CDS features for protein extraction.

    Args:
        annotation_path: Path to annotation file.

    Raises:
        CompGeneError: If no CDS features found.
    """
    has_cds = False
    for feature in parse_annotation(annotation_path):
        if feature.type == "CDS":
            has_cds = True
            break

    if not has_cds:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"No CDS features found in annotation: {annotation_path}",
            details="Annotation must contain CDS features for protein extraction"
        )


# =============================================================================
# File Standardization
# =============================================================================

def _open_input_file(path: Path, mode: str = "rt"):
    """Open a file, auto-detecting gzip compression."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def standardize_genome(
    input_path: Path,
    output_path: Path
) -> int:
    """
    Standardize a genome FASTA file to the standard output path.

    Copies and compresses (if not already compressed) the genome file.
    Uses atomic write pattern for checkpoint safety.

    Args:
        input_path: Path to input genome FASTA.
        output_path: Path to output genome.fa.gz.

    Returns:
        Number of sequences processed.

    Raises:
        CompGeneError: If input file not found.
    """
    if not input_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Genome file not found: {input_path}",
        )

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Count sequences and copy with compression
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    count = 0

    try:
        with _open_input_file(input_path) as fin:
            with gzip.open(temp_path, "wt", encoding="utf-8") as fout:
                for line in fin:
                    if line.startswith(">"):
                        count += 1
                    fout.write(line)

        # Atomic rename
        os.rename(temp_path, output_path)
    except Exception:
        if temp_path.exists():
            temp_path.unlink()
        raise

    return count


def standardize_annotation(
    input_path: Path,
    output_path: Path
) -> int:
    """
    Standardize an annotation GFF3/GTF file to the standard output path.

    Copies and compresses the annotation file, converting to GFF3 format.
    Uses atomic write pattern for checkpoint safety.

    Args:
        input_path: Path to input annotation file.
        output_path: Path to output annotation.gff3.gz.

    Returns:
        Number of feature lines processed.

    Raises:
        CompGeneError: If input file not found.
    """
    if not input_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Annotation file not found: {input_path}",
        )

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Copy with compression
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    count = 0

    try:
        with _open_input_file(input_path) as fin:
            with gzip.open(temp_path, "wt", encoding="utf-8") as fout:
                for line in fin:
                    fout.write(line)
                    if not line.startswith("#") and line.strip():
                        count += 1

        # Atomic rename
        os.rename(temp_path, output_path)
    except Exception:
        if temp_path.exists():
            temp_path.unlink()
        raise

    return count


# =============================================================================
# DNA/Protein Translation
# =============================================================================

def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq.upper()))


def translate_cds(cds_sequence: str, include_stop: bool = False) -> str:
    """
    Translate a CDS nucleotide sequence to amino acids.

    Uses the standard genetic code (NCBI Table 1).

    Args:
        cds_sequence: CDS nucleotide sequence (already oriented 5' to 3').
        include_stop: Whether to include stop codon (*) in output.

    Returns:
        Amino acid sequence.
    """
    seq = cds_sequence.upper()
    protein = []

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon, 'X')  # X for unknown codons
        if aa == '*' and not include_stop:
            break
        protein.append(aa)

    return ''.join(protein)


def extract_cds_sequence(
    transcript: TranscriptModel,
    genome_sequences: dict[str, str]
) -> Optional[str]:
    """
    Extract and concatenate CDS sequence for a transcript.

    Handles both positive and negative strand transcripts.

    Args:
        transcript: TranscriptModel with CDS features.
        genome_sequences: Dictionary mapping seqid to sequence.

    Returns:
        Concatenated CDS sequence (5' to 3'), or None if no CDS.
    """
    if not transcript.cds_features:
        return None

    # Get strand from first CDS
    strand = transcript.cds_features[0].strand

    # Sort CDS by position
    sorted_cds = sorted(transcript.cds_features, key=lambda c: c.start)

    # Extract sequences
    cds_parts = []
    for cds in sorted_cds:
        seqid = cds.seqid
        if seqid not in genome_sequences:
            continue

        genome_seq = genome_sequences[seqid]
        # GFF coordinates are 1-based, Python is 0-based
        start = cds.start - 1
        end = cds.end
        cds_seq = genome_seq[start:end]
        cds_parts.append(cds_seq)

    if not cds_parts:
        return None

    # Concatenate
    full_cds = ''.join(cds_parts)

    # Reverse complement for negative strand
    if strand == '-':
        full_cds = reverse_complement(full_cds)

    return full_cds


# =============================================================================
# Representative Transcript Selection
# =============================================================================

def select_representative_transcripts(
    genes: dict[str, GeneModel]
) -> list[RepresentativeTranscript]:
    """
    Select representative transcripts for each gene.

    Uses the GeneModel.get_representative() method which selects based on:
    1. Longest CDS total length
    2. Longest exon total length (if no CDS)
    3. Lexicographically smallest ID (tiebreaker)

    Args:
        genes: Dictionary of gene_id to GeneModel.

    Returns:
        List of RepresentativeTranscript records.
    """
    representatives = []

    for gene_id, gene in genes.items():
        rep_transcript = gene.get_representative()
        if rep_transcript:
            representatives.append(RepresentativeTranscript(
                gene_id=gene_id,
                transcript_id=rep_transcript.transcript_id,
                cds_length=rep_transcript.cds_length,
                exon_length=rep_transcript.exon_length,
                is_representative=True
            ))

    return representatives


def write_representative_transcripts_tsv(
    representatives: list[RepresentativeTranscript],
    output_path: Path
) -> int:
    """
    Write representative transcripts mapping to TSV file.

    Args:
        representatives: List of RepresentativeTranscript records.
        output_path: Path to output TSV file.

    Returns:
        Number of records written.
    """
    lines = ["gene_id\ttranscript_id\tcds_length\texon_length\tis_representative"]
    for rep in representatives:
        lines.append(
            f"{rep.gene_id}\t{rep.transcript_id}\t{rep.cds_length}\t"
            f"{rep.exon_length}\t{1 if rep.is_representative else 0}"
        )

    content = "\n".join(lines) + "\n"
    atomic_write(output_path, content)
    return len(representatives)


# =============================================================================
# Protein Extraction
# =============================================================================

def extract_proteins(
    species_id: str,
    annotation_path: Path,
    genome_path: Path,
    proteins_output: Path,
    mapping_output: Path,
    min_protein_length: int = 10
) -> tuple[int, int]:
    """
    Extract representative protein sequences from annotation and genome.

    This is the main function that:
    1. Parses annotation and builds gene hierarchy
    2. Selects representative transcripts
    3. Extracts CDS sequences from genome
    4. Translates to protein
    5. Writes output files

    Args:
        species_id: Species identifier for protein IDs.
        annotation_path: Path to annotation GFF3/GTF file.
        genome_path: Path to genome FASTA file.
        proteins_output: Path to output proteins FASTA (.gz).
        mapping_output: Path to output representative transcripts TSV.
        min_protein_length: Minimum protein length to include.

    Returns:
        Tuple of (genes_count, proteins_count).

    Raises:
        CompGeneError: If input files missing or invalid.
    """
    # Validate inputs
    if not annotation_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Annotation file not found: {annotation_path}",
        )
    if not genome_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Genome file not found: {genome_path}",
        )

    # Ensure output directories exist
    proteins_output.parent.mkdir(parents=True, exist_ok=True)
    mapping_output.parent.mkdir(parents=True, exist_ok=True)

    # Load genome sequences into memory
    genome_sequences = {}
    for record in parse_fasta(genome_path):
        genome_sequences[record.id] = record.sequence

    # Parse annotation and build gene hierarchy
    genes = build_gene_hierarchy(parse_annotation(annotation_path))

    # Select representative transcripts
    representatives = select_representative_transcripts(genes)

    # Extract proteins
    protein_records = []
    for gene_id, gene in genes.items():
        rep_transcript = gene.get_representative()
        if not rep_transcript:
            continue

        # Extract CDS sequence
        cds_seq = extract_cds_sequence(rep_transcript, genome_sequences)
        if not cds_seq:
            continue

        # Translate to protein
        protein_seq = translate_cds(cds_seq, include_stop=False)

        # Skip short proteins
        if len(protein_seq) < min_protein_length:
            continue

        # Create protein record with OrthoFinder-compatible ID
        protein_id = f"{species_id}|{gene_id}"
        protein_records.append(FastaRecord(
            id=protein_id,
            description=f"{protein_id} {rep_transcript.transcript_id}",
            sequence=protein_seq
        ))

    # Write outputs
    write_fasta_gzip(protein_records, proteins_output)
    write_representative_transcripts_tsv(representatives, mapping_output)

    return len(genes), len(protein_records)


# =============================================================================
# High-Level Standardization
# =============================================================================

def standardize_species(species: SpeciesData) -> dict:
    """
    Fully standardize a species' data to the output directory.

    Creates:
    - {output_dir}/genome.fa.gz
    - {output_dir}/annotation.gff3.gz
    - {output_dir}/proteins.longest.fa.gz
    - {output_dir}/representative_transcripts.tsv

    Args:
        species: SpeciesData configuration.

    Returns:
        Dictionary with statistics:
        - genome_sequences: Number of genome sequences
        - annotation_features: Number of annotation features
        - genes: Number of genes
        - proteins: Number of proteins extracted

    Raises:
        CompGeneError: If validation or processing fails.
    """
    # Validate inputs
    validate_species_inputs(species)

    # Define output paths
    genome_output = species.output_dir / "genome.fa.gz"
    annotation_output = species.output_dir / "annotation.gff3.gz"
    proteins_output = species.output_dir / "proteins.longest.fa.gz"
    mapping_output = species.output_dir / "representative_transcripts.tsv"

    # Standardize genome
    genome_count = standardize_genome(
        species.genome_path,
        genome_output
    )

    # Standardize annotation
    annotation_count = standardize_annotation(
        species.annotation_path,
        annotation_output
    )

    # Extract proteins (using standardized files)
    genes_count, proteins_count = extract_proteins(
        species_id=species.species_id,
        annotation_path=annotation_output,
        genome_path=genome_output,
        proteins_output=proteins_output,
        mapping_output=mapping_output
    )

    return {
        "genome_sequences": genome_count,
        "annotation_features": annotation_count,
        "genes": genes_count,
        "proteins": proteins_count,
    }
