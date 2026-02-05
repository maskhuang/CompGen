#!/bin/bash
#$ -N liftoff_mmu
#$ -q eichler-short.q
#$ -l h_rt=4:00:00
#$ -l h_vmem=4G
#$ -pe serial 16
#$ -cwd
#$ -o logs/qsub_liftoff_mmu.stdout
#$ -e logs/qsub_liftoff_mmu.stderr
#$ -V

set -euo pipefail

module load minimap2/2.26

WORKDIR=/net/eichler/vol28/projects/nhp_rearrangements/nobackups/lemur_annotations/compgene
cd "$WORKDIR"

REF_GFF=results/liftoff/ncbi_mmu_to_lab_mmu/decompressed/annotation.gff3
REF_FA=results/liftoff/ncbi_mmu_to_lab_mmu/decompressed/genome.fa
TARGET_FA=/net/eichler/vol28/projects/nhp_rearrangements/nobackups/nhp_assemble/Inina_4.0/assembly_haphiC/fasta/Inina4_hap1.chr.fasta
OUTDIR=results/liftoff/ncbi_mmu_to_lab_mmu
CHROMS=$OUTDIR/chroms.txt

echo "=== Liftoff MMU (Microcebus murinus) ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Threads: $NSLOTS"
free -g

mkdir -p "$OUTDIR/intermediate" logs

~/.local/bin/liftoff \
  -g "$REF_GFF" \
  -o "$OUTDIR/lifted_annotation.gff3" \
  -u "$OUTDIR/unmapped_features.txt" \
  -dir "$OUTDIR/intermediate" \
  -chroms "$CHROMS" \
  -p "$NSLOTS" -a 0.5 -s 0.5 \
  "$TARGET_FA" "$REF_FA"

echo "Exit code: $?"
echo "End: $(date)"

# Quick stats
if [ -f "$OUTDIR/lifted_annotation.gff3" ]; then
  GENES=$(grep -c $'\tgene\t' "$OUTDIR/lifted_annotation.gff3" || true)
  UNMAPPED=$(wc -l < "$OUTDIR/unmapped_features.txt" || echo 0)
  echo "Lifted genes: $GENES"
  echo "Unmapped features: $UNMAPPED"
fi
