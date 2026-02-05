#!/bin/bash
#$ -N liftoff_lct
#$ -q eichler-short.q
#$ -l h_rt=4:00:00
#$ -l h_vmem=4G
#$ -pe serial 16
#$ -cwd
#$ -o logs/qsub_liftoff_lct.stdout
#$ -e logs/qsub_liftoff_lct.stderr
#$ -V

set -euo pipefail

module load minimap2/2.26

WORKDIR=/net/eichler/vol28/projects/nhp_rearrangements/nobackups/lemur_annotations/compgene
cd "$WORKDIR"

REF_GFF=results/liftoff/ncbi_lct_to_lab_lct/decompressed/annotation.gff3
REF_FA=results/liftoff/ncbi_lct_to_lab_lct/decompressed/genome.fa
TARGET_FA=/net/eichler/vol28/projects/nhp_rearrangements/nobackups/nhp_assemble/Onyx_LCT/parental_verkko/fasta_rag/new_Onyx_LCT_hap1.fasta
OUTDIR=results/liftoff/ncbi_lct_to_lab_lct

echo "=== Liftoff LCT (Lemur catta) ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Threads: $NSLOTS"
free -g

mkdir -p "$OUTDIR/intermediate" logs

# Decompress NCBI files if not already done
if [ ! -f "$REF_GFF" ]; then
  mkdir -p "$OUTDIR/decompressed"
  echo "Decompressing annotation..."
  zcat /net/eichler/vol28/projects/nhp_rearrangements/nobackups/lemur_annotations/ncbi_downloads/GCF_020740605.2/annotation.gff.gz > "$REF_GFF"
fi
if [ ! -f "$REF_FA" ]; then
  echo "Decompressing genome..."
  zcat /net/eichler/vol28/projects/nhp_rearrangements/nobackups/lemur_annotations/ncbi_downloads/GCF_020740605.2/genome.fna.gz > "$REF_FA"
fi

# Build chroms mapping if not exists
CHROMS=$OUTDIR/chroms.txt
if [ ! -f "$CHROMS" ]; then
  echo "Building chromosome mapping..."
  python3 -c "
import re
ncbi = {}
with open('$REF_FA') as f:
    for line in f:
        if line.startswith('>'):
            acc = line.split()[0][1:]
            m = re.search(r'chromosome (\w+)', line)
            if m: ncbi[acc] = m.group(1)
import subprocess
r = subprocess.run(['grep', '^>', '$TARGET_FA'], capture_output=True, text=True)
lab = {}
for line in r.stdout.strip().split('\n'):
    name = line[1:].split()[0]
    # LCT uses NC_*_RagTag_hap1 naming
    m = re.match(r'(NC_\d+\.\d+)_RagTag', name)
    if m: lab[m.group(1)] = name
with open('$CHROMS', 'w') as f:
    for acc in ncbi:
        if acc in lab:
            f.write(f'{acc},{lab[acc]}\n')
print(f'Mapped {sum(1 for a in ncbi if a in lab)} chromosomes')
"
fi

# Remove old output
rm -f "$OUTDIR/lifted_annotation.gff3" "$OUTDIR/unmapped_features.txt"

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

if [ -f "$OUTDIR/lifted_annotation.gff3" ]; then
  GENES=$(grep -c $'\tgene\t' "$OUTDIR/lifted_annotation.gff3" || true)
  UNMAPPED=$(wc -l < "$OUTDIR/unmapped_features.txt" || echo 0)
  echo "Lifted genes: $GENES"
  echo "Unmapped features: $UNMAPPED"
fi
