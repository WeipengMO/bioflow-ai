#!/usr/bin/env bash
# Download E. coli K-12 MG1655 reference FASTA from NCBI and build a Bowtie2 index.
# Usage: bash scripts/build_spikein_index.sh [--output-dir data/spikein] [--threads 8]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

OUTPUT_DIR="${1:-$WORKFLOW_DIR/data/spikein}"
THREADS="${2:-8}"

ACCESSION="GCF_000005845.2"
ASSEMBLY="ASM584v2"
FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/${ACCESSION}_${ASSEMBLY}/${ACCESSION}_${ASSEMBLY}_genomic.fna.gz"
FASTA_FILE="$OUTPUT_DIR/ecoli_mg1655.fna"
INDEX_PREFIX="$OUTPUT_DIR/ecoli_mg1655"

mkdir -p "$OUTPUT_DIR"

# Check required tools
for cmd in bowtie2-build curl gunzip; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found. Please install it first." >&2
        exit 1
    fi
done

# Download FASTA if not present
COMPRESSED="$OUTPUT_DIR/${ACCESSION}_${ASSEMBLY}_genomic.fna.gz"
if [[ -f "$FASTA_FILE" ]]; then
    echo "FASTA already exists: $FASTA_FILE"
elif [[ -f "$COMPRESSED" ]]; then
    echo "Decompressing $COMPRESSED ..."
    gunzip -c "$COMPRESSED" > "$FASTA_FILE"
else
    echo "Downloading E. coli K-12 MG1655 from NCBI ..."
    curl -fSL -o "$COMPRESSED" "$FASTA_URL"
    echo "Decompressing ..."
    gunzip -c "$COMPRESSED" > "$FASTA_FILE"
fi

# Build Bowtie2 index if not present
if ls "$INDEX_PREFIX".*.bt2 1>/dev/null 2>&1 || ls "$INDEX_PREFIX".*.bt2l 1>/dev/null 2>&1; then
    echo "Bowtie2 index already exists at $INDEX_PREFIX"
else
    echo "Building Bowtie2 index with $THREADS threads ..."
    bowtie2-build --threads "$THREADS" "$FASTA_FILE" "$INDEX_PREFIX"
    echo "Bowtie2 index built at $INDEX_PREFIX"
fi

echo ""
echo "Done. Add this to your config/config.yml:"
echo "  spikein_genome: $INDEX_PREFIX"
