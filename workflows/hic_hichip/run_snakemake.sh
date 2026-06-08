#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

MODE="${1:-dry-run}"
CORES="${2:-8}"
SNAKEFILE="Snakefile"

case "$MODE" in
  dry-run|dryrun|np)
    snakemake -s "$SNAKEFILE" -np --cores "$CORES"
    ;;
  run)
    snakemake -s "$SNAKEFILE" --use-conda --cores "$CORES" --rerun-incomplete --printshellcmds
    ;;
  run-singularity)
    # HiC-Pro itself is run through the configured Singularity image.
    # This flag is only for rules that declare Snakemake containers in future extensions.
    snakemake -s "$SNAKEFILE" --use-conda --use-singularity --cores "$CORES" --rerun-incomplete --printshellcmds
    ;;
  unlock)
    snakemake -s "$SNAKEFILE" --unlock
    ;;
  clean-metadata)
    snakemake -s "$SNAKEFILE" --cleanup-metadata results || true
    ;;
  *)
    echo "Usage: $0 {dry-run|run|run-singularity|unlock|clean-metadata} [cores]" >&2
    exit 2
    ;;
esac
