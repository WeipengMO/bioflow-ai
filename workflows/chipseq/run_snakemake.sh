#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="${SNAKEFILE:-$SCRIPT_DIR/Snakefile}"
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-110}"
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-/data/user/mowp/snakemake_conda_envs}"
SNAKEMAKE_CACHE_DIR="${SNAKEMAKE_CACHE_DIR:-$SCRIPT_DIR/.snakemake/cache}"
SNAKEMAKE_TMPDIR="${SNAKEMAKE_TMPDIR:-$SNAKEMAKE_CACHE_DIR/tmp}"

mkdir -p "$SNAKEMAKE_CONDA_PREFIX/pkgs" "$SNAKEMAKE_CACHE_DIR" "$SNAKEMAKE_TMPDIR"

export XDG_CACHE_HOME="$SNAKEMAKE_CACHE_DIR"
export TMPDIR="$SNAKEMAKE_TMPDIR"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$SNAKEMAKE_CONDA_PREFIX/pkgs}"

# Optional proxy support. Define SNAKEMAKE_PROXY before running if needed.
if [[ -n "${SNAKEMAKE_PROXY:-}" ]]; then
    export http_proxy="${http_proxy:-$SNAKEMAKE_PROXY}"
    export https_proxy="${https_proxy:-$SNAKEMAKE_PROXY}"
fi

snakemake \
    -s "$SNAKEFILE" \
    --use-conda \
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
    -j "$SNAKEMAKE_CORES" \
    --rerun-triggers mtime \
    "$@"
