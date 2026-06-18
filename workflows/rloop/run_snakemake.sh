#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEMAKE_RUNTIME_CONFIG="${SNAKEMAKE_RUNTIME_CONFIG:-$SCRIPT_DIR/config/run_snakemake.env}"

if [[ -f "$SNAKEMAKE_RUNTIME_CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$SNAKEMAKE_RUNTIME_CONFIG"
fi

SNAKEFILE="${SNAKEFILE:-$SCRIPT_DIR/Snakefile}"
: "${SNAKEMAKE_CORES:?Set SNAKEMAKE_CORES in config/run_snakemake.env or the environment.}"
: "${SNAKEMAKE_CONDA_PREFIX:?Set SNAKEMAKE_CONDA_PREFIX in config/run_snakemake.env or the environment.}"
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
    "$@" \
    --rerun-triggers mtime
