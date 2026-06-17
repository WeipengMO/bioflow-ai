#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEMAKE_RUNTIME_CONFIG="${SNAKEMAKE_RUNTIME_CONFIG:-$SCRIPT_DIR/config/run_snakemake.env}"

if [[ -f "$SNAKEMAKE_RUNTIME_CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$SNAKEMAKE_RUNTIME_CONFIG"
fi

SNAKEFILE="${SNAKEFILE:-$SCRIPT_DIR/Snakefile}"
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-}"
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$SCRIPT_DIR/.snakemake/conda}"
SNAKEMAKE_CACHE_DIR="${SNAKEMAKE_CACHE_DIR:-$SCRIPT_DIR/.snakemake/cache}"
SNAKEMAKE_TMPDIR="${SNAKEMAKE_TMPDIR:-$SNAKEMAKE_CACHE_DIR/tmp}"

mkdir -p "$SNAKEMAKE_CONDA_PREFIX/pkgs" "$SNAKEMAKE_CACHE_DIR" "$SNAKEMAKE_TMPDIR"

export XDG_CACHE_HOME="$SNAKEMAKE_CACHE_DIR"
export TMPDIR="$SNAKEMAKE_TMPDIR"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$SNAKEMAKE_CONDA_PREFIX/pkgs}"

if [[ -n "${SNAKEMAKE_PROXY:-}" ]]; then
    export http_proxy="${http_proxy:-$SNAKEMAKE_PROXY}"
    export https_proxy="${https_proxy:-$SNAKEMAKE_PROXY}"
fi

if [[ -n "${FITHICHIP_HOME:-}" && -z "${FITHICHIP_SCRIPT:-}" ]]; then
    export FITHICHIP_SCRIPT="$FITHICHIP_HOME/FitHiChIP_HiCPro.sh"
fi

has_cores_arg=0
for arg in "$@"; do
    if [[ "$arg" == "--cores" || "$arg" == "-j" || "$arg" == --cores=* || "$arg" == -j* ]]; then
        has_cores_arg=1
        break
    fi
done

cores_args=()
if [[ "$has_cores_arg" == "0" ]]; then
    cores_args=(-j "${SNAKEMAKE_CORES:-1}")
fi

snakemake \
    -s "$SNAKEFILE" \
    --use-conda \
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
    "${cores_args[@]}" \
    --rerun-triggers mtime \
    "$@"
