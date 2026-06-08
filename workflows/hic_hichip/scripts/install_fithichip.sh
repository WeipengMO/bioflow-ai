#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  install_fithichip.sh [--prefix DIR] [--env-name NAME] [--skip-env] [--force]

Install or update FitHiChIP source code and prepare its conda dependencies.

Defaults:
  --prefix     <workflow_dir>/external/FitHiChIP
  --env-name   ${FITHICHIP_ENV_NAME:-bioflow-hic-hichip-fithichip}

Examples:
  scripts/install_fithichip.sh
  ./scripts/install_fithichip.sh --prefix tools/FitHiChIP
  ./scripts/install_fithichip.sh --prefix tools/FitHiChIP --env-name bioflow-hic-hichip-fithichip

The installer clones https://github.com/ay-lab/FitHiChIP.git and creates or
updates the conda environment described by envs/fithichip.yml.
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PREFIX="$WORKFLOW_DIR/external/FitHiChIP"
ENV_NAME="${FITHICHIP_ENV_NAME:-bioflow-hic-hichip-fithichip}"
SKIP_ENV=0
FORCE=0
REPO_URL="https://github.com/ay-lab/FitHiChIP.git"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --prefix)
            PREFIX="${2:?Missing value for --prefix}"
            shift 2
            ;;
        --env-name)
            ENV_NAME="${2:?Missing value for --env-name}"
            shift 2
            ;;
        --skip-env)
            SKIP_ENV=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            echo "Unexpected argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

find_conda_frontend() {
    if command -v mamba >/dev/null 2>&1; then
        echo mamba
    elif command -v conda >/dev/null 2>&1; then
        echo conda
    elif command -v micromamba >/dev/null 2>&1; then
        echo micromamba
    else
        echo ""
    fi
}

run_in_env() {
    local frontend="$1"
    shift
    if [[ "$frontend" == "micromamba" ]]; then
        micromamba run -n "$ENV_NAME" "$@"
    else
        "$frontend" run -n "$ENV_NAME" "$@"
    fi
}

create_or_update_env() {
    local frontend="$1"
    local env_file="$WORKFLOW_DIR/envs/fithichip.yml"
    if [[ ! -f "$env_file" ]]; then
        echo "Missing env file: $env_file" >&2
        exit 1
    fi
    if [[ "$frontend" == "micromamba" ]]; then
        micromamba env update -n "$ENV_NAME" -f "$env_file"
    else
        "$frontend" env update -n "$ENV_NAME" -f "$env_file"
    fi
}

PREFIX="$(mkdir -p "$(dirname "$PREFIX")" && cd "$(dirname "$PREFIX")" && pwd)/$(basename "$PREFIX")"

if [[ -e "$PREFIX" && ! -d "$PREFIX/.git" ]]; then
    if [[ "$FORCE" != "1" ]]; then
        echo "Prefix exists but is not a FitHiChIP git checkout: $PREFIX" >&2
        echo "Use --force to replace it." >&2
        exit 1
    fi
    rm -rf "$PREFIX"
fi

if [[ -d "$PREFIX/.git" ]]; then
    echo "Updating FitHiChIP in $PREFIX"
    git -C "$PREFIX" fetch --tags origin
    git -C "$PREFIX" pull --ff-only
else
    echo "Cloning FitHiChIP into $PREFIX"
    git clone "$REPO_URL" "$PREFIX"
fi

CONDA_FRONTEND="$(find_conda_frontend)"
if [[ "$SKIP_ENV" != "1" ]]; then
    if [[ -z "$CONDA_FRONTEND" ]]; then
        echo "No conda-compatible frontend found. Install mamba/conda/micromamba or rerun with --skip-env." >&2
        exit 1
    fi
    echo "Creating/updating conda env: $ENV_NAME"
    create_or_update_env "$CONDA_FRONTEND"
fi

SCRIPT="$PREFIX/FitHiChIP_HiCPro.sh"
if [[ ! -s "$SCRIPT" ]]; then
    echo "FitHiChIP runner not found after install: $SCRIPT" >&2
    exit 1
fi

if [[ -z "$CONDA_FRONTEND" || "$SKIP_ENV" == "1" ]]; then
    echo "Skipping dependency checks because conda env creation was skipped."
else
    echo "Checking FitHiChIP runtime dependencies in env: $ENV_NAME"
    run_in_env "$CONDA_FRONTEND" Rscript -e 'pkgs <- c("optparse","ggplot2","data.table","fdrtool","plyr","dplyr","GenomicRanges","edgeR"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop(paste("Missing R packages:", paste(missing, collapse=", ")))'
    run_in_env "$CONDA_FRONTEND" python -c 'import numpy, networkx, cooler, pandas, hicstraw'
    for exe in bedtools samtools bgzip tabix bowtie2 macs2; do
        run_in_env "$CONDA_FRONTEND" bash -lc "command -v $exe >/dev/null"
    done
fi

cat <<EOF

FitHiChIP installation ready:
  FITHICHIP_HOME=$PREFIX
  FITHICHIP_SCRIPT=$SCRIPT
  FITHICHIP_ENV_NAME=$ENV_NAME

For deployed projects, add these values to config/run_snakemake.env or set
loop_calling.fithichip_script in config/config.yml.
EOF
