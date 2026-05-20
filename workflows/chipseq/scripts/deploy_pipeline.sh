#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  deploy_pipeline.sh <project_dir> [--force]

Create a ChIP-seq/CUT&Tag Snakemake project directory.

The script copies editable example files:
    config/config.example.yml     -> <project_dir>/config/config.yml
    config/samples.example.yml    -> <project_dir>/config/samples.yml
    config/run_snakemake.env.example -> <project_dir>/config/run_snakemake.env

The script links reusable pipeline files:
    Snakefile, rules/, envs/, scripts/, run_snakemake.sh

Options:
    --force    overwrite existing copied config files and links
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_DIR=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
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
            if [[ -n "$PROJECT_DIR" ]]; then
                echo "Only one project directory can be specified." >&2
                usage >&2
                exit 2
            fi
            PROJECT_DIR="$1"
            shift
            ;;
    esac
done

if [[ -z "$PROJECT_DIR" ]]; then
    usage >&2
    exit 2
fi

mkdir -p "$PROJECT_DIR/config" "$PROJECT_DIR/raw_data"

copy_example() {
    local source_name="$1"
    local target_name="$2"
    local source_path="$WORKFLOW_DIR/$source_name"
    local target_path="$PROJECT_DIR/$target_name"

    if [[ ! -f "$source_path" ]]; then
        echo "Missing template file: $source_path" >&2
        exit 1
    fi

    if [[ -e "$target_path" && "$FORCE" != "1" ]]; then
        echo "Keep existing $target_path"
        return
    fi

    cp "$source_path" "$target_path"
    echo "Copied $target_name"
}

link_item() {
    local name="$1"
    local source_path="$WORKFLOW_DIR/$name"
    local target_path="$PROJECT_DIR/$name"

    if [[ ! -e "$source_path" ]]; then
        echo "Missing pipeline item: $source_path" >&2
        exit 1
    fi

    if [[ -e "$target_path" || -L "$target_path" ]]; then
        if [[ "$FORCE" != "1" ]]; then
            echo "Keep existing $target_path"
            return
        fi
        rm -rf "$target_path"
    fi

    ln -s "$source_path" "$target_path"
    echo "Linked $name"
}

copy_example "config/config.example.yml" "config/config.yml"
copy_example "config/samples.example.yml" "config/samples.yml"
copy_example "config/run_snakemake.env.example" "config/run_snakemake.env"

link_item "Snakefile"
link_item "rules"
link_item "envs"
link_item "scripts"
link_item "run_snakemake.sh"

cat <<EOF

Deployment complete:
  $PROJECT_DIR

Next steps:
  1. Put FASTQ files in $PROJECT_DIR/raw_data
  2. Edit $PROJECT_DIR/config/config.yml and $PROJECT_DIR/config/samples.yml
  3. Optional: edit $PROJECT_DIR/config/run_snakemake.env for local runtime settings
  4. Run: cd "$PROJECT_DIR" && ./run_snakemake.sh -n --cores 1
EOF
