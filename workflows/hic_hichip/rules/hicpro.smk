rule make_hicpro_config:
    input:
        template=lambda wildcards: CTX.hicpro_template,
        samples=lambda wildcards: config.get("samples", "config/samples.tsv")
    output:
        PATHS.hicpro_generated_config
    params:
        config_yml=lambda wildcards: config.get("config_yml", "config/config.yml")
    log:
        PATHS.log("hicpro/make_hicpro_config")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
python scripts/make_hicpro_config.py \
    --config {params.config_yml:q} \
    --samples {input.samples:q} \
    --template {input.template:q} \
    --output {output:q} \
    &> {log:q}
test -s {output:q}
        """

rule run_hicpro:
    input:
        config=PATHS.hicpro_generated_config,
        sample_inputs=expand(PATHS.hicpro_input_done("{sample}"), sample=CTX.sample_names)
    output:
        done=PATHS.hicpro_done()
    params:
        input_dir=PATHS.hicpro_input,
        output_dir=PATHS.hicpro,
        container=lambda wildcards: CTX.hicpro_container,
        bind_paths=lambda wildcards: ",".join(CTX.hicpro_bind_paths),
        singularity_args=lambda wildcards: CTX.hicpro_singularity_args,
        extra=lambda wildcards: (config.get("hicpro", {}) or {}).get("extra", ""),
        benchmark_dir=PATHS.benchmarks,
        hic_results_data=f"{PATHS.hicpro}/hic_results/data"
    log:
        PATHS.log("hicpro/run_hicpro")
    benchmark:
        PATHS.benchmark("hicpro/run_hicpro")
    threads:
        workflow_threads("hicpro", 32)
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {params.output_dir:q}) $(dirname {log:q}) {params.benchmark_dir:q}
if [[ ! -s {params.container:q} ]]; then
    echo "ERROR: HiC-Pro Singularity image not found: {params.container}" >&2
    exit 1
fi
rmdir {params.output_dir:q} 2>/dev/null || true
singularity exec --cleanenv \
    -B {params.bind_paths:q} \
    {params.singularity_args} \
    {params.container:q} \
    HiC-Pro \
        -i {params.input_dir:q} \
        -o {params.output_dir:q} \
        -c {input.config:q} \
        {params.extra} \
        &> {log:q}
test -d {params.hic_results_data:q}
date > {output.done:q}
test -s {output.done:q}
        """

rule link_valid_pairs:
    input:
        done=PATHS.hicpro_done()
    output:
        valid_pairs=PATHS.sample_valid_pairs("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    log:
        PATHS.log("hicpro/link_valid_pairs/{sample}")
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.valid_pairs:q}) $(dirname {log:q})
base={PATHS.hicpro:q}
sample={wildcards.sample:q}
mapfile -t candidates < <(find "$base" -path "*/hic_results/data/${{sample}}/*allValidPairs*" -type f | sort)
if [[ ${{#candidates[@]}} -eq 0 ]]; then
    echo "ERROR: Could not find HiC-Pro allValidPairs for sample $sample under $base" | tee {log:q}
    exit 1
fi
if [[ ${{#candidates[@]}} -gt 1 ]]; then
    echo "ERROR: multiple HiC-Pro allValidPairs candidates found for $sample:" | tee {log:q}
    printf '%s\n' "${{candidates[@]}}" | tee -a {log:q}
    exit 1
fi
ln -sf $(realpath "${{candidates[0]}}") {output.valid_pairs:q}
test -s {output.valid_pairs:q}
        """

rule link_hicpro_matrices:
    input:
        done=PATHS.hicpro_done()
    output:
        raw_matrix=PATHS.raw_matrix("{sample}", "{resolution}"),
        raw_bins=PATHS.raw_bins("{sample}", "{resolution}"),
        iced_matrix=PATHS.iced_matrix("{sample}", "{resolution}"),
        iced_bins=PATHS.iced_bins("{sample}", "{resolution}")
    wildcard_constraints:
        sample=CTX.sample_pattern,
        resolution="|".join(CTX.resolutions)
    log:
        PATHS.log("hicpro/link_matrices/{sample}.{resolution}")
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.raw_matrix:q}) $(dirname {output.iced_matrix:q}) $(dirname {log:q})
base={PATHS.hicpro:q}
sample={wildcards.sample:q}
res={wildcards.resolution:q}
mapfile -t raw_dirs < <(find "$base" -path "*/hic_results/matrix/${{sample}}/raw/${{res}}" -type d | sort)
mapfile -t iced_dirs < <(find "$base" -path "*/hic_results/matrix/${{sample}}/iced/${{res}}" -type d | sort)
if [[ ${{#raw_dirs[@]}} -eq 0 ]]; then
    echo "ERROR: raw matrix directory not found for $sample resolution $res" | tee {log:q}
    exit 1
fi
if [[ ${{#raw_dirs[@]}} -gt 1 ]]; then
    echo "ERROR: multiple raw matrix directories found for $sample resolution $res:" | tee {log:q}
    printf '%s\n' "${{raw_dirs[@]}}" | tee -a {log:q}
    exit 1
fi
if [[ ${{#iced_dirs[@]}} -eq 0 ]]; then
    echo "ERROR: iced matrix directory not found for $sample resolution $res" | tee {log:q}
    exit 1
fi
if [[ ${{#iced_dirs[@]}} -gt 1 ]]; then
    echo "ERROR: multiple iced matrix directories found for $sample resolution $res:" | tee {log:q}
    printf '%s\n' "${{iced_dirs[@]}}" | tee -a {log:q}
    exit 1
fi
raw_dir="${{raw_dirs[0]}}"
iced_dir="${{iced_dirs[0]}}"
mapfile -t raw_matrices < <(find "$raw_dir" -name "*.matrix" -type f | sort)
mapfile -t raw_bins_files < <(find "$raw_dir" -name "*abs.bed" -type f | sort)
mapfile -t iced_matrices < <(find "$iced_dir" -name "*.matrix" -type f | sort)
mapfile -t iced_bins_files < <(find "$iced_dir" -name "*abs.bed" -type f | sort)
if [[ ${{#raw_matrices[@]}} -ne 1 || ${{#raw_bins_files[@]}} -ne 1 || ${{#iced_matrices[@]}} -ne 1 ]]; then
    echo "ERROR: expected one raw matrix, one raw abs bed, and one iced matrix for $sample resolution $res" | tee {log:q}
    printf 'raw_matrices:\n%s\nraw_bins:\n%s\niced_matrices:\n%s\n' "${{raw_matrices[*]}}" "${{raw_bins_files[*]}}" "${{iced_matrices[*]}}" | tee -a {log:q}
    exit 1
fi
raw_matrix="${{raw_matrices[0]}}"
raw_bins="${{raw_bins_files[0]}}"
iced_matrix="${{iced_matrices[0]}}"
if [[ ${{#iced_bins_files[@]}} -eq 0 ]]; then
    iced_bins="$raw_bins"
elif [[ ${{#iced_bins_files[@]}} -eq 1 ]]; then
    iced_bins="${{iced_bins_files[0]}}"
else
    echo "ERROR: multiple iced abs bed candidates found for $sample resolution $res:" | tee {log:q}
    printf '%s\n' "${{iced_bins_files[@]}}" | tee -a {log:q}
    exit 1
fi
if [[ -z "$raw_matrix" || -z "$raw_bins" || -z "$iced_matrix" || -z "$iced_bins" ]]; then
    echo "ERROR: matrix files are incomplete for $sample resolution $res" | tee {log:q}
    exit 1
fi
ln -sf $(realpath "$raw_matrix") {output.raw_matrix:q}
ln -sf $(realpath "$raw_bins") {output.raw_bins:q}
ln -sf $(realpath "$iced_matrix") {output.iced_matrix:q}
ln -sf $(realpath "$iced_bins") {output.iced_bins:q}
test -s {output.raw_matrix:q}
test -s {output.raw_bins:q}
test -s {output.iced_matrix:q}
test -s {output.iced_bins:q}
        """
