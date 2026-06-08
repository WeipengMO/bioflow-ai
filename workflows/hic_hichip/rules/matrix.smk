rule hicpro_matrix_to_cool:
    input:
        matrix=PATHS.iced_matrix("{sample}", "{resolution}"),
        bins=PATHS.iced_bins("{sample}", "{resolution}")
    output:
        cool=PATHS.cool_file("{sample}", "{resolution}")
    wildcard_constraints:
        sample=CTX.sample_pattern,
        resolution="|".join(CTX.resolutions)
    params:
        chrom_sizes=lambda wildcards: (config.get("genome", {}) or {}).get("chrom_sizes", ""),
        balance=lambda wildcards: str((config.get("matrix", {}) or {}).get("balance", True)).lower()
    log:
        PATHS.log("matrix/hicpro_to_cool/{sample}.{resolution}")
    benchmark:
        PATHS.benchmark("matrix/hicpro_to_cool/{sample}.{resolution}")
    threads:
        workflow_threads("matrix_convert", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.cool:q}) $(dirname {log:q}) $(dirname {benchmark:q})
python scripts/convert_hicpro_matrix.py \
    --matrix {input.matrix:q} \
    --bins {input.bins:q} \
    --chrom-sizes {params.chrom_sizes:q} \
    --resolution {wildcards.resolution:q} \
    --output {output.cool:q} \
    --balance {params.balance:q} \
    &> {log:q}
test -s {output.cool:q}
        """

rule make_mcool:
    input:
        cool=lambda wildcards: PATHS.cool_file(wildcards.sample, CTX.primary_resolution)
    output:
        mcool=PATHS.mcool_file("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    params:
        resolutions=lambda wildcards: ",".join(CTX.resolutions)
    log:
        PATHS.log("matrix/mcool/{sample}")
    benchmark:
        PATHS.benchmark("matrix/mcool/{sample}")
    threads:
        workflow_threads("matrix_convert", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.mcool:q}) $(dirname {log:q}) $(dirname {benchmark:q})
rm -f {output.mcool:q}
cooler zoomify \
    -n {threads} \
    -r {params.resolutions:q} \
    -o {output.mcool:q} \
    {input.cool:q} \
    &> {log:q}
test -s {output.mcool:q}
        """
