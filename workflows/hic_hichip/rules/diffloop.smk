rule quantify_loop_counts:
    input:
        universe=PATHS.loop_universe(),
        valid_pairs=expand(PATHS.sample_valid_pairs("{sample}"), sample=CTX.sample_names),
        cools=expand(PATHS.cool_file("{sample}", CTX.primary_resolution), sample=CTX.sample_names)
    output:
        counts=PATHS.loop_counts(),
        loop_metadata=PATHS.loop_metadata(),
        sample_metadata=PATHS.sample_metadata()
    params:
        samples=",".join(CTX.sample_names),
        sample_table=lambda wildcards: config.get("samples", "config/samples.tsv"),
        source=lambda wildcards: CTX.quantification_source,
        resolution=lambda wildcards: CTX.primary_resolution
    log:
        PATHS.log("diffloop/quantify_loop_counts")
    benchmark:
        PATHS.benchmark("diffloop/quantify_loop_counts")
    threads:
        workflow_threads("diffloop", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.counts:q}) $(dirname {log:q}) $(dirname {benchmark:q})
if [[ {params.source:q} == "validpairs" ]]; then
    python scripts/quantify_loops_from_validpairs.py \
        --universe {input.universe:q} \
        --samples {params.samples:q} \
        --sample-table {params.sample_table:q} \
        --validpairs {input.valid_pairs:q} \
        --counts {output.counts:q} \
        --loop-metadata {output.loop_metadata:q} \
        --sample-metadata {output.sample_metadata:q} \
        --threads {threads} \
        &> {log:q}
else
    python scripts/quantify_loops_from_cool.py \
        --universe {input.universe:q} \
        --samples {params.samples:q} \
        --sample-table {params.sample_table:q} \
        --cools {input.cools:q} \
        --resolution {params.resolution:q} \
        --counts {output.counts:q} \
        --loop-metadata {output.loop_metadata:q} \
        --sample-metadata {output.sample_metadata:q} \
        --threads {threads} \
        &> {log:q}
fi
test -s {output.counts:q}
test -s {output.loop_metadata:q}
test -s {output.sample_metadata:q}
        """

rule differential_loops:
    input:
        counts=PATHS.loop_counts(),
        loop_metadata=PATHS.loop_metadata(),
        sample_metadata=PATHS.sample_metadata()
    output:
        table=PATHS.diffloops("{comparison}"),
        significant=PATHS.diffloops_significant("{comparison}"),
        up=PATHS.diffloops_up("{comparison}"),
        down=PATHS.diffloops_down("{comparison}"),
        volcano=PATHS.diffloop_volcano("{comparison}"),
        ma=PATHS.diffloop_ma("{comparison}"),
        heatmap=PATHS.diffloop_heatmap("{comparison}")
    wildcard_constraints:
        comparison=CTX.comparison_pattern
    params:
        group1=comparison_group1,
        group2=comparison_group2,
        fdr=lambda wildcards: CTX.diff_fdr,
        lfc_cutoff=lambda wildcards: CTX.diff_lfc_cutoff,
        min_count=lambda wildcards: (config.get("diffloop", {}) or {}).get("min_count", 10)
    log:
        PATHS.log("diffloop/{comparison}")
    benchmark:
        PATHS.benchmark("diffloop/{comparison}")
    conda:
        DIFF_ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.table:q}) $(dirname {log:q}) $(dirname {benchmark:q})
Rscript scripts/diffloop_deseq2.R \
    --counts {input.counts:q} \
    --loop-metadata {input.loop_metadata:q} \
    --sample-metadata {input.sample_metadata:q} \
    --comparison {wildcards.comparison:q} \
    --group1 {params.group1:q} \
    --group2 {params.group2:q} \
    --fdr {params.fdr:q} \
    --lfc-cutoff {params.lfc_cutoff:q} \
    --min-count {params.min_count:q} \
    --out-table {output.table:q} \
    --out-significant {output.significant:q} \
    --out-up {output.up:q} \
    --out-down {output.down:q} \
    --volcano {output.volcano:q} \
    --ma-plot {output.ma:q} \
    --heatmap {output.heatmap:q} \
    &> {log:q}
test -s {output.table:q}
test -s {output.significant:q}
test -s {output.up:q}
test -s {output.down:q}
test -s {output.volcano:q}
test -s {output.ma:q}
test -s {output.heatmap:q}
        """
