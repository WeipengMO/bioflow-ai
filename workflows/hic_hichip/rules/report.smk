rule annotate_loops:
    input:
        universe=PATHS.loop_universe()
    output:
        loops_to_genes=PATHS.loops_to_genes(),
        promoter_enhancer=PATHS.promoter_enhancer_loops()
    params:
        promoters=lambda wildcards: (config.get("annotation", {}) or {}).get("promoters", ""),
        enhancers=lambda wildcards: (config.get("annotation", {}) or {}).get("enhancers", ""),
        gtf=lambda wildcards: (config.get("annotation", {}) or {}).get("gtf", ""),
        promoter_upstream=lambda wildcards: (config.get("annotation", {}) or {}).get("promoter_upstream", 2000),
        promoter_downstream=lambda wildcards: (config.get("annotation", {}) or {}).get("promoter_downstream", 2000)
    log:
        PATHS.log("annotation/annotate_loops")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.loops_to_genes:q}) $(dirname {log:q})
python scripts/annotate_loops.py \
    --loops {input.universe:q} \
    --promoters {params.promoters:q} \
    --enhancers {params.enhancers:q} \
    --gtf {params.gtf:q} \
    --promoter-upstream {params.promoter_upstream:q} \
    --promoter-downstream {params.promoter_downstream:q} \
    --loops-to-genes {output.loops_to_genes:q} \
    --promoter-enhancer {output.promoter_enhancer:q} \
    &> {log:q}
test -s {output.loops_to_genes:q}
test -s {output.promoter_enhancer:q}
        """

rule workflow_report:
    input:
        qc=PATHS.hicpro_qc_summary(),
        peak_sources=PATHS.peak_source_table() if CTX.enable_loops else [],
        loops=expand(PATHS.sample_loops_bedpe("{sample}"), sample=CTX.sample_names) if CTX.enable_loops else []
    output:
        html=PATHS.workflow_report()
    params:
        outdir=PATHS.outdir,
        samples=",".join(CTX.sample_names),
        comparisons=",".join(CTX.comparison_names)
    log:
        PATHS.log("report/workflow_report")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.html:q}) $(dirname {log:q})
python scripts/make_report.py \
    --qc {input.qc:q} \
    --peak-sources {input.peak_sources:q} \
    --loops {input.loops:q} \
    --outdir {params.outdir:q} \
    --samples {params.samples:q} \
    --comparisons "{params.comparisons}" \
    --output {output.html:q} \
    &> {log:q}
test -s {output.html:q}
        """
