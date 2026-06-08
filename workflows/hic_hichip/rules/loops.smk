rule call_loops:
    input:
        valid_pairs=PATHS.sample_valid_pairs("{sample}"),
        cool=lambda wildcards: PATHS.cool_file(wildcards.sample, CTX.primary_resolution),
        mcool=PATHS.mcool_file("{sample}")
    output:
        bedpe=PATHS.sample_loops_bedpe("{sample}"),
        tsv=PATHS.sample_loops_tsv("{sample}"),
        anchors=PATHS.sample_anchors_bed("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    params:
        sample="{sample}",
        assay=lambda wildcards: CTX.samples[wildcards.sample]["assay"],
        target=lambda wildcards: CTX.samples[wildcards.sample].get("target", ""),
        peaks=sample_peak_bed,
        caller=lambda wildcards: CTX.loop_caller,
        resolution=lambda wildcards: CTX.primary_resolution,
        fdr=lambda wildcards: CTX.loop_fdr,
        config_yml="config/config.yml"
    log:
        PATHS.log("loops/call_loops/{sample}")
    benchmark:
        PATHS.benchmark("loops/call_loops/{sample}")
    threads:
        workflow_threads("loop_calling", 16)
    conda:
        LOOP_ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bedpe:q}) $(dirname {log:q}) $(dirname {benchmark:q})
python scripts/call_loops.py \
    --config {params.config_yml:q} \
    --sample {wildcards.sample:q} \
    --assay {params.assay:q} \
    --target {params.target:q} \
    --caller {params.caller:q} \
    --resolution {params.resolution:q} \
    --fdr {params.fdr:q} \
    --valid-pairs {input.valid_pairs:q} \
    --cool {input.cool:q} \
    --mcool {input.mcool:q} \
    --peaks {params.peaks:q} \
    --threads {threads} \
    --bedpe {output.bedpe:q} \
    --tsv {output.tsv:q} \
    --anchors {output.anchors:q} \
    &> {log:q}
test -s {output.bedpe:q}
test -s {output.tsv:q}
test -s {output.anchors:q}
        """

rule make_loop_universe:
    input:
        loops=expand(PATHS.sample_loops_bedpe("{sample}"), sample=CTX.sample_names)
    output:
        universe=PATHS.loop_universe(),
        group_consensus=PATHS.group_consensus_loops(),
        annotation=PATHS.loop_annotation()
    params:
        samples=",".join(CTX.sample_names),
        sample_table=lambda wildcards: config.get("samples", "config/samples.tsv"),
        mode=lambda wildcards: CTX.loop_universe_mode,
        anchor_slop=lambda wildcards: CTX.anchor_slop
    log:
        PATHS.log("loops/make_loop_universe")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.universe:q}) $(dirname {log:q})
python scripts/make_loop_universe.py \
    --loops {input.loops:q} \
    --samples {params.samples:q} \
    --sample-table {params.sample_table:q} \
    --mode {params.mode:q} \
    --anchor-slop {params.anchor_slop:q} \
    --universe {output.universe:q} \
    --group-consensus {output.group_consensus:q} \
    --annotation {output.annotation:q} \
    &> {log:q}
test -s {output.universe:q}
test -s {output.group_consensus:q}
test -s {output.annotation:q}
        """
