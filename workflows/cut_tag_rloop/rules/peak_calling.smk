rule macs3_callpeak:
    input:
        bam=peak_input_bam,
        bai=peak_input_bai,
        control=macs3_control_input
    output:
        PATHS.peak("{sample}")
    wildcard_constraints:
        sample=TREATMENT_PATTERN
    params:
        outdir=lambda wildcards: PATHS.peaks,
        name=lambda wildcards: wildcards.sample + ".cut_tag_rloop",
        genome_size=lambda wildcards: config["gsize"],
        fmt=lambda wildcards: macs3_format(),
        control_arg=macs3_control_arg,
        extra=lambda wildcards: macs3_extra()
    log:
        PATHS.log("{sample}.macs3")
    threads:
        workflow_threads("macs3", 1)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p {params.outdir:q} $(dirname {log:q})
macs3 callpeak \
    -t {input.bam:q} \
    {params.control_arg} \
    -f {params.fmt} \
    -g {params.genome_size:q} \
    -n {params.name:q} \
    --keep-dup all \
    --outdir {params.outdir:q} \
    {params.extra} \
    &> {log:q}
test -s {output:q}
        """


if CTX.count_matrix_enabled:
    rule peak_count_matrix:
        input:
            peaks=count_matrix_peak_inputs,
            bams=count_matrix_bam_inputs,
            normalization=PATHS.normalization_metrics() if CTX.has_spikein else []
        output:
            universe=PATHS.peak_universe(),
            saf=PATHS.peak_universe_saf(),
            featurecounts=PATHS.peak_counts_featurecounts(),
            raw=PATHS.peak_counts_raw(),
            cpm=PATHS.peak_counts_cpm(),
            spikein=PATHS.peak_counts_spikein_normalized() if CTX.has_spikein else [],
            annotation=PATHS.peak_annotation_input()
        params:
            script="scripts/peak_count_matrix.py",
            samples=",".join(SAMPLES),
            mode=lambda wildcards: MODE,
            universe_mode=lambda wildcards: CTX.count_matrix_peak_universe,
            min_support=lambda wildcards: CTX.consensus_min_support,
            min_width=lambda wildcards: CTX.count_matrix_min_peak_width,
            merge_distance=lambda wildcards: int(CTX.count_matrix_merge_distance),
            min_mapq=lambda wildcards: int(CTX.count_matrix_min_mapq),
            featurecounts_extra=lambda wildcards: CTX.count_matrix_featurecounts_extra,
            normalization_metrics=lambda wildcards: PATHS.normalization_metrics() if CTX.has_spikein else "",
            output_spikein=lambda wildcards: PATHS.peak_counts_spikein_normalized() if CTX.has_spikein else ""
        log:
            PATHS.log("peak_count_matrix")
        threads:
            workflow_threads("count_matrix", 4)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.universe:q}) $(dirname {log:q})
python {params.script:q} \
    --samples {params.samples:q} \
    --mode {params.mode:q} \
    --peak-files {input.peaks:q} \
    --bams {input.bams:q} \
    --universe-mode {params.universe_mode:q} \
    --min-support {params.min_support} \
    --min-peak-width {params.min_width} \
    --merge-distance {params.merge_distance} \
    --min-mapq {params.min_mapq} \
    --featurecounts-extra {params.featurecounts_extra:q} \
    --normalization-metrics {params.normalization_metrics:q} \
    --output-universe {output.universe:q} \
    --output-saf {output.saf:q} \
    --output-featurecounts {output.featurecounts:q} \
    --output-raw {output.raw:q} \
    --output-cpm {output.cpm:q} \
    --output-spikein {params.output_spikein:q} \
    --output-annotation-bed {output.annotation:q} \
    --threads {threads} \
    &> {log:q}
test -s {output.universe:q}
test -s {output.raw:q}
test -s {output.cpm:q}
test -s {output.annotation:q}
            """


if CTX.rnaseh_enabled and CTX.rnaseh_mode in {"ratio", "both"}:
    rule rnaseh_sensitive_ratio:
        input:
            treatment_peaks=PATHS.peak("{sample}"),
            treatment_bw=rnaseh_sensitive_treatment_bigwig,
            rnaseh_bw=rnaseh_sensitive_control_bigwig
        output:
            table=PATHS.rnaseh_ratio_table("{sample}"),
            bed=PATHS.rnaseh_ratio_bed("{sample}"),
            summary=PATHS.rnaseh_ratio_summary("{sample}")
        wildcard_constraints:
            sample=TREATMENT_PATTERN
        params:
            script="scripts/rnaseh_signal.py",
            rnaseh_sample=lambda wildcards: CTX.rnaseh_controls[wildcards.sample],
            min_fc=lambda wildcards: CTX.rnaseh_signal_min_fold_change,
            min_treatment_signal=lambda wildcards: CTX.rnaseh_signal_min_treatment_signal,
            min_abs_signal_diff=lambda wildcards: CTX.rnaseh_min_abs_signal_diff,
            pseudocount=lambda wildcards: CTX.rnaseh_pseudocount,
            scale_method=lambda wildcards: "matched_ref_spikein" if CTX.has_spikein and "matched_ref_spikein" in CTX.scale_methods else "cpm"
        log:
            PATHS.log("{sample}.rnaseh_sensitive_ratio")
        threads:
            workflow_threads("rnaseh_sensitive", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.table:q}) $(dirname {output.bed:q}) $(dirname {output.summary:q}) $(dirname {log:q})
python {params.script:q} \
    --sample {wildcards.sample:q} \
    --rnaseh-sample {params.rnaseh_sample:q} \
    --treatment-peaks {input.treatment_peaks:q} \
    --treatment-bw {input.treatment_bw:q} \
    --rnaseh-bw {input.rnaseh_bw:q} \
    --output-tsv {output.table:q} \
    --output-bed {output.bed:q} \
    --summary {output.summary:q} \
    --min-fold-change {params.min_fc} \
    --min-treatment-signal {params.min_treatment_signal} \
    --min-abs-signal-diff {params.min_abs_signal_diff} \
    --pseudocount {params.pseudocount} \
    --scale-method {params.scale_method:q} \
    &> {log:q}
test -s {output.table:q}
test -s {output.summary:q}
            """


if CTX.rnaseh_enabled and CTX.rnaseh_mode in {"deseq2", "both"} and CTX.count_matrix_enabled:
    rule rnaseh_depleted_deseq2:
        input:
            raw=PATHS.peak_counts_raw(),
            universe=PATHS.peak_universe(),
            spikein=PATHS.normalization_metrics() if CTX.has_spikein else []
        output:
            tsv=PATHS.rnaseh_deseq2_tsv("{contrast}"),
            bed=PATHS.rnaseh_deseq2_bed("{contrast}")
        wildcard_constraints:
            contrast=GROUP_PATTERN
        params:
            script="scripts/rnaseh_depleted_deseq2.R",
            design=lambda wildcards: rnaseh_contrast_pairs_lines(wildcards),
            spikein_metrics=lambda wildcards: PATHS.normalization_metrics() if CTX.has_spikein else "",
            fdr=lambda wildcards: CTX.rnaseh_fdr_threshold,
            log2fc=lambda wildcards: CTX.rnaseh_log2fc_threshold
        log:
            PATHS.log("{contrast}.rnaseh_depleted_deseq2")
        threads:
            workflow_threads("rnaseh_deseq2", 1)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.tsv:q}) $(dirname {log:q})
design=$(mktemp)
trap 'rm -f "$design"' EXIT
printf '%s\n' {params.design:q} > "$design"
Rscript {params.script:q} \
    --counts {input.raw:q} \
    --universe {input.universe:q} \
    --design "$design" \
    --contrast {wildcards.contrast:q} \
    --spikein-metrics {params.spikein_metrics:q} \
    --fdr {params.fdr} \
    --log2fc {params.log2fc} \
    --output-tsv {output.tsv:q} \
    --output-bed {output.bed:q} \
    &> {log:q}
test -s {output.tsv:q}
            """


rule intersect_peaks:
    input:
        peaks=replicate_peak_inputs
    output:
        PATHS.intersect_peak("{group}")
    wildcard_constraints:
        group=GROUP_PATTERN
    log:
        PATHS.log("{group}.intersect_peaks")
    threads:
        workflow_threads("intersect_peaks", 2)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
tmpbase="$(dirname {output:q})/.tmp"
mkdir -p "$tmpbase"
tmpdir=$(mktemp -d "$tmpbase/{wildcards.group}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

{{
    i=0
    for f in {input.peaks:q}; do
        i=$((i + 1))
        awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' "$f" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > "$tmpdir/rep_${{i}}.bed"
    done
    if [ "$i" -lt 2 ]; then
        echo "At least two replicate peak files are required." >&2
        exit 1
    fi
    cp "$tmpdir/rep_1.bed" "$tmpdir/current.bed"
    for rep in "$tmpdir"/rep_*.bed; do
        if [[ "$rep" == "$tmpdir/rep_1.bed" ]]; then
            continue
        fi
        bedtools intersect -a "$tmpdir/current.bed" -b "$rep" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > "$tmpdir/next.bed"
        mv "$tmpdir/next.bed" "$tmpdir/current.bed"
    done
    awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "intersect_peak_"NR}}' "$tmpdir/current.bed" > {output:q}
    echo "Intersect peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """


rule consensus_peaks:
    input:
        peaks=replicate_peak_inputs
    output:
        PATHS.consensus_peak("{group}")
    wildcard_constraints:
        group=GROUP_PATTERN
    log:
        PATHS.log("{group}.consensus_peaks")
    params:
        min_support=consensus_min_support
    threads:
        workflow_threads("consensus_peaks", 2)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
tmpbase="$(dirname {output:q})/.tmp"
mkdir -p "$tmpbase"
tmpdir=$(mktemp -d "$tmpbase/{wildcards.group}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

{{
    i=0
    for f in {input.peaks:q}; do
        i=$((i + 1))
        awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' "$f" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > "$tmpdir/rep_${{i}}.bed"
    done
    if [ "$i" -lt 2 ]; then
        echo "At least two replicate peak files are required." >&2
        exit 1
    fi
    cat "$tmpdir"/rep_*.bed \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "candidate_"NR}}' \
        > "$tmpdir/candidates.bed"
    : > "$tmpdir/support_ids.txt"
    for rep in "$tmpdir"/rep_*.bed; do
        bedtools intersect -u -a "$tmpdir/candidates.bed" -b "$rep" | cut -f4 >> "$tmpdir/support_ids.txt"
    done
    awk 'BEGIN{{OFS="\t"}} {{support[$1]++}} END{{for (id in support) print id, support[id]}}' "$tmpdir/support_ids.txt" > "$tmpdir/support.tsv"
    awk -v min_support="{params.min_support}" 'BEGIN{{OFS="\t"}}
        FNR==NR {{support[$1]=$2; next}}
        {{
            s = support[$4] + 0
            if (s >= min_support) print $1, $2, $3, "consensus_peak_"++n, s, "."
        }}' "$tmpdir/support.tsv" "$tmpdir/candidates.bed" > {output:q}
    echo "Consensus peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
