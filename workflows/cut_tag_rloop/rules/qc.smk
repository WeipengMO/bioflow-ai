if MODE == "pe" and ENABLE_INSERT_SIZE_QC:
    rule collect_insert_size_metrics:
        input:
            bam=PATHS.signal_bam("{sample}"),
            bai=PATHS.signal_bai("{sample}")
        output:
            metrics=PATHS.insert_size_metrics("{sample}"),
            histogram=PATHS.insert_size_histogram("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            picard=lambda wildcards: config.get("picard_path", "picard"),
            extra=lambda wildcards: config.get("insert_size_extra", "")
        log:
            PATHS.log("{sample}.insert_size")
        threads:
            workflow_threads("insert_size", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.metrics:q}) $(dirname {log:q})
if [[ "{params.picard}" == *.jar ]]; then
  PICARD_CMD="java -jar {params.picard}"
else
  PICARD_CMD="{params.picard}"
fi
$PICARD_CMD CollectInsertSizeMetrics \
    I={input.bam:q} \
    O={output.metrics:q} \
    H={output.histogram:q} \
    M=0.5 \
    {params.extra} \
    &> {log:q}
test -s {output.metrics:q}
test -s {output.histogram:q}
        """


rule frip_score:
    input:
        bam=PATHS.signal_bam("{sample}"),
        bai=PATHS.signal_bai("{sample}"),
        peaks=lambda wildcards: PATHS.peak(wildcards.sample)
    output:
        PATHS.frip("{sample}")
    wildcard_constraints:
        sample=TREATMENT_PATTERN
    log:
        PATHS.log("{sample}.frip")
    threads:
        workflow_threads("frip", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
total=$(samtools view -@ {threads} -c {input.bam:q} 2> {log:q})
if [[ -s {input.peaks:q} ]]; then
    in_peaks=$(bedtools intersect -u -abam {input.bam:q} -b {input.peaks:q} 2>> {log:q} | samtools view -@ {threads} -c - 2>> {log:q})
else
    in_peaks=0
    echo "Peak file is empty; FRiP alignments_in_peaks set to 0: {input.peaks}" >> {log:q}
fi
awk -v sample="{wildcards.sample}" -v total="$total" -v in_peaks="$in_peaks" 'BEGIN{{
    frip = total > 0 ? in_peaks / total : 0
    print "sample\ttotal_alignments\talignments_in_peaks\tfrip"
    printf "%s\t%d\t%d\t%.6f\n", sample, total, in_peaks, frip
}}' > {output:q}
test -s {output:q}
        """


rule rloop_qc_summary:
    input:
        counts=PATHS.peak_counts_raw() if CTX.count_matrix_enabled else [],
        peaks=expand(PATHS.peak("{sample}"), sample=TREATMENTS),
        frip=expand(PATHS.frip("{sample}"), sample=TREATMENTS),
        spikein=PATHS.normalization_metrics() if CTX.has_spikein else [],
        rnaseh_ratio=expand(PATHS.rnaseh_ratio_summary("{sample}"), sample=rnaseh_treatments(CTX)) if CTX.rnaseh_enabled and CTX.rnaseh_mode in {"ratio", "both"} else [],
        rnaseh_deseq2=expand(PATHS.rnaseh_deseq2_tsv("{contrast}"), contrast=rnaseh_contrasts(CTX)) if CTX.rnaseh_enabled and CTX.rnaseh_mode in {"deseq2", "both"} and CTX.count_matrix_enabled else []
    output:
        replicate_correlation=PATHS.replicate_correlation(),
        replicate_correlation_heatmap=PATHS.replicate_correlation_heatmap(),
        sample_pca=PATHS.sample_pca(),
        peak_width_distribution=PATHS.peak_width_distribution(),
        rnaseh_depletion_summary=PATHS.rnaseh_depletion_summary(),
        spikein_summary=PATHS.spikein_summary(),
        blacklist_mito_summary=PATHS.blacklist_mito_summary(),
        frip_fragment_level=PATHS.frip_fragment_level()
    params:
        script="scripts/rloop_qc_summary.py",
        spikein=lambda wildcards: PATHS.normalization_metrics() if CTX.has_spikein else ""
    log:
        PATHS.log("rloop_qc_summary")
    threads:
        workflow_threads("rloop_qc", 1)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.replicate_correlation:q}) $(dirname {log:q})
python {params.script:q} \
    --counts {input.counts:q} \
    --peaks {input.peaks:q} \
    --frip {input.frip:q} \
    --spikein {params.spikein:q} \
    --rnaseh-ratio-summaries {input.rnaseh_ratio:q} \
    --rnaseh-deseq2 {input.rnaseh_deseq2:q} \
    --replicate-correlation {output.replicate_correlation:q} \
    --replicate-correlation-heatmap {output.replicate_correlation_heatmap:q} \
    --sample-pca {output.sample_pca:q} \
    --peak-width-distribution {output.peak_width_distribution:q} \
    --rnaseh-depletion-summary {output.rnaseh_depletion_summary:q} \
    --spikein-summary {output.spikein_summary:q} \
    --blacklist-mito-summary {output.blacklist_mito_summary:q} \
    --frip-fragment-level {output.frip_fragment_level:q} \
    &> {log:q}
test -s {output.replicate_correlation:q}
test -s {output.replicate_correlation_heatmap:q}
test -s {output.sample_pca:q}
        """


if ENABLE_MULTIQC:
    rule multiqc_report:
        input:
            fastp=expand(PATHS.fastp_json("{sample}"), sample=SAMPLES),
            markdup=expand(PATHS.mark_duplicates_metrics("{sample}"), sample=SAMPLES),
            insert_size=expand(PATHS.insert_size_metrics("{sample}"), sample=SAMPLES) if MODE == "pe" and ENABLE_INSERT_SIZE_QC else [],
            frip=expand(PATHS.frip("{sample}"), sample=TREATMENTS),
            bigwig_scale=expand(PATHS.bigwig_header_qc("{sample}"), sample=SAMPLES) if "CPM" in SCALE_METHODS else [],
            normalization=PATHS.normalization_metrics() if CTX.has_spikein else [],
            rloop_qc=PATHS.replicate_correlation() if CTX.count_matrix_enabled else [],
            rnaseh=expand(PATHS.rnaseh_ratio_summary("{sample}"), sample=rnaseh_treatments(CTX)) if CTX.rnaseh_enabled and CTX.rnaseh_mode in {"ratio", "both"} else []
        output:
            PATHS.multiqc_report()
        params:
            outdir=PATHS.reports,
            extra=lambda wildcards: config.get("multiqc_extra", "")
        log:
            PATHS.log("multiqc")
        threads:
            workflow_threads("multiqc", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p {params.outdir:q} $(dirname {log:q})
multiqc \
    {OUTDIR:q} \
    -o {params.outdir:q} \
    -n multiqc_report.html \
    --force \
    {params.extra} \
    &> {log:q}
test -s {output:q}
            """
