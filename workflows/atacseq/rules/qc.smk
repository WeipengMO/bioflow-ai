rule collect_insert_size_metrics:
    input:
        bam=PATHS.filtered_bam("{sample}"),
        bai=PATHS.filtered_bai("{sample}")
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
        bam=PATHS.filtered_bam("{sample}"),
        bai=PATHS.filtered_bai("{sample}"),
        peaks=PATHS.peak("{sample}")
    output:
        PATHS.frip("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
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
in_peaks=$(bedtools intersect -u -abam {input.bam:q} -b {input.peaks:q} 2>> {log:q} | samtools view -@ {threads} -c - 2>> {log:q})
awk -v sample="{wildcards.sample}" -v total="$total" -v in_peaks="$in_peaks" 'BEGIN{{
    frip = total > 0 ? in_peaks / total : 0
    print "sample\ttotal_alignments\talignments_in_peaks\tfrip"
    printf "%s\t%d\t%d\t%.6f\n", sample, total, in_peaks, frip
}}' > {output:q}
test -s {output:q}
        """


rule tss_enrichment:
    input:
        bw=PATHS.bigwig_track("{sample}")
    output:
        matrix=temp(PATHS.tss_matrix("{sample}")),
        png=PATHS.tss_profile("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        tss=lambda wildcards: config["tss_bed"],
        matrix_extra=lambda wildcards: config.get("compute_matrix_tss_extra", "reference-point --referencePoint TSS -b 2000 -a 2000 --skipZeros"),
        plot_extra=lambda wildcards: config.get("plot_profile_extra", "--perGroup")
    log:
        PATHS.log("{sample}.tss_enrichment")
    threads:
        workflow_threads("tss_enrichment", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.matrix:q}) $(dirname {log:q})
computeMatrix \
    {params.matrix_extra} \
    -R {params.tss:q} \
    -S {input.bw:q} \
    -o {output.matrix:q} \
    -p {threads} \
    &> {log:q}
plotProfile \
    -m {output.matrix:q} \
    -out {output.png:q} \
    {params.plot_extra} \
    2>> {log:q}
test -s {output.matrix:q}
test -s {output.png:q}
        """


if ENABLE_MULTIQC:
    rule multiqc_report:
        input:
            fastp=expand(PATHS.fastp_json("{sample}"), sample=SAMPLES),
            markdup=expand(PATHS.mark_duplicates_metrics("{sample}"), sample=SAMPLES),
            insert_size=expand(PATHS.insert_size_metrics("{sample}"), sample=SAMPLES) if MODE == "pe" else [],
            frip=expand(PATHS.frip("{sample}"), sample=SAMPLES)
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
    {params.extra} \
    &> {log:q}
test -s {output:q}
            """
