if MODE == "pe":
    rule fastp_pe:
        input:
            r1=raw_read1,
            r2=raw_read2
        output:
            r1=temp(PATHS.clean_r1("{sample}", "pe")),
            r2=temp(PATHS.clean_r2("{sample}")),
            html=PATHS.fastp_html("{sample}"),
            json=PATHS.fastp_json("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            extra=lambda wildcards: " ".join(
                part for part in [config.get("fastp_extra", ""), config.get("fastp_pe_extra", "--detect_adapter_for_pe")] if part
            )
        log:
            PATHS.log("{sample}.fastp")
        threads:
            workflow_threads("fastp", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {output.html:q}) $(dirname {log:q})
fastp \
    -i {input.r1:q} \
    -I {input.r2:q} \
    -o {output.r1:q} \
    -O {output.r2:q} \
    -w {threads} \
    -h {output.html:q} \
    -j {output.json:q} \
    {params.extra} \
    &> {log:q}
test -s {output.r1:q}
test -s {output.r2:q}
test -s {output.html:q}
test -s {output.json:q}
        """


if MODE == "se":
    rule fastp_se:
        input:
            r1=raw_read1
        output:
            r1=temp(PATHS.clean_r1("{sample}", "se")),
            html=PATHS.fastp_html("{sample}"),
            json=PATHS.fastp_json("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            extra=lambda wildcards: " ".join(
                part for part in [config.get("fastp_extra", ""), config.get("fastp_se_extra", "")] if part
            )
        log:
            PATHS.log("{sample}.fastp")
        threads:
            workflow_threads("fastp", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {output.html:q}) $(dirname {log:q})
fastp \
    -i {input.r1:q} \
    -o {output.r1:q} \
    -w {threads} \
    -h {output.html:q} \
    -j {output.json:q} \
    {params.extra} \
    &> {log:q}
test -s {output.r1:q}
test -s {output.html:q}
test -s {output.json:q}
        """


rule align_reads:
    input:
        reads=clean_reads
    output:
        bam=temp(PATHS.sorted_bam("{sample}")),
        bai=temp(PATHS.sorted_bai("{sample}"))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        genome=lambda wildcards: config["genome"],
        reads_arg=bowtie2_reads_arg,
        read_group_arg=bowtie2_read_group_arg,
        extra=lambda wildcards: config.get("bowtie2_extra", "--very-sensitive")
    log:
        PATHS.log("{sample}.bowtie2")
    threads:
        workflow_threads("bowtie2", 16)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
bowtie2 \
    -t \
    -p {threads} \
    {params.extra} \
    {params.read_group_arg} \
    -x {params.genome:q} \
    {params.reads_arg} \
    2> {log:q} \
    | samtools sort -@ {threads} -O bam -o {output.bam:q} - 2>> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule filter_aligned:
    input:
        bam=PATHS.sorted_bam("{sample}"),
        bai=PATHS.sorted_bai("{sample}")
    output:
        bam=temp(PATHS.alignment_filtered_bam("{sample}")),
        bai=temp(PATHS.alignment_filtered_bai("{sample}"))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        min_mapq=lambda wildcards: int(config.get("alignment_filter_min_mapq", 30)),
        view_extra=lambda wildcards: config.get("alignment_filter_view_extra") or ("-F 3852" if MODE == "pe" else "-F 3844"),
        mito_chromosomes=lambda wildcards: ",".join(as_list(config.get("mitochondrial_chromosomes", ["chrM", "MT"])))
    log:
        PATHS.log("{sample}.alignment_filter")
    threads:
        workflow_threads("alignment_filter", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})

tmpbase="$(dirname {output.bam:q})/.tmp"
mkdir -p "$tmpbase"
tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

samtools view \
    -@ {threads} \
    -b \
    -q {params.min_mapq} \
    {params.view_extra} \
    {input.bam:q} \
    -o "$tmpdir/{wildcards.sample}.mapq.bam" \
    &> {log:q}

mito_chromosomes={params.mito_chromosomes:q}
if [[ -n "$mito_chromosomes" ]]; then
    samtools view -@ {threads} -h "$tmpdir/{wildcards.sample}.mapq.bam" \
        | awk -v mito="$mito_chromosomes" 'BEGIN{{
        split(mito, names, ",")
        for (i in names) skip[names[i]] = 1
    }}
    /^@/ {{print; next}}
    !($3 in skip) {{print}}' \
        | samtools view -@ {threads} -b -o {output.bam:q} - \
        2>> {log:q}
else
    mv "$tmpdir/{wildcards.sample}.mapq.bam" {output.bam:q}
fi

samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule mark_duplicates:
    input:
        bam=PATHS.alignment_filtered_bam("{sample}"),
        bai=PATHS.alignment_filtered_bai("{sample}")
    output:
        bam=temp(PATHS.markdup_bam("{sample}")),
        bai=temp(PATHS.markdup_bai("{sample}")),
        metrics=PATHS.mark_duplicates_metrics("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        picard=lambda wildcards: config.get("picard_path", "picard"),
        extra=lambda wildcards: CTX.mark_duplicates_extra
    log:
        PATHS.log("{sample}.mark_duplicates")
    threads:
        workflow_threads("mark_duplicates", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {output.metrics:q}) $(dirname {log:q})

if [[ "{params.picard}" == *.jar ]]; then
  PICARD_CMD="java -jar {params.picard}"
else
  PICARD_CMD="{params.picard}"
fi
$PICARD_CMD MarkDuplicates \
    I={input.bam:q} \
    O={output.bam:q} \
    M={output.metrics:q} \
    REMOVE_DUPLICATES=false \
    {params.extra} \
    &> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
test -s {output.metrics:q}
        """


rule remove_duplicates_for_peaks:
    input:
        bam=PATHS.markdup_bam("{sample}"),
        bai=PATHS.markdup_bai("{sample}")
    output:
        bam=temp(PATHS.dedup_bam("{sample}")),
        bai=temp(PATHS.dedup_bai("{sample}"))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        mode=lambda wildcards: CTX.peak_duplicate_mode
    log:
        PATHS.log("{sample}.remove_duplicates_for_peaks")
    threads:
        workflow_threads("remove_duplicates", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
if [[ "{params.mode}" == "remove" ]]; then
    samtools view -@ {threads} -b -F 1024 {input.bam:q} -o {output.bam:q} &> {log:q}
else
    samtools view -@ {threads} -b {input.bam:q} -o {output.bam:q} &> {log:q}
fi
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule filter_blacklist_signal:
    input:
        bam=PATHS.markdup_bam("{sample}"),
        bai=PATHS.markdup_bai("{sample}")
    output:
        bam=PATHS.signal_bam("{sample}"),
        bai=PATHS.signal_bai("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        blacklist=lambda wildcards: str(config.get("blacklist", "") or "")
    log:
        PATHS.log("{sample}.blacklist_filter.signal")
    threads:
        workflow_threads("blacklist_filter", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
blacklist={params.blacklist:q}
if [[ -n "$blacklist" ]]; then
    tmpbase="$(dirname {output.bam:q})/.tmp"
    mkdir -p "$tmpbase"
    tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.signal.blacklist.XXXXXX")
    trap 'rm -rf "$tmpdir"' EXIT
    samtools view -H {input.bam:q} > "$tmpdir/header.sam" 2> {log:q}
    bedtools intersect -v -abam {input.bam:q} -b "$blacklist" > "$tmpdir/filtered.bam" 2>> {log:q}
    samtools reheader "$tmpdir/header.sam" "$tmpdir/filtered.bam" > {output.bam:q} 2>> {log:q}
else
    samtools view -@ {threads} -b {input.bam:q} -o {output.bam:q} &> {log:q}
fi
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule filter_blacklist_peak:
    input:
        bam=PATHS.dedup_bam("{sample}"),
        bai=PATHS.dedup_bai("{sample}")
    output:
        bam=PATHS.filtered_bam("{sample}"),
        bai=PATHS.filtered_bai("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        blacklist=lambda wildcards: str(config.get("blacklist", "") or "")
    log:
        PATHS.log("{sample}.blacklist_filter.peak")
    threads:
        workflow_threads("blacklist_filter", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
blacklist={params.blacklist:q}
if [[ -n "$blacklist" ]]; then
    tmpbase="$(dirname {output.bam:q})/.tmp"
    mkdir -p "$tmpbase"
    tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.peak.blacklist.XXXXXX")
    trap 'rm -rf "$tmpdir"' EXIT
    samtools view -H {input.bam:q} > "$tmpdir/header.sam" 2> {log:q}
    bedtools intersect -v -abam {input.bam:q} -b "$blacklist" > "$tmpdir/filtered.bam" 2>> {log:q}
    samtools reheader "$tmpdir/header.sam" "$tmpdir/filtered.bam" > {output.bam:q} 2>> {log:q}
else
    samtools view -@ {threads} -b {input.bam:q} -o {output.bam:q} &> {log:q}
fi
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


if CTX.signal_scale_factor_method == "spikein":
    rule align_spikein:
        input:
            reads=clean_reads
        output:
            bam=temp(PATHS.spikein_bam("{sample}")),
            bai=temp(PATHS.spikein_bai("{sample}"))
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            genome=lambda wildcards: config["spikein_genome"],
            reads_arg=bowtie2_reads_arg,
            extra=lambda wildcards: config.get("bowtie2_extra", "--very-sensitive")
        log:
            PATHS.log("{sample}.bowtie2_spikein")
        threads:
            workflow_threads("bowtie2", 16)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
bowtie2 \
    -t \
    -p {threads} \
    {params.extra} \
    --no-unal \
    -x {params.genome:q} \
    {params.reads_arg} \
    2> {log:q} \
    | samtools sort -@ {threads} -O bam -o {output.bam:q} - 2>> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
            """

    rule count_spikein_alignments:
        input:
            bam=PATHS.spikein_bam("{sample}")
        output:
            PATHS.spikein_counts("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        log:
            PATHS.log("{sample}.spikein_count")
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
count=$(samtools view -c -F 4 {input.bam:q})
printf "sample\tspikein_aligned\n%s\t%s\n" "{wildcards.sample}" "$count" > {output:q}
test -s {output:q}
            """


rule signal_scale_factors:
    input:
        bams=expand(PATHS.signal_bam("{sample}"), sample=SAMPLES),
        source=lambda wildcards: config.get("signal_scale_factors_tsv", "") if CTX.signal_scale_factor_method == "tsv" else [],
        spikein_counts=lambda wildcards: expand(PATHS.spikein_counts("{sample}"), sample=SAMPLES) if CTX.signal_scale_factor_method == "spikein" else []
    output:
        PATHS.signal_scale_factors()
    params:
        script="scripts/compute_scale_factors.py",
        samples=",".join(SAMPLES),
        method=lambda wildcards: CTX.signal_scale_factor_method,
        mode=lambda wildcards: MODE,
        source=lambda wildcards: config.get("signal_scale_factors_tsv", "") if CTX.signal_scale_factor_method == "tsv" else "",
        spikein_method=lambda wildcards: CTX.spikein_method if CTX.signal_scale_factor_method == "spikein" else "",
        spikein_ref=lambda wildcards: CTX.spikein_reference_sample if CTX.signal_scale_factor_method == "spikein" else ""
    log:
        PATHS.log("signal_scale_factors")
    threads:
        workflow_threads("scale_factors", 2)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
source_arg=""
if [[ -n "{params.source}" ]]; then
    source_arg="--source-tsv {params.source:q}"
fi
spikein_arg=""
if [[ -n "{params.spikein_method}" ]]; then
    merged=$(mktemp)
    first=1
    for f in {input.spikein_counts:q}; do
        if [[ "$first" == "1" ]]; then
            cat "$f" > "$merged"
            first=0
        else
            tail -n +2 "$f" >> "$merged"
        fi
    done
    spikein_arg="--spikein-counts $merged --spikein-method {params.spikein_method} --spikein-reference-sample {params.spikein_ref}"
    trap 'rm -f "$merged"' EXIT
fi
python {params.script:q} \
    --samples {params.samples:q} \
    --bams {input.bams:q} \
    --mode {params.mode:q} \
    --method {params.method:q} \
    $source_arg \
    $spikein_arg \
    --output {output:q} \
    &> {log:q}
test -s {output:q}
        """


rule dedup_scale_factors:
    input:
        bams=expand(PATHS.filtered_bam("{sample}"), sample=SAMPLES)
    output:
        PATHS.dedup_scale_factors()
    params:
        script="scripts/compute_scale_factors.py",
        samples=",".join(SAMPLES),
        mode=lambda wildcards: MODE
    log:
        PATHS.log("dedup_scale_factors")
    threads:
        workflow_threads("scale_factors", 2)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
python {params.script:q} \
    --samples {params.samples:q} \
    --bams {input.bams:q} \
    --mode {params.mode:q} \
    --method effective_fragments \
    --output {output:q} \
    &> {log:q}
test -s {output:q}
        """


rule bam_coverage_cpm:
    input:
        bam=PATHS.signal_bam("{sample}"),
        bai=PATHS.signal_bai("{sample}")
    output:
        PATHS.bigwig_track("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        extra=lambda wildcards: config.get("bam_coverage_extra", "--binSize 10 --normalizeUsing CPM --skipNonCoveredRegions")
    log:
        PATHS.log("{sample}.bamCoverage.CPM")
    threads:
        workflow_threads("bam_coverage", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.extra} \
    &> {log:q}
test -s {output:q}
        """


if CTX.enable_common_scale_bigwig:
    rule bam_coverage_common_scale:
        input:
            bam=PATHS.signal_bam("{sample}"),
            bai=PATHS.signal_bai("{sample}"),
            scales=PATHS.signal_scale_factors()
        output:
            PATHS.common_scale_bigwig_track("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            args=lambda wildcards: common_scale_base_args()
        log:
            PATHS.log("{sample}.bamCoverage.common_scale")
        threads:
            workflow_threads("bam_coverage_common_scale", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{next}} $1==sample {{print $3; found=1; exit}} END {{if (!found) exit 2}}' {input.scales:q})
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.args} \
    --scaleFactor "$scale_factor" \
    &> {log:q}
test -s {output:q}
            """


if CTX.enable_bigwig_debug_tracks:
    rule bam_coverage_debug_raw_scale:
        input:
            bam=PATHS.signal_bam("{sample}"),
            bai=PATHS.signal_bai("{sample}")
        output:
            PATHS.debug_raw_bigwig_track("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            extra=lambda wildcards: config.get("debug_raw_bam_coverage_extra", "--binSize 10 --normalizeUsing None")
        log:
            PATHS.log("{sample}.bamCoverage.debug_raw_scale")
        threads:
            workflow_threads("bam_coverage_debug", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.extra} \
    --scaleFactor 1 \
    &> {log:q}
test -s {output:q}
            """

    rule bam_coverage_debug_dedup_scale:
        input:
            bam=PATHS.filtered_bam("{sample}"),
            bai=PATHS.filtered_bai("{sample}"),
            scales=PATHS.dedup_scale_factors()
        output:
            PATHS.debug_dedup_bigwig_track("{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            args=lambda wildcards: common_scale_base_args()
        log:
            PATHS.log("{sample}.bamCoverage.debug_dedup_scale")
        threads:
            workflow_threads("bam_coverage_debug", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{next}} $1==sample {{print $3; found=1; exit}} END {{if (!found) exit 2}}' {input.scales:q})
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.args} \
    --scaleFactor "$scale_factor" \
    &> {log:q}
test -s {output:q}
            """


rule bigwig_header_qc:
    input:
        cpm=PATHS.bigwig_track("{sample}"),
        common=lambda wildcards: PATHS.common_scale_bigwig_track(wildcards.sample) if CTX.enable_common_scale_bigwig else [],
        scales=lambda wildcards: PATHS.signal_scale_factors() if CTX.enable_common_scale_bigwig else [],
        debug_raw=lambda wildcards: PATHS.debug_raw_bigwig_track(wildcards.sample) if CTX.enable_bigwig_debug_tracks else [],
        debug_dedup=lambda wildcards: PATHS.debug_dedup_bigwig_track(wildcards.sample) if CTX.enable_bigwig_debug_tracks else [],
        dedup_scales=lambda wildcards: PATHS.dedup_scale_factors() if CTX.enable_bigwig_debug_tracks else []
    output:
        PATHS.bigwig_header_qc("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        script="scripts/bigwig_header_qc.py",
        common_scale_enabled=lambda wildcards: "true" if CTX.enable_common_scale_bigwig else "false",
        debug_enabled=lambda wildcards: "true" if CTX.enable_bigwig_debug_tracks else "false",
        common_scale_normalization=lambda wildcards: CTX.common_scale_normalization if CTX.enable_common_scale_bigwig else "NA"
    log:
        PATHS.log("{sample}.bigwig_header_qc")
    threads:
        workflow_threads("bigwig_header_qc", 1)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
common_arg=""
if [[ "{params.common_scale_enabled}" == "true" ]]; then
    scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{next}} $1==sample {{print $3; found=1; exit}} END {{if (!found) exit 2}}' {input.scales:q})
    common_arg="--common-scale {input.common:q} --scale-factor $scale_factor --common-scale-normalization {params.common_scale_normalization}"
fi
debug_arg=""
if [[ "{params.debug_enabled}" == "true" ]]; then
    dedup_scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{next}} $1==sample {{print $3; found=1; exit}} END {{if (!found) exit 2}}' {input.dedup_scales:q})
    debug_arg="--debug-raw {input.debug_raw:q} --debug-dedup {input.debug_dedup:q} --dedup-scale-factor $dedup_scale_factor"
fi
python {params.script:q} \
    --sample {wildcards.sample:q} \
    --cpm {input.cpm:q} \
    $common_arg \
    $debug_arg \
    --output {output:q} \
    &> {log:q}
test -s {output:q}
        """
