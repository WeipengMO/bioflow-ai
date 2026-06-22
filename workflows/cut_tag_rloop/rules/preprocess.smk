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


if CTX.has_spikein:
    rule align_spikein:
        input:
            reads=clean_reads
        output:
            bam=temp(PATHS.spikein_bam("{sample}")),
            bai=temp(PATHS.spikein_bai("{sample}"))
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            genome=lambda wildcards: CTX.spikein_index,
            reads_arg=bowtie2_reads_arg,
            extra=lambda wildcards: spikein_bowtie2_extra()
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

    rule normalization_metrics:
        input:
            human_bams=expand(PATHS.signal_bam("{sample}"), sample=SAMPLES),
            human_unique_bams=expand(PATHS.filtered_bam("{sample}"), sample=SAMPLES),
            spikein_bams=expand(PATHS.spikein_bam("{sample}"), sample=SAMPLES)
        output:
            metrics=PATHS.normalization_metrics(),
            warning_tsv=PATHS.spikein_warning_tsv(),
            warning_txt=PATHS.spikein_warning_txt()
        params:
            script="scripts/normalization_metrics.py",
            samples=",".join(SAMPLES),
            mode=lambda wildcards: MODE,
            fastq_lines=lambda wildcards: fastq_metrics_lines(),
            spikein_group_lines=lambda wildcards: spikein_group_lines(),
            spikein_genome=lambda wildcards: CTX.spikein_genome,
            spikein_counting_mode=lambda wildcards: CTX.spikein_counting_mode,
            min_mapped=lambda wildcards: CTX.spikein_min_mapped_reads,
            min_fraction=lambda wildcards: CTX.spikein_min_fraction
        log:
            PATHS.log("normalization_metrics")
        threads:
            workflow_threads("normalization", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.metrics:q}) $(dirname {output.warning_tsv:q}) $(dirname {log:q})
manifest=$(mktemp)
groups=$(mktemp)
trap 'rm -f "$manifest" "$groups"' EXIT
printf '%s\n' {params.fastq_lines:q} > "$manifest"
printf '%s\n' {params.spikein_group_lines:q} > "$groups"
python {params.script:q} \
    --samples {params.samples:q} \
    --mode {params.mode:q} \
    --fastq-manifest "$manifest" \
    --spikein-groups "$groups" \
    --human-bams {input.human_bams:q} \
    --human-unique-bams {input.human_unique_bams:q} \
    --spikein-bams {input.spikein_bams:q} \
    --spikein-genome {params.spikein_genome:q} \
    --spikein-counting-mode {params.spikein_counting_mode:q} \
    --min-spikein-reads {params.min_mapped} \
    --warn-low-fraction {params.min_fraction} \
    --output {output.metrics:q} \
    --warning-tsv {output.warning_tsv:q} \
    --warning-txt {output.warning_txt:q} \
    &> {log:q}
test -s {output.metrics:q}
test -s {output.warning_tsv:q}
test -s {output.warning_txt:q}
            """


rule bam_coverage_cpm:
    input:
        bam=PATHS.signal_bam("{sample}"),
        bai=PATHS.signal_bai("{sample}")
    output:
        PATHS.bigwig_track("cpm", "{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        extra=lambda wildcards: config.get("bam_coverage_extra", "--binSize 10 --normalizeUsing CPM --skipNonCoveredRegions")
    log:
        PATHS.log("{sample}.bamCoverage.cpm")
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


if CTX.has_spikein:
    rule bam_coverage_matched_ref_spikein:
        input:
            bam=PATHS.signal_bam("{sample}"),
            bai=PATHS.signal_bai("{sample}"),
            metrics=PATHS.normalization_metrics()
        output:
            PATHS.bigwig_track("matched_ref_spikein", "{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            extra=lambda wildcards: config.get("spikein_bam_coverage_extra", "--binSize 10 --normalizeUsing None")
        log:
            PATHS.log("{sample}.bamCoverage.matched_ref_spikein")
        threads:
            workflow_threads("bam_coverage_scaled", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{for (i=1; i<=NF; i++) if ($i=="matched_anchor_spikein_scale_factor") col=i; next}} $1==sample {{print $col; found=1; exit}} END {{if (!found || !col) exit 2}}' {input.metrics:q})
if [[ "$scale_factor" == "NA" || -z "$scale_factor" ]]; then
    echo "Missing matched_anchor_spikein_scale_factor for {wildcards.sample}; see normalization metrics and warnings." > {log:q}
    exit 1
fi
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.extra} \
    --scaleFactor "$scale_factor" \
    &> {log:q}
test -s {output:q}
            """


    rule bam_coverage_absolute_spikein:
        input:
            bam=PATHS.signal_bam("{sample}"),
            bai=PATHS.signal_bai("{sample}"),
            metrics=PATHS.normalization_metrics()
        output:
            PATHS.bigwig_track("absolute_spikein", "{sample}")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            extra=lambda wildcards: config.get("spikein_bam_coverage_extra", "--binSize 10 --normalizeUsing None")
        log:
            PATHS.log("{sample}.bamCoverage.absolute_spikein")
        threads:
            workflow_threads("bam_coverage_scaled", 8)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
scale_factor=$(awk -F '\t' -v sample={wildcards.sample:q} 'NR==1 {{for (i=1; i<=NF; i++) if ($i=="absolute_spikein_scale_factor") col=i; next}} $1==sample {{print $col; found=1; exit}} END {{if (!found || !col) exit 2}}' {input.metrics:q})
if [[ "$scale_factor" == "NA" || -z "$scale_factor" ]]; then
    echo "Missing absolute_spikein_scale_factor for {wildcards.sample}; see normalization metrics and warnings." > {log:q}
    exit 1
fi
bamCoverage \
    --bam {input.bam:q} \
    -o {output:q} \
    --numberOfProcessors {threads} \
    {params.extra} \
    --scaleFactor "$scale_factor" \
    &> {log:q}
test -s {output:q}
            """


rule bigwig_header_qc:
    input:
        cpm=PATHS.bigwig_track("cpm", "{sample}")
    output:
        PATHS.bigwig_header_qc("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        script="scripts/bigwig_header_qc.py"
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
python {params.script:q} \
    --sample {wildcards.sample:q} \
    --cpm {input.cpm:q} \
    --output {output:q} \
    &> {log:q}
test -s {output:q}
        """
