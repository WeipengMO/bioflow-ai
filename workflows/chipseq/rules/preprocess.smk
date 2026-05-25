ruleorder: fastp_pe > fastp_se


rule fastp_pe:
    input:
        r1=raw_read1,
        r2=raw_read2
    output:
        r1=temp(PATHS.clean_r1("{sample}", "pe")),
        r2=temp(PATHS.clean_r2("{sample}"))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        html=PATHS.fastp_html("{sample}"),
        json=PATHS.fastp_json("{sample}"),
        extra=lambda wildcards: config.get("fastp_extra", "")
    log:
        PATHS.log("{sample}.fastp")
    threads:
        workflow_threads("fastp", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {params.html:q}) $(dirname {log:q})
fastp \
    -i {input.r1:q} \
    -I {input.r2:q} \
    -o {output.r1:q} \
    -O {output.r2:q} \
    -w {threads} \
    -h {params.html:q} \
    -j {params.json:q} \
    {params.extra} \
    &> {log:q}
test -s {output.r1:q}
test -s {output.r2:q}
        """


rule fastp_se:
    input:
        r1=raw_read1
    output:
        r1=temp(PATHS.clean_r1("{sample}", "se"))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        html=PATHS.fastp_html("{sample}"),
        json=PATHS.fastp_json("{sample}"),
        extra=lambda wildcards: config.get("fastp_extra", "")
    log:
        PATHS.log("{sample}.fastp")
    threads:
        workflow_threads("fastp", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {params.html:q}) $(dirname {log:q})
fastp \
    -i {input.r1:q} \
    -o {output.r1:q} \
    -w {threads} \
    -h {params.html:q} \
    -j {params.json:q} \
    {params.extra} \
    &> {log:q}
test -s {output.r1:q}
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
        extra=lambda wildcards: config.get("bowtie2_extra", "--dovetail -X 1000" if MODE == "pe" else "")
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
        view_extra=alignment_filter_view_extra
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

samtools view \
    -@ {threads} \
    -b \
    -q {params.min_mapq} \
    {params.view_extra} \
    {input.bam:q} \
    -o {output.bam:q} \
    &> {log:q}

samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule mark_duplicates:
    input:
        bam=mark_duplicates_input_bam,
        bai=mark_duplicates_input_bai
    output:
        bam=temp(PATHS.dedup_bam("{sample}")) if CTX.enable_blacklist_filter else PATHS.dedup_bam("{sample}"),
        bai=temp(PATHS.dedup_bai("{sample}")) if CTX.enable_blacklist_filter else PATHS.dedup_bai("{sample}"),
        metrics=PATHS.mark_duplicates_metrics("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        picard=lambda wildcards: config.get("picard_path", "picard"),
        extra=lambda wildcards: config.get("mark_duplicates_extra", "REMOVE_DUPLICATES=true")
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
    {params.extra} \
    &> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
test -s {output.metrics:q}
        """


rule filter_blacklist:
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
        PATHS.log("{sample}.blacklist_filter")
    threads:
        workflow_threads("blacklist_filter", 4)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})

if [[ -n {params.blacklist:q} ]]; then
    tmpbase="$(dirname {output.bam:q})/.tmp"
    mkdir -p "$tmpbase"
    tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.blacklist.XXXXXX")
    trap 'rm -rf "$tmpdir"' EXIT

    samtools view -H {input.bam:q} > "$tmpdir/header.sam" 2> {log:q}
    bedtools intersect \
        -v \
        -abam {input.bam:q} \
        -b {params.blacklist:q} \
        > "$tmpdir/filtered.bam" \
        2>> {log:q}
    samtools reheader "$tmpdir/header.sam" "$tmpdir/filtered.bam" > {output.bam:q} 2>> {log:q}
else
    samtools view -@ {threads} -b {input.bam:q} -o {output.bam:q} &> {log:q}
fi

samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
        """


rule bam_coverage:
    input:
        bam=lambda wildcards: analysis_bam(wildcards.sample),
        bai=lambda wildcards: analysis_bai(wildcards.sample)
    output:
        PATHS.bigwig_track("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        extra=lambda wildcards: config.get("bam_coverage_extra", "--binSize 10 --normalizeUsing CPM --skipNonCoveredRegions")
    log:
        PATHS.log("{sample}.bamCoverage")
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


rule compute_matrix_profile:
    input:
        bw=PATHS.bigwig_track("{sample}")
    output:
        matrix=temp(PATHS.matrix("{sample}")),
        png=PATHS.profile_png("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        regions=lambda wildcards: config["regions_bed"],
        matrix_extra=lambda wildcards: config.get("compute_matrix_extra", "scale-regions -b 1000 -a 1000 --skipZeros"),
        plot_extra=lambda wildcards: config.get("plot_profile_extra", "")
    log:
        PATHS.log("{sample}.computeMatrix")
    threads:
        workflow_threads("compute_matrix", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.matrix:q}) $(dirname {log:q})
computeMatrix \
    {params.matrix_extra} \
    -R {params.regions:q} \
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
