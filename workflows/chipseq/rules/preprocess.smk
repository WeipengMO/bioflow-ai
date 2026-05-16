if MODE == "pe":
    rule fastp:
        input:
            r1=raw_read1,
            r2=raw_read2
        output:
            r1=temp("clean_data/{sample}.R1.clean.fq.gz"),
            r2=temp("clean_data/{sample}.R2.clean.fq.gz")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            html="qc/fastp/{sample}.html",
            json="qc/fastp/{sample}.json",
            extra=lambda wildcards: config.get("fastp_extra", "")
        log:
            "logs/{sample}.fastp.log"
        threads:
            workflow_threads("fastp", 8)
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p clean_data qc/fastp logs
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


    rule align_reads:
        input:
            r1="clean_data/{sample}.R1.clean.fq.gz",
            r2="clean_data/{sample}.R2.clean.fq.gz"
        output:
            bam=temp("aligned_data/{sample}.sorted.bam"),
            bai=temp("aligned_data/{sample}.sorted.bam.bai")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            genome=lambda wildcards: config["genome"],
            extra=lambda wildcards: config.get("bowtie2_extra", "--dovetail -X 1000")
        log:
            "logs/{sample}.bowtie2.log"
        threads:
            workflow_threads("bowtie2", 16)
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p aligned_data logs
bowtie2 \
    -t \
    -p {threads} \
    {params.extra} \
    -x {params.genome:q} \
    -1 {input.r1:q} \
    -2 {input.r2:q} \
    2> {log:q} \
    | samtools sort -@ {threads} -O bam -o {output.bam:q} - 2>> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
            """
else:
    rule fastp:
        input:
            r1=raw_read1
        output:
            r1=temp("clean_data/{sample}.clean.fq.gz")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            html="qc/fastp/{sample}.html",
            json="qc/fastp/{sample}.json",
            extra=lambda wildcards: config.get("fastp_extra", "")
        log:
            "logs/{sample}.fastp.log"
        threads:
            workflow_threads("fastp", 8)
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p clean_data qc/fastp logs
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
            r1="clean_data/{sample}.clean.fq.gz"
        output:
            bam=temp("aligned_data/{sample}.sorted.bam"),
            bai=temp("aligned_data/{sample}.sorted.bam.bai")
        wildcard_constraints:
            sample=SAMPLE_PATTERN
        params:
            genome=lambda wildcards: config["genome"],
            extra=lambda wildcards: config.get("bowtie2_extra", "")
        log:
            "logs/{sample}.bowtie2.log"
        threads:
            workflow_threads("bowtie2", 16)
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p aligned_data logs
bowtie2 \
    -t \
    -p {threads} \
    {params.extra} \
    -x {params.genome:q} \
    -U {input.r1:q} \
    2> {log:q} \
    | samtools sort -@ {threads} -O bam -o {output.bam:q} - 2>> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
            """


rule mark_duplicates:
    input:
        bam="aligned_data/{sample}.sorted.bam",
        bai="aligned_data/{sample}.sorted.bam.bai"
    output:
        bam="aligned_data/{sample}.sorted.rmdup.bam",
        bai="aligned_data/{sample}.sorted.rmdup.bam.bai",
        metrics="qc/mark_duplicates/{sample}.metrics.txt"
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        picard=lambda wildcards: config.get("picard_path", "picard"),
        extra=lambda wildcards: config.get("mark_duplicates_extra", "REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01")
    log:
        "logs/{sample}.mark_duplicates.log"
    threads:
        workflow_threads("mark_duplicates", 4)
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p aligned_data qc/mark_duplicates logs
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


rule bam_coverage:
    input:
        bam="aligned_data/{sample}.sorted.rmdup.bam",
        bai="aligned_data/{sample}.sorted.rmdup.bam.bai"
    output:
        "tracks/{sample}.sorted.rmdup.CPM.bw"
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        extra=lambda wildcards: config.get("bam_coverage_extra", "--binSize 10 --normalizeUsing CPM --skipNonCoveredRegions")
    log:
        "logs/{sample}.bamCoverage.log"
    threads:
        workflow_threads("bam_coverage", 8)
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p tracks logs
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
        bw="tracks/{sample}.sorted.rmdup.CPM.bw"
    output:
        matrix=temp("deeptools_profile/{sample}.matrix.gz"),
        png="deeptools_profile/{sample}.scale.png"
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        regions=lambda wildcards: config["regions_bed"],
        matrix_extra=lambda wildcards: config.get("compute_matrix_extra", "scale-regions -b 1000 -a 1000 --skipZeros"),
        plot_extra=lambda wildcards: config.get("plot_profile_extra", "")
    log:
        "logs/{sample}.computeMatrix.log"
    threads:
        workflow_threads("compute_matrix", 8)
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p deeptools_profile logs
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
