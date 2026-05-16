import shlex


def macs3_input(wildcards):
    peak_mode, peak_type = parse_peak_dir(wildcards.peak_dir)
    expected_ext = PEAK_EXTENSIONS[peak_type]
    if wildcards.peak_ext != expected_ext:
        raise ValueError(
            f"Inconsistent MACS3 output extension for {wildcards.peak_dir}: "
            f"expected {expected_ext}, got {wildcards.peak_ext}."
        )

    inputs = {
        "treatment": f"aligned_data/{wildcards.treatment}.sorted.rmdup.bam",
    }
    if peak_mode == "with_control":
        if wildcards.treatment not in TREATMENT_TO_CONTROL:
            raise ValueError(f"No control configured for treatment: {wildcards.treatment}")
        control = TREATMENT_TO_CONTROL[wildcards.treatment]
        inputs["control"] = f"aligned_data/{control}.sorted.rmdup.bam"
    return inputs


def macs3_control_arg(wildcards, input):
    peak_mode, _ = parse_peak_dir(wildcards.peak_dir)
    if peak_mode == "with_control":
        return "-c " + shlex.quote(str(input.control))
    return ""


def macs3_broad_arg(wildcards):
    _, peak_type = parse_peak_dir(wildcards.peak_dir)
    if peak_type == "broad":
        broad_cutoff = config.get("macs3_broad_cutoff", 0.1)
        return f"--broad --broad-cutoff {broad_cutoff}"
    return ""


def pooled_control_inputs(wildcards):
    return expand(
        "aligned_data/{sample}.sorted.rmdup.bam",
        sample=POOLED_CONTROL_GROUPS[str(wildcards.pool)],
    )


if USE_POOLED_CONTROL:
    rule merge_pooled_control:
        input:
            pooled_control_inputs
        output:
            bam="aligned_data/{pool}.sorted.rmdup.bam",
            bai="aligned_data/{pool}.sorted.rmdup.bam.bai"
        wildcard_constraints:
            pool=POOLED_CONTROL_PATTERN
        log:
            "logs/{pool}.merge_control.log"
        threads:
            workflow_threads("merge_control", 4)
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p aligned_data logs
samtools merge -@ {threads} -f {output.bam:q} {input:q} &> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
            """


rule macs3_callpeak:
    input:
        unpack(macs3_input)
    output:
        "macs3_results/{peak_dir}/{treatment}_peaks.{peak_ext}"
    wildcard_constraints:
        peak_dir="narrow|broad|narrow_no_control|broad_no_control",
        peak_ext="narrowPeak|broadPeak",
        treatment=TREATMENT_PATTERN
    params:
        outdir=lambda wildcards: f"macs3_results/{wildcards.peak_dir}",
        name=lambda wildcards: wildcards.treatment,
        genome_size=lambda wildcards: config["gsize"],
        fmt=lambda wildcards: macs3_format(),
        control_arg=macs3_control_arg,
        broad_arg=macs3_broad_arg,
        extra=lambda wildcards: config.get("macs3_extra", "")
    log:
        "logs/{treatment}.macs3.{peak_dir}.{peak_ext}.log"
    threads:
        workflow_threads("macs3", 1)
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p {params.outdir:q} logs
macs3 callpeak \
    -t {input.treatment:q} \
    {params.control_arg} \
    -f {params.fmt} \
    -g {params.genome_size:q} \
    -n {params.name:q} \
    --outdir {params.outdir:q} \
    {params.broad_arg} \
    {params.extra} \
    &> {log:q}
test -s {output:q}
        """


rule replicate_intersect:
    input:
        peaks=replicate_peak_inputs
    output:
        "replicate_intersect/{group}_intersect.bed"
    params:
        min_support=replicate_min_support
    log:
        "logs/{group}.replicate_intersect.log"
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p replicate_intersect logs

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

{{
    echo "Input peak files:"
    printf '%s\n' {input.peaks:q}
    echo "Minimum replicate support: {params.min_support}"

    i=0
    for f in {input.peaks:q}; do
        i=$((i + 1))
        awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' "$f" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > "$tmpdir/rep_${{i}}.bed"
        test -s "$tmpdir/rep_${{i}}.bed"
    done

    bedtools multiinter -i "$tmpdir"/rep_*.bed \
        | awk -v k={params.min_support} 'BEGIN{{OFS="\t"}} $4 >= k {{print $1, $2, $3, "support_"$4, $4, "."}}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - -c 5 -o max \
        | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "consensus_"NR, $4, "."}}' \
        > {output:q}

    test -s {output:q}
    echo "Consensus peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
