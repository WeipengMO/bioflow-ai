if USE_POOLED_CONTROL:
    rule merge_pooled_control:
        input:
            pooled_control_inputs
        output:
            bam=analysis_bam("{pool}"),
            bai=analysis_bai("{pool}")
        wildcard_constraints:
            pool=POOLED_CONTROL_PATTERN
        log:
            PATHS.log("{pool}.merge_control")
        threads:
            workflow_threads("merge_control", 4)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
samtools merge -@ {threads} -f {output.bam:q} {input:q} &> {log:q}
samtools index -@ {threads} {output.bam:q} 2>> {log:q}
test -s {output.bam:q}
test -s {output.bai:q}
            """


rule macs3_callpeak:
    input:
        unpack(macs3_input)
    output:
        f"{PATHS.macs3_results}/{{peak_dir}}/{{treatment}}_peaks.{{peak_ext}}"
    wildcard_constraints:
        peak_dir="narrow|broad|narrow_no_control|broad_no_control",
        peak_ext="narrowPeak|broadPeak",
        treatment=TREATMENT_PATTERN
    params:
        outdir=lambda wildcards: f"{PATHS.macs3_results}/{wildcards.peak_dir}",
        name=lambda wildcards: wildcards.treatment,
        genome_size=lambda wildcards: config["gsize"],
        fmt=lambda wildcards: macs3_format(MODE),
        control_arg=macs3_control_arg,
        broad_arg=macs3_broad_arg,
        extra=lambda wildcards: config.get("macs3_extra", "")
    log:
        PATHS.log("{treatment}.macs3.{peak_dir}.{peak_ext}")
    threads:
        workflow_threads("macs3", 1)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p {params.outdir:q} $(dirname {log:q})
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
        PATHS.replicate_consensus("{group}")
    log:
        PATHS.log("{group}.replicate_intersect")
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
    echo "Input peak files:"
    printf '%s\n' {input.peaks:q}

    i=0
    for f in {input.peaks:q}; do
        i=$((i + 1))
        awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' "$f" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > "$tmpdir/rep_${{i}}.bed"
        test -s "$tmpdir/rep_${{i}}.bed"
    done

    if [ "$i" -lt 2 ]; then
        echo "At least two replicate peak files are required." >&2
        exit 1
    fi
    echo "Replicate support required: $i of $i"

    bedtools multiinter -i "$tmpdir"/rep_*.bed \
        | awk -v k="$i" 'BEGIN{{OFS="\t"}} $4 >= k {{print $1, $2, $3, "support_"$4, $4, "."}}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - -c 5 -o max \
        | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "consensus_"NR, $4, "."}}' \
        > {output:q}

    test -s {output:q}
    echo "Consensus peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
