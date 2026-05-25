rule macs3_callpeak:
    input:
        bam=PATHS.filtered_bam("{sample}"),
        bai=PATHS.filtered_bai("{sample}"),
        control=macs3_control_input
    output:
        PATHS.peak("{sample}")
    wildcard_constraints:
        sample=PEAK_SAMPLE_PATTERN
    params:
        outdir=lambda wildcards: f"{PATHS.macs3_results}/{CTX.peak_type}",
        name=lambda wildcards: wildcards.sample,
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
    --outdir {params.outdir:q} \
    {params.extra} \
    &> {log:q}
test -s {output:q}
        """


if ENABLE_RNASEH_SUBTRACTION:
    rule rnaseh_subtract_peaks:
        input:
            treatment=PATHS.peak("{sample}"),
            rnaseh=rnaseh_control_peak
        output:
            PATHS.rnaseh_sensitive_peak("{sample}")
        wildcard_constraints:
            sample=TREATMENT_WITH_RNASEH_PATTERN
        params:
            blacklist=lambda wildcards: str(config.get("blacklist", "") or "")
        log:
            PATHS.log("{sample}.rnaseh_subtract")
        threads:
            workflow_threads("rnaseh_subtract", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})

tmpbase="$(dirname {output:q})/.tmp"
mkdir -p "$tmpbase"
tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print}}' {input.treatment:q} \
    | sort -k1,1 -k2,2n \
    > "$tmpdir/treatment.sorted.bed"
awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' {input.rnaseh:q} \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$tmpdir/rnaseh.merged.bed"

if [[ -s "$tmpdir/rnaseh.merged.bed" ]]; then
    bedtools intersect -v -a "$tmpdir/treatment.sorted.bed" -b "$tmpdir/rnaseh.merged.bed" > "$tmpdir/rnaseh_sensitive.bed" 2> {log:q}
else
    cp "$tmpdir/treatment.sorted.bed" "$tmpdir/rnaseh_sensitive.bed"
    echo "RNase H control peak file is empty after BED cleanup; copied treatment peaks." > {log:q}
fi

if [[ -n {params.blacklist:q} ]]; then
    bedtools intersect -v -a "$tmpdir/rnaseh_sensitive.bed" -b {params.blacklist:q} > {output:q} 2>> {log:q}
else
    cp "$tmpdir/rnaseh_sensitive.bed" {output:q}
fi

test -s {output:q}
            """


rule replicate_consensus:
    input:
        peaks=replicate_peak_inputs
    output:
        PATHS.replicate_consensus_peak("{group}")
    wildcard_constraints:
        group=GROUP_PATTERN
    log:
        PATHS.log("{group}.replicate_consensus")
    params:
        min_support=consensus_min_support
    threads:
        workflow_threads("replicate_consensus", 2)
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

    nrep="$i"
    if [ "$nrep" -lt 2 ]; then
        echo "At least two replicate peak files are required." >&2
        exit 1
    fi

    cat "$tmpdir"/rep_*.bed \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "candidate_"NR}}' \
        > "$tmpdir/candidates.bed"
    test -s "$tmpdir/candidates.bed"

    : > "$tmpdir/support_ids.txt"
    for rep in "$tmpdir"/rep_*.bed; do
        bedtools intersect -u -a "$tmpdir/candidates.bed" -b "$rep" \
            | cut -f4 \
            >> "$tmpdir/support_ids.txt"
    done

    awk 'BEGIN{{OFS="\t"}} {{support[$1]++}} END{{for (id in support) print id, support[id]}}' \
        "$tmpdir/support_ids.txt" \
        > "$tmpdir/support.tsv"

    awk -v min_support="{params.min_support}" 'BEGIN{{OFS="\t"}}
        FNR==NR {{support[$1]=$2; next}}
        {{
            s = support[$4] + 0
            if (s >= min_support) print $1, $2, $3, "consensus_"++n, s, "."
        }}' "$tmpdir/support.tsv" "$tmpdir/candidates.bed" \
        > {output:q}

    test -s {output:q}
    echo "Candidate peak universe:"
    wc -l "$tmpdir/candidates.bed"
    echo "Consensus peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
