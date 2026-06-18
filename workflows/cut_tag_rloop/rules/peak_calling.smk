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


rule rnaseh_sensitive_regions:
    input:
        treatment_peaks=PATHS.peak("{sample}"),
        treatment_bw=rnaseh_sensitive_treatment_bigwig,
        rnaseh_bw=rnaseh_sensitive_control_bigwig
    output:
        table=PATHS.rnaseh_sensitive_signal_table("{scale_method}", "{sample}"),
        bed=PATHS.rnaseh_sensitive_regions("{scale_method}", "{sample}"),
        summary=PATHS.rnaseh_sensitive_summary("{scale_method}", "{sample}")
    wildcard_constraints:
        sample=TREATMENT_PATTERN,
        scale_method=SCALE_METHOD_PATTERN
    params:
        script="scripts/rnaseh_signal.py",
        rnaseh_sample=lambda wildcards: CTX.rnaseh_controls[wildcards.sample],
        min_fc=lambda wildcards: CTX.rnaseh_signal_min_fold_change,
        min_treatment_signal=lambda wildcards: CTX.rnaseh_signal_min_treatment_signal
    log:
        PATHS.log("{sample}.{scale_method}.rnaseh_sensitive")
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
    --scale-method {wildcards.scale_method:q} \
    &> {log:q}
test -s {output.table:q}
test -s {output.summary:q}
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


rule rnaseh_sensitive_consensus:
    input:
        peaks=replicate_sensitive_inputs
    output:
        PATHS.rnaseh_sensitive_consensus_peak("{scale_method}", "{group}")
    wildcard_constraints:
        group=GROUP_PATTERN,
        scale_method=SCALE_METHOD_PATTERN
    log:
        PATHS.log("{group}.{scale_method}.rnaseh_sensitive_consensus")
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
tmpdir=$(mktemp -d "$tmpbase/{wildcards.group}.{wildcards.scale_method}.XXXXXX")
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
        : > {output:q}
        echo "Fewer than two RNaseH-sensitive replicate files for {wildcards.group}; wrote empty consensus."
        exit 0
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
            if (s >= min_support) print $1, $2, $3, "rnaseh_sensitive_consensus_"++n, s, "."
        }}' "$tmpdir/support.tsv" "$tmpdir/candidates.bed" > {output:q}
    echo "RNaseH-sensitive consensus regions:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
