rule macs3_callpeak:
    input:
        bam=peak_input_bam,
        bai=peak_input_bai,
        control=macs3_control_input
    output:
        PATHS.peak("{sample}")
    wildcard_constraints:
        sample=SAMPLE_PATTERN
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
    --keep-dup all \
    --outdir {params.outdir:q} \
    {params.extra} \
    &> {log:q}
test -s {output:q}
        """


if ENABLE_RNASEH_SUBTRACTION:
    rule rnaseh_no_overlap_peaks:
        input:
            treatment=PATHS.peak("{sample}"),
            rnaseh=rnaseh_control_peak
        output:
            no_overlap=PATHS.rnaseh_no_overlap_peak("{sample}")
        wildcard_constraints:
            sample=TREATMENT_PATTERN
        params:
            blacklist=lambda wildcards: str(config.get("blacklist", "") or "")
        log:
            PATHS.log("{sample}.rnaseh_no_overlap")
        threads:
            workflow_threads("rnaseh_subtract", 2)
        conda:
            ENV
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output.no_overlap:q}) $(dirname {log:q})

tmpbase="$(dirname {output.no_overlap:q})/.tmp"
mkdir -p "$tmpbase"
tmpdir=$(mktemp -d "$tmpbase/{wildcards.sample}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

{{
    echo "RNase H no-overlap filtering for {wildcards.sample}."
    echo "Method: bedtools intersect -v against merged RNase H control peaks."
    echo "Interpretation: this is overlap filtering only, not a quantitative RNaseH-depletion test."
}} > {log:q}

awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print}}' {input.treatment:q} \
    | sort -k1,1 -k2,2n \
    > "$tmpdir/treatment.sorted.bed"
awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $1, $2, $3}}' {input.rnaseh:q} \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$tmpdir/rnaseh.merged.bed"

if [[ -s "$tmpdir/rnaseh.merged.bed" ]]; then
    bedtools intersect -v -a "$tmpdir/treatment.sorted.bed" -b "$tmpdir/rnaseh.merged.bed" > "$tmpdir/rnaseh_no_overlap.bed" 2>> {log:q}
else
    cp "$tmpdir/treatment.sorted.bed" "$tmpdir/rnaseh_no_overlap.bed"
    echo "RNase H control peak file is empty after BED cleanup; copied treatment peaks." >> {log:q}
fi

blacklist={params.blacklist:q}
if [[ -n "$blacklist" ]]; then
    bedtools intersect -v -a "$tmpdir/rnaseh_no_overlap.bed" -b "$blacklist" > {output.no_overlap:q} 2>> {log:q}
else
    cp "$tmpdir/rnaseh_no_overlap.bed" {output.no_overlap:q}
fi

test -e {output.no_overlap:q}
            """

    if CTX.write_deprecated_rnaseh_sensitive_alias:
        rule rnaseh_sensitive_alias:
            input:
                PATHS.rnaseh_no_overlap_peak("{sample}")
            output:
                PATHS.rnaseh_sensitive_peak("{sample}")
            wildcard_constraints:
                sample=TREATMENT_PATTERN
            log:
                PATHS.log("{sample}.rnaseh_sensitive_alias")
            conda:
                ENV
            shell:
                r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})
cp {input:q} {output:q}
echo "Deprecated compatibility alias for RNase H no-overlap peaks; not quantitative RNaseH sensitivity." > {log:q}
test -e {output:q}
                """

    rule rnaseh_signal:
        input:
            treatment_peaks=PATHS.peak("{sample}"),
            no_overlap=PATHS.rnaseh_no_overlap_peak("{sample}"),
            treatment_bw=rnaseh_signal_treatment_bigwig,
            rnaseh_bw=rnaseh_signal_control_bigwig
        output:
            table=PATHS.rnaseh_signal_table("{sample}"),
            bed=PATHS.rnaseh_depleted_peak("{sample}"),
            summary=PATHS.rnaseh_sensitivity_summary("{sample}")
        wildcard_constraints:
            sample=TREATMENT_PATTERN
        params:
            script="scripts/rnaseh_signal.py",
            rnaseh_sample=lambda wildcards: CTX.rnaseh_controls[wildcards.sample],
            min_fc=lambda wildcards: CTX.rnaseh_signal_min_fold_change,
            min_treatment_signal=lambda wildcards: CTX.rnaseh_signal_min_treatment_signal,
            signal_track_type=lambda wildcards: rnaseh_signal_track_type()
        log:
            PATHS.log("{sample}.rnaseh_signal")
        threads:
            workflow_threads("rnaseh_signal", 2)
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
    --no-overlap-peaks {input.no_overlap:q} \
    --treatment-bw {input.treatment_bw:q} \
    --rnaseh-bw {input.rnaseh_bw:q} \
    --output-tsv {output.table:q} \
    --output-bed {output.bed:q} \
    --summary {output.summary:q} \
    --min-fold-change {params.min_fc} \
    --min-treatment-signal {params.min_treatment_signal} \
    --signal-track-type {params.signal_track_type:q} \
    &> {log:q}
test -s {output.table:q}
test -s {output.summary:q}
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

    if [[ ! -s "$tmpdir/candidates.bed" ]]; then
        : > {output:q}
        echo "No candidate peaks remained after replicate peak cleanup."
        exit 0
    fi

    : > "$tmpdir/support_ids.txt"
    for rep in "$tmpdir"/rep_*.bed; do
        if [[ -s "$rep" ]]; then
            bedtools intersect -u -a "$tmpdir/candidates.bed" -b "$rep" \
                | cut -f4 \
                >> "$tmpdir/support_ids.txt"
        fi
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

    echo "Candidate peak universe:"
    wc -l "$tmpdir/candidates.bed"
    echo "Consensus peaks:"
    wc -l {output:q}
}} > {log:q} 2>&1
        """
