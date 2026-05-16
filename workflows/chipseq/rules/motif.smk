HOMER_REPORT_DIR = config.get("homer_report_dir", "reports/homer")
HOMER_INPUT_DIR = config.get("homer_input_dir", "homer_inputs")
HOMER_GENOME = config.get("homer_genome", "")
HOMER_SIZE = config.get("homer_size", 200 if HOMER_PEAK_TYPE == "narrow" else "given")
HOMER_EXTRA = config.get("homer_extra", "-mask")
HOMER_USE_SUMMIT = bool(config.get("homer_use_summit", HOMER_PEAK_TYPE == "narrow"))
HOMER_BLACKLIST = config.get("homer_blacklist", "")
HOMER_MIN_PEAKS = int(config.get("homer_min_peaks", 50))


def homer_raw_peak_input(wildcards):
    if HOMER_INPUT_SOURCE == "replicate_intersect":
        return f"replicate_intersect/{wildcards.target}_intersect.bed"
    return peak_path(wildcards.target, HOMER_PEAK_MODE, HOMER_PEAK_TYPE)


rule homer_prepare_peaks:
    input:
        peaks=homer_raw_peak_input
    output:
        temp(f"{HOMER_INPUT_DIR}/{HOMER_RESULT_TAG}/{{target}}.bed")
    params:
        use_summit=HOMER_USE_SUMMIT,
        peak_type=HOMER_PEAK_TYPE,
        blacklist=HOMER_BLACKLIST
    log:
        f"logs/{{target}}.homer.prepare.{HOMER_RESULT_TAG}.log"
    conda:
        "../envs/chipseq.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) logs

tmp1=$(mktemp)
tmp2=$(mktemp)
trap 'rm -f "$tmp1" "$tmp2"' EXIT

# If the input is a narrowPeak file with a valid summit column, use the MACS summit.
# Otherwise, use the interval center so consensus BED files are also supported.
if [[ "{params.peak_type}" == "narrow" && "{params.use_summit}" == "True" ]]; then
    awk 'BEGIN{{OFS="\t"}}
        NF >= 10 && $10 >= 0 {{
            s = $2 + $10;
            if (s < 0) s = 0;
            print $1, s, s + 1, $4, $5, $6, $7, $8, $9, $10;
            next;
        }}
        NF >= 3 && $2 < $3 {{
            s = int(($2 + $3) / 2);
            if (s < 0) s = 0;
            name = (NF >= 4 ? $4 : ".");
            score = (NF >= 5 ? $5 : 0);
            strand = (NF >= 6 ? $6 : ".");
            print $1, s, s + 1, name, score, strand;
        }}
    ' {input.peaks:q} > "$tmp1"
else
    awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $0}}' {input.peaks:q} > "$tmp1"
fi

blacklist={params.blacklist:q}
if [[ -n "$blacklist" ]]; then
    bedtools intersect -v -a "$tmp1" -b "$blacklist" > "$tmp2"
else
    cp "$tmp1" "$tmp2"
fi

sort -k1,1 -k2,2n "$tmp2" > {output:q}

before_count=$(awk 'NF >= 3 && $2 < $3' "$tmp1" | wc -l)
after_count=$(awk 'NF >= 3 && $2 < $3' {output:q} | wc -l)

{{
    echo "Input peak file: {input.peaks}"
    echo "Output peak file: {output}"
    echo "Blacklist: {params.blacklist}"
    echo "Peaks before blacklist: $before_count"
    echo "Peaks after blacklist: $after_count"
}} > {log:q}

test -f {output:q}
        """


rule homer_find_motifs:
    input:
        peaks=rules.homer_prepare_peaks.output
    output:
        touch(f"{HOMER_REPORT_DIR}/{HOMER_RESULT_TAG}/{{target}}/.homer_complete")
    params:
        genome=HOMER_GENOME,
        outdir=f"{HOMER_REPORT_DIR}/{HOMER_RESULT_TAG}/{{target}}",
        size=HOMER_SIZE,
        extra=HOMER_EXTRA,
        min_peaks=HOMER_MIN_PEAKS
    log:
        f"logs/{{target}}.homer.{HOMER_RESULT_TAG}.log"
    threads:
        workflow_threads("homer", 8)
    conda:
        "../envs/homer.yml"
    shell:
        r"""
set -euo pipefail
rm -rf {params.outdir:q}
mkdir -p {params.outdir:q} logs

n_peaks=$(awk 'NF >= 3 && $2 < $3' {input.peaks:q} | wc -l)

{{
    echo "Input peaks for HOMER: $n_peaks"
    echo "Minimum peaks required: {params.min_peaks}"
}} > {log:q}

if [[ "$n_peaks" -lt "{params.min_peaks}" ]]; then
    echo "Skip HOMER: too few peaks after filtering." >> {log:q}
    echo "Skipped HOMER because only $n_peaks peaks remained after filtering." \
        > {params.outdir:q}/SKIPPED.too_few_peaks.txt
    exit 0
fi

findMotifsGenome.pl \
    {input.peaks:q} \
    {params.genome:q} \
    {params.outdir:q} \
    -size {params.size} \
    -p {threads} \
    {params.extra} \
    >> {log:q} 2>&1

if grep -q '^!!!!Genome .* not found' {log:q}; then
    exit 1
fi

test -s {params.outdir:q}/knownResults.txt
test -s {params.outdir:q}/knownResults.html
test -s {params.outdir:q}/homerResults.html
        """
