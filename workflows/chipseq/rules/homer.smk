HOMER_REPORT_DIR = str(config.get("homer_report_dir") or PATHS.default_homer_report)
HOMER_INPUT_DIR = str(config.get("homer_input_dir") or PATHS.default_homer_input)
HOMER_GENOME = config.get("homer_genome", "")
HOMER_GTF = str(config.get("gtf", "") or "")
HOMER_SIZE = config.get("homer_size", 200 if HOMER_PEAK_TYPE == "narrow" else "given")
HOMER_EXTRA = config.get("homer_extra", "-mask")
HOMER_ANNOTATE_EXTRA = config.get("homer_annotate_extra", "")
HOMER_USE_SUMMIT = config_bool("homer_use_summit", HOMER_PEAK_TYPE == "narrow")
BLACKLIST = config.get("blacklist", "")
HOMER_MIN_PEAKS = int(config.get("homer_min_peaks", 50))


rule homer_prepare_motif_peaks:
    input:
        peaks=homer_raw_peak_input
    output:
        temp(f"{HOMER_INPUT_DIR}/{HOMER_RESULT_TAG}/{{target}}.motif.bed")
    params:
        use_summit=HOMER_USE_SUMMIT,
        peak_type=HOMER_PEAK_TYPE,
        blacklist=BLACKLIST
    log:
        PATHS.log(f"{{target}}.homer.prepare_motif.{HOMER_RESULT_TAG}")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})

tmpbase="$(dirname {output:q})/.tmp"
mkdir -p "$tmpbase"
tmp1=$(mktemp "$tmpbase/{wildcards.target}.motif.tmp1.XXXXXX")
tmp2=$(mktemp "$tmpbase/{wildcards.target}.motif.tmp2.XXXXXX")
trap 'rm -f "$tmp1" "$tmp2"' EXIT

# Motif discovery can use MACS summits for narrowPeak files; consensus BED inputs fall back to centers.
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
    echo "Output motif peak file: {output}"
    echo "Use summit: {params.use_summit}"
    echo "Blacklist: {params.blacklist}"
    echo "Peaks before blacklist: $before_count"
    echo "Peaks after blacklist: $after_count"
}} > {log:q}

test -f {output:q}
        """


rule homer_prepare_annotation_peaks:
    input:
        peaks=homer_raw_peak_input
    output:
        temp(f"{HOMER_INPUT_DIR}/{HOMER_RESULT_TAG}/{{target}}.annotation.bed")
    params:
        blacklist=BLACKLIST
    log:
        PATHS.log(f"{{target}}.homer.prepare_annotation.{HOMER_RESULT_TAG}")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) $(dirname {log:q})

tmpbase="$(dirname {output:q})/.tmp"
mkdir -p "$tmpbase"
tmp1=$(mktemp "$tmpbase/{wildcards.target}.annotation.tmp1.XXXXXX")
tmp2=$(mktemp "$tmpbase/{wildcards.target}.annotation.tmp2.XXXXXX")
trap 'rm -f "$tmp1" "$tmp2"' EXIT

awk 'BEGIN{{OFS="\t"}} NF >= 3 && $2 < $3 {{print $0}}' {input.peaks:q} > "$tmp1"

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
    echo "Output annotation peak file: {output}"
    echo "Blacklist: {params.blacklist}"
    echo "Peaks before blacklist: $before_count"
    echo "Peaks after blacklist: $after_count"
}} > {log:q}

test -f {output:q}
        """


rule homer_find_motifs:
    input:
        peaks=rules.homer_prepare_motif_peaks.output
    output:
        touch(f"{HOMER_REPORT_DIR}/{HOMER_RESULT_TAG}/{{target}}/.motifs_complete")
    params:
        genome=HOMER_GENOME,
        outdir=f"{HOMER_REPORT_DIR}/{HOMER_RESULT_TAG}/{{target}}",
        size=HOMER_SIZE,
        extra=HOMER_EXTRA,
        min_peaks=HOMER_MIN_PEAKS
    log:
        PATHS.log(f"{{target}}.homer.{HOMER_RESULT_TAG}")
    threads:
        workflow_threads("homer", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p {params.outdir:q} $(dirname {log:q})

rm -rf \
    {params.outdir:q}/knownResults \
    {params.outdir:q}/homerResults \
    {params.outdir:q}/knownResults.txt \
    {params.outdir:q}/knownResults.html \
    {params.outdir:q}/homerResults.html \
    {params.outdir:q}/homerMotifs.all.motifs \
    {params.outdir:q}/motifFindingParameters.txt \
    {params.outdir:q}/SKIPPED.too_few_peaks.txt \
    {output:q}

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


rule homer_annotate_peaks:
    input:
        peaks=rules.homer_prepare_annotation_peaks.output
    output:
        annotation=f"{HOMER_REPORT_DIR}/{HOMER_RESULT_TAG}/{{target}}/annotatePeaks.txt"
    params:
        genome=HOMER_GENOME,
        gtf=HOMER_GTF,
        extra=HOMER_ANNOTATE_EXTRA
    log:
        PATHS.log(f"{{target}}.homer.annotate.{HOMER_RESULT_TAG}")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.annotation:q}) $(dirname {log:q})
rm -f {output.annotation:q}

annotatePeaks.pl \
    {input.peaks:q} \
    {params.genome:q} \
    -gtf {params.gtf:q} \
    {params.extra} \
    > {output.annotation:q} \
    2> {log:q}

test -s {output.annotation:q}
        """
