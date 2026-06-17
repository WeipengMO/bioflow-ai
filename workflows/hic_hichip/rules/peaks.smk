rule link_hicpro_peak_bam:
    input:
        done=PATHS.hicpro_done()
    output:
        bam=PATHS.hicpro_mapped_bam("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    params:
        assembly=lambda wildcards: (config.get("genome", {}) or {}).get("assembly", "")
    log:
        PATHS.log("peaks/link_hicpro_bam/{sample}")
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
base={PATHS.hicpro:q}
sample={wildcards.sample:q}
assembly={params.assembly:q}
mapfile -t candidates < <(
    find "$base" -path "*/bowtie_results/bwt2/${{sample}}/${{sample}}__*.bwt2pairs.bam" -type f | sort
)
if [[ ${{#candidates[@]}} -eq 0 && -n "$assembly" ]]; then
    mapfile -t candidates < <(
        find "$base" -path "*/bowtie_results/bwt2/${{sample}}/${{sample}}__${{assembly}}.bwt2pairs.bam" -type f | sort
    )
fi
if [[ ${{#candidates[@]}} -eq 0 ]]; then
    echo "ERROR: no HiC-Pro paired BAM found for $sample under $base" | tee {log:q}
    exit 1
fi
if [[ ${{#candidates[@]}} -gt 1 ]]; then
    echo "ERROR: multiple HiC-Pro paired BAM candidates found for $sample:" | tee {log:q}
    printf '%s\n' "${{candidates[@]}}" | tee -a {log:q}
    exit 1
fi
ln -sf $(realpath "${{candidates[0]}}") {output.bam:q}
test -s {output.bam:q}
        """

rule call_hichip_peaks:
    input:
        bam=PATHS.hicpro_mapped_bam("{sample}")
    output:
        narrowpeak=PATHS.auto_peak_narrowpeak("{sample}"),
        bed=PATHS.auto_peak_bed("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    params:
        caller=lambda wildcards: CTX.peak_caller,
        genome_size=lambda wildcards: CTX.peak_genome_size,
        qvalue=lambda wildcards: CTX.peak_qvalue,
        mode=lambda wildcards: CTX.peak_mode_macs,
        extra=lambda wildcards: CTX.peak_extra,
        allow_empty=lambda wildcards: "true" if CTX.peak_allow_empty else "false",
        outdir=lambda wildcards: f"{PATHS.peaks}/{wildcards.sample}",
        prefix=lambda wildcards: f"{wildcards.sample}.auto_peaks",
    log:
        PATHS.log("peaks/call_hichip_peaks/{sample}")
    conda:
        PEAK_ENV
    shell:
        r"""
set -euo pipefail
mkdir -p {params.outdir:q} $(dirname {log:q})
{params.caller:q} callpeak \
    -t {input.bam:q} \
    -f {params.mode:q} \
    -g {params.genome_size:q} \
    -n {params.prefix:q} \
    --outdir {params.outdir:q} \
    -q {params.qvalue:q} \
    {params.extra} \
    &> {log:q}
mapfile -t peaks < <(find {params.outdir:q} -maxdepth 1 -name "{params.prefix}_peaks.narrowPeak" -type f | sort)
if [[ ${{#peaks[@]}} -eq 0 ]]; then
    if [[ {params.allow_empty:q} == "true" ]]; then
        echo "WARNING: no MACS narrowPeak produced for {wildcards.sample}; writing empty auto peak outputs." | tee -a {log:q}
        : > {output.narrowpeak:q}
        : > {output.bed:q}
        exit 0
    fi
    echo "ERROR: MACS did not produce {params.prefix}_peaks.narrowPeak in {params.outdir}" | tee -a {log:q}
    exit 1
fi
if [[ ${{#peaks[@]}} -gt 1 ]]; then
    echo "ERROR: multiple narrowPeak candidates found:" | tee -a {log:q}
    printf '%s\n' "${{peaks[@]}}" | tee -a {log:q}
    exit 1
fi
cp "${{peaks[0]}}" {output.narrowpeak:q}
cut -f1-6 {output.narrowpeak:q} > {output.bed:q}
if [[ {params.allow_empty:q} != "true" ]]; then
    test -s {output.narrowpeak:q}
    test -s {output.bed:q}
fi
        """

rule write_peak_sources:
    output:
        table=PATHS.peak_source_table()
    run:
        from pathlib import Path
        Path(output.table).parent.mkdir(parents=True, exist_ok=True)
        with open(output.table, "w") as handle:
            handle.write("sample\tpeak_source\tpeak_file\n")
            for sample in CTX.sample_names:
                handle.write(f"{sample}\t{sample_peak_source(sample)}\t{sample_peak_path(sample)}\n")
