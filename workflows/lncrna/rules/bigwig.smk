rule bam_to_bigwig:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai"
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}.bw"
    params:
        bin_size=config.get("bigwig", {}).get("bin_size", 10),
        normalize_using=config.get("bigwig", {}).get("normalize_using", "CPM"),
        effective_genome_size=config.get("bigwig", {}).get("effective_genome_size", None),
        extend_reads=config.get("bigwig", {}).get("extend_reads", False),
        extra=config.get("bigwig", {}).get("extra", "")
    log:
        f"{OUTDIR}/logs/bigwig/{{sample}}.log"
    threads:
        config.get("bigwig", {}).get("threads", 4)
    conda:
        ENV_RNASEQ
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/bigwig {OUTDIR}/logs/bigwig
        EGS_ARG=""
        if [ "{params.effective_genome_size}" != "None" ] && [ "{params.effective_genome_size}" != "null" ]; then
          EGS_ARG="--effectiveGenomeSize {params.effective_genome_size}"
        fi
        EXTEND_ARG=""
        if [ "{params.extend_reads}" = "True" ] || [ "{params.extend_reads}" = "true" ]; then
          EXTEND_ARG="--extendReads"
        fi
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize {params.bin_size} \
          --normalizeUsing {params.normalize_using} \
          -p {threads} \
          $EGS_ARG \
          $EXTEND_ARG \
          {params.extra} \
          &> {log}
        test -s {output.bw}
        """
