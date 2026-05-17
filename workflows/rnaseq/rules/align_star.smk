rule star_align:
    input:
        reads=clean_reads
    output:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai",
        log_final=f"{OUTDIR}/qc/star/{{sample}}.Log.final.out"
    params:
        index=config["genome"]["star_index"],
        prefix=lambda wc: f"{OUTDIR}/star/{wc.sample}/",
        extra=config.get("star", {}).get("extra", "")
    log:
        f"{OUTDIR}/logs/star/{{sample}}.log"
    threads:
        config.get("star", {}).get("threads", 12)
    conda:
        ENV
    shell:
        r"""
        set -euo pipefail
        ulimit -n 65535 || true
        mkdir -p {OUTDIR}/star/{wildcards.sample} {OUTDIR}/bam {OUTDIR}/qc/star {OUTDIR}/logs/star
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.reads} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.prefix} \
          --outSAMtype BAM SortedByCoordinate \
          {params.extra} \
          &> {log}
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index -@ {threads} {output.bam}
        cp {params.prefix}Log.final.out {output.log_final}
        test -s {output.bam}
        test -s {output.bai}
        test -s {output.log_final}
        """
