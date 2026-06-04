rule featurecounts_lncrna_sample:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai",
        gtf=f"{OUTDIR}/lncrna/high_confidence_lncRNA.gtf"
    output:
        counts=f"{OUTDIR}/lncrna_counts/{{sample}}.featureCounts.txt",
        summary=f"{OUTDIR}/lncrna_counts/{{sample}}.featureCounts.txt.summary"
    params:
        strand=strand_flag_featurecounts(),
        feature_type=config.get("featurecounts", {}).get("feature_type", "exon"),
        attribute_type=config.get("featurecounts", {}).get("attribute_type", "gene_id"),
        extra=featurecounts_extra
    log:
        f"{OUTDIR}/logs/featurecounts_lncrna/{{sample}}.log"
    threads:
        config.get("featurecounts", {}).get("threads", 8)
    conda:
        ENV_RNASEQ
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/lncrna_counts {OUTDIR}/logs/featurecounts_lncrna
        featureCounts \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          -t {params.feature_type} \
          -g {params.attribute_type} \
          -s {params.strand} \
          {params.extra} \
          {input.bam} \
          &> {log}
        test -s {output.counts}
        test -s {output.summary}
        """

rule featurecounts_gene_sample:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai"
    output:
        counts=f"{OUTDIR}/counts/{{sample}}.featureCounts.txt",
        summary=f"{OUTDIR}/counts/{{sample}}.featureCounts.txt.summary"
    params:
        gtf=config["genome"]["gtf"],
        strand=strand_flag_featurecounts(),
        feature_type=config.get("featurecounts", {}).get("feature_type", "exon"),
        attribute_type=config.get("featurecounts", {}).get("attribute_type", "gene_id"),
        extra=featurecounts_extra
    log:
        f"{OUTDIR}/logs/featurecounts/{{sample}}.log"
    threads:
        config.get("featurecounts", {}).get("threads", 8)
    conda:
        ENV_RNASEQ
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/counts {OUTDIR}/logs/featurecounts
        featureCounts \
          -T {threads} \
          -a {params.gtf} \
          -o {output.counts} \
          -t {params.feature_type} \
          -g {params.attribute_type} \
          -s {params.strand} \
          {params.extra} \
          {input.bam} \
          &> {log}
        test -s {output.counts}
        test -s {output.summary}
        """
