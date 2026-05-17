rule all:
    default_target: True
    input:
        f"{OUTDIR}/matrix/gene_counts.tsv",
        f"{OUTDIR}/matrix/gene_fpkm.tsv",
        f"{OUTDIR}/matrix/gene_tpm.tsv",
        f"{OUTDIR}/matrix/gene_length.tsv",
        expand(f"{OUTDIR}/bam/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        expand(f"{OUTDIR}/bigwig/{{sample}}.bw", sample=SAMPLES),
        f"{OUTDIR}/report/multiqc_report.html"
