rule gene_length:
    input:
        gtf=config["genome"]["gtf"]
    output:
        temp(f"{OUTDIR}/matrix/gene_length.tsv")
    log:
        f"{OUTDIR}/logs/matrix/gene_length.log"
    conda:
        ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/matrix {OUTDIR}/logs/matrix
        python scripts/build_gene_length_table.py \
          --gtf {input.gtf} \
          --output {output} \
          &> {log}
        test -s {output}
        """

rule expression_matrices:
    input:
        counts=expand(f"{OUTDIR}/counts/{{sample}}.featureCounts.txt", sample=SAMPLES),
        gene_length=f"{OUTDIR}/matrix/gene_length.tsv"
    output:
        raw=f"{OUTDIR}/matrix/gene_counts.tsv",
        fpkm=f"{OUTDIR}/matrix/gene_fpkm.tsv",
        tpm=f"{OUTDIR}/matrix/gene_tpm.tsv"
    params:
        sample_sheet=SAMPLES_TSV
    log:
        f"{OUTDIR}/logs/matrix/expression_matrices.log"
    conda:
        ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/matrix {OUTDIR}/logs/matrix
        python scripts/counts_to_fpkm_tpm.py \
          --featurecounts {input.counts} \
          --gene-length {input.gene_length} \
          --sample-sheet {params.sample_sheet} \
          --out-counts {output.raw} \
          --out-fpkm {output.fpkm} \
          --out-tpm {output.tpm} \
          &> {log}
        test -s {output.raw}
        test -s {output.fpkm}
        test -s {output.tpm}
        """
