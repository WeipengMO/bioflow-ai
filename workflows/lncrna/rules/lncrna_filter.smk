rule filter_candidate_lncrna:
    input:
        annotated=f"{OUTDIR}/gffcompare/merged.annotated.gtf",
        fasta=config["genome"]["fasta"]
    output:
        gtf=f"{OUTDIR}/lncrna/candidate_lncRNA.gtf",
        summary=f"{OUTDIR}/lncrna/candidate_lncRNA.summary.tsv"
    params:
        min_length=config.get("lncrna_filter", {}).get("min_length", 200),
        min_exons=config.get("lncrna_filter", {}).get("min_exons", 2),
        class_codes=",".join(config.get("lncrna_filter", {}).get("class_codes", ["u", "i", "x", "o", "e"]))
    log:
        f"{OUTDIR}/logs/lncrna/filter_candidate_lncRNA.log"
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/lncrna {OUTDIR}/logs/lncrna
        python scripts/filter_lncRNA_candidates.py \
          --gtf {input.annotated} \
          --output-gtf {output.gtf} \
          --summary {output.summary} \
          --min-length {params.min_length} \
          --min-exons {params.min_exons} \
          --class-codes {params.class_codes} \
          --fasta {input.fasta} \
          &> {log}
        test -s {output.gtf}
        test -s {output.summary}
        """

rule candidate_lncrna_fasta:
    input:
        gtf=f"{OUTDIR}/lncrna/candidate_lncRNA.gtf",
        fasta=config["genome"]["fasta"]
    output:
        fasta=f"{OUTDIR}/lncrna/candidate_lncRNA.fa"
    log:
        f"{OUTDIR}/logs/lncrna/candidate_lncRNA_fasta.log"
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/lncrna {OUTDIR}/logs/lncrna
        gffread {input.gtf} \
          -g {input.fasta} \
          -w {output.fasta} \
          &> {log}
        test -s {output.fasta}
        """
