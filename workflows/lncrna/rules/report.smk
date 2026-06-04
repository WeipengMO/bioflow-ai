rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
        expand(f"{OUTDIR}/qc/star/{{sample}}.Log.final.out", sample=SAMPLES),
        expand(f"{OUTDIR}/counts/{{sample}}.featureCounts.txt.summary", sample=SAMPLES),
        expand(f"{OUTDIR}/lncrna_counts/{{sample}}.featureCounts.txt.summary", sample=SAMPLES),
        f"{OUTDIR}/gffcompare/merged.stats",
        f"{OUTDIR}/lncrna/candidate_lncRNA.summary.tsv",
        f"{OUTDIR}/lncrna/high_confidence_lncRNA.summary.tsv",
        f"{OUTDIR}/lncrna/coding_filter/orf_metrics.tsv"
    output:
        html=f"{OUTDIR}/report/multiqc_report.html"
    log:
        f"{OUTDIR}/logs/multiqc/multiqc.log"
    conda:
        ENV_RNASEQ
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/report {OUTDIR}/logs/multiqc
        multiqc {OUTDIR} -o {OUTDIR}/report -n multiqc_report.html --force &> {log}
        test -s {output.html}
        """
