rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
        expand(f"{OUTDIR}/qc/star/{{sample}}.Log.final.out", sample=SAMPLES),
        expand(f"{OUTDIR}/counts/{{sample}}.featureCounts.txt.summary", sample=SAMPLES)
    output:
        html=f"{OUTDIR}/report/multiqc_report.html"
    log:
        f"{OUTDIR}/logs/multiqc/multiqc.log"
    conda:
        ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/report {OUTDIR}/logs/multiqc
        multiqc {OUTDIR} -o {OUTDIR}/report -n multiqc_report.html &> {log}
        test -s {output.html}
        """
