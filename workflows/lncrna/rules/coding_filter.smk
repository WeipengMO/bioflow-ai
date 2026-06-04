rule high_confidence_lncrna:
    input:
        gtf=f"{OUTDIR}/lncrna/candidate_lncRNA.gtf",
        fasta=f"{OUTDIR}/lncrna/candidate_lncRNA.fa",
        summary=f"{OUTDIR}/lncrna/candidate_lncRNA.summary.tsv"
    output:
        gtf=f"{OUTDIR}/lncrna/high_confidence_lncRNA.gtf",
        fasta=f"{OUTDIR}/lncrna/high_confidence_lncRNA.fa",
        summary=f"{OUTDIR}/lncrna/high_confidence_lncRNA.summary.tsv",
        orf=f"{OUTDIR}/lncrna/coding_filter/orf_metrics.tsv"
    params:
        max_orf_aa=config.get("coding_filter", {}).get("max_orf_aa", 100),
        exclude=config.get("coding_filter", {}).get("exclude_transcripts", "")
    log:
        f"{OUTDIR}/logs/lncrna/high_confidence_lncRNA.log"
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/lncrna/coding_filter {OUTDIR}/logs/lncrna
        EXCLUDE_ARG=""
        if [ -n "{params.exclude}" ] && [ "{params.exclude}" != "None" ] && [ "{params.exclude}" != "null" ]; then
          EXCLUDE_ARG="--exclude-transcripts {params.exclude}"
        fi
        python scripts/filter_high_confidence_lncRNA.py \
          --candidate-gtf {input.gtf} \
          --candidate-fasta {input.fasta} \
          --candidate-summary {input.summary} \
          --output-gtf {output.gtf} \
          --output-fasta {output.fasta} \
          --output-summary {output.summary} \
          --orf-metrics {output.orf} \
          --max-orf-aa {params.max_orf_aa} \
          $EXCLUDE_ARG \
          &> {log}
        test -s {output.gtf}
        test -s {output.fasta}
        test -s {output.summary}
        test -s {output.orf}
        """
