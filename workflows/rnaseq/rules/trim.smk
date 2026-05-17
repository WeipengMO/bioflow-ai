ruleorder: fastp_pe > fastp_se

rule fastp_pe:
    input:
        fq1=raw_fq1,
        fq2=raw_fq2
    output:
        fq1=temp(f"{OUTDIR}/clean_fastq/{{sample}}.R1.clean.fastq.gz"),
        fq2=temp(f"{OUTDIR}/clean_fastq/{{sample}}.R2.clean.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json"
    log:
        f"{OUTDIR}/logs/fastp/{{sample}}.log"
    threads:
        config.get("fastp", {}).get("threads", 4)
    conda:
        ENV
    params:
        q=config.get("fastp", {}).get("qualified_quality_phred", 20),
        length=config.get("fastp", {}).get("length_required", 30),
        extra=config.get("fastp", {}).get("extra", "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/clean_fastq {OUTDIR}/qc/fastp {OUTDIR}/logs/fastp
        fastp \
          -i {input.fq1} -I {input.fq2} \
          -o {output.fq1} -O {output.fq2} \
          --thread {threads} \
          --qualified_quality_phred {params.q} \
          --length_required {params.length} \
          --detect_adapter_for_pe \
          --html {output.html} \
          --json {output.json} \
          {params.extra} \
          &> {log}
        test -s {output.fq1}
        test -s {output.fq2}
        """

rule fastp_se:
    input:
        fq1=raw_fq1
    output:
        fq1=temp(f"{OUTDIR}/clean_fastq/{{sample}}.clean.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}.fastp.json"
    log:
        f"{OUTDIR}/logs/fastp/{{sample}}.log"
    threads:
        config.get("fastp", {}).get("threads", 4)
    conda:
        ENV
    params:
        q=config.get("fastp", {}).get("qualified_quality_phred", 20),
        length=config.get("fastp", {}).get("length_required", 30),
        extra=config.get("fastp", {}).get("extra", "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/clean_fastq {OUTDIR}/qc/fastp {OUTDIR}/logs/fastp
        fastp \
          -i {input.fq1} \
          -o {output.fq1} \
          --thread {threads} \
          --qualified_quality_phred {params.q} \
          --length_required {params.length} \
          --html {output.html} \
          --json {output.json} \
          {params.extra} \
          &> {log}
        test -s {output.fq1}
        """
