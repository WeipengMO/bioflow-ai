rule stringtie_assemble:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai"
    output:
        gtf=f"{OUTDIR}/stringtie/{{sample}}.gtf"
    params:
        ref=config["genome"]["gtf"],
        extra=config.get("stringtie", {}).get("assemble_extra", "")
    log:
        f"{OUTDIR}/logs/stringtie/{{sample}}.assemble.log"
    threads:
        config.get("stringtie", {}).get("threads", 8)
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/stringtie {OUTDIR}/logs/stringtie
        stringtie {input.bam} \
          -p {threads} \
          -G {params.ref} \
          -o {output.gtf} \
          {params.extra} \
          &> {log}
        test -s {output.gtf}
        """

rule stringtie_merge_list:
    input:
        expand(f"{OUTDIR}/stringtie/{{sample}}.gtf", sample=SAMPLES)
    output:
        list=f"{OUTDIR}/stringtie/mergelist.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/stringtie
        printf "%s\n" {input} > {output.list}
        test -s {output.list}
        """

rule stringtie_merge:
    input:
        list=f"{OUTDIR}/stringtie/mergelist.txt"
    output:
        gtf=f"{OUTDIR}/assembly/stringtie_merged.gtf"
    params:
        ref=config["genome"]["gtf"],
        extra=config.get("stringtie", {}).get("merge_extra", "")
    log:
        f"{OUTDIR}/logs/stringtie/merge.log"
    threads:
        config.get("stringtie", {}).get("threads", 8)
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/assembly {OUTDIR}/logs/stringtie
        stringtie --merge \
          -p {threads} \
          -G {params.ref} \
          -o {output.gtf} \
          {params.extra} \
          {input.list} \
          &> {log}
        test -s {output.gtf}
        """

rule gffcompare_merged:
    input:
        gtf=f"{OUTDIR}/assembly/stringtie_merged.gtf"
    output:
        annotated=f"{OUTDIR}/gffcompare/merged.annotated.gtf",
        tracking=f"{OUTDIR}/gffcompare/merged.tracking",
        stats=f"{OUTDIR}/gffcompare/merged.stats"
    params:
        ref=config["genome"]["gtf"],
        prefix=f"{OUTDIR}/gffcompare/merged"
    log:
        f"{OUTDIR}/logs/gffcompare/merged.log"
    conda:
        ENV_LNCRNA
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/gffcompare {OUTDIR}/logs/gffcompare
        gffcompare \
          -r {params.ref} \
          -o {params.prefix} \
          {input.gtf} \
          &> {log}
        test -s {output.annotated}
        test -s {output.tracking}
        test -s {output.stats}
        """
