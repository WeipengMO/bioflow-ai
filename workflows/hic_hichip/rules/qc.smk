rule collect_hicpro_qc:
    input:
        done=PATHS.hicpro_done(),
        valid_pairs=expand(PATHS.sample_valid_pairs("{sample}"), sample=CTX.sample_names)
    output:
        summary=PATHS.hicpro_qc_summary()
    params:
        samples=",".join(CTX.sample_names),
        hicpro_dir=PATHS.hicpro
    log:
        PATHS.log("qc/collect_hicpro_qc")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.summary:q}) $(dirname {log:q})
python scripts/collect_hicpro_qc.py \
    --hicpro-dir {params.hicpro_dir:q} \
    --samples {params.samples:q} \
    --validpairs {input.valid_pairs:q} \
    --output {output.summary:q} \
    &> {log:q}
test -s {output.summary:q}
        """

rule plot_hic_qc:
    input:
        summary=PATHS.hicpro_qc_summary(),
        mcools=expand(PATHS.mcool_file("{sample}"), sample=CTX.sample_names)
    output:
        distance_pdf=PATHS.distance_decay_pdf(),
        distance_tsv=PATHS.distance_decay_tsv(),
        correlation_pdf=PATHS.sample_correlation_pdf()
    params:
        samples=",".join(CTX.sample_names),
        resolution=CTX.primary_resolution
    log:
        PATHS.log("qc/plot_hic_qc")
    conda:
        DIFF_ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.distance_pdf:q}) $(dirname {log:q})
Rscript scripts/plot_hic_qc.R \
    --summary {input.summary:q} \
    --mcools "{input.mcools}" \
    --samples {params.samples:q} \
    --resolution {params.resolution:q} \
    --distance-pdf {output.distance_pdf:q} \
    --distance-tsv {output.distance_tsv:q} \
    --correlation-pdf {output.correlation_pdf:q} \
    &> {log:q}
test -s {output.distance_pdf:q}
test -s {output.distance_tsv:q}
test -s {output.correlation_pdf:q}
        """

rule multiqc_report:
    input:
        summary=PATHS.hicpro_qc_summary()
    output:
        html=PATHS.multiqc_report()
    params:
        search_dir=lambda wildcards: PATHS.outdir
    log:
        PATHS.log("qc/multiqc")
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.html:q}) $(dirname {log:q})
multiqc {params.search_dir:q} \
    -o $(dirname {output.html:q}) \
    -n $(basename {output.html:q}) \
    --force \
    &> {log:q}
test -s {output.html:q}
        """
