CHIPQC_REPORT_DIR = config.get("chipqc_report_dir", "reports/chipqc")
CHIPQC_USER_SAMPLE_SHEET = config.get("chipqc_sample_sheet", "")
CHIPQC_SAMPLE_SHEET = CHIPQC_USER_SAMPLE_SHEET or f"{CHIPQC_REPORT_DIR}/chipqc_sample_sheet.generated.tsv"
CHIPQC_DONE = f"{CHIPQC_REPORT_DIR}/.chipqc_complete"
CHIPQC_ANNOTATION = config.get("chipqc_annotation", "")
CHIPQC_ANNOTATION_ARG = CHIPQC_ANNOTATION or "NA"
CHIPQC_REPORT_FACET = "TRUE" if config.get("chipqc_report_facet", False) else "FALSE"


if ENABLE_CHIPQC and not CHIPQC_USER_SAMPLE_SHEET:
    rule generate_chipqc_sample_sheet:
        input:
            samples=SAMPLES_CONFIG,
            replicates=REPLICATE_CONFIG if Path(REPLICATE_CONFIG).exists() else []
        output:
            CHIPQC_SAMPLE_SHEET
        params:
            peak_mode=CHIPQC_PEAK_MODE,
            peak_type=CHIPQC_PEAK_TYPE,
            bam_dir="aligned_data",
            no_control_value="NA"
        log:
            "logs/generate_chipqc_sample_sheet.log"
        conda:
            "../envs/chipseq.yml"
        shell:
            r"""
set -euo pipefail
mkdir -p $(dirname {output:q}) logs
python scripts/generate_chipqc_sample_sheet.py \
    --samples {input.samples:q} \
    --replicates {REPLICATE_CONFIG:q} \
    --output {output:q} \
    --peak-mode {params.peak_mode:q} \
    --peak-type {params.peak_type:q} \
    --bam-dir {params.bam_dir:q} \
    --no-control-value {params.no_control_value:q} \
    &> {log:q}
test -s {output:q}
            """


rule chipqc_report:
    input:
        sample_sheet=CHIPQC_SAMPLE_SHEET,
        bams=expand("aligned_data/{sample}.sorted.rmdup.bam", sample=PEAK_TREATMENTS),
        bais=expand("aligned_data/{sample}.sorted.rmdup.bam.bai", sample=PEAK_TREATMENTS),
        controls=expand("aligned_data/{control}.sorted.rmdup.bam", control=EFFECTIVE_CONTROLS) if HAS_EFFECTIVE_CONTROL else [],
        control_bais=expand("aligned_data/{control}.sorted.rmdup.bam.bai", control=EFFECTIVE_CONTROLS) if HAS_EFFECTIVE_CONTROL else [],
        peaks=peak_selection_outputs(PEAK_TREATMENTS, CHIPQC_PEAK_MODE, CHIPQC_PEAK_TYPE)
    output:
        touch(CHIPQC_DONE)
    params:
        report_dir=CHIPQC_REPORT_DIR,
        annotation=CHIPQC_ANNOTATION_ARG,
        facet=CHIPQC_REPORT_FACET
    log:
        "logs/chipqc.log"
    threads:
        workflow_threads("chipqc", 8)
    conda:
        "../envs/chipqc.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p {params.report_dir:q} logs
Rscript - {input.sample_sheet:q} {params.report_dir:q} {params.annotation:q} {params.facet:q} {threads} <<'RSCRIPT' &> {log:q}
args <- commandArgs(trailingOnly = TRUE)
sample_sheet <- args[[1]]
report_dir <- args[[2]]
annotation <- args[[3]]
if (annotation == "NA") annotation <- ""
facet_report <- as.logical(args[[4]])
workers <- as.integer(args[[5]])

suppressPackageStartupMessages({{
    library(BiocParallel)
    library(ChIPQC)
}})

chipqc_bpparam <- MulticoreParam(workers = workers)
register(chipqc_bpparam, default = TRUE)

samples <- read.delim(
    sample_sheet,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(0)
)
args <- list(experiment = samples)
if (nzchar(annotation)) args$annotation <- annotation
experiment <- do.call(ChIPQC, args)
ChIPQCreport(
    experiment,
    reportName = "ChIPQC",
    reportFolder = report_dir,
    facet = facet_report
)
RSCRIPT
        """
