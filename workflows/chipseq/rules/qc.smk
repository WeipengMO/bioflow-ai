CHIPQC_REPORT_DIR = str(config.get("chipqc_report_dir") or PATHS.default_chipqc_report)
CHIPQC_USER_SAMPLE_SHEET = config.get("chipqc_sample_sheet", "")
CHIPQC_SAMPLE_SHEET = CHIPQC_USER_SAMPLE_SHEET or f"{CHIPQC_REPORT_DIR}/chipqc_sample_sheet.generated.tsv"
CHIPQC_DONE = f"{CHIPQC_REPORT_DIR}/.chipqc_complete"
CHIPQC_ANNOTATION_ARG = config.get("chipqc_annotation", "") or "NA"
CHIPQC_REPORT_FACET = "TRUE" if config_bool("chipqc_report_facet", False) else "FALSE"


if ENABLE_CHIPQC and not CHIPQC_USER_SAMPLE_SHEET:
    rule generate_chipqc_sample_sheet:
        input:
            config_file="config/config.yml",
            sample_config=str(config.get("sample_config", "config/samples.yml"))
        output:
            CHIPQC_SAMPLE_SHEET
        log:
            PATHS.log("generate_chipqc_sample_sheet")
        run:
            Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
            rows = chipqc_sample_sheet_rows(CTX)
            write_chipqc_sample_sheet(output[0], rows)
            with Path(log[0]).open("w") as handle:
                handle.write(f"Wrote {len(rows)} rows to {output[0]}\n")


rule chipqc_report:
    input:
        sample_sheet=CHIPQC_SAMPLE_SHEET,
        bams=expand(analysis_bam("{sample}"), sample=TREATMENTS),
        bais=expand(analysis_bai("{sample}"), sample=TREATMENTS),
        controls=expand(analysis_bam("{control}"), control=EFFECTIVE_CONTROLS) if HAS_EFFECTIVE_CONTROL else [],
        control_bais=expand(analysis_bai("{control}"), control=EFFECTIVE_CONTROLS) if HAS_EFFECTIVE_CONTROL else [],
        peaks=peak_selection_outputs(TREATMENTS, CHIPQC_PEAK_MODE, CHIPQC_PEAK_TYPE)
    output:
        touch(CHIPQC_DONE)
    params:
        report_dir=CHIPQC_REPORT_DIR,
        annotation=CHIPQC_ANNOTATION_ARG,
        facet=CHIPQC_REPORT_FACET
    log:
        PATHS.log("chipqc")
    threads:
        workflow_threads("chipqc", 8)
    conda:
        "../envs/chipqc.yml"
    shell:
        r"""
set -euo pipefail
mkdir -p {params.report_dir:q} $(dirname {log:q})
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
    library(DiffBind)
}})

chipqc_bpparam <- SerialParam()
register(chipqc_bpparam, default = TRUE)

samples <- read.delim(
    sample_sheet,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(0)
)
for (column in intersect(c("ControlID", "bamControl"), colnames(samples))) {{
    samples[[column]][samples[[column]] %in% c("", "NA")] <- NA_character_
}}

chipqc_dba <- DiffBind::dba(sampleSheet = samples, peakCaller = "bed")
if (is.null(rownames(chipqc_dba$class))) {{
    chipqc_class_fields <- c(
        "ID",
        "Tissue",
        "Factor",
        "Condition",
        "Consensus",
        "Caller",
        "Control",
        "Reads",
        "Replicate",
        "bamRead",
        "bamControl",
        "Treatment",
        "Peaks"
    )
    if (nrow(chipqc_dba$class) != length(chipqc_class_fields)) {{
        stop(
            "Unexpected DiffBind class layout: expected ",
            length(chipqc_class_fields),
            " rows, found ",
            nrow(chipqc_dba$class),
            "."
        )
    }}
    rownames(chipqc_dba$class) <- chipqc_class_fields
}}

args <- list(experiment = chipqc_dba)
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
