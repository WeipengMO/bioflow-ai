#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  args[[idx + 1]]
}

counts_path <- get_arg("--counts")
universe_path <- get_arg("--universe")
design_path <- get_arg("--design")
contrast_name <- get_arg("--contrast", "contrast")
spikein_metrics_path <- get_arg("--spikein-metrics", "")
fdr_threshold <- as.numeric(get_arg("--fdr", "0.05"))
log2fc_threshold <- as.numeric(get_arg("--log2fc", "1.0"))
output_tsv <- get_arg("--output-tsv")
output_bed <- get_arg("--output-bed")

counts <- read.delim(counts_path, check.names = FALSE)
design <- read.delim(design_path, check.names = FALSE)
universe <- read.delim(universe_path, header = FALSE, check.names = FALSE)
colnames(universe)[1:4] <- c("chrom", "start", "end", "peak_id")

sample_cols <- intersect(design$sample, colnames(counts))
if (length(sample_cols) < 4) {
  stop("DESeq2 RNaseH mode requires at least two treatment/RNaseH pairs.")
}

count_matrix <- as.matrix(counts[, sample_cols, drop = FALSE])
rownames(count_matrix) <- counts$peak_id
storage.mode(count_matrix) <- "integer"

coldata <- design[match(sample_cols, design$sample), , drop = FALSE]
rownames(coldata) <- coldata$sample
coldata$condition <- relevel(factor(coldata$condition), ref = "RNaseH")
coldata$pair <- factor(coldata$pair)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ pair + condition
)
dds <- dds[rowSums(counts(dds)) > 0, ]
if (!is.null(spikein_metrics_path) && nzchar(spikein_metrics_path) && file.exists(spikein_metrics_path)) {
  spikein <- read.delim(spikein_metrics_path, check.names = FALSE)
  if (all(c("sample", "spikein_scale_factor") %in% colnames(spikein))) {
    factors <- spikein$spikein_scale_factor[match(colnames(dds), spikein$sample)]
    factors <- suppressWarnings(as.numeric(factors))
    if (all(is.finite(factors)) && all(factors > 0)) {
      sf <- 1 / factors
      sf <- sf / exp(mean(log(sf)))
      sizeFactors(dds) <- sf
    }
  }
}
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "NoRNaseH", "RNaseH"))
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df <- merge(counts[, c("peak_id", "chrom", "start", "end")], res_df, by = "peak_id", all.x = TRUE)
res_df$contrast <- contrast_name
res_df$classification <- ifelse(
  !is.na(res_df$padj) & res_df$padj <= fdr_threshold & !is.na(res_df$log2FoldChange) & res_df$log2FoldChange >= log2fc_threshold,
  "rnaseh_depleted",
  "not_significant"
)

dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(res_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

passed <- res_df[res_df$classification == "rnaseh_depleted", , drop = FALSE]
bed <- data.frame(
  chrom = passed$chrom,
  start = passed$start,
  end = passed$end,
  name = passed$peak_id,
  score = ifelse(is.na(passed$padj), 0, pmin(1000, round(-log10(passed$padj) * 100))),
  strand = ".",
  log2FoldChange = passed$log2FoldChange,
  padj = passed$padj
)
write.table(bed, output_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
