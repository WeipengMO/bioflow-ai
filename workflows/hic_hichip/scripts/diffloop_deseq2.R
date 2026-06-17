suppressPackageStartupMessages({
  has_optparse <- requireNamespace("optparse", quietly = TRUE)
})

if (!has_optparse) {
  stop("R package 'optparse' is required")
}

option_list <- list(
  optparse::make_option("--counts", type="character"),
  optparse::make_option(c("--loop-metadata"), dest="loop_metadata", type="character"),
  optparse::make_option(c("--sample-metadata"), dest="sample_metadata", type="character"),
  optparse::make_option("--comparison", type="character"),
  optparse::make_option("--group1", type="character"),
  optparse::make_option("--group2", type="character"),
  optparse::make_option("--fdr", type="double", default=0.05),
  optparse::make_option(c("--lfc-cutoff"), dest="lfc_cutoff", type="double", default=0),
  optparse::make_option(c("--min-count"), dest="min_count", type="integer", default=10),
  optparse::make_option(c("--design-formula"), dest="design_formula", type="character", default="~ group"),
  optparse::make_option(c("--out-table"), dest="out_table", type="character"),
  optparse::make_option(c("--out-significant"), dest="out_significant", type="character"),
  optparse::make_option(c("--out-up"), dest="out_up", type="character"),
  optparse::make_option(c("--out-down"), dest="out_down", type="character"),
  optparse::make_option("--volcano", type="character"),
  optparse::make_option(c("--ma-plot"), dest="ma_plot", type="character"),
  optparse::make_option("--heatmap", type="character")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

counts <- read.delim(opt$counts, check.names=FALSE)
loop_meta <- read.delim(opt$loop_metadata, check.names=FALSE)
sample_meta <- read.delim(opt$sample_metadata, check.names=FALSE)
stopifnot("loop_id" %in% colnames(counts), "sample" %in% colnames(sample_meta), "group" %in% colnames(sample_meta))

samples <- sample_meta$sample[sample_meta$group %in% c(opt$group1, opt$group2)]
count_mat <- counts[, c("loop_id", samples), drop=FALSE]
rownames(count_mat) <- count_mat$loop_id
count_mat$loop_id <- NULL
count_mat <- round(as.matrix(count_mat))
coldata <- sample_meta[match(colnames(count_mat), sample_meta$sample), , drop=FALSE]
coldata$group <- factor(coldata$group, levels=c(opt$group2, opt$group1))
keep <- rowSums(count_mat) >= opt$min_count
count_mat_f <- count_mat[keep, , drop=FALSE]

run_deseq <- requireNamespace("DESeq2", quietly=TRUE) &&
  all(table(coldata$group) >= 2) &&
  nrow(count_mat_f) >= 2

if (run_deseq) {
  design_formula <- as.formula(opt$design_formula)
  needed_terms <- all.vars(design_formula)
  missing_terms <- setdiff(needed_terms, colnames(coldata))
  if (length(missing_terms) > 0) {
    stop(paste("Design formula references missing sample metadata columns:", paste(missing_terms, collapse=", ")))
  }
  if (!("group" %in% needed_terms)) {
    warning("design_formula does not contain group; DESeq2 contrast still uses group and may be invalid.")
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=count_mat_f, colData=coldata, design=design_formula)
  dds <- DESeq2::DESeq(dds, quiet=TRUE)
  res <- as.data.frame(DESeq2::results(dds, contrast=c("group", opt$group1, opt$group2)))
  res$loop_id <- rownames(res)
  res <- res[, c("loop_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
} else {
  warning("DESeq2 unavailable, too few replicates, or too few loops after filtering. Writing descriptive fold-change results.")
  group1_samples <- coldata$sample[coldata$group == opt$group1]
  group2_samples <- coldata$sample[coldata$group == opt$group2]
  mat <- count_mat_f
  g1 <- rowMeans(mat[, group1_samples, drop=FALSE] + 1)
  g2 <- rowMeans(mat[, group2_samples, drop=FALSE] + 1)
  res <- data.frame(
    loop_id=rownames(mat),
    baseMean=rowMeans(mat),
    log2FoldChange=log2(g1/g2),
    lfcSE=NA_real_, stat=NA_real_, pvalue=NA_real_, padj=NA_real_
  )
}

# Add group means for all loops in result.
group1_samples <- coldata$sample[coldata$group == opt$group1]
group2_samples <- coldata$sample[coldata$group == opt$group2]
mat_all <- count_mat[res$loop_id, , drop=FALSE]
res$group1_mean <- rowMeans(mat_all[, group1_samples, drop=FALSE])
res$group2_mean <- rowMeans(mat_all[, group2_samples, drop=FALSE])
res$direction <- ifelse(res$log2FoldChange > 0, paste0(opt$group1, "_up"), ifelse(res$log2FoldChange < 0, paste0(opt$group2, "_up"), "unchanged"))
res$comparison <- opt$comparison

merged <- merge(loop_meta, res, by="loop_id", all.y=TRUE)
merged <- merged[order(merged$padj, -abs(merged$log2FoldChange), na.last=TRUE), ]

sig <- merged[!is.na(merged$padj) & merged$padj <= opt$fdr & abs(merged$log2FoldChange) >= opt$lfc_cutoff, , drop=FALSE]
up <- sig[sig$log2FoldChange > 0, , drop=FALSE]
down <- sig[sig$log2FoldChange < 0, , drop=FALSE]

bedpe_cols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "loop_id", "log2FoldChange", "padj")
for (path in c(opt$out_table, opt$out_significant, opt$out_up, opt$out_down, opt$volcano, opt$ma_plot, opt$heatmap)) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
}
write.table(merged, opt$out_table, sep="\t", quote=FALSE, row.names=FALSE)
write.table(sig, opt$out_significant, sep="\t", quote=FALSE, row.names=FALSE)
write.table(up[, intersect(bedpe_cols, colnames(up)), drop=FALSE], opt$out_up, sep="\t", quote=FALSE, row.names=FALSE)
write.table(down[, intersect(bedpe_cols, colnames(down)), drop=FALSE], opt$out_down, sep="\t", quote=FALSE, row.names=FALSE)

pdf(opt$volcano)
if (nrow(merged) == 0) {
  plot.new(); text(0.5, 0.5, "No loops passed filtering")
} else if (all(is.na(merged$padj))) {
  plot(merged$log2FoldChange, log10(merged$baseMean + 1), pch=16, cex=0.6,
       xlab="log2FC", ylab="log10(baseMean + 1)",
       main=paste("Descriptive fold-change only", opt$comparison))
  abline(v=0, lty=2)
  mtext("No formal DESeq2 statistics; p-values/FDR are NA", side=3, line=0.2, cex=0.75)
} else {
  plot(merged$log2FoldChange, -log10(merged$padj), pch=16, cex=0.6, xlab="log2FC", ylab="-log10(FDR)", main=paste("Differential loops", opt$comparison))
  abline(v=c(-opt$lfc_cutoff, opt$lfc_cutoff), lty=2)
  abline(h=-log10(opt$fdr), lty=2)
}
dev.off()

pdf(opt$ma_plot)
if (nrow(merged) == 0) {
  plot.new(); text(0.5, 0.5, "No loops passed filtering")
} else {
  plot(log10(merged$baseMean + 1), merged$log2FoldChange, pch=16, cex=0.6, xlab="log10(baseMean + 1)", ylab="log2FC", main=paste("MA plot", opt$comparison))
  abline(h=0, lty=2)
}
dev.off()

pdf(opt$heatmap)
top_ids <- head(merged$loop_id[order(abs(merged$log2FoldChange), decreasing=TRUE)], 50)
if (length(top_ids) > 1) {
  hm <- log2(count_mat[top_ids, , drop=FALSE] + 1)
  heatmap(hm, scale="row", margins=c(8, 8), main=paste("Top loops", opt$comparison))
} else {
  plot.new(); text(0.5, 0.5, "Not enough loops for heatmap")
}
dev.off()
