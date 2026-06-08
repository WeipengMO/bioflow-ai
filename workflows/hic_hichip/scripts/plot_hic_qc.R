suppressPackageStartupMessages({
  has_optparse <- requireNamespace("optparse", quietly = TRUE)
})
if (!has_optparse) stop("R package 'optparse' is required")
option_list <- list(
  optparse::make_option("--summary", type="character"),
  optparse::make_option("--mcools", type="character", action="append"),
  optparse::make_option("--samples", type="character"),
  optparse::make_option("--resolution", type="character"),
  optparse::make_option(c("--distance-pdf"), dest="distance_pdf", type="character"),
  optparse::make_option(c("--distance-tsv"), dest="distance_tsv", type="character"),
  optparse::make_option(c("--correlation-pdf"), dest="correlation_pdf", type="character")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
summary <- read.delim(opt$summary, check.names=FALSE)
samples <- strsplit(opt$samples, ",")[[1]]
for (p in c(opt$distance_pdf, opt$distance_tsv, opt$correlation_pdf)) dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)

# Placeholder distance-decay summary derived from counted valid pairs. For production, replace with cooltools expected-cis.
df <- data.frame(sample=samples, distance_bin="all", contacts=summary$valid_pairs_counted[match(samples, summary$sample)])
write.table(df, opt$distance_tsv, sep="\t", quote=FALSE, row.names=FALSE)

pdf(opt$distance_pdf)
barplot(as.numeric(df$contacts), names.arg=df$sample, las=2, ylab="valid pairs", main="Valid interaction pairs")
dev.off()

pdf(opt$correlation_pdf)
if (length(samples) > 1 && all(!is.na(as.numeric(df$contacts)))) {
  m <- matrix(as.numeric(df$contacts), nrow=1)
  colnames(m) <- samples
  plot.new(); text(0.5, 0.5, "Sample correlation requires per-bin contact vectors.\nThis placeholder confirms QC plotting completed.")
} else {
  plot.new(); text(0.5, 0.5, "Not enough samples for correlation")
}
dev.off()
