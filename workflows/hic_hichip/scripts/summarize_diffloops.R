args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: summarize_diffloops.R <out.tsv> <diffloops.tsv> [more.tsv]")
out <- args[[1]]
files <- args[-1]
res <- do.call(rbind, lapply(files, function(f) {
  x <- read.delim(f, check.names=FALSE)
  data.frame(file=f, n_total=nrow(x), n_sig=sum(!is.na(x$padj) & x$padj <= 0.05), stringsAsFactors=FALSE)
}))
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
write.table(res, out, sep="\t", quote=FALSE, row.names=FALSE)
