#!/usr/bin/env Rscript
# Minimal report generator for ChIP-seq QC
# Usage: Rscript compile_chipseq_qc_report.R qc_results/

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide qc_results directory, e.g. Rscript compile_chipseq_qc_report.R qc_results/")
}
qc_dir <- args[1]
if (!dir.exists(qc_dir)) stop(paste("Directory not found:", qc_dir))

# Helpers
read_first_number <- function(path, default=NA_real_) {
  if (!file.exists(path)) return(default)
  x <- tryCatch(readLines(path, n=1, warn=FALSE), error=function(e) NULL)
  if (is.null(x)) return(default)
  as.numeric(gsub("[^0-9eE+\\.-]", "", x))
}

# Collect samples from *.spp.txt or *.flagstat.txt
spp_files <- list.files(qc_dir, pattern="spp\\.txt$", full.names=TRUE)
flagstat_files <- list.files(qc_dir, pattern="flagstat\\.txt$", full.names=TRUE)
samples <- unique(gsub("_spp\\.txt$|\\.flagstat\\.txt$", "", basename(c(spp_files, flagstat_files))))

# Parse simple metrics from known files
parse_spp <- function(f){
  # phantompeakqualtools outputs tab-delimited; we will try to read and search for NSC/RSC
  res <- list(NSC=NA_real_, RSC=NA_real_, FragLen=NA_integer_)
  if (!file.exists(f)) return(res)
  txt <- tryCatch(readLines(f, warn=FALSE), error=function(e) character())
  # heuristics
  nsc_line <- grep("NSC", txt, value=TRUE)
  rsc_line <- grep("RSC", txt, value=TRUE)
  frag_line <- grep("peak strand cross-correlation at fragment length|estimated fragment length", txt, value=TRUE, ignore.case=TRUE)
  num <- function(s) { as.numeric(gsub(".*?([0-9]+\\.?[0-9eE\\+-]*).*", "\\1", s)) }
  if (length(nsc_line)>0) res$NSC <- suppressWarnings(num(nsc_line[1]))
  if (length(rsc_line)>0) res$RSC <- suppressWarnings(num(rsc_line[1]))
  if (length(frag_line)>0) res$FragLen <- suppressWarnings(as.integer(num(frag_line[1])))
  res
}

parse_frip <- function(sample){
  frip_path <- file.path(qc_dir, paste0(sample, ".frip.txt"))
  if (!file.exists(frip_path)) return(NA_real_)
  val <- suppressWarnings(as.numeric(readLines(frip_path, warn=FALSE)[1]))
  val
}

# Build a data frame
NSC <- RSC <- FragLen <- FRiP <- rep(NA_real_, length(samples))
for (i in seq_along(samples)) {
  s <- samples[i]
  spp_path <- file.path(qc_dir, paste0(s, "_spp.txt"))
  if (!file.exists(spp_path)) spp_path <- file.path(qc_dir, paste0(s, "spp.txt"))
  met <- parse_spp(spp_path)
  NSC[i] <- met$NSC
  RSC[i] <- met$RSC
  FragLen[i] <- met$FragLen
  FRiP[i] <- parse_frip(s)
}

df <- data.frame(
  sample=samples,
  NSC=round(NSC,3),
  RSC=round(RSC,3),
  FragLen=FragLen,
  FRiP=signif(FRiP,3),
  stringsAsFactors=FALSE
)

# Output PDF
pdf(file.path(qc_dir, "ChIPseq_QC_Report.pdf"), width=8.5, height=11)
op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)

# Title Page
plot.new(); title("ChIP-seq QC Report", cex.main=1.6)
mtext(paste("Generated:", Sys.time()), side=3, line=-2, cex=0.8)
mtext(paste("QC directory:", normalizePath(qc_dir)), side=1, line=2, cex=0.7)

# Summary Table
plot.new()
title("Summary Metrics (NSC, RSC, FragLen, FRiP)")
if (nrow(df)>0) {
  # Draw as table
  grid <- function(){
    par(mar=c(1,1,1,1)); plot.new()
  }
  # Render table using base
  text(0.5, 0.95, "Per-sample summary", cex=1.2)
  y <- 0.85
  headers <- c("Sample","NSC","RSC","FragLen","FRiP")
  text(x=c(0.2,0.45,0.6,0.75,0.9), y=y, labels=headers, font=2)
  y <- y - 0.05
  for (i in seq_len(nrow(df))) {
    text(0.2, y, df$sample[i], adj=0)
    text(0.45, y, ifelse(is.na(df$NSC[i]), "NA", df$NSC[i]))
    text(0.6,  y, ifelse(is.na(df$RSC[i]), "NA", df$RSC[i]))
    text(0.75, y, ifelse(is.na(df$FragLen[i]), "NA", df$FragLen[i]))
    text(0.9,  y, ifelse(is.na(df$FRiP[i]), "NA", df$FRiP[i]))
    y <- y - 0.04
    if (y < 0.1 && i < nrow(df)) { plot.new(); y <- 0.9 }
  }
} else {
  text(0.5, 0.5, "No metrics found.", cex=1.2)
}

# Guidance Page
plot.new(); title("Interpretation Guidelines")
txt <- c(
  "NSC > 1.05 (ok), > 1.10 (good)",
  "RSC > 0.8 (ok), > 1.0 (good)",
  "NRF > 0.8, PBC > 0.8 suggest good library complexity",
  "FRiP: TFs > 0.05 (good), Histone > 0.1 (good)",
  "IDR < 0.05 indicates strong replicate reproducibility",
  "Correlation between replicates r > 0.8 is a good sign"
)
y <- 0.8
for (line in txt) { text(0.1, y, line, adj=0); y <- y - 0.08 }

dev.off()

# Also write a TSV for programmatic consumption
utils::write.table(df, file=file.path(qc_dir, "ChIPseq_QC_Summary.tsv"),
                   sep="\t", quote=FALSE, row.names=FALSE)
cat("Report written to", file.path(qc_dir, "ChIPseq_QC_Report.pdf"), "\\n")