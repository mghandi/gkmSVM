#!/usr/bin/env Rscript
# Runs the human (hg38) part of tutorials/gkmsvm-tutorial.md unattended, in a work directory:
#   Rscript tutorials/run_tutorial.R [workdir]     (default: ./tutorial_run)
# Downloads the input files from beerlab.org if they are not present, then runs genNullSeqs ->
# gkmsvm_kernel -> gkmsvm_trainCV -> gkmsvm_classify (default and v0.80-compatible) -> gkmsvm_delta,
# and prints the AUROC/AUPRC, the top 12 10-mers in both modes and a summary of the delta scores.
# Needs the packages listed under "Installation" in the tutorial (including
# BSgenome.Hsapiens.UCSC.hg38.masked); it takes a few minutes on a laptop.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[1] else "tutorial_run"
dir.create(workdir, showWarnings = FALSE, recursive = TRUE); setwd(workdir)
suppressPackageStartupMessages({ library(gkmSVM); library(BSgenome.Hsapiens.UCSC.hg38.masked) })
set.seed(1)   # the negative set and the CV folds are sampled; fix them so a re-run is comparable

for (f in c("CTCF_GM12878_hg38_top5k.bed", "nr10mers.fa", "ref.fa", "alt.fa")) {
  if (!file.exists(f)) download.file(paste0("https://www.beerlab.org/gkmsvm/downloads/", f), f, quiet = TRUE)
}
step <- function(msg) cat("\n==", msg, "\n")

step("1. genNullSeqs")
t <- system.time(genNullSeqs('CTCF_GM12878_hg38_top5k.bed', nMaxTrials=10, xfold=1,
                             genome=BSgenome.Hsapiens.UCSC.hg38.masked,
                             outputPosFastaFN='CTCF_GM12878_hg38_top5k.fa',
                             outputBedFN='neg1x_CTCF_GM12878_hg38_top5k.bed',
                             outputNegFastaFN='neg1x_CTCF_GM12878_hg38_top5k.fa'))
npos <- length(grep("^>", readLines("CTCF_GM12878_hg38_top5k.fa"))); nneg <- length(grep("^>", readLines("neg1x_CTCF_GM12878_hg38_top5k.fa")))
cat(sprintf("   %d positives, %d negatives, %.0f s\n", npos, nneg, t[["elapsed"]]))

step("2. gkmsvm_kernel")
t <- system.time(invisible(capture.output(gkmsvm_kernel('CTCF_GM12878_hg38_top5k.fa', 'neg1x_CTCF_GM12878_hg38_top5k.fa', 'CTCF_GM12878_vs_neg1x_kernel.out'))))
cat(sprintf("   %.0f s, %s\n", t[["elapsed"]], format(structure(file.size('CTCF_GM12878_vs_neg1x_kernel.out'), class = "object_size"), units = "MB")))

step("3. gkmsvm_trainCV")
t <- system.time(res <- gkmsvm_trainCV('CTCF_GM12878_vs_neg1x_kernel.out', 'CTCF_GM12878_hg38_top5k.fa',
                                       'neg1x_CTCF_GM12878_hg38_top5k.fa', svmfnprfx='CTCF_GM12878_vs_neg1x',
                                       outputCVpredfn='CTCF_GM12878_vs_neg1x_cvpred.out',
                                       outputROCfn='CTCF_GM12878_vs_neg1x_roc.out', showPlots = FALSE))
cat(sprintf("   AUROC %.3f  AUPRC %.3f  (tutorial: 0.965 / 0.964), %.0f s\n", res$aucss[["avgauc"]], res$aucss[["avgprc"]], t[["elapsed"]]))

step("4. gkmsvm_classify: 10-mer weights, default (consistent norm, bias applied)")
t <- system.time(invisible(capture.output(gkmsvm_classify('nr10mers.fa', svmfnprfx='CTCF_GM12878_vs_neg1x', 'CTCF_GM12878_vs_neg1x_weights.out'))))
wts <- read.table('CTCF_GM12878_vs_neg1x_weights.out'); cat(sprintf("   %.0f s\n", t[["elapsed"]]))
print(wts[order(wts$V2, decreasing = TRUE)[1:12], ], row.names = FALSE)

step("4b. gkmsvm_classify: v0.80-compatible (legacyNorm = TRUE, legacy model pair)")
t <- system.time(invisible(capture.output(gkmsvm_classify('nr10mers.fa', svmfnprfx='CTCF_GM12878_vs_neg1x', 'CTCF_GM12878_vs_neg1x_weights_v080.out',
                                                          legacyNorm = TRUE, svseqfile = 'CTCF_GM12878_vs_neg1x_svseq.fa', alphafile = 'CTCF_GM12878_vs_neg1x_svalpha.out'))))
wts0 <- read.table('CTCF_GM12878_vs_neg1x_weights_v080.out'); cat(sprintf("   %.0f s\n", t[["elapsed"]]))
print(wts0[order(wts0$V2, decreasing = TRUE)[1:12], ], row.names = FALSE)
cat(sprintf("   rank correlation between the two modes: %.6f; same top-12 set: %s\n",
            cor(wts$V2, wts0$V2, method = "spearman"), setequal(wts$V1[order(-wts$V2)[1:12]], wts0$V1[order(-wts0$V2)[1:12]])))

step("5. gkmsvm_delta (legacyNorm = TRUE)")
d <- gkmsvm_delta('ref.fa', 'alt.fa', svmfnprfx='CTCF_GM12878_vs_neg1x', 'dsvm_CTCF_GM12878_vs_neg1x.out', legacyNorm = TRUE)
cat(sprintf("   %d variants scored; deltaSVM range %.3f .. %.3f\n", length(d), min(d), max(d)))
cat("\nDone. Outputs are in", normalizePath("."), "\n")
