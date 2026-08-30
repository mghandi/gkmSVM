#!/usr/bin/env Rscript
# Runs tutorials/gkmsvm-multitrack-tutorial.md unattended, in a work directory:
#   Rscript tutorials/run_multitrack_tutorial.R [workdir]     (default: ./multitrack_run)
# Generates a synthetic DNA + methylation data set (no download), then runs gkmsvm_kernel ->
# gkmsvm_trainCV -> gkmsvm_train -> gkmsvm_classify with alphabets = c("dna", "01"), and the same
# with the DNA track alone and the methylation track alone, and prints the cross-validated AUROC of
# each. Needs gkmSVM, kernlab and ROCR; it takes well under a minute.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[1] else "multitrack_run"
dir.create(workdir, showWarnings = FALSE, recursive = TRUE); setwd(workdir)
suppressPackageStartupMessages(library(gkmSVM))
set.seed(1)

## ---- 1. a synthetic data set: two motifs, methylated or not ------------------------------------
## Positives carry TGACTCA unmethylated and CACGTG methylated; negatives carry the same two motifs
## with the methylation states swapped. Neither track alone separates the classes; the joint
## (base, methylation) pattern does. Methylation is a 0/1 flag per base ("meth2s": both bases of a
## CpG are flagged, so reverse complements stay consistent).
rand_dna <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
make_record <- function(len, meth_tgactca, meth_cacgtg) {
  dna <- strsplit(rand_dna(len), "")[[1]]
  ann <- as.character(rbinom(len, 1, 0.3))
  repeat {
    p1 <- sample(len - 7 + 1, 1); p2 <- sample(len - 6 + 1, 1)
    if (p1 + 7 <= p2 || p2 + 6 <= p1) break
  }
  dna[p1:(p1 + 6)] <- strsplit("TGACTCA", "")[[1]]; ann[p1:(p1 + 6)] <- meth_tgactca
  dna[p2:(p2 + 5)] <- strsplit("CACGTG", "")[[1]];  ann[p2:(p2 + 5)] <- meth_cacgtg
  c(paste(dna, collapse = ""), paste(ann, collapse = ""))
}
n <- 60; len <- 40
pos <- lapply(seq_len(n), function(i) make_record(len, "0", "1")); names(pos) <- paste0("pos", seq_len(n))
neg <- lapply(seq_len(n), function(i) make_record(len, "1", "0")); names(neg) <- paste0("neg", seq_len(n))
write_mfa(pos, "pos.mfa"); write_mfa(neg, "neg.mfa")
cat("wrote pos.mfa / neg.mfa:", n, "+", n, "records of", len, "bp, two tracks each\n")
cat(readLines("pos.mfa", 3), sep = "\n")

## ---- 2. kernel, cross-validation and a model with both tracks -----------------------------------
alphabets <- c("dna", "01")
step <- function(msg) cat("\n==", msg, "\n")
step("gkmsvm_kernel, alphabets = c('dna', '01')")
invisible(capture.output(gkmsvm_kernel("pos.mfa", "neg.mfa", "kernel_joint.txt", L = 6, K = 4, maxnmm = -1, useTgkm = 0, alphabets = alphabets)))
step("gkmsvm_trainCV (5-fold)")
cv <- gkmsvm_trainCV("kernel_joint.txt", "pos.mfa", "neg.mfa", svmfnprfx = "model_joint", nCV = 5, showPlots = FALSE, alphabets = alphabets)
cat(sprintf("   joint (DNA + methylation): AUROC %.3f  AUPRC %.3f\n", cv$aucss[2], cv$aucss[3]))   # aucss = (C, AUROC, AUPRC) of the best C
step("gkmsvm_train + gkmsvm_classify")
invisible(capture.output(gkmsvm_train("kernel_joint.txt", "pos.mfa", "neg.mfa", "model_joint", alphabets = alphabets)))
cat("   model_joint.gkmmodel records its alphabets:", grep("^#alphabets", readLines("model_joint.gkmmodel"), value = TRUE), "\n")
invisible(capture.output(gkmsvm_classify("pos.mfa", "model_joint", "scores_pos.txt", L = 6, K = 4, maxnmm = -1, useTgkm = 0)))
sp <- read.delim("scores_pos.txt", header = FALSE)[[2]]
cat(sprintf("   training positives scored by the model: min %.3f, mean %.3f\n", min(sp), mean(sp)))

## ---- 3. each track alone, for comparison --------------------------------------------------------
one_track <- function(recs, i) { r <- lapply(recs, function(x) x[i]); names(r) <- names(recs); r }
write_mfa(one_track(pos, 1), "pos_dna.mfa"); write_mfa(one_track(neg, 1), "neg_dna.mfa")     # plain FASTA
write_mfa(one_track(pos, 2), "pos_meth.mfa"); write_mfa(one_track(neg, 2), "neg_meth.mfa")
step("DNA track only (the 2014 gkm-SVM)")
invisible(capture.output(gkmsvm_kernel("pos_dna.mfa", "neg_dna.mfa", "kernel_dna.txt", L = 6, K = 4, maxnmm = -1, useTgkm = 0)))
cv <- gkmsvm_trainCV("kernel_dna.txt", "pos_dna.mfa", "neg_dna.mfa", nCV = 5, showPlots = FALSE)
cat(sprintf("   DNA only:                  AUROC %.3f\n", cv$aucss[2]))
step("methylation track only")
invisible(capture.output(gkmsvm_kernel("pos_meth.mfa", "neg_meth.mfa", "kernel_meth.txt", L = 6, K = 4, maxnmm = -1, useTgkm = 0, addRC = FALSE, alphabets = "01")))
cv <- gkmsvm_trainCV("kernel_meth.txt", "pos_meth.mfa", "neg_meth.mfa", nCV = 5, showPlots = FALSE)
cat(sprintf("   methylation only:          AUROC %.3f\n", cv$aucss[2]))
