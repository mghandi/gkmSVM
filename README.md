# gkmSVM

R package (with standalone command-line tools) for the **gapped k-mer support vector machine**:
sequence classifiers that learn regulatory sequence features, e.g. from ChIP-seq or DNase-seq peaks,
and score any sequence or variant with them. The kernel counts gapped k-mers shared between two
sequences (a fast trie-based algorithm with a mismatch bound), the SVM is trained on the kernel
matrix, and the resulting model gives per-sequence scores, per-k-mer weights and variant effect
scores (deltaSVM).

If you use it, please cite

* M. Ghandi, D. Lee, M. Mohammad-Noori, M. A. Beer. *Enhanced regulatory sequence prediction using
  gapped k-mer features.* PLoS Computational Biology 10(7): e1003711 (2014).
  <https://doi.org/10.1371/journal.pcbi.1003711>
* M. Ghandi, M. Mohammad-Noori, N. Ghareghani, D. Lee, L. Garraway, M. A. Beer. *gkmSVM: an R
  package for gapped-kmer SVM.* Bioinformatics 32(14): 2205–2207 (2016).
  <https://doi.org/10.1093/bioinformatics/btw203>
* For variant scoring (deltaSVM): D. Lee, D. U. Gorkin, M. Baker, B. J. Strober, A. L. Asoni,
  A. S. McCallion, M. A. Beer. *A method to predict the impact of regulatory variants from DNA
  sequence.* Nature Genetics 47: 955–961 (2015). <https://doi.org/10.1038/ng.3331>

## Getting started

* **Tutorial:** [`tutorials/gkmsvm-tutorial.md`](tutorials/gkmsvm-tutorial.md) — installation and a
  complete CTCF example (matched negative set, kernel, cross-validated training, 10-mer weights,
  variant scoring); `tutorials/run_tutorial.R` runs it end to end.
* **Version used in the papers:** the code as of 2018, version 0.80, is kept as the
  [v0.80 release](https://github.com/mghandi/gkmSVM/releases/tag/v0.80) and is what the original
  tutorial at <https://www.beerlab.org/gkmsvm/gkmsvm-tutorial.htm> describes. The current version
  produces identical kernels; the classifier scores are normalised slightly differently, and
  `gkmsvm_classify(..., legacyNorm = TRUE)` reproduces the v0.80 scores exactly (see the tutorial).
* **Several tracks per sequence (DNA + methylation, protein + structure, sensor channels, …):**
  [`tutorials/gkmsvm-multitrack-tutorial.md`](tutorials/gkmsvm-multitrack-tutorial.md) — the
  gapped k-mer kernel over words whose positions come from alphabets of different sizes
  (M. Mohammad-Noori, N. Ghareghani, M. Ghandi, *Generalized Gapped k-mer Filters for Robust
  Frequency Estimation*), `gkmsvm_kernel(..., alphabets = c("dna", "01"))`; the design is in
  [`dev/PHASE7_PLAN.md`](dev/PHASE7_PLAN.md).
* **What changed since 0.80:** [`NEWS.md`](NEWS.md); the refactoring plan and its execution log are
  in [`dev/REFACTORING_PLAN.md`](dev/REFACTORING_PLAN.md).

Install from a clone with `R CMD INSTALL gkmSVM` (the Bioconductor genome packages are only needed
for `genNullSeqs`); `make` in the same tree builds the `gkmsvm_kernel`, `gkmsvm_train` and
`gkmsvm_classify` binaries. The CRAN release is at <https://cran.r-project.org/package=gkmSVM>.

## License

GPL (≥ 2). The embedded LIBSVM (`src/libsvm`) is BSD-3-Clause.
