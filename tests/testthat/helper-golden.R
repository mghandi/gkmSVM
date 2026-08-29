# Shared helpers for the golden tests: locate the corpus, read cases.tsv, and map one case to a call
# of the R wrapper. The corpus lives in tests/golden (kept out of the built package by .Rbuildignore),
# so these tests only run from a source checkout: GKMSVM_GOLDEN_DIR or ../golden relative to this dir.

golden_dir <- function() {
  d <- Sys.getenv("GKMSVM_GOLDEN_DIR", unset = "")
  if (nzchar(d)) return(d)
  for (cand in c(file.path(testthat::test_path(), "..", "golden"),
                 file.path(getwd(), "tests", "golden"), file.path(getwd(), "..", "golden"))) {
    if (file.exists(file.path(cand, "cases.tsv"))) return(normalizePath(cand))
  }
  NULL
}

skip_if_no_golden <- function() {
  if (is.null(golden_dir())) testthat::skip("golden corpus (tests/golden) not available")
}

golden_cases <- function(tool = NULL, tag = NULL) {
  cases <- utils::read.delim(file.path(golden_dir(), "cases.tsv"), stringsAsFactors = FALSE,
                             colClasses = "character", na.strings = character(0))
  if (!is.null(tool)) cases <- cases[cases$tool == tool, ]
  if (!is.null(tag)) cases <- cases[vapply(strsplit(cases$tags, " "), function(t) tag %in% t, logical(1)), ]
  cases
}

num_or <- function(x, default) if (nzchar(x)) as.numeric(x) else default

# Run one case through the R wrapper (same parameter mapping as tests/golden/golden.py build_argv),
# returning the path of the produced output file.
run_golden_case_R <- function(case, outfile) {
  fx <- function(f) file.path(golden_dir(), "fixtures", f)
  alphabetFN <- if (nzchar(case$A)) fx(case$A) else "NULL"
  common <- list(L = num_or(case$L, 10), K = num_or(case$K, 6), maxnmm = num_or(case$d, 3),
                 maxseqlen = num_or(case$m, 10000), maxnumseq = num_or(case$n, 1000000),
                 useTgkm = num_or(case$t, 1), alg = num_or(case$alg, 2), addRC = case$RC != "0",
                 usePseudocnt = case$p == "1", wildcardLambda = num_or(case$lam, 1.0),
                 wildcardMismatchM = num_or(case$M, 2), alphabetFN = alphabetFN)
  if (case$tool == "kernel") {
    args <- c(list(posfile = fx(case$pos), negfile = fx(case$neg), outfile = outfile), common,
              list(nmaxThreads = num_or(case$T, 1000), mergeByName = identical(case$N, "1"),
                   format = if (identical(case$bin, "1")) "binary" else "text"))
    do.call(gkmSVM::gkmsvm_kernel, args)
  } else {
    args <- c(list(seqfile = fx(case$seq), svmfnprfx = NA, outfile = outfile), common,
              list(batchSize = num_or(case$batch, 100000), svseqfile = fx(case$svseq), alphafile = fx(case$svalpha)))
    do.call(gkmSVM::gkmsvm_classify, args)
  }
  outfile
}

read_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)

has_tag <- function(case, tag) tag %in% strsplit(case$tags, " ")[[1]]
