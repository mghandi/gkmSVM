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
                 useTgkm = num_or(case$t, 1), alg = num_or(case$alg, 2), addRC = case$RC != "0",
                 usePseudocnt = case$p == "1", wildcardLambda = num_or(case$lam, 1.0),
                 wildcardMismatchM = num_or(case$M, 2), alphabetFN = alphabetFN)
  if (case$tool == "kernel") {
    args <- c(list(posfile = fx(case$pos), negfile = fx(case$neg), outfile = outfile), common,
              list(nmaxThreads = num_or(case$T, 1000)))
    do.call(gkmSVM::gkmsvm_kernel, args)
  } else {
    args <- c(list(seqfile = fx(case$seq), svmfnprfx = NA, outfile = outfile), common,
              list(batchSize = num_or(case$batch, 100000), svseqfile = fx(case$svseq), alphafile = fx(case$svalpha)))
    do.call(gkmSVM::gkmsvm_classify, args)
  }
  outfile
}

read_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)

# Run one case in a fresh R process. Until Phase 1 lands, ~38 consecutive gkmsvm_kernel() calls in
# one R session abort the process (heap corruption from the known reader bugs: seqBaseId[-1] is
# written at every EOF, names >= 100 chars overflow, ...; measured 2026-08-29 on macOS/R 4.5.3,
# the same call succeeds on its own). Phase 1 switches the golden tests back to in-process calls.
run_golden_case_subprocess <- function(case, outfile) {
  script <- tempfile(fileext = ".R")
  rds <- tempfile(fileext = ".rds")
  saveRDS(case, rds)
  writeLines(sprintf('
    suppressPackageStartupMessages(library(gkmSVM))
    Sys.setenv(GKMSVM_GOLDEN_DIR = "%s")
    source("%s")
    case <- readRDS("%s")
    invisible(capture.output(run_golden_case_R(case, "%s")))
  ', golden_dir(), file.path(golden_dir(), "..", "testthat", "helper-golden.R"), rds, outfile), script)
  rc <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script), stdout = FALSE, stderr = FALSE)
  if (rc != 0) stop(sprintf("case %s: Rscript exited with status %d", case$name, rc))
  outfile
}
