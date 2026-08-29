# The R wrapper gkmsvm_kernel must reproduce the frozen CLI outputs byte-for-byte: same C++ core,
# same file format. Cases whose fixture the R layer itself rejects (duplicate names -> stop()) or
# that crash today (xfail-*) are skipped here; they are covered by tests/golden/golden.py.

test_that("gkmsvm_kernel reproduces the golden corpus", {
  skip_if_no_golden()
  cases <- golden_cases("kernel")
  # longnames: 120-char names overflow new char[100] and abort R (exit 134); the CLI survives by luck. Phase 1.
  skip_tags <- c("dupnames", "longnames", "xfail-phase1", "xfail-phase5")
  keep <- !vapply(strsplit(cases$tags, " "), function(t) any(t %in% skip_tags), logical(1))
  cases <- cases[keep, ]
  expect_gt(nrow(cases), 40)
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    out <- tempfile(fileext = ".txt")
    run_golden_case_subprocess(case, out)
    expected <- file.path(golden_dir(), "expected", paste0(case$name, ".out"))
    expect_true(file.exists(out), info = case$name)
    expect_identical(read_bytes(out), read_bytes(expected), info = case$name)
    unlink(out)
  }
})

test_that("gkmsvm_kernel rejects duplicated sequence names at the R level", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  expect_error(suppressWarnings(capture.output(
    gkmsvm_kernel(file.path(fx, "dupnames_pos.fa"), file.path(fx, "dna_small_neg.fa"), tempfile()))),
    "duplicated sequence ID")
})

test_that("the multithreaded kernel is identical to the single-threaded one", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  o1 <- tempfile(); o4 <- tempfile()
  capture.output(gkmsvm_kernel(file.path(fx, "dna_medium_pos.fa"), file.path(fx, "dna_medium_neg.fa"), o1, nmaxThreads = 1))
  capture.output(gkmsvm_kernel(file.path(fx, "dna_medium_pos.fa"), file.path(fx, "dna_medium_neg.fa"), o4, nmaxThreads = 4))
  expect_identical(read_bytes(o1), read_bytes(o4))
})

test_that("many gkmsvm_kernel calls in one R session do not crash", {
  skip_if_no_golden()
  skip("known (Phase 1): ~38 consecutive calls abort the R session with heap corruption; see helper-golden.R")
  fx <- file.path(golden_dir(), "fixtures")
  for (i in 1:60) capture.output(gkmsvm_kernel(file.path(fx, "dna_small_pos.fa"), file.path(fx, "dna_small_neg.fa"), tempfile()))
  succeed()
})
