# The R wrapper gkmsvm_kernel must reproduce the frozen CLI outputs byte-for-byte: same C++ core,
# same file format. All cases run in this one R process (Phase 1 fixed the heap corruption that
# used to abort R after ~38 calls). Cases whose fixture the R layer itself rejects (duplicate names
# -> stop()) or that still crash (xfail-*) are skipped; they are covered by tests/golden/golden.py.
# expect-error cases must surface as R errors.

test_that("gkmsvm_kernel reproduces the golden corpus", {
  skip_if_no_golden()
  cases <- golden_cases("kernel")
  skip_tags <- c("xfail-phase1", "xfail-phase5")
  keep <- !vapply(strsplit(cases$tags, " "), function(t) any(t %in% skip_tags), logical(1))
  cases <- cases[keep, ]
  expect_gt(nrow(cases), 40)
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    out <- tempfile(fileext = ".txt")
    if (has_tag(case, "expect-error")) {
      expect_error(capture.output(run_golden_case_R(case, out)), "gkmsvm_kernel failed", info = case$name)
      expect_false(file.exists(out), info = case$name)
      next
    }
    capture.output(run_golden_case_R(case, out))
    expected <- file.path(golden_dir(), "expected", paste0(case$name, ".out"))
    expect_true(file.exists(out), info = case$name)
    expect_identical(read_bytes(out), read_bytes(expected), info = case$name)
    unlink(out)
  }
})


test_that("the multithreaded kernel is identical to the single-threaded one", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  o1 <- tempfile(); o4 <- tempfile()
  capture.output(gkmsvm_kernel(file.path(fx, "dna_medium_pos.fa"), file.path(fx, "dna_medium_neg.fa"), o1, nmaxThreads = 1))
  capture.output(gkmsvm_kernel(file.path(fx, "dna_medium_pos.fa"), file.path(fx, "dna_medium_neg.fa"), o4, nmaxThreads = 4))
  expect_identical(read_bytes(o1), read_bytes(o4))
})

test_that("many gkmsvm_kernel calls in one R session do not crash and give identical output", {
  # before Phase 1, ~38 consecutive calls aborted R (SIGABRT) from accumulated heap corruption
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  expected <- read_bytes(file.path(golden_dir(), "expected", "k_small_t1_a2.out"))
  for (i in 1:60) {
    out <- tempfile()
    capture.output(gkmsvm_kernel(file.path(fx, "dna_small_pos.fa"), file.path(fx, "dna_small_neg.fa"), out))
    expect_identical(read_bytes(out), expected, info = paste("call", i))
    unlink(out)
  }
})

test_that("an alphabet file does not leak into the next call", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  o1 <- tempfile(); o2 <- tempfile()
  capture.output(gkmsvm_kernel(file.path(fx, "b2_pos.fa"), file.path(fx, "b2_neg.fa"), o1, L = 8, K = 5, maxnmm = -1,
                               alphabetFN = file.path(fx, "alphabet_b2.txt")))
  capture.output(gkmsvm_kernel(file.path(fx, "dna_small_pos.fa"), file.path(fx, "dna_small_neg.fa"), o2))
  expect_identical(read_bytes(o1), read_bytes(file.path(golden_dir(), "expected", "k_b2_t1.out")))
  expect_identical(read_bytes(o2), read_bytes(file.path(golden_dir(), "expected", "k_small_t1_a2.out")))
})

test_that("an alphabet larger than the supported maximum (32) is an R error, not a crash", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  big <- tempfile(fileext = ".txt")
  writeLines(c(LETTERS, letters[1:7]), big)   # 33 symbols
  expect_error(capture.output(gkmsvm_kernel(file.path(fx, "b5_pos.fa"), file.path(fx, "b5_neg.fa"), tempfile(),
                                            alphabetFN = big)), "gkmsvm_kernel failed")
})
