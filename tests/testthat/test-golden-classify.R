test_that("gkmsvm_classify reproduces the golden corpus", {
  skip_if_no_golden()
  cases <- golden_cases("classify")
  expect_gt(nrow(cases), 10)
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    out <- tempfile(fileext = ".txt")
    if (has_tag(case, "expect-error")) {
      expect_error(capture.output(run_golden_case_R(case, out)), "gkmsvm_classify failed", info = case$name)
      next
    }
    capture.output(run_golden_case_R(case, out))
    expected <- file.path(golden_dir(), "expected", paste0(case$name, ".out"))
    expect_identical(read_bytes(out), read_bytes(expected), info = case$name)
    unlink(out)
  }
})

test_that("gkmsvm_delta subtracts scores by row position", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  out <- tempfile()
  capture.output(res <- gkmsvm_delta(file.path(fx, "test_seqs.fa"), file.path(fx, "test_seqs.fa"), svmfnprfx = NA,
                                     outfile = out, svseqfile = file.path(fx, "sv_svseq.fa"),
                                     alphafile = file.path(fx, "sv_svalpha.out")))
  expect_length(res, 10)
  expect_true(all(res == 0))
  tab <- utils::read.table(out, stringsAsFactors = FALSE)
  expect_equal(tab[[1]], paste0("t", 0:9))
})

test_that("classify output is complete inside a long-lived R session (alg=1 used to never fclose)", {
  skip_if_no_golden()
  fx <- file.path(golden_dir(), "fixtures")
  out <- tempfile()
  capture.output(gkmsvm_classify(file.path(fx, "test_seqs.fa"), svmfnprfx = NA, outfile = out, alg = 1,
                                 svseqfile = file.path(fx, "sv_svseq.fa"), alphafile = file.path(fx, "sv_svalpha.out")))
  expect_identical(read_bytes(out), read_bytes(file.path(golden_dir(), "expected", "c_t1_a1.out")))
})
