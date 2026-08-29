# Phase 7: multi-track input (gkm-SVM over heterogeneous alphabets) through the R interface.

fx <- function(f) file.path(golden_dir(), "fixtures", f)
gold <- function(f) file.path(golden_dir(), "expected", f)

test_that("a DNA + methylation two-track kernel matches the frozen output (alphabets = c('dna','01'))", {
  skip_if_no_golden()
  o1 <- tempfile(); o2 <- tempfile()
  capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), o1, L = 6, K = 4, maxnmm = -1, useTgkm = 0, alphabets = c("dna", "01")))
  expect_identical(read_bytes(o1), read_bytes(gold("k_mt_meth_t0.out")))
  # the same through the explicit spec string, and with the gkm counts (t=2)
  capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), o2, L = 6, K = 4, maxnmm = -1, useTgkm = 2, alphabetFN = "dna,=01"))
  expect_identical(read_bytes(o2), read_bytes(gold("k_mt_meth_t2.out")))
  idx <- attr(read_gkm_kernel(o1), "index")
  expect_equal(nrow(idx), 40)
  expect_equal(idx$nlmers[1], 2 * (40 - 6 + 1))   # both strands
})

test_that("binary output and a bounded mismatch count agree with the text/golden results", {
  skip_if_no_golden()
  ob <- tempfile(fileext = ".gkmk"); ot <- tempfile()
  capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), ob, L = 6, K = 4, maxnmm = -1, useTgkm = 0, alphabets = c("dna", "01"), format = "binary"))
  kb <- read_gkm_kernel(ob); kt <- read_gkm_kernel(gold("k_mt_meth_t0.out"))
  plain <- function(k) matrix(as.numeric(k), nrow(k))
  expect_equal(plain(kb), plain(kt), tolerance = 1e-6)   # float32 vs "%e": ~7 significant digits
  expect_equal(attr(kb, "params")$alphabet, "dna,=01")
  capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), ot, L = 6, K = 4, maxnmm = 4, useTgkm = 0, addRC = FALSE, alphabets = c("dna", "01")))
  expect_identical(read_bytes(ot), read_bytes(gold("k_mt_meth_t0_d4_noRC.out")))
})

test_that("invalid multi-track requests are errors, not crashes", {
  skip_if_no_golden()
  o <- tempfile()
  expect_error(capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), o, L = 6, K = 4, alphabets = c("dna", "01", "abc"))))
  expect_error(capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("meth_neg.mfa"), o, L = 6, K = 4, alphabets = c("dna", "01"))))
  expect_error(capture.output(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), o, L = 6, K = 4, useTgkm = 3, alphabets = c("dna", "01"))))
  expect_error(gkmsvm_kernel(fx("meth_pos.mfa"), fx("meth_neg.mfa"), o, alphabets = c("dna", "01"), alphabetFN = "dna"))
  expect_false(file.exists(o))
})
