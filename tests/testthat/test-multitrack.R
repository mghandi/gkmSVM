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

test_that("multi-track classify: the .gkmmodel (alphabets recorded) and the legacy pair match the frozen scores", {
  skip_if_no_golden()
  o <- tempfile()
  capture.output(gkmsvm_classify(fx("test_meth.mfa"), NA, o, L = 6, K = 4, maxnmm = -1, useTgkm = 0, svseqfile = fx("sv_meth.gkmmodel"), alphafile = fx("sv_meth.gkmmodel")))
  expect_identical(read_bytes(o), read_bytes(gold("c_mt_model_t0.out")))
  capture.output(gkmsvm_classify(fx("test_meth.mfa"), NA, o, L = 6, K = 4, maxnmm = -1, useTgkm = 0, addRC = FALSE, svseqfile = fx("sv_meth_svseq.mfa"), alphafile = fx("sv_meth_svalpha.out"), alphabets = c("dna", "01")))
  expect_identical(read_bytes(o), read_bytes(gold("c_mt_pair_t0_noRC.out")))
  expect_error(capture.output(gkmsvm_classify(fx("test_meth.mfa"), NA, o, L = 6, K = 4, svseqfile = fx("sv_meth.gkmmodel"), alphafile = fx("sv_meth.gkmmodel"), alphabets = c("dna", "NUM"))))
})

test_that("multi-track train (kernlab and libsvm) writes multi-track models that classify consistently", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  kern <- gold("k_mt_meth_t0.out")
  pk <- tempfile("kl"); pl <- tempfile("ls")
  capture.output(gkmsvm_train(kern, fx("meth_pos.mfa"), fx("meth_neg.mfa"), pk, C = 1, shrinking = FALSE, alphabets = c("dna", "01")))
  capture.output(gkmsvm_train(kern, fx("meth_pos.mfa"), fx("meth_neg.mfa"), pl, C = 1, shrinking = FALSE, backend = "libsvm", alphabets = c("dna", "01")))
  for (p in c(pk, pl)) {
    m <- readLines(paste0(p, ".gkmmodel"))
    expect_true("#alphabets dna,=01" %in% m)
    sv <- read_mfa(paste0(p, "_svseq.fa"), 2)
    expect_true(length(sv) > 0)
    expect_true(all(nchar(sapply(sv, `[`, 1)) == 40))
  }
  score <- function(prefix, f) { o <- tempfile(); capture.output(gkmsvm_classify(fx(f), prefix, o, L = 6, K = 4, maxnmm = -1, useTgkm = 0)); utils::read.delim(o, header = FALSE)[[2]] }
  sp <- score(pk, "meth_pos.mfa"); sn <- score(pk, "meth_neg.mfa")
  expect_gt(min(sp), max(sn))     # the synthetic example is separable
  expect_equal(score(pl, "meth_pos.mfa"), sp, tolerance = 1e-3)
  expect_equal(score(pl, "meth_neg.mfa"), sn, tolerance = 1e-3)
  # the R and the C++ readers agree on the multi-track fixtures
  r <- read_mfa(fx("meth_pos.mfa"), 2)
  expect_equal(length(r), 20); expect_equal(names(r)[1], "pos1"); expect_equal(nchar(r[[1]][2]), 40)
  f <- tempfile(); write_mfa(r, f); expect_identical(read_mfa(f, 2), r)
})
