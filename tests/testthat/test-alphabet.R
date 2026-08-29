# Phase 5: alphabets other than DNA (b up to 32) through the R interface.

fx <- function(f) file.path(golden_dir(), "fixtures", f)
gold <- function(f) file.path(golden_dir(), "expected", f)

test_that("protein (b=20) kernels run and match the frozen output; the alphabet does not leak", {
  skip_if_no_golden()
  o1 <- tempfile(); o2 <- tempfile(); o3 <- tempfile()
  capture.output(gkmsvm_kernel(fx("protein_pos.fa"), fx("protein_neg.fa"), o1, L = 6, K = 4, maxnmm = -1, alphabetFN = "protein"))
  capture.output(gkmsvm_kernel(fx("protein_pos.fa"), fx("protein_neg.fa"), o2, L = 6, K = 4, maxnmm = -1, alphabetFN = fx("alphabet_protein.txt")))
  capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), o3))
  expect_identical(read_bytes(o1), read_bytes(gold("k_protein_t1.out")))
  expect_identical(read_bytes(o2), read_bytes(gold("k_protein_t1.out")))
  expect_identical(read_bytes(o3), read_bytes(gold("k_small_t1_a2.out")))
  expect_equal(attr(read_gkm_kernel(o1), "index")$nlmers[1], 60 - 6 + 1)   # no reverse complement for protein
})

test_that("RNA: complement pairs from a file or the keyword give the DNA result", {
  skip_if_no_golden()
  o1 <- tempfile(); o2 <- tempfile()
  capture.output(gkmsvm_kernel(fx("rna_pos.fa"), fx("rna_neg.fa"), o1, alphabetFN = "rna"))
  capture.output(gkmsvm_kernel(fx("rna_pos.fa"), fx("rna_neg.fa"), o2, alphabetFN = fx("alphabet_rna_pairs.txt")))
  expect_identical(read_bytes(o1), read_bytes(o2))
  dna <- read_gkm_kernel(gold("k_small_t1_a2.out")); rna <- read_gkm_kernel(o1)
  expect_equal(matrix(as.numeric(rna), 12)[1:6, 1:6], matrix(as.numeric(dna), 24)[1:6, 1:6])
  expect_equal(matrix(as.numeric(rna), 12)[7:12, 7:12], matrix(as.numeric(dna), 24)[13:18, 13:18])
})

test_that("a 5-letter alphabet runs (it used to segfault) and partial complement pairs are rejected", {
  skip_if_no_golden()
  out <- tempfile()
  capture.output(gkmsvm_kernel(fx("b5_pos.fa"), fx("b5_neg.fa"), out, L = 8, K = 5, maxnmm = -1, alphabetFN = fx("alphabet_b5.txt")))
  expect_identical(read_bytes(out), read_bytes(gold("k_b5_t1.out")))
  expect_error(capture.output(gkmsvm_kernel(fx("b5_pos.fa"), fx("b5_neg.fa"), tempfile(), L = 8, K = 5, maxnmm = -1,
                                            alphabetFN = fx("alphabet_b5_pairs_bad.txt"))), "failed")
})

test_that("classify works with a protein model", {
  skip_if_no_golden()
  out <- tempfile()
  capture.output(gkmsvm_classify(fx("protein_pos.fa"), svmfnprfx = NA, outfile = out, L = 6, K = 4, maxnmm = -1,
                                 alphabetFN = "protein", svseqfile = fx("protein_neg.fa"), alphafile = fx("sv_protein_svalpha.out")))
  expect_identical(read_bytes(out), read_bytes(gold("c_protein_t1.out")))
})
