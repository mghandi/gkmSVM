# Phase 6: tiled mismatch profile (bounded memory) gives byte-identical kernels.

fx <- function(f) file.path(golden_dir(), "fixtures", f)
gold <- function(f) file.path(golden_dir(), "expected", f)

test_that("tileRows does not change the kernel", {
  skip_if_no_golden()
  for (r in c(1, 5, 24)) {
    out <- tempfile()
    capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), out, tileRows = r))
    expect_identical(read_bytes(out), read_bytes(gold("k_small_t1_a2.out")), info = paste("tileRows", r))
  }
  out <- tempfile()
  capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), out, usePseudocnt = TRUE, tileRows = 7))
  expect_identical(read_bytes(out), read_bytes(gold("k_small_pseudo_a2.out")))
  out <- tempfile()   # a tiny memory budget forces many tiles
  capture.output(gkmsvm_kernel(fx("dna_medium_pos.fa"), fx("dna_medium_neg.fa"), out, tileMemoryMB = 1))
  expect_identical(read_bytes(out), read_bytes(gold("k_medium_T1.out")))
})
