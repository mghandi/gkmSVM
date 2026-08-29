# legacyNorm = TRUE reproduces the scores of versions up to 0.80 exactly (the c_*_legacy golden
# files are byte-identical to the Phase 0 freeze of the unmodified code).

test_that("legacyNorm reproduces pre-0.90 classify scores and defaults to off", {
  skip_if_no_golden()
  fx <- function(f) file.path(golden_dir(), "fixtures", f)
  gold <- function(f) file.path(golden_dir(), "expected", f)
  o_new <- tempfile(); o_old <- tempfile()
  capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = o_new, svseqfile = fx("sv_svseq.fa"), alphafile = fx("sv_svalpha.out")))
  capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = o_old, svseqfile = fx("sv_svseq.fa"), alphafile = fx("sv_svalpha.out"), legacyNorm = TRUE))
  expect_identical(read_bytes(o_new), read_bytes(gold("c_t1_a2.out")))
  expect_identical(read_bytes(o_old), read_bytes(gold("c_t1_a2_legacy.out")))
  expect_false(identical(read_bytes(o_new), read_bytes(o_old)))
  d <- tempfile()   # gkmsvm_delta passes it through
  capture.output(gkmsvm_delta(fx("test_seqs.fa"), fx("test_seqs.fa"), svmfnprfx = NA, outfile = d, svseqfile = fx("sv_svseq.fa"), alphafile = fx("sv_svalpha.out"), legacyNorm = TRUE))
  expect_true(all(utils::read.table(d)[[2]] == 0))
})
