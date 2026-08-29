# gkmsvm_train (kernlab C-SVC on a frozen kernel) and the round trip through gkmsvm_classify.
# Alphas are not compared to a golden file (dual degeneracy, kernlab version drift): the test pins
# the output *format* and checks that the model ranks the motif-carrying test sequences higher.

test_that("gkmsvm_train writes the svalpha/svseq pair in the documented format", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  fx <- file.path(golden_dir(), "fixtures")
  kern <- file.path(golden_dir(), "expected", "k_small_t1_a2.out")
  prefix <- tempfile("model")
  capture.output(gkmsvm_train(kern, file.path(fx, "dna_small_pos.fa"), file.path(fx, "dna_small_neg.fa"), prefix, C = 1))
  alpha <- utils::read.delim(paste0(prefix, "_svalpha.out"), header = FALSE, stringsAsFactors = FALSE)
  expect_equal(ncol(alpha), 2)
  expect_true(all(grepl("^(pos|neg)[0-9]+$", alpha[[1]])))
  expect_true(all(alpha[[2]][grepl("^neg", alpha[[1]])] < 0))
  expect_true(all(alpha[[2]][grepl("^pos", alpha[[1]])] > 0))
  sv <- seqinr::read.fasta(paste0(prefix, "_svseq.fa"))
  expect_equal(names(sv), alpha[[1]])
  # classify the training sequences with the trained model: positives must score higher on average
  out <- tempfile()
  capture.output(gkmsvm_classify(file.path(fx, "dna_small_pos.fa"), prefix, out))
  sp <- utils::read.delim(out, header = FALSE)[[2]]
  capture.output(gkmsvm_classify(file.path(fx, "dna_small_neg.fa"), prefix, out))
  sn <- utils::read.delim(out, header = FALSE)[[2]]
  expect_gt(mean(sp), mean(sn))
})
