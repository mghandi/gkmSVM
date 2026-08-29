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

# Phase 4b: the libsvm backend (vendored LIBSVM, precomputed kernel) vs kernlab, and vs the CLI.

test_that("backend='libsvm' agrees with kernlab on decision values (same solver family)", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  fx <- function(f) file.path(golden_dir(), "fixtures", f)
  kern <- file.path(golden_dir(), "expected", "k_medium_T1.out")
  pk <- tempfile("kl"); pl <- tempfile("ls")
  capture.output(gkmsvm_train(kern, fx("dna_medium_pos.fa"), fx("dna_medium_neg.fa"), pk, C = 1, shrinking = FALSE))
  capture.output(gkmsvm_train(kern, fx("dna_medium_pos.fa"), fx("dna_medium_neg.fa"), pl, C = 1, backend = "libsvm", shrinking = FALSE))
  expect_true(file.exists(paste0(pl, ".gkmmodel")))
  # score the same test sequences with both models (legacy pair: bias-free scores on both sides)
  score <- function(prefix) {
    out <- tempfile()
    capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = out,
                                   svseqfile = paste0(prefix, "_svseq.fa"), alphafile = paste0(prefix, "_svalpha.out")))
    utils::read.delim(out, header = FALSE)[[2]]
  }
  sk <- score(pk); sl <- score(pl)
  expect_gt(cor(sk, sl), 0.99999)
  expect_lt(max(abs(sk - sl)), 1e-3 * max(abs(sk)))
  # and the biases agree (kernlab's b is libsvm's rho)
  rho <- as.numeric(sub("^#rho ", "", readLines(paste0(pl, ".gkmmodel"), n = 2)[2]))
  # (kernlab's b is not exposed in the legacy files; compare it through a fresh fit)
  mat <- read_gkm_kernel(kern); y <- c(rep(1, 40), rep(0, 40))
  svp <- kernlab::ksvm(kernlab::as.kernelMatrix(mat), y, type = "C-svc", C = 1, shrinking = FALSE)
  expect_equal(rho, svp@b, tolerance = 1e-2)
})

test_that("the R libsvm backend and the CLI produce the same model", {
  skip_if_no_golden()
  fx <- function(f) file.path(golden_dir(), "fixtures", f)
  kern <- file.path(golden_dir(), "expected", "k_small_t1_a2.out")
  prefix <- tempfile("r")
  capture.output(gkmsvm_train(kern, fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), prefix, C = 1, backend = "libsvm"))
  got <- utils::read.delim(paste0(prefix, "_svalpha.out"), header = FALSE, stringsAsFactors = FALSE)
  exp <- utils::read.delim(file.path(golden_dir(), "expected", "t_small.out"), header = FALSE, stringsAsFactors = FALSE)
  expect_equal(got[[1]], exp[[1]])
  expect_equal(got[[2]], exp[[2]], tolerance = 1e-5)
})
