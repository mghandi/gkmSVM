# Phase 3: sequence identity is positional; names are metadata (dev/REFACTORING_PLAN.md, Phase 3).

fx <- function(f) file.path(golden_dir(), "fixtures", f)

test_that("one FASTA record is one kernel row by default; mergeByName=TRUE restores merging", {
  skip_if_no_golden()
  o1 <- tempfile(); o2 <- tempfile()
  capture.output(gkmsvm_kernel(fx("dupnames_pos.fa"), fx("dna_small_neg.fa"), o1))
  capture.output(gkmsvm_kernel(fx("dupnames_pos.fa"), fx("dna_small_neg.fa"), o2, mergeByName = TRUE))
  expect_equal(length(readLines(o1)), 24)
  expect_equal(length(readLines(o2)), 23)
  expect_identical(read_bytes(o1), read_bytes(file.path(golden_dir(), "expected", "k_dupnames.out")))
  expect_identical(read_bytes(o2), read_bytes(file.path(golden_dir(), "expected", "k_dupnames_merge.out")))
})

test_that("the .index sidecar carries row identity", {
  skip_if_no_golden()
  out <- tempfile()
  capture.output(gkmsvm_kernel(fx("names_edge_pos.fa"), fx("dna_small_neg.fa"), out))
  idx <- gkmSVM:::.gkm_read_index(out)
  expect_equal(nrow(idx), 20)
  expect_equal(idx$row, 0:19)
  expect_equal(idx$id, c(paste0("pos", 1:8), paste0("neg", 1:12)))
  expect_equal(idx$label, c(rep(1, 8), rep(-1, 12)))
  expect_equal(idx$name[1:6], c("", "a", "séq_üñ", strrep("X", 200), "neg0", "neg0"))
  expect_true(all(idx$length == 100))
  expect_true(all(idx$nlmers == 2 * (100 - 10 + 1)))
})

test_that("duplicate names across the positive and negative files are accepted", {
  skip_if_no_golden()
  out <- tempfile()
  expect_message(capture.output(gkmsvm_kernel(fx("names_edge_pos.fa"), fx("dna_small_neg.fa"), out)), "not unique")
  expect_equal(length(readLines(out)), 20)
})

test_that("gkmsvm_train uses the sidecar and falls back to row ids when names are not unique", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  kern <- tempfile()
  capture.output(gkmsvm_kernel(fx("dupnames_pos.fa"), fx("dna_small_neg.fa"), kern))
  prefix <- tempfile("dupmodel")
  expect_message(capture.output(gkmsvm_train(kern, fx("dupnames_pos.fa"), fx("dna_small_neg.fa"), prefix, C = 1)), "row ids")
  alpha <- utils::read.delim(paste0(prefix, "_svalpha.out"), header = FALSE, stringsAsFactors = FALSE)
  expect_true(all(grepl("^(pos|neg)[0-9]+$", alpha[[1]])))
  expect_false(anyDuplicated(alpha[[1]]) > 0)
  sv <- seqinr::read.fasta(paste0(prefix, "_svseq.fa"))
  expect_equal(names(sv), alpha[[1]])
  # and the model is usable by gkmsvm_classify (name matching on the unique ids)
  out <- tempfile()
  capture.output(gkmsvm_classify(fx("test_seqs.fa"), prefix, out))
  expect_equal(nrow(utils::read.delim(out, header = FALSE)), 10)
})

test_that("gkmsvm_train still works on a kernel without a sidecar (older versions)", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  kern <- tempfile()
  file.copy(file.path(golden_dir(), "expected", "k_small_t1_a2.out"), kern)
  expect_false(file.exists(paste0(kern, ".index")))
  prefix <- tempfile("nosidecar")
  capture.output(gkmsvm_train(kern, fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), prefix, C = 1))
  expect_true(file.exists(paste0(prefix, "_svalpha.out")))
})

test_that("a model whose alpha file lists a missing or duplicated support vector is rejected", {
  skip_if_no_golden()
  out <- tempfile()
  expect_error(capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = out,
                                              svseqfile = fx("sv_svseq.fa"), alphafile = fx("sv_missing_svalpha.out"))), "failed")
  expect_error(capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = out,
                                              svseqfile = fx("sv_svseq.fa"), alphafile = fx("sv_dupalpha_svalpha.out"))), "failed")
})
