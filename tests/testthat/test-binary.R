# Phase 4: binary kernel (.gkmk) and single-file model (.gkmmodel).

fx <- function(f) file.path(golden_dir(), "fixtures", f)
gold <- function(f) file.path(golden_dir(), "expected", f)
plain <- function(m) matrix(as.numeric(m), nrow(m))   # drop dimnames and the npos/index/params attributes

test_that("read_gkm_kernel auto-detects text and binary and returns the same matrix", {
  skip_if_no_golden()
  txt <- read_gkm_kernel(gold("k_small_t1_a2.out"))      # frozen text kernel (+ .out.index sidecar)
  bin <- read_gkm_kernel(gold("k_small_binary.out"))     # the same kernel written with -b
  expect_false(is_gkm_binary(gold("k_small_t1_a2.out")))
  expect_true(is_gkm_binary(gold("k_small_binary.out")))
  expect_equal(dim(txt), c(24, 24)); expect_equal(dim(bin), c(24, 24))
  expect_true(isSymmetric(unname(txt))); expect_true(isSymmetric(unname(bin)))
  expect_equal(unname(diag(bin)), rep(1, 24))
  expect_equal(plain(bin), plain(txt), tolerance = 1e-6)   # float32 vs "%e": ~7 significant digits
  expect_equal(attr(bin, "npos"), 12)
  expect_equal(attr(bin, "params")$L, 10); expect_equal(attr(bin, "params")$alphabet, "ACGT")
  expect_equal(attr(bin, "index")$id, c(paste0("pos", 1:12), paste0("neg", 1:12)))
  expect_equal(attr(bin, "index")$name, attr(txt, "index")$name)
})

test_that("write_gkm_kernel round-trips (binary float32/float64 and text)", {
  skip_if_no_golden()
  bin <- read_gkm_kernel(gold("k_small_binary.out"))
  f32 <- tempfile(); f64 <- tempfile(); ft <- tempfile()
  write_gkm_kernel(bin, f32, "binary", "float32")
  write_gkm_kernel(bin, f64, "binary", "float64")
  write_gkm_kernel(bin, ft, "text")
  # the R writer reproduces the C++ writer byte for byte (same header, table, payload, crc)
  expect_identical(read_bytes(f32), read_bytes(gold("k_small_binary.out")))
  expect_equal(plain(read_gkm_kernel(f64)), plain(bin))
  expect_equal(plain(read_gkm_kernel(ft)), plain(bin), tolerance = 1e-6)
  expect_true(file.exists(paste0(ft, ".index")))
  expect_equal(attr(read_gkm_kernel(ft), "index")$name, attr(bin, "index")$name)
})

test_that("a corrupted binary payload is detected", {
  skip_if_no_golden()
  b <- read_bytes(gold("k_small_binary.out"))
  b[length(b) - 10] <- xor(b[length(b) - 10], as.raw(0xFF))
  f <- tempfile(); writeBin(b, f)
  expect_error(read_gkm_kernel(f), "checksum")
  expect_silent(read_gkm_kernel(f, verify = FALSE))
})

test_that("gkmsvm_kernel(format='binary') and gkmsvm_train on it give the same model as text", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  kt <- tempfile(); kb <- tempfile()
  capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), kt))
  capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), kb, format = "binary"))
  expect_identical(read_bytes(kb), read_bytes(gold("k_small_binary.out")))
  pt <- tempfile("mt"); pb <- tempfile("mb")
  capture.output(gkmsvm_train(kt, fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), pt, C = 1))
  capture.output(gkmsvm_train(kb, fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), pb, C = 1))
  at <- utils::read.delim(paste0(pt, "_svalpha.out"), header = FALSE); ab <- utils::read.delim(paste0(pb, "_svalpha.out"), header = FALSE)
  expect_equal(at[[1]], ab[[1]])
  expect_equal(at[[2]], ab[[2]], tolerance = 1e-4)   # float32 kernel -> alphas agree to the solver tolerance
})

test_that("the .gkmmodel single-file model scores like the legacy pair minus rho", {
  skip_if_no_golden()
  skip_if_not_installed("kernlab")
  kt <- tempfile(); prefix <- tempfile("m")
  capture.output(gkmsvm_kernel(fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), kt))
  capture.output(gkmsvm_train(kt, fx("dna_small_pos.fa"), fx("dna_small_neg.fa"), prefix, C = 1))
  model <- paste0(prefix, ".gkmmodel")
  expect_true(file.exists(model))
  hdr <- readLines(model, n = 2)
  expect_equal(hdr[1], "#gkmmodel 1")
  rho <- as.numeric(sub("^#rho ", "", hdr[2]))
  o_new <- tempfile(); o_old <- tempfile()
  capture.output(gkmsvm_classify(fx("test_seqs.fa"), prefix, o_new))                        # picks the .gkmmodel
  capture.output(gkmsvm_classify(fx("test_seqs.fa"), svmfnprfx = NA, outfile = o_old,
                                 svseqfile = paste0(prefix, "_svseq.fa"), alphafile = paste0(prefix, "_svalpha.out")))
  s_new <- utils::read.delim(o_new, header = FALSE); s_old <- utils::read.delim(o_old, header = FALSE)
  expect_equal(s_new[[1]], s_old[[1]])
  expect_equal(s_new[[2]], s_old[[2]] - rho, tolerance = 2e-6)
})
