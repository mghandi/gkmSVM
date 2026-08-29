# Kernel-file I/O (Phase 4 of dev/REFACTORING_PLAN.md).
#
# Two on-disk formats, auto-detected by read_gkm_kernel():
#   * text  : the historical format - lower triangle including the diagonal, tab separated, one row
#             per sequence, "%e" values, "1.0" on the diagonal; row identity in <file>.index (Phase 3)
#   * binary: ".gkmk" v1 (magic "GKMK"), documented in src/KernelFile.h - float32 by default,
#             carries n, npos, the parameters, the alphabet and the row-identity table, crc32-checked.
# Both return the full symmetric n x n matrix with attributes "npos" and "index" (a data.frame with
# columns row, id, name, label, length, nlmers, or NULL when unknown) and, for binary files, "params".

is_gkm_binary <- function(file) {
  con <- file(file, "rb"); on.exit(close(con))
  identical(readBin(con, "raw", 4), charToRaw("GKMK"))
}

read_gkm_kernel <- function(file, verify = TRUE) {
  if (is_gkm_binary(file)) .gkm_read_binary(file, verify) else .gkm_read_text(file)
}

.gkm_read_text <- function(file) {
  idx <- .gkm_read_index(file)
  n <- if (!is.null(idx)) nrow(idx) else length(readLines(file, warn = FALSE))
  mat <- data.matrix(utils::read.table(file, fill = TRUE, col.names = paste("V", 1:n)))
  if (nrow(mat) != n) stop(sprintf("%s: %d rows read but %d expected", file, nrow(mat), n))
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  rownames(mat) <- colnames(mat)
  attr(mat, "npos") <- if (!is.null(idx)) sum(idx$label > 0) else NA_integer_
  attr(mat, "index") <- idx
  mat
}

.gkm_read_binary <- function(file, verify = TRUE) {
  con <- file(file, "rb"); on.exit(close(con))
  rd <- function(what, n = 1L, size = NA_integer_, signed = TRUE) readBin(con, what, n = n, size = size, endian = "little", signed = signed)
  rdstr <- function() { l <- rd("integer", size = 4); if (l == 0) "" else rawToChar(rd("raw", l)) }
  if (!identical(rd("raw", 4), charToRaw("GKMK"))) stop(file, ": not a gkmSVM binary kernel file")
  version <- rd("integer", size = 1, signed = FALSE)
  if (version != 1) stop(file, ": unsupported .gkmk version ", version)
  dtype <- rd("integer", size = 1, signed = FALSE); layout <- rd("integer", size = 1, signed = FALSE); flags <- rd("integer", size = 1, signed = FALSE)
  if (layout != 0) stop(file, ": unsupported layout ", layout)
  n <- rd("integer", size = 4); npos <- rd("integer", size = 4)
  prov <- rd("integer", n = 7, size = 4)
  params <- list(L = prov[1], K = prov[2], maxnmm = prov[3], useTgkm = prov[4], b = prov[5], addRC = prov[6] == 1, usePseudocnt = prov[7] == 1,
                 alphabet = rdstr())
  idx <- NULL
  if (bitwAnd(flags, 1L) == 1L) {
    idx <- data.frame(row = integer(n), id = character(n), name = character(n), label = integer(n), length = numeric(n), nlmers = numeric(n),
                      stringsAsFactors = FALSE)
    for (i in seq_len(n)) {
      idx$row[i] <- i - 1L; idx$id[i] <- rdstr(); idx$name[i] <- rdstr()
      idx$label[i] <- rd("integer", size = 1)
      idx$length[i] <- .gkm_read_int64(con); idx$nlmers[i] <- .gkm_read_int64(con)
    }
  }
  m <- n * (n + 1) / 2
  size <- if (dtype == 0) 4L else 8L
  payload <- rd("raw", m * size)
  crc <- rd("integer", size = 4)
  if (verify) {
    got <- .gkm_crc32_cpp(payload)
    if (got != crc) stop(file, ": payload checksum mismatch (file corrupted or truncated)")
  }
  vals <- readBin(payload, "numeric", n = m, size = size, endian = "little")
  mat <- matrix(0, n, n)
  mat[upper.tri(mat, diag = TRUE)] <- vals   # row-major lower triangle == column-major upper triangle
  mat <- t(mat)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  dimnames(mat) <- list(paste("V", 1:n), paste("V", 1:n))
  attr(mat, "npos") <- npos
  attr(mat, "index") <- idx
  attr(mat, "params") <- params
  mat
}

# int64 as two little-endian uint32 halves (R has no 64-bit integer; values here fit in a double)
.gkm_read_int64 <- function(con) {
  lo <- readBin(con, "integer", size = 4, endian = "little"); hi <- readBin(con, "integer", size = 4, endian = "little")
  if (lo < 0) lo <- lo + 2^32
  hi * 2^32 + lo
}

# Write a (symmetric) kernel matrix. format "binary" writes .gkmk v1 (dtype float32 or float64);
# format "text" writes the historical text format (and <file>.index when an index is supplied).
write_gkm_kernel <- function(x, file, format = c("binary", "text"), dtype = c("float32", "float64"),
                             npos = attr(x, "npos"), index = attr(x, "index"), params = attr(x, "params")) {
  format <- match.arg(format); dtype <- match.arg(dtype)
  n <- nrow(x)
  if (is.null(n) || n != ncol(x)) stop("x must be a square matrix")
  if (is.null(npos) || is.na(npos)) npos <- if (!is.null(index)) sum(index$label > 0) else 0L
  vals <- t(x)[upper.tri(x, diag = TRUE)]           # row-major lower triangle
  if (format == "text") {
    con <- file(file, "w"); on.exit(close(con))
    k <- 1
    for (i in seq_len(n)) {
      row <- vals[k:(k + i - 1)]; k <- k + i
      writeLines(paste0(paste(c(sprintf("%e", row[-i]), "1.0"), collapse = "\t"), "\t"), con)
    }
    if (!is.null(index)) utils::write.table(index, paste0(file, ".index"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    return(invisible(file))
  }
  if (is.null(params)) params <- list(L = 0L, K = 0L, maxnmm = 0L, useTgkm = 0L, b = 0L, addRC = FALSE, usePseudocnt = FALSE, alphabet = "")
  con <- file(file, "wb"); on.exit(close(con))
  wr <- function(what, size) writeBin(what, con, size = size, endian = "little")
  wrstr <- function(s) { r <- charToRaw(enc2utf8(s)); wr(length(r), 4L); if (length(r)) writeBin(r, con) }
  writeBin(charToRaw("GKMK"), con)
  writeBin(as.raw(c(1L, if (dtype == "float32") 0L else 1L, 0L, if (is.null(index)) 0L else 1L)), con)
  wr(as.integer(c(n, npos)), 4L)
  wr(as.integer(c(params$L, params$K, params$maxnmm, params$useTgkm, params$b, as.integer(params$addRC), as.integer(params$usePseudocnt))), 4L)
  wrstr(params$alphabet)
  if (!is.null(index)) {
    for (i in seq_len(n)) {
      wrstr(index$id[i]); wrstr(index$name[i]); writeBin(as.raw(bitwAnd(as.integer(index$label[i]), 255L)), con)
      for (v in c(index$length[i], index$nlmers[i])) { wr(as.integer(v %% 2^32 - (v %% 2^32 >= 2^31) * 2^32), 4L); wr(as.integer(v %/% 2^32), 4L) }
    }
  }
  payload <- writeBin(as.numeric(vals), raw(), size = if (dtype == "float32") 4L else 8L, endian = "little")
  writeBin(payload, con)
  wr(.gkm_crc32_cpp(payload), 4L)
  invisible(file)
}
