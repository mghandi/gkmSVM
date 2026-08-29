# Phase 7: multi-track FASTA (".mfa") -- a header line followed by one line per track, all of the
# same length (e.g. a DNA line and a 0/1 methylation line). This is the input format of
# gkmsvm_kernel / gkmsvm_classify / gkmsvm_train when `alphabets` has two or more entries.

# Reads a multi-track FASTA file into a named list: one character vector of `ntracks` strings per
# record (names = the first token of the header). Blank lines and lines starting with "#" or ";"
# are skipped.
read_mfa <- function(file, ntracks) {
  ntracks <- as.integer(ntracks)
  if (length(ntracks) != 1 || is.na(ntracks) || ntracks < 1) stop("'ntracks' must be a positive integer")
  ln <- readLines(file, warn = FALSE)
  ln <- sub("[\r[:space:]]+$", "", ln)
  ln <- ln[nzchar(ln) & !startsWith(ln, "#") & !startsWith(ln, ";")]
  hdr <- which(startsWith(ln, ">"))
  if (length(hdr) == 0) return(structure(list(), names = character(0)))
  ends <- c(hdr[-1] - 1, length(ln))
  recs <- vector("list", length(hdr))
  nm <- character(length(hdr))
  for (i in seq_along(hdr)) {
    body <- if (ends[i] >= hdr[i] + 1) ln[(hdr[i] + 1):ends[i]] else character(0)
    nm[i] <- sub("[[:space:]].*$", "", substring(ln[hdr[i]], 2))
    if (length(body) != ntracks) stop(sprintf("record '%s' has %d track line(s); expected %d", nm[i], length(body), ntracks))
    if (length(unique(nchar(body))) != 1) stop(sprintf("record '%s': the tracks have different lengths", nm[i]))
    recs[[i]] <- body
  }
  names(recs) <- nm
  recs
}

# Writes records (a list of character vectors, one string per track) as multi-track FASTA.
write_mfa <- function(records, file, names = NULL) {
  if (is.null(names)) names <- names(records)
  if (is.null(names) || length(names) != length(records)) stop("one name per record is required")
  con <- file(file, "w"); on.exit(close(con))
  for (i in seq_along(records)) writeLines(c(paste0(">", names[i]), as.character(records[[i]])), con)
  invisible(file)
}
