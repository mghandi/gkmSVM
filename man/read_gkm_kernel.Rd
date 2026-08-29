\name{read_gkm_kernel}
\alias{read_gkm_kernel}
\alias{write_gkm_kernel}
\alias{is_gkm_binary}
\title{Read and write gkmSVM kernel matrices}
\description{
  \code{read_gkm_kernel} reads a kernel matrix written by \code{\link{gkmsvm_kernel}} in either the
  text format or the binary \code{.gkmk} format (auto-detected) and returns the full symmetric matrix.
  \code{write_gkm_kernel} writes a matrix in either format. \code{is_gkm_binary} tests the format.
}
\usage{
read_gkm_kernel(file, verify = TRUE)
write_gkm_kernel(x, file, format = c("binary", "text"), dtype = c("float32", "float64"),
                 npos = attr(x, "npos"), index = attr(x, "index"), params = attr(x, "params"))
is_gkm_binary(file)
}
\arguments{
\item{file}{kernel file name}
\item{verify}{check the CRC-32 of a binary payload (default TRUE)}
\item{x}{a square (symmetric) kernel matrix}
\item{format}{\code{"binary"} (\code{.gkmk} v1) or \code{"text"} (the historical format)}
\item{dtype}{storage type of the binary payload; \code{"float32"} keeps the ~7 significant digits the text format has}
\item{npos}{number of positive sequences (rows come first)}
\item{index}{optional row-identity table (columns \code{row, id, name, label, length, nlmers}), as written by \code{gkmsvm_kernel} to \code{<file>.index}}
\item{params}{optional list of the kernel parameters stored in a binary file (\code{L, K, maxnmm, useTgkm, b, addRC, usePseudocnt, alphabet})}
}
\value{
  \code{read_gkm_kernel} returns an n x n numeric matrix with attributes \code{"npos"}, \code{"index"}
  (a data.frame, or NULL when no identity information is available) and, for binary files, \code{"params"}.
}
\details{
  The binary layout is documented in \code{src/KernelFile.h} of the package sources: a \code{"GKMK"}
  magic, a version byte, the data type, the number of sequences, the parameters, the alphabet, the
  row-identity table, the lower triangle of the matrix (row-major, including the diagonal) and a
  CRC-32 of the payload. All integers are little-endian regardless of the host.
}
\author{Mahmoud Ghandi}
