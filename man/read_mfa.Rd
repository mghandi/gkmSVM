\name{read_mfa}
\alias{read_mfa}
\alias{write_mfa}
\title{Read and write multi-track FASTA files}
\description{Multi-track FASTA (".mfa") is the input format of \code{\link{gkmsvm_kernel}}, \code{\link{gkmsvm_classify}} and \code{\link{gkmsvm_train}} when \code{alphabets} lists two or more alphabets: every record is a header line (\code{>name}) followed by one line per track, all of the same length -- for example a DNA line and a 0/1 methylation line. \code{read_mfa} reads such a file into a named list with one character vector (one string per track) per record; \code{write_mfa} writes one.}
\usage{read_mfa(file, ntracks)
write_mfa(records, file, names = NULL)}
\arguments{
  \item{file}{file name}
  \item{ntracks}{number of tracks per record (the length of \code{alphabets})}
  \item{records}{a list of character vectors, one string per track}
  \item{names}{record names; by default \code{names(records)}}
}
\value{\code{read_mfa}: a named list of character vectors. \code{write_mfa}: the file name, invisibly.}
\examples{
f <- tempfile(fileext = ".mfa")
write_mfa(list(s1 = c("ACGTACGT", "01100101"), s2 = c("TTGACTCA", "00000000")), f)
read_mfa(f, 2)
}
