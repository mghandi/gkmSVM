#include <Rcpp.h>
#include "crc32.h"

// CRC-32 of a raw vector, as a signed 32-bit integer (R has no unsigned type); used by the .gkmk
// reader/writer in R/gkm_kernel_io.R to match the checksum written by src/KernelFile.h.
// [[Rcpp::export(name = ".gkm_crc32_cpp")]]
int gkm_crc32_cpp(Rcpp::RawVector x) {
  uint32_t c = gkm_crc32(x.size() ? &x[0] : NULL, (size_t)x.size());
  return (int)c;
}
