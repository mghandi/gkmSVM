#include <Rcpp.h>
#include "gkmOptions.h"
using namespace Rcpp;

// [[Rcpp::export]]
void gkmsvm_kernel(SEXP params){
  Rcpp::List rparam(params); // Get parameters in params.

  OptsGkmKernel opt;
  opt.L = Rcpp::as<int>(rparam["L"]);
  opt.K = Rcpp::as<int>(rparam["K"]);
  opt.maxnmm = Rcpp::as<int>(rparam["maxnmm"]);
  opt.maxseqlen = Rcpp::as<int>(rparam["maxseqlen"]);
  opt.maxnumseq = Rcpp::as<int>(rparam["maxnumseq"]);
  opt.useTgkm = Rcpp::as<int>(rparam["useTgkm"]);
  opt.alg = Rcpp::as<int>(rparam["alg"]);
  opt.addRC = Rcpp::as<bool>(rparam["addRC"]);
  opt.usePseudocnt = Rcpp::as<bool>(rparam["usePseudocnt"]);
  opt.OutputBinary = Rcpp::as<bool>(rparam["OutputBinary"]);
  opt.posfile = Rcpp::as<std::string>(rparam["posfile"]);
  opt.negfile = Rcpp::as<std::string>(rparam["negfile"]);
  opt.outfile = Rcpp::as<std::string>(rparam["outfile"]);
  opt.wildcardLambda = Rcpp::as<double>(rparam["wildcardLambda"]);
  opt.wildcardMismatchM = Rcpp::as<int>(rparam["wildcardMismatchM"]);
  std::string alphabetFN = Rcpp::as<std::string>(rparam["alphabetFN"]);
  if (alphabetFN != "NULL") opt.alphabetFN = alphabetFN; // the R default is the string "NULL"
  opt.maxnThread = Rcpp::as<int>(rparam["nmaxThreads"]);
  opt.mergeByName = Rcpp::as<bool>(rparam["mergeByName"]);

  if (gkmKernelRun(opt) != 0) Rcpp::stop("gkmsvm_kernel failed (see the messages above)");
}
