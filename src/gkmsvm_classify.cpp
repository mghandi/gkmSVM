#include <Rcpp.h>
#include "gkmOptions.h"
using namespace Rcpp;

// [[Rcpp::export(name = ".gkm_classify_cpp")]]
void gkm_classify_cpp(SEXP params){
  Rcpp::List rparam(params); // Get parameters in params.

  OptsSVMClassify opt;
  opt.L = Rcpp::as<int>(rparam["L"]);
  opt.K = Rcpp::as<int>(rparam["K"]);
  opt.maxnmm = Rcpp::as<int>(rparam["maxnmm"]);
  opt.maxseqlen = Rcpp::as<int>(rparam["maxseqlen"]);
  opt.maxnumseq = Rcpp::as<int>(rparam["maxnumseq"]);
  opt.batchSize = Rcpp::as<int>(rparam["batchSize"]);
  opt.useTgkm = Rcpp::as<int>(rparam["useTgkm"]);
  opt.alg = Rcpp::as<int>(rparam["alg"]);
  opt.addRC = Rcpp::as<bool>(rparam["addRC"]);
  opt.usePseudocnt = Rcpp::as<bool>(rparam["usePseudocnt"]);
  opt.seqfile = Rcpp::as<std::string>(rparam["seqfile"]);
  opt.svseqfile = Rcpp::as<std::string>(rparam["svseqfile"]);
  opt.alphafile = Rcpp::as<std::string>(rparam["alphafile"]);
  opt.outfile = Rcpp::as<std::string>(rparam["outfile"]);
  opt.wildcardLambda = Rcpp::as<double>(rparam["wildcardLambda"]);
  opt.wildcardMismatchM = Rcpp::as<int>(rparam["wildcardMismatchM"]);
  std::string alphabetFN = Rcpp::as<std::string>(rparam["alphabetFN"]);
  if (alphabetFN != "NULL") opt.alphabetFN = alphabetFN;

  if (svmClassifyRun(opt) != 0) Rcpp::stop("gkmsvm_classify failed (see the messages above)");
}
