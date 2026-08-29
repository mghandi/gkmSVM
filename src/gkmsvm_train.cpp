#include <Rcpp.h>
#include "gkmOptions.h"
using namespace Rcpp;

// C-SVC training on a precomputed kernel (LIBSVM backend); see mainSVMtrain.cpp.
// [[Rcpp::export(name = ".gkm_train_cpp")]]
void gkm_train_cpp(SEXP params){
  Rcpp::List rparam(params);
  OptsGkmTrain opt;
  opt.kernelfile = Rcpp::as<std::string>(rparam["kernelfn"]);
  opt.posfile = Rcpp::as<std::string>(rparam["posfn"]);
  opt.negfile = Rcpp::as<std::string>(rparam["negfn"]);
  opt.outprefix = Rcpp::as<std::string>(rparam["svmfnprfx"]);
  opt.C = Rcpp::as<double>(rparam["C"]);
  opt.posWeight = Rcpp::as<double>(rparam["posWeight"]);
  opt.eps = Rcpp::as<double>(rparam["eps"]);
  opt.shrinking = Rcpp::as<bool>(rparam["shrinking"]);
  opt.nfold = Rcpp::as<int>(rparam["nfold"]);
  opt.quiet = Rcpp::as<bool>(rparam["quiet"]);
  if (rparam.containsElementNamed("alphabetFN")) { std::string a = Rcpp::as<std::string>(rparam["alphabetFN"]); if (a != "NULL") opt.alphabetFN = a; }
  if (gkmTrainRun(opt) != 0) Rcpp::stop("gkmsvm_train failed (see the messages above)");
}
