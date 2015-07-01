#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
int mainSVMclassify(int argc, char** argv); 
  
// [[Rcpp::export]]
void gkmsvm_classify(SEXP params){
  Rcpp::List rparam(params); // Get parameters in params.

  int L = Rcpp::as<int>(rparam["L"]);
  int K = Rcpp::as<int>(rparam["K"]);
  int maxnmm= Rcpp::as<int>(rparam["maxnmm"]);
  int maxseqlen= Rcpp::as<int>(rparam["maxseqlen"]);
  int maxnumseq= Rcpp::as<int>(rparam["maxnumseq"]);
  int batchSize= Rcpp::as<int>(rparam["batchSize"]);
  int useTgkm= Rcpp::as<int>(rparam["useTgkm"]);
  int alg= Rcpp::as<int>(rparam["alg"]);
  bool addRC= Rcpp::as<bool>(rparam["addRC"]);
  bool usePseudocnt= Rcpp::as<bool>(rparam["usePseudocnt"]);
  std::string seqfile = Rcpp::as<std::string>(rparam["seqfile"]);
  std::string svseqfile = Rcpp::as<std::string>(rparam["svseqfile"]);
  std::string alphafile = Rcpp::as<std::string>(rparam["alphafile"]);
  std::string outfile = Rcpp::as<std::string>(rparam["outfile"]);
  double wildcardLambda= Rcpp::as<double>(rparam["wildcardLambda"]);
  int wildcardMismatchM= Rcpp::as<int>(rparam["wildcardMismatchM"]);
  std::string alphabetFN = Rcpp::as<std::string>(rparam["alphabetFN"]);
 
   
    int argc=0; 
    char** argv = new char *[30]; 
    for(int i=0;i<30;i++){
      argv[i]=new char[5000]; 
    }
    
    sprintf(argv[argc++], "gkmsvm_classify"); 
    
    sprintf(argv[argc++], "-l"); 
    sprintf(argv[argc++], "%d", L); 
    
    sprintf(argv[argc++], "-k"); 
    sprintf(argv[argc++], "%d", K); 
    
    sprintf(argv[argc++], "-d"); 
    sprintf(argv[argc++], "%d", maxnmm); 
    
    sprintf(argv[argc++], "-m"); 
    sprintf(argv[argc++], "%d", maxseqlen); 
    
    sprintf(argv[argc++], "-n"); 
    sprintf(argv[argc++], "%d", maxnumseq); 
    
    sprintf(argv[argc++], "-t"); 
    sprintf(argv[argc++], "%d", useTgkm); 
    
    sprintf(argv[argc++], "-a"); 
    sprintf(argv[argc++], "%d", alg); 
    
    sprintf(argv[argc++], "-b"); 
    sprintf(argv[argc++], "%d", batchSize); 

    sprintf(argv[argc++], "-M"); 
    sprintf(argv[argc++], "%d", wildcardMismatchM); 
    
    sprintf(argv[argc++], "-L"); 
    sprintf(argv[argc++], "%lf", wildcardLambda); 
    
    if(strcmp(alphabetFN.c_str(), "NULL")!=0){
      sprintf(argv[argc++], "-A"); 
      sprintf(argv[argc++], "%s", alphabetFN.c_str()); 
    }

    if(addRC==false){
      sprintf(argv[argc++], "-R"); 
    }
    if(usePseudocnt==true){
      sprintf(argv[argc++], "-p"); 
    }

    sprintf(argv[argc++], "%s", seqfile.c_str()); 
    sprintf(argv[argc++], "%s", svseqfile.c_str()); 
    sprintf(argv[argc++], "%s", alphafile.c_str()); 
    sprintf(argv[argc++], "%s", outfile.c_str()); 

    mainSVMclassify( argc, argv); 
   
   for(int i=0;i<30;i++){
     delete []argv[i]; 
   }
   delete[]argv;   
   
}
  
