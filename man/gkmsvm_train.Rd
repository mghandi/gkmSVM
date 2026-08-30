\name{gkmsvm_train}
\alias{gkmsvm_train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Training the SVM model}
\description{Using the kernel matrix created by 'gkmsvm_kernel', this function trains the SVM classifier, either through the 'kernlab' package (default) or through the LIBSVM solver embedded in the package (\code{backend="libsvm"}).}
\usage{gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx,  Type="C-svc", C=1, shrinking=FALSE,
             backend=c("kernlab", "libsvm"), posWeight=1, eps=0.001, nfold=0,
             quiet=TRUE, alphabets=NULL, ...)}
\arguments{
  \item{kernelfn}{kernel matrix file name}
  \item{posfn}{positive sequences file name}
  \item{negfn}{negative sequences file name}
  \item{svmfnprfx}{output SVM model file name prefix}
  \item{Type}{optional: SVM type (default='C-svc'), see 'kernlab' documentation for more details.}
  \item{C}{optional: SVM parameter C (default=1), see 'kernlab' documentation for more details.}
  \item{shrinking}{optional: shrinking parameter for kernlab (default=FALSE), see 'kernlab' documentation for more details.}
  \item{backend}{"kernlab" (default) or "libsvm": the C-SVC solver of LIBSVM 3.37 embedded in the package, the same solver family that kernlab uses; it reads text and binary kernels directly and is what the standalone \code{gkmsvm_train} binary runs.}
  \item{posWeight}{libsvm backend only: weight of the positive class (the negative class has weight 1), default 1}
  \item{eps}{libsvm backend only: stopping tolerance, default 0.001}
  \item{nfold}{libsvm backend only: if > 1, the nfold cross-validation accuracy is printed before training}
  \item{quiet}{libsvm backend only: suppress the solver's progress output (default TRUE)}
  \item{alphabets}{one alphabet per track when the kernel was computed from multi-track input (see \code{\link{gkmsvm_kernel}}): the sequence files are then multi-track FASTA, the support vectors are written as multi-track records and the \code{.gkmmodel} file records the alphabets for \code{\link{gkmsvm_classify}}.}
  \item{...}{optional: additional SVM parameters, see 'kernlab' documentation for more details.}
}
\details{Trains SVM classifier and generates three files: [svmfnprfx]_svalpha.out for SVM alphas, [svmfnprfx]_svseq.fa for the corresponding SV sequences, and the single-file model [svmfnprfx].gkmmodel (support vectors with their alphas and the bias rho), which \code{gkmsvm_classify} prefers when present. The kernel file may be in the text or the binary (\code{.gkmk}) format.}
\author{Mahmoud Ghandi}
\examples{
  #Input file names:  
  posfn= 'test_positives.fa'   #positive set (FASTA format)
  negfn= 'test_negatives.fa'   #negative set (FASTA format)
  testfn= 'test_testset.fa'    #test set (FASTA format)
  
  #Output file names:  
  kernelfn= 'test_kernel.txt' #kernel matrix
  svmfnprfx= 'test_svmtrain'  #SVM files 
  outfn =   'output.txt'      #output scores for sequences in the test set       

#  gkmsvm_kernel(posfn, negfn, kernelfn);                #computes kernel 
#  gkmsvm_train(kernelfn,posfn, negfn, svmfnprfx);       #trains SVM
#  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{gkmsvm_train}
