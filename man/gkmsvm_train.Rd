\name{gkmsvm_train}
\alias{gkmsvm_train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Training the SVM model}
\description{Using the kernel matrix created by 'gkmsvm_kernel', this function trains the SVM classifier. Here we rely on the 'kernlab' package, and merely provide a wrapper function.}
\usage{gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx,  Type="C-svc", C=1, shrinking=FALSE, ...)}
\arguments{
  \item{kernelfn}{kernel matrix file name}
  \item{posfn}{positive sequences file name}
  \item{negfn}{negative sequences file name}
  \item{svmfnprfx}{output SVM model file name prefix}
  \item{Type}{optional: SVM type (default='C-svc'), see 'kernlab' documentation for more details.}
  \item{C}{optional: SVM parameter C (default=1), see 'kernlab' documentation for more details.}
  \item{shrinking}{optional: shrinking parameter for kernlab (default=FALSE), see 'kernlab' documentation for more details.}
  \item{...}{optional: additional SVM parameters, see 'kernlab' documentation for more details.}
}
\details{Trains SVM classifier and generates two files: [svmfnprfx]_svalpha.out for SVM alphas and the other for the corresponding SV sequences ([svmfnprfx]_svseq.fa) }
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
