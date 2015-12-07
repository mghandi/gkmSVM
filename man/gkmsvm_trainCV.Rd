\name{gkmsvm_trainCV}
\alias{gkmsvm_trainCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Training the SVM model, using repeated CV to tune parameter C and plot ROC curves}
\description{Using the kernel matrix created by 'gkmsvm_kernel', this function trains the SVM classifier. 
It uses repeated CV to find optimum SVM parameter C. Also generates ROC and PRC curves.}

\usage{gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx=NA, 
  nCV=5, nrepeat=1, cv=NA, Type="C-svc", C=1, shrinking=FALSE, 
  showPlots=TRUE, outputPDFfn=NA,  outputCVpredfn=NA, outputROCfn=NA, ...)}
  
\arguments{
  \item{kernelfn}{kernel matrix file name}
  \item{posfn}{positive sequences file name}
  \item{negfn}{negative sequences file name}
  \item{svmfnprfx}{(optional) output SVM model file name prefix }
  \item{nCV}{(optional) number of CV folds}
  \item{nrepeat}{(optional) number of repeated CVs}
  \item{cv}{(optional) CV group label. An array of length (npos+nneg), containing CV group number (between 1 an nCV) for each sequence}
  \item{Type}{(optional) SVM type (default='C-svc'), see 'kernlab' documentation for more details.}
  \item{C}{(optional)a vector of all values of C (SVM parameter) to be tested. (default=1), see 'kernlab' documentation for more details.}
  \item{shrinking}{optional: shrinking parameter for kernlab (default=FALSE), see 'kernlab' documentation for more details.}
  \item{showPlots}{generate plots (default==TRUE)}
  \item{outputPDFfn}{filename for output PDF, default=NA (no PDF output)}
  \item{outputCVpredfn}{filename for output cvpred (predicted CV values), default=NA (no output)}
  \item{outputROCfn}{filename for output auROC (Area Under an ROC Curve) and auPRC (Area Under the Precision Recall Curve) values, default=NA (no output)}
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
#  cvres = gkmsvm_trainCV(kernelfn,posfn, negfn, svmfnprfx, 
#      outputPDFfn='ROC.pdf', outputCVpredfn='cvpred.out');    
#      #trains SVM, plots ROC and PRC curves, and outputs model predictions.
#  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{gkmsvm_train}
