\name{gkmSVM-package}
\alias{gkmSVM-package}
\alias{gkmSVM}
\docType{package}
\title{
Gapped-Kmer Support Vector Machine
}
\description{
Imports the 'gkmSVM' v2.0 functionalities into R <http://www.beerlab.org/gkmsvm/> . It also uses the 'kernlab' library (separate R package by different authors) for various SVM algorithms.   
}
\details{

  The gkm-SVM provides implementation of a new SVM kernel method using gapped
  k-mers as features for DNA or Protein sequences.

  There are three main functions in the gkmSVM package: 

  gkmsvm_kernel: computes the kernel matrix
  
  gkmsvm_train:  computes the SVM coefficients
  
  gkmsvm_classify: scores new sequences using the SVM model  


Tutorial

========

  We introduce the users to the basic workflow of our gkmSVM step-by-step.
  Please refer to help messages for more detailed information of each function.

  1) making a kernel matrix
  
  First of all, we should calculate a full kernel matrix before training SVM
  classifiers. In this tutorial, we are going to use test_positives.fa
  as a positive set, and test_negatives.fa as a negative set. 

  #Input file names:  
  
  posfn= 'test_positives.fa'   #positive set (FASTA format)
  
  negfn= 'test_negatives.fa'   #negative set (FASTA format)
  
  testfn= 'test_testset.fa'    #test set (FASTA format)
  
  Alternatively if the negative set is not available, and positive set is provided as a bed file, genNullSeqs function could be used to generate the negative set and positive set sequences. 
  
  #Output file names:  
  
  kernelfn= 'test_kernel.txt' #kernel matrix
  
  svmfnprfx= 'test_svmtrain'  #SVM files 
  
  outfn =   'output.txt'      #output scores for sequences in the test set       

  gkmsvm_kernel(posfn, negfn, kernelfn);                #computes kernel 

  2) training SVM
  
  We can now train a SVM classifier using the kernel matrix generated above. For that we use gkmsvm_train function It takes four arguments; kernel file, positive sequences file, negative sequences file, and prefix of output file names for the svm model.
  
  gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx);        #trains SVM

  It will generate two files, test_svmtrain_svalpha.out and
  test_svmtrain_svseq.fa, which will then be used for classification/scoring
  of test sequences as described below.

  3) classification using SVM
  
  gkmsvm_classify can be used to score any set of sequences. Here, we will 
  score the test sequences which are given in test_testset.fa. Note that the same
  set of parameters used in the gkmsvm_kernel should always be specified for
  optimal classification (here we used default parameters).
  
  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences 


  In a more advanced example, we set the word length L=18, and the number of non-gapped positions K=7, and maximum number of mismatches maxnmm=4:
  
  gkmsvm_kernel(posfn, negfn, kernelfn, L=18, K=7, maxnmm=4);     #computes kernel 
  
  gkmsvm_train(kernelfn,posfn, negfn, svmfnprfx);                 #trains SVM
  
  gkmsvm_classify(testfn, svmfnprfx, outfn, L=18, K=7, maxnmm=4); #scores test sequences 

  In another example, we run a 5-fold cross validation to plot the ROC curves: 
  
  gkmsvm_kernel(posfn, negfn, kernelfn);     #computes kernel 

  cvres = gkmsvm_trainCV(kernelfn, posfn, negfn, svmfnprfx, outputPDFfn='ROC.pdf', outputCVpredfn='cvpred.out');       #trains SVM, plots ROC and PRC curves, and outputs model predictions.

}
\author{
Mahmoud Ghandi

Maintainer: Mike Beer <mbeer@jhu.edu>
}
\references{
Ghandi M, Lee D, Mohammad-Noori M, Beer MA. 2014. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. PLoS Comput Biol 10: e1003711. 

Ghandi M, Mohammad-Noori M, Ghareghani N, Lee D, Garraway LA, and Beer MA. 2016. gkmSVM an R package for gapped-kmer SVM, Bioinformatics 32 (14), 2205-2207.
}
\keyword{package}
\keyword{gkmSVM}
\keyword{kernel}
\keyword{SVM}
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
#  gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx);       #trains SVM
#  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences 

#  using L=18, K=7, maxnmm=4

#  gkmsvm_kernel(posfn, negfn, kernelfn, L=18, K=7, maxnmm=4);     #computes kernel 
#  gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx);                 #trains SVM
#  gkmsvm_classify(testfn, svmfnprfx, outfn, L=18, K=7, maxnmm=4); #scores test sequences 

}
