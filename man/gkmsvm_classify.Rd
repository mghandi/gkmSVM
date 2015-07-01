\name{gkmsvm_classify}
\alias{gkmsvm_classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classifying(/scoring) new sequences using the gkmSVM model}
\description{Given support vectors SVs and corresponding coefficients alphas and a set of 
  sequences, calculates the SVM scores for the sequences.}
\usage{gkmsvm_classify(seqfile, svmfnprfx, outfile, L=10, K=6, maxnmm=3, 
maxseqlen=10000, maxnumseq=1000000, useTgkm=1, alg=0, addRC=TRUE, usePseudocnt=FALSE, 
batchSize=100000, wildcardLambda=1.0, wildcardMismatchM=2, alphabetFN="NULL", 
svseqfile=NA, alphafile=NA)}
\arguments{
\item{seqfile}{input sequences file name (FASTA format)}
\item{svmfnprfx}{SVM model file name prefix}
\item{outfile}{output file name}
\item{L}{word length, default=10}
\item{K}{number of informative columns, default=6}
\item{maxnmm}{maximum number of mismatches to consider, default=3}
\item{maxseqlen}{maximum sequence length in the sequence files, default=10000}
\item{maxnumseq}{maximum number of sequences in the sequence files, default=1000000}
\item{useTgkm}{filter type: 0(use full filter), 1(use truncated filter: this gaurantees non-negative counts for all L-mers), 2(use h[m], gkm count vector), 3(wildcard), 4(mismatch), default=1}
\item{alg}{algorithm type: 0(auto), 1(XOR Hashtable), 2(tree), default=0}
\item{addRC}{adds reverse complement sequences, default=TRUE}
\item{usePseudocnt}{adds a constant to count estimates, default=FALSE}
\item{batchSize}{number of sequences to compute scores for in batch, default=100000}
\item{wildcardLambda}{lambda for wildcard kernel, defaul=0.9}
\item{wildcardMismatchM}{max mismatch for Mismatch kernel or wildcard kernel, default=2}
\item{alphabetFN}{alphabets file name, if not specified, it is assumed the inputs are DNA sequences}
\item{svseqfile}{SVM support vectors sequence file name (not needed if svmfnprfx is provided)}
\item{alphafile}{SVM support vectors weights file name (not needed if svmfnprfx is provided)}

}
\details{ classification using SVM:
  gkmsvm_classify can be used to score any set of sequences. Note that the same
  set of parameters (L, K, maxnmm) used in the gkmsvm_kernel should be specified for
  optimal classification.
  
  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences  }
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
%\keyword{gkmsvm_classify}
