\name{gkmsvm_delta}
\alias{gkmsvm_delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating deltaSVM scores 
}
\description{
Given support vectors SVs and corresponding coefficients alphas and a pair of file test sequence files (one for reference allele, and one for alternate allele), calculates the deltaSVM scores for the sequences.}
\usage{
gkmsvm_delta(seqfile_allele1, seqfile_allele2, svmfnprfx, outfile, L = 10, 
K = 6, maxnmm = 3, maxseqlen = 10000, maxnumseq = 1e+06, useTgkm = 1, alg = 2, 
addRC = TRUE, usePseudocnt = FALSE, batchSize = 1e+05, wildcardLambda = 1, 
wildcardMismatchM = 2, alphabetFN = "NULL", svseqfile = NA,alphafile = NA, 
outfile_allele1 = NA, outfile_allele2 = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqfile_allele1}{
fasta file containing the test sequences (reference allele)}
  \item{seqfile_allele2}{
fasta file containing the test sequences (alternate allele). The sequences in this file should be in the exact same order as in seqfile_allele1.}
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
\item{outfile_allele1}{
output filename for gkmSVM scores for the reference sequences (optional)
}
  \item{outfile_allele2}{
output filename for gkmSVM scores for the alternate sequences (optional)
}
}
\details{
  predicting the effect of variants using gkmSVM model:
  gkmsvm_delta can be used to predict the effect of sequence variants. The sequences corresponding to reference allele and alternate alleles are given in two separate files. gkmSVM is used to score each set of sequences, and the difference in the gkmSVM score for the reference and alternate allele is reported. 
  Note that the same set of parameters (L, K, maxnmm) used in the gkmsvm_kernel should be specified for
  optimal scoring.
  
  gkmsvm_kernel(seqfile_allele1, seqfile_allele2, svmfnprfx, outfn);            #scores test sequences }
\value{
  deltaSVM scores
}
\references{
Ghandi M, Mohammad-Noori M, Ghareghani N, Lee D, Garraway LA, and Beer MA. gkmSVM: an R package for gapped-kmer SVM, Bioinformatics 2016.

Ghandi M, Lee D, Mohammad-Noori M, Beer MA. 2014. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. PLoS Comput Biol 10: e1003711. 

Lee D, Gorkin DU, Baker M, Strober BJ, Asoni AL, McCallion AS, and Beer MA. A method to predict the impact of regulatory variants from DNA sequence. Nature Genetics 2015. 
}
\author{
Mahmoud Ghandi
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

  #Input file names:  
  posfn= 'test_positives.fa'   #positive set (FASTA format)
  negfn= 'test_negatives.fa'   #negative set (FASTA format)
  testfn_ref= 'test_testsetRef.fa'    #test set (reference allele) (FASTA format)
  testfn_alt= 'test_testsetAlt.fa'    #test set (alternate allele) (FASTA format)
  
  #Output file names:  
  kernelfn= 'test_kernel.txt' #kernel matrix
  svmfnprfx= 'test_svmtrain'  #SVM files 
  outfn =   'output.txt'      #output delta svm scores for sequences in the test set       

#  gkmsvm_kernel(posfn, negfn, kernelfn);                #computes kernel 
#  gkmsvm_train(kernelfn,posfn, negfn, svmfnprfx);       #trains SVM
#  gkmsvm_delta(testfn_ref,testfn_alt, svmfnprfx, outfn);            #scores test sequences 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mutation impact}
\keyword{variant impact}% __ONLY ONE__ keyword per line
