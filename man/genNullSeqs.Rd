
\name{genNullSeqs}
\alias{genNullSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating GC/repeat matched randomly selected genomic sequences for the negative set
}
\description{
Generates null sequences (negative set) with matching repeat and GC content as the input bed file for positive set regions.
}
\usage{
genNullSeqs(
  inputBedFN, 
  genomeVersion='hg19', 
  outputBedFN = 'negSet.bed', 
  outputPosFastaFN = 'posSet.fa',
  outputNegFastaFN = 'negSet.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputBedFN}{
  positive set regions 
}
  \item{genomeVersion}{
  genome version: 'hg19' and 'hg18' are supported. Default='hg19'. For other genomes, provide the BSgenome object using parameter 'genome'
}
  \item{outputBedFN}{
  output file name for the null sequences genomic regions. Default='negSet.bed'
}
  \item{outputPosFastaFN}{
     output file name for the positive set sequences. Default='posSet.fa'
     }
  \item{outputNegFastaFN}{
     output file name for the negative set sequences. Default='negSet.fa'
     }
  \item{xfold}{
     controls the desired number of sequences in the negative set. Default=1 (same number as in positive set)
     }
  \item{repeat_match_tol}{
     tolerance for difference in repeat ratio. Default=0.02 (repeat content difference of 0.02 or less is acceptable)      }
  \item{GC_match_tol}{
     tolerance for difference in GC content. Default=0.02
     }

  \item{length_match_tol}{
     tolerance for difference in relative sequence length. Default=0.02
     }
  \item{batchsize}{
     number of candidate random sequences tested in each trial. Default=5000
     }
  \item{nMaxTrials}{
     maximum number of trials. Default=20. 
     }
  \item{genome}{
     BSgenome object. Default=NULL. If this parameter is used, parameter genomeVersion is ignored. 
     }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Writes the null sequences to files with the provided filenames. Outputs the filename for the output negative sequences file. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Mahmoud Ghandi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# Example 1: 
# genNullSeqs('ctcfpos.bed' );

#Example 2:
# genNullSeqs('ctcfpos.bed', nMaxTrials=3, xfold=2, genomeVersion = 'hg18' );

#Example 3:
# genNullSeqs('ctcfpos.bed', xfold=2, genomeVersion = 'hg18', outputBedFN = 'ctcf_negSet.bed',
# outputPosFastaFN = 'ctcf_posSet.fa',outputNegFastaFN = 'ctcf_negSet.fa' );

#Example 4:
# Input file names:  
  
  posBedFN = 'test_positives.bed' # positive set genomic ranges (bed format)
  genomeVer = 'hg19' #genome version 
  testfn= 'test_testset.fa'    #test set (FASTA format)
  
# output file names:  
  posfn= 'test_positives.fa'   #positive set (FASTA format)
  negfn= 'test_negatives.fa'   #negative set (FASTA format)
  kernelfn= 'test_kernel.txt' #kernel matrix
  svmfnprfx= 'test_svmtrain'  #SVM files 
  outfn =   'output.txt'      #output scores for sequences in the test set       

#  genNullSeqs(posBedFN, genomeVersion = genomeVer, 
#    outputPosFastaFN = posfn, outputNegFastaFN = negfn );

#  gkmsvm_kernel(posfn, negfn, kernelfn);                #computes kernel 
#  gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx);       #trains SVM
#  gkmsvm_classify(testfn, svmfnprfx, outfn);            #scores test sequences 

#  using L=18, K=7, maxnmm=4

#  gkmsvm_kernel(posfn, negfn, kernelfn, L=18, K=7, maxnmm=4);     #computes kernel 
#  gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx);                 #trains SVM
#  gkmsvm_classify(testfn, svmfnprfx, outfn, L=18, K=7, maxnmm=4); #scores test sequences 


}