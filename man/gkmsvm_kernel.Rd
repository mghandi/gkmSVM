\name{gkmsvm_kernel}
\alias{gkmsvm_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computing the kernel matrix}
\description{Generates a lower triangle of kernel matrix (i.e. pairwise similarities)
  between the sequences.}
\usage{gkmsvm_kernel(posfile, negfile, outfile, L=10, K=6, maxnmm=3, maxseqlen=10000, 
maxnumseq=1000000, useTgkm=1, alg=2, addRC=TRUE, usePseudocnt=FALSE, wildcardLambda=1.0,
wildcardMismatchM=2, alphabetFN="NULL", nmaxThreads=1000, mergeByName=FALSE,
format=c("text", "binary"), tileRows=0, tileMemoryMB=1024, alphabets=NULL, passDesign=0,
accMode=0, accMemoryMB=4096)}

\arguments{
\item{posfile}{positive sequences file name (FASTA format)}
\item{negfile}{negative sequences file name (FASTA format)}
\item{outfile}{output file name}
\item{L}{word length, default=10}
\item{K}{number of informative columns, default=6}
\item{maxnmm}{maximum number of mismatches to consider, default=3}
\item{maxseqlen}{maximum sequence length in the sequence files, default=10000}
\item{maxnumseq}{maximum number of sequences in the sequence files, default=1000000}
\item{useTgkm}{filter type: 0(use full filter), 1(use truncated filter: this gaurantees non-negative counts for all L-mers), 2(use h[m], gkm count vector), 3(wildcard), 4(mismatch), default=1}
\item{alg}{algorithm type: 0(auto), 1(XOR Hashtable), 2(tree), default=2}
\item{addRC}{adds reverse complement sequences, default=TRUE}
\item{usePseudocnt}{adds a constant to count estimates, default=FALSE}
\item{wildcardLambda}{lambda for wildcard kernel, defaul=0.9}
\item{wildcardMismatchM}{max mismatch for Mismatch kernel or wildcard kernel, default=2}
\item{alphabetFN}{alphabet: one of the keywords "dna" (the default), "rna", "protein", or the name of a file listing one symbol per line. A second symbol on a line declares the complement of the first (e.g. "A U"); reverse complements (\code{addRC}) are only added when every symbol has a complement (built in for dna/rna and for any 4-letter file, which is paired 1st-4th, 2nd-3rd). Alphabets of up to 32 symbols are supported; DNA keeps a dedicated fast path.}
\item{nmaxThreads}{maximum number of threads, defaul=2*L}
\item{format}{"text" (default) writes the historical tab-separated lower-triangle text file; "binary" writes the compact \code{.gkmk} format (float32, about 3x smaller and much faster to load; see \code{\link{read_gkm_kernel}}). \code{gkmsvm_train} and \code{gkmsvm_trainCV} read both.}
\item{tileRows}{number of kernel rows computed at once. The mismatch profile of a tile needs about \code{tileRows * (maxnmm+1) * 4 * n} bytes; the default 0 chooses the largest tile that fits \code{tileMemoryMB}. Results do not depend on the tiling (integer counts); more tiles cost more time.}
\item{tileMemoryMB}{memory budget (MB) for the mismatch profile when \code{tileRows = 0}, default 1024}
\item{passDesign}{how the mismatch patterns are split into passes (advanced; the result is the same): 0 = automatic (default), 1 = the greedy iDL design of the 2016 paper (needs a table of every pattern with at most \code{maxnmm} mismatches), 2 = prefix-split passes, which need no table and are used automatically for long multi-track words and, with 4 or more threads, for capped \\code{maxnmm} on words of 16 or more positions; 4 = the greedy design with its 2016 cost model (1 now measures per-position match probabilities and per-order trie widths on the actual l-mers).}
\item{accMode}{how the kernel is accumulated (advanced; the result is the same): 0 = automatic (default), 1 = the mismatch-class profile counters (tiled to \code{tileMemoryMB}), 2 = kernel mode: the coefficient of each mismatch class, scaled to an integer, is added directly into one 64-bit integer lower triangle per thread, so no per-class counters, no tiles and no repeated traversal. Kernel mode applies to the full filter (\code{useTgkm = 0}) and the gapped k-mer counts (\code{useTgkm = 2}) without pseudocounts, and is chosen automatically there; the truncated filter and the wildcard/mismatch kernels keep the profile.}
\item{accMemoryMB}{memory budget (MB) for the kernel-mode accumulators, default 4096: one private triangle per thread while they fit (plain adds), fewer shared triangles with atomic adds otherwise. Results do not depend on it.}
\item{alphabets}{one alphabet per \emph{track}, for multi-track input (gkm-SVM over heterogeneous alphabets): a character vector whose entries are a keyword ("dna", "rna", "protein"), the path of an alphabet file (as for \code{alphabetFN}), or a string of symbols (e.g. "01" for a binary annotation, "NUM" for a 3-state one; a string that is also the name of an existing file is read as a file). With two or more entries the input files are multi-track FASTA: each record is a header line followed by one line per track, all of the same length, e.g. a DNA line and a 0/1 methylation line. A window of length \code{L} over the T tracks is a word of T*L positions whose mismatches are counted per alphabet size ("block"); \code{maxnmm} bounds the total number of mismatches (-1 = all). A symbol outside its track's alphabet (e.g. "N") skips every window that covers it. Only \code{useTgkm} 0, 1, 2 are available for multi-track input; reverse complements complement track 1 and reverse the others. Two tracks are supported in this version. Cannot be combined with \code{alphabetFN}.}
\item{mergeByName}{if TRUE, records that share a name within one FASTA file are merged into a single kernel row (the behaviour of versions before 0.90). Default FALSE: one record = one row. Row identity (id, name, label, length) is written to \code{<outfile>.index}.}

}
\details{  It calculates the full kernel matrix that can be then used to train an SVM classifier. 
gkmsvm_kernel(posfn, negfn, kernelfn);

}
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
