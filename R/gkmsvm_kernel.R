
gkmsvm_kernel <- function( posfile, 
                           negfile,
                           outfile,
                           L=10, 
                           K=6, 
                           maxnmm=3, 
                           maxseqlen=10000,
                           maxnumseq=1000000, 
                           useTgkm=1,
                           alg=2, 
                           addRC=TRUE, 
                           usePseudocnt=FALSE, 
                           wildcardLambda=1.0, 
                           wildcardMismatchM=2,
                           alphabetFN="NULL",
                           nmaxThreads=1000){
  
  params = list(L=L, 
                K=K, 
                maxnmm=maxnmm, 
                maxseqlen=maxseqlen,
                maxnumseq=maxnumseq, 
                useTgkm=useTgkm,
                alg=alg, 
                addRC=addRC, 
                usePseudocnt=usePseudocnt, 
                OutputBinary=FALSE, 
                posfile=normalizePath(posfile, mustWork = TRUE), 
                negfile=normalizePath(negfile, mustWork = TRUE),
                outfile=normalizePath(outfile, mustWork = FALSE),
                wildcardLambda=wildcardLambda, 
                wildcardMismatchM=wildcardMismatchM,
                alphabetFN=alphabetFN, 
                nmaxThreads=nmaxThreads
                ); 
 # print(params)
  
  
  ## test for duplicate id 
  if (requireNamespace("seqinr", quietly = TRUE)){
    
    posfn = posfile; 
    negfn = negfile;
    
    pos = seqinr::read.fasta(posfn)
    neg = seqinr::read.fasta(negfn)
    
    if(length(which(duplicated(names(pos))))>0){
      print(paste("Error: duplicated sequence ID in", posfn))
      print(names(pos)[which(duplicated(names(pos)))])
      stop("Error: duplicated sequence ID");
    }
    if(length(which(duplicated(names(neg))))>0){
      print(paste("Error: duplicated sequence ID in", negfn))
      print(names(neg)[which(duplicated(names(neg)))])
      stop("Error: duplicated sequence ID");
    }
    if(length(which(duplicated(c(names(pos),names(neg)))))>0){
      print(paste("Error: Same sequence ID found in positive and negative sets:", posfn, negfn))
      print(c(names(pos),names(neg))[which(duplicated(c(names(pos),names(neg))))])
      stop("Error: duplicated sequence ID");
    }
  }
  
  
 invisible(.Call('gkmSVM_gkmsvm_kernel', PACKAGE = 'gkmSVM', params))
}
