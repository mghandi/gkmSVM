gkmsvm_delta<- function( seqfile_allele1, 
                         seqfile_allele2,
                             svmfnprfx,
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
                             batchSize=100000, 
                             wildcardLambda=1.0, 
                             wildcardMismatchM=2,
                             alphabetFN="NULL",
                             svseqfile=NA,
                             alphafile=NA, 
                             outfile_allele1 = NA, 
                             outfile_allele2 = NA 
                         ){
  
  
  if(is.na(outfile_allele1)){
    outfile_allele1 = tempfile()
  }
  if(is.na(outfile_allele2)){
    outfile_allele2 = tempfile()
  }
  
  gkmsvm_classify(seqfile = seqfile_allele1, svmfnprfx = svmfnprfx, outfile = outfile_allele1, 
                  L=L ,K=K, maxnmm=maxnmm,maxseqlen=maxseqlen,maxnumseq=maxnumseq,useTgkm=useTgkm,alg=alg,addRC=addRC,
                  usePseudocnt=usePseudocnt,batchSize=batchSize,wildcardLambda=wildcardLambda,wildcardMismatchM=wildcardMismatchM,
                  alphabetFN=alphabetFN,svseqfile=svseqfile,alphafile=alphafile)
  gkmsvm_classify(seqfile = seqfile_allele2, svmfnprfx = svmfnprfx, outfile = outfile_allele2, 
                  L=L ,K=K, maxnmm=maxnmm,maxseqlen=maxseqlen,maxnumseq=maxnumseq,useTgkm=useTgkm,alg=alg,addRC=addRC,
                  usePseudocnt=usePseudocnt,batchSize=batchSize,wildcardLambda=wildcardLambda,wildcardMismatchM=wildcardMismatchM,
                  alphabetFN=alphabetFN,svseqfile=svseqfile,alphafile=alphafile)
  
  s1 = utils::read.delim(outfile_allele1, header=FALSE)
  s2 = utils::read.delim(outfile_allele2, header=FALSE)
  res =  as.numeric(s2[,2])-as.numeric(s1[,2]); names(res)=as.character(s1[,1]); 
  utils::write.table(cbind(as.character(s1[,1]), res), file = outfile, quote = FALSE, row.names=FALSE, col.names = FALSE) 
  res=res; 
}
