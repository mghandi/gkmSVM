
gkmsvm_classify <- function( seqfile, 
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
                           alphafile=NA){

                             if(is.na(svseqfile)){
                               svseqfile= paste(svmfnprfx, 'svseq.fa', sep='_')
                               alphafile= paste(svmfnprfx, 'svalpha.out', sep='_')
                             }    
                             
                             params = list(seqfile=normalizePath(seqfile, mustWork =  TRUE), 
                                           svseqfile=normalizePath(svseqfile, mustWork = TRUE),
                                           alphafile=normalizePath(alphafile, mustWork = TRUE),
                                           outfile=normalizePath(outfile, mustWork = FALSE),
                                           L=L, 
                                           K=K, 
                                           maxnmm=maxnmm, 
                                           maxseqlen=maxseqlen,
                                           maxnumseq=maxnumseq, 
                                           useTgkm=useTgkm,
                                           alg=alg, 
                                           addRC=addRC, 
                                           usePseudocnt=usePseudocnt, 
                                           batchSize=batchSize,
                                           wildcardLambda=wildcardLambda, 
                                           wildcardMismatchM=wildcardMismatchM,
                                           alphabetFN=alphabetFN
                             ); 
                             # print(params)
                             
                             invisible(.Call('gkmSVM_gkmsvm_classify', PACKAGE = 'gkmSVM', params))
                           }





