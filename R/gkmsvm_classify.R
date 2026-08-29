
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
                           alphafile=NA,
                           legacyNorm=FALSE){

                             if(is.na(svseqfile)){
                               # a single-file model (Phase 4) is preferred over the legacy pair; it
                               # carries the bias rho, which is applied to the scores
                               model = paste0(svmfnprfx, '.gkmmodel')
                               if (file.exists(model)) {
                                 svseqfile = model; alphafile = model
                               } else {
                                 svseqfile= paste(svmfnprfx, 'svseq.fa', sep='_')
                                 alphafile= paste(svmfnprfx, 'svalpha.out', sep='_')
                               }
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
                                           alphabetFN=alphabetFN,
                                           legacyNorm=isTRUE(legacyNorm)
                             ); 
                             # print(params)
                             
                             invisible(.gkm_classify_cpp(params))
                           }





