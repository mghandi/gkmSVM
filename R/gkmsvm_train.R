
gkmsvm_train = function (kernelfn, posfn, negfn, svmfnprfx,  Type="C-svc", C=1, shrinking=FALSE, ...){
#TODO: add bootstrapping and cv capabilities -- also autyomatic choise of C  . check if kernlab does that 
  

  #  library(seqinr)
  #  library(kernlab)
  #  library(utils)
  if (requireNamespace("seqinr", quietly = TRUE)&
      requireNamespace("utils", quietly = TRUE)&
      requireNamespace("kernlab", quietly = TRUE)){
        
      
    
  #  negfn= '/Users/mghandi/gkmsvm/test/testneg9.fa'
  #  posfn= '/Users/mghandi/gkmsvm/test/testpos9.fa'
  #  kernelfn= '/Users/mghandi/gkmsvm/test/test9kernel.txt'
    
    pos = seqinr::read.fasta(posfn)
    neg = seqinr::read.fasta(negfn)
    idx = .gkm_read_index(kernelfn)
    if (!is.null(idx)) {
      # rows are identified by position (Phase 3): the sidecar is authoritative
      npos = sum(idx$label > 0); nneg = sum(idx$label < 0); nseq = nrow(idx)
      seqnames = idx$name
    } else {
      # kernel written by an older version: it merged same-named records within a file
      npos = length(unique(names(pos))); nneg = length(unique(names(neg))); nseq = npos+nneg
      seqnames = c(unique(names(pos)), unique(names(neg)))
    }
    
    mat <- data.matrix( utils::read.table(file=kernelfn, fill=TRUE, col.names=paste("V", 1:nseq)))
    mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
    rownames(mat)=colnames(mat)
    K <- kernlab::as.kernelMatrix(mat)
    y = c(rep(1, npos), rep(0, nneg)); names(y)=rownames(mat)
  
  #  svp <- ksvm(K, y, type="C-svc", C=1)
    svp <- kernlab::ksvm(K, y, type=Type, C=C, shrinking=shrinking, ...)
    
    if(svp@nSV>0){
      alpha = unlist(svp@alpha )
      ii = unlist(svp@SVindex)
      jj = which(ii>npos); 
      alpha[jj]= -alpha[jj];
      
      # The legacy svalpha/svseq pair identifies support vectors by name, so the names written
      # must be unique: fall back to the positional ids (pos1.., neg1..) when they are not.
      svnames = seqnames
      if (anyDuplicated(seqnames) && !is.null(idx)) {
        message("Note: sequence names are not unique; the model files use the row ids from ", kernelfn, ".index instead")
        svnames = idx$id
      }
      utils::write.table(cbind(svnames[ii], sprintf("%11.6e",alpha)),
                  file = paste(svmfnprfx, 'svalpha.out', sep='_'),
                  col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
      
      seqs = c(pos, neg)
      if (length(seqs) == nseq) {
        svseqs = seqs[ii]                                  # one record per row
      } else {
        svseqs = seqs[match(seqnames[ii], names(seqs))]    # merged rows: first record with that name
      }
      seqinr::write.fasta(svseqs, svnames[ii],  file.out= paste(svmfnprfx, 'svseq.fa', sep='_'))
    }
  }
}
