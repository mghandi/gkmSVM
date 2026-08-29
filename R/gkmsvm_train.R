
gkmsvm_train = function (kernelfn, posfn, negfn, svmfnprfx,  Type="C-svc", C=1, shrinking=FALSE,
                         backend=c("kernlab", "libsvm"), posWeight=1, eps=0.001, nfold=0, quiet=TRUE,
                         alphabets=NULL, ...){
  backend <- match.arg(backend)
  # Phase 7: multi-track inputs (two or more alphabets): records are header + one line per track
  spec <- .gkm_alphabet_spec(alphabets)
  multitrack <- !is.null(alphabets) && length(alphabets) > 1
  if (backend == "libsvm") {
    # Phase 4b: the C++ trainer (vendored LIBSVM, precomputed kernel). Same solver family as
    # kernlab's C-svc; writes the legacy pair and <prefix>.gkmmodel (with rho) itself.
    if (Type != "C-svc") stop("backend='libsvm' supports Type='C-svc' only")
    params <- list(kernelfn = normalizePath(kernelfn, mustWork = TRUE), posfn = normalizePath(posfn, mustWork = TRUE),
                   negfn = normalizePath(negfn, mustWork = TRUE), svmfnprfx = svmfnprfx,
                   C = as.numeric(C), posWeight = as.numeric(posWeight), eps = as.numeric(eps),
                   shrinking = isTRUE(shrinking), nfold = as.integer(nfold), quiet = isTRUE(quiet), alphabetFN = spec)
    return(invisible(.gkm_train_cpp(params)))
  }

  #  library(seqinr)
  #  library(kernlab)
  #  library(utils)
  if (requireNamespace("seqinr", quietly = TRUE)&
      requireNamespace("utils", quietly = TRUE)&
      requireNamespace("kernlab", quietly = TRUE)){
        
      
    
  #  negfn= '/Users/mghandi/gkmsvm/test/testneg9.fa'
  #  posfn= '/Users/mghandi/gkmsvm/test/testpos9.fa'
  #  kernelfn= '/Users/mghandi/gkmsvm/test/test9kernel.txt'
    
    if (multitrack) {
      pos = read_mfa(posfn, length(alphabets)); neg = read_mfa(negfn, length(alphabets))
    } else {
      pos = seqinr::read.fasta(posfn)
      neg = seqinr::read.fasta(negfn)
    }
    mat = read_gkm_kernel(kernelfn)   # text or binary (.gkmk), auto-detected
    idx = attr(mat, "index")
    if (!is.null(idx)) {
      # rows are identified by position (Phase 3): the sidecar is authoritative
      npos = sum(idx$label > 0); nneg = sum(idx$label < 0); nseq = nrow(idx)
      seqnames = idx$name
    } else {
      # kernel written by an older version: it merged same-named records within a file
      npos = length(unique(names(pos))); nneg = length(unique(names(neg))); nseq = npos+nneg
      seqnames = c(unique(names(pos)), unique(names(neg)))
    }
    
    if (nrow(mat) != nseq) stop(sprintf("kernel %s has %d rows but the sequence files have %d sequences", kernelfn, nrow(mat), nseq))
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
      if (multitrack) {
        write_mfa(svseqs, paste(svmfnprfx, 'svseq.fa', sep='_'), names = svnames[ii])
      } else {
        seqinr::write.fasta(svseqs, svnames[ii],  file.out= paste(svmfnprfx, 'svseq.fa', sep='_'))
      }
      
      # single-file model (Phase 4): FASTA with ">id<TAB>alpha" headers, "#" header lines carry the
      # bias rho = svp@b (kernlab's b), which gkmsvm_classify subtracts from the scores
      con = file(paste0(svmfnprfx, '.gkmmodel'), "w")
      writeLines(c("#gkmmodel 1", sprintf("#rho %.10e", svp@b), sprintf("#nsv %d", length(ii)), sprintf("#npos %d", npos), sprintf("#nneg %d", nneg)), con)
      if (multitrack) writeLines(sprintf("#alphabets %s", spec), con)
      for (s in seq_along(ii)) {
        if (multitrack) writeLines(c(sprintf(">%s\t%11.6e", svnames[ii[s]], alpha[s]), svseqs[[s]]), con)
        else writeLines(c(sprintf(">%s\t%11.6e", svnames[ii[s]], alpha[s]), toupper(paste(svseqs[[s]], collapse = ""))), con)
      }
      close(con)
    }
  }
}
