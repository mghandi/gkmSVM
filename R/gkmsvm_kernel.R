
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
                           nmaxThreads=1000,
                           mergeByName=FALSE,
                           format=c("text", "binary"),
                           tileRows=0,
                           tileMemoryMB=1024,
                           alphabets=NULL,
                           passDesign=0,
                           accMode=0,
                           accMemoryMB=4096){
  format <- match.arg(format)
  alphabetFN <- .gkm_alphabet_spec(alphabets, alphabetFN)
  
  params = list(L=L, 
                K=K, 
                maxnmm=maxnmm, 
                maxseqlen=maxseqlen,
                maxnumseq=maxnumseq, 
                useTgkm=useTgkm,
                alg=alg, 
                addRC=addRC, 
                usePseudocnt=usePseudocnt, 
                OutputBinary=(format == "binary"), 
                posfile=normalizePath(posfile, mustWork = TRUE), 
                negfile=normalizePath(negfile, mustWork = TRUE),
                outfile=normalizePath(outfile, mustWork = FALSE),
                wildcardLambda=wildcardLambda, 
                wildcardMismatchM=wildcardMismatchM,
                alphabetFN=alphabetFN, 
                nmaxThreads=nmaxThreads,
                mergeByName=mergeByName,
                tileRows=as.integer(tileRows),
                tileMemoryMB=as.integer(tileMemoryMB),
                passDesign=as.integer(passDesign),
                accMode=as.integer(accMode),
                accMemoryMB=as.integer(accMemoryMB)
                ); 
 # print(params)
  
  
  ## Identity is positional since Phase 3 (one FASTA record = one kernel row, written to
  ## <outfile>.index); duplicated names are allowed and only reported for information.
  dups <- unique(c(intersect(.gkm_fasta_names(posfile), .gkm_fasta_names(negfile)),
                   .gkm_fasta_names(posfile)[duplicated(.gkm_fasta_names(posfile))],
                   .gkm_fasta_names(negfile)[duplicated(.gkm_fasta_names(negfile))]))
  if (length(dups) > 0) {
    message(sprintf("Note: %d sequence name(s) are not unique (e.g. '%s'); rows are identified by position (see %s.index)%s",
                    length(dups), dups[1], outfile, if (mergeByName) "; records with the same name within a file are merged" else ""))
  }

 invisible(.gkm_kernel_cpp(params))
}

# First whitespace-delimited token of every FASTA header (what the C++ reader uses as the name).
.gkm_fasta_names <- function(fn) {
  h <- grep("^>", readLines(fn, warn = FALSE), value = TRUE)
  sub("[[:space:]].*$", "", substring(h, 2))
}

# Row identity for a kernel file: the <kernelfn>.index sidecar written by gkmsvm_kernel (Phase 3+),
# or NULL for kernels produced by older versions.
.gkm_read_index <- function(kernelfn) {
  fn <- paste0(kernelfn, ".index")
  if (!file.exists(fn)) return(NULL)
  utils::read.delim(fn, header = TRUE, stringsAsFactors = FALSE, colClasses = c("integer", "character", "character", "integer", "numeric", "numeric"),
                    quote = "", comment.char = "", na.strings = character(0))
}

# Phase 7: per-track alphabets. `alphabets` is a character vector with one entry per track: a keyword
# ("dna", "rna", "protein"), the path of an alphabet file, or a string of symbols (anything else, e.g.
# "01"); it is turned into the C++ spec list ("dna,=01"). Two or more entries mean the input files are
# multi-track FASTA (a header line followed by one line per track). Returns the string passed as
# `alphabetFN` to the C++ code ("NULL" = built-in DNA).
.gkm_alphabet_spec <- function(alphabets, alphabetFN = "NULL") {
  if (is.null(alphabets)) return(alphabetFN)
  if (!identical(alphabetFN, "NULL")) stop("give either 'alphabets' or 'alphabetFN', not both")
  alphabets <- as.character(alphabets)
  if (length(alphabets) < 1 || any(!nzchar(alphabets))) stop("'alphabets' must contain at least one non-empty entry")
  spec <- vapply(alphabets, function(a) {
    if (tolower(a) %in% c("dna", "rna", "protein")) return(tolower(a))
    if (startsWith(a, "=")) return(a)
    if (grepl(",", a, fixed = TRUE)) stop("an alphabet entry must not contain a comma: ", a)
    if (file.exists(a) && !dir.exists(a)) return(normalizePath(a, mustWork = TRUE))
    paste0("=", a)
  }, character(1), USE.NAMES = FALSE)
  paste(spec, collapse = ",")
}
