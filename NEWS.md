# gkmSVM NEWS

## Unreleased (development version, refactoring in progress)

See `dev/REFACTORING_PLAN.md` for the plan and its execution log.

### Phase 1 — latent-bug fixes

Two of the fixes change numerical results; everything else is behaviour-preserving (the golden
corpus in `tests/golden/` is byte-identical for all kernel outputs).

* **`gkmsvm_classify` (tree algorithm, `alg=2`, the default): norms are now truncated at `maxnmm`
  like the scores.** The score sums mismatch levels `m <= maxnmm` (`-d`), but for DNA the
  normalisation summed all `m <= L`; the two were inconsistent whenever `maxnmm` was smaller than
  the support of the kernel weights (the non-DNA path and the kernel matrix were already consistent).
  Scores therefore change slightly when `maxnmm` (default 3) is small. Measured on the test corpus
  (10 test sequences, a fixed 10-SV model, L=10, K=6, d=3): max |Δscore| = 0.00016 (t=1, default),
  0.00003 (t=0), 0.0008 (t=2), 0.0023 (t=4); score ratios new/old within 0.98–1.02 for t=1; no sign
  flips. `alg=1` and `maxnmm=-1` (all levels) are unchanged. The new scores are verified against an
  exact-arithmetic reference (`tests/oracle/oracle_check.py`, "classify scores").
* **`gkmsvm_classify` with a non-DNA alphabet and `useTgkm=3` (wildcard) or `4` (mismatch)** used
  the DNA alphabet size (a literal 4) when computing the weights; it now uses the actual alphabet
  size, as `gkmsvm_kernel` already did. DNA results are unchanged.

Crashes and memory errors fixed (all reproduced before the fix, all now ASAN/UBSAN-clean):

* An empty FASTA record (or, silently, every end-of-file) wrote one element before the start of a
  heap buffer; ~38 consecutive `gkmsvm_kernel()` calls in one R session then aborted R.
* FASTA names of 100 characters or more overflowed a fixed buffer (aborted R immediately).
* Sequences longer than `maxseqlen` overflowed four buffers: now an error naming the sequence.
* More sequences than `maxnumseq` overran several arrays: now an error.
* An alphabet file with more than `MAX_ALPHABET_SIZE` (4) symbols crashed with SIGSEGV: now an
  error. A missing alphabet file is an error too (it was a NULL dereference).
* The alphabet loaded with `alphabetFN` persisted into the next call in the same R session.
* `gkmsvm_classify` crashed when the number of test sequences was an exact multiple of
  `batchSize` (an empty final batch wrote into a zero-length array).
* `gkmsvm_classify(alg=1)` never closed its output file, so in a long-lived R session the scores
  stayed in the stdio buffer.
* Out-of-bounds read of the `wm[]` table on every kernel run; `h[]` partly uninitialised;
  `delete`/`delete[]` mismatches; several leaks per call.
* Errors in the C++ core are now reported to R as errors (`stop()`); previously the R functions
  returned normally with no or partial output.
* The mismatch kernel (`useTgkm=4`) with a 2-letter alphabet produced NaN weights (`0^negative`
  in the weight formula); it now produces finite ones. No other alphabet is affected.
* Two Linux-only defects: `getopt` was not re-initialised between calls from R on glibc (a later
  call could misparse its options and silently write nothing), and the reverse-complement code
  copied a sequence name onto itself (undefined; flagged by AddressSanitizer on Linux).

### Phase 2a — core extraction (behaviour-preserving)

* Invalid parameters (`K > L`, `maxnmm > L`, `useTgkm` outside 0–4, `alg` outside 0–2, `batchSize < 1`)
  are now errors (`stop()` in R, exit status 1 on the command line). They used to print the usage
  text and return success, leaving no output file.
* The R bindings pass their arguments to the C++ core directly instead of re-serialising them into a
  command line; file paths are no longer silently truncated at 5000 characters.
* Thirteen unreachable source files were moved to `src/legacy/` and are no longer compiled.
* Kernel and classify outputs are byte-identical to the previous version.

### Phase 3 — sequence identity (one intentional behaviour change)

* **One FASTA record is now one kernel row.** Previously, records that shared a name *within* one
  file were silently merged into a single row (a 60-record file with one repeated name produced a
  59-row kernel). That behaviour is still available as `gkmsvm_kernel(..., mergeByName = TRUE)` /
  `gkmsvm_kernel -N` and reproduces the old output byte-for-byte; the default is off.
* `gkmsvm_kernel` writes `<outfile>.index` next to the kernel: `row id name label length nlmers`,
  one line per row. Ids (`pos1..`, `neg1..`) are unique by construction; names may be empty,
  repeated, non-ASCII or arbitrarily long. `gkmsvm_train` and `gkmsvm_trainCV` use the sidecar
  when it exists and fall back to the old `read.fasta` counting for kernels made by older versions.
* The R functions no longer `stop()` on duplicated sequence names (a duplicate is reported with a
  `message()`), and no longer read both FASTA files a second time just to check.
* When names are not unique, `gkmsvm_train` writes the row ids instead of the names into
  `_svalpha.out` / `_svseq.fa`, so the legacy model files stay unambiguous.
* `gkmsvm_classify` matches support vectors to sequences with a hash map instead of a linear scan
  and **fails** when a name listed in `_svalpha.out` is duplicated there, missing from `_svseq.fa`,
  or present twice in `_svseq.fa`. It used to print "the sequences for only N out of M ... were
  found" and score with a partial model.

### Phase 4 — binary formats (additive)

* **Binary kernel files.** `gkmsvm_kernel(..., format = "binary")` (CLI `-b`, which was parsed
  but did nothing since 2.0) writes the `.gkmk` v1 format: magic `GKMK`, parameters, alphabet, the
  row-identity table, the lower triangle as float32 and a CRC-32. About 3x smaller than the text
  file and far faster to load. `gkmsvm_train` / `gkmsvm_trainCV` read both formats transparently
  (detected by the magic bytes). New helpers `read_gkm_kernel()`, `write_gkm_kernel()`,
  `is_gkm_binary()`. The default stays `"text"` for this release.
* **Single-file model.** `gkmsvm_train` now also writes `<prefix>.gkmmodel` — a FASTA file with
  `>id<TAB>alpha` headers and `#rho <bias>` in its header — next to the legacy
  `_svalpha.out`/`_svseq.fa` pair, which is still written. `gkmsvm_classify(seqfile, svmfnprfx)`
  prefers the `.gkmmodel` when it exists and then **applies the SVM bias** (`score − rho`, kernlab's
  `b`); the legacy pair stays bias-free as before, so rankings and `gkmsvm_delta` are unaffected.
  The C++ classifier accepts a `.gkmmodel` as both `<sv_seqfile>` and `<sv_alphafile>`.
* FASTA lines starting with `#` are now treated as comments (like `;`).
* The R/C++ bridge uses registered routines (`Rcpp::compileAttributes()` output).

### Phase 4b — SVM training in C++ (LIBSVM)

* The standalone `gkmsvm_train` binary is now a real C-SVC solver: it embeds **LIBSVM 3.37**
  (`src/libsvm/`, BSD-3-Clause, `COPYRIGHT` retained) with `kernel_type = PRECOMPUTED`. It reads
  text or binary kernels, takes `-c C`, `-w` (positive-class weight), `-e eps`, `-v nfold`
  (cross-validation accuracy), `-S` (no shrinking), `-q`; `-n niter20` is accepted and ignored.
  It writes the legacy `_svalpha.out`/`_svseq.fa` pair and `<prefix>.gkmmodel` (with `rho`). The
  previous solver (`CSVMtrain`, a fixed-iteration heuristic with no `C` and no stopping criterion)
  is deleted.
* `gkmsvm_train(..., backend = "libsvm")` calls the same code from R; the default `backend =
  "kernlab"` is unchanged for this release.

### Phase 5 — alphabets other than DNA

* Alphabets of **up to 32 symbols** work (protein, DNA+N, …). Previously anything larger than 4
  symbols crashed with SIGSEGV. DNA keeps a dedicated compiled fast path (identical code and node
  layout to before), so DNA speed and results are unchanged; larger alphabets use a wider trie.
* `alphabetFN` accepts the keywords `"dna"`, `"rna"` and `"protein"` in addition to a file name.
* An alphabet file may declare complements (`A U` on a line); reverse complements are then added
  for that alphabet. Without declared pairs, `addRC` is turned off for non-4-letter alphabets as
  before (with a clearer message). Partial or asymmetric pair declarations are errors.
* Protein and 5-letter kernels are verified against the exact-arithmetic reference
  (`tests/oracle`), and frozen as golden cases.

Build: C++17 (`SystemRequirements: C++17`, `src/Makevars`). The R package routes all console output
through `Rprintf` and random numbers through R's RNG (`-DRPACKAGE`); the standalone CLI is unchanged.
