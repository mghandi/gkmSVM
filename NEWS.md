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

### Phase 2a — core extraction (behaviour-preserving)

* Invalid parameters (`K > L`, `maxnmm > L`, `useTgkm` outside 0–4, `alg` outside 0–2, `batchSize < 1`)
  are now errors (`stop()` in R, exit status 1 on the command line). They used to print the usage
  text and return success, leaving no output file.
* The R bindings pass their arguments to the C++ core directly instead of re-serialising them into a
  command line; file paths are no longer silently truncated at 5000 characters.
* Thirteen unreachable source files were moved to `src/legacy/` and are no longer compiled.
* Kernel and classify outputs are byte-identical to the previous version.

Build: C++17 (`SystemRequirements: C++17`, `src/Makevars`). The R package routes all console output
through `Rprintf` and random numbers through R's RNG (`-DRPACKAGE`); the standalone CLI is unchanged.
