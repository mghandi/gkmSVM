# gkmSVM refactoring plan

**Status: approved (§5). Execution log at the bottom (§7).**
Author: drafted 2026-08-29 against `master` @ `222cc50` (package version 0.80).
Everything marked *(measured)* was produced by `dev/baseline.sh` on this machine
(Apple M3 Ultra, clang 17, R 4.5.3); the script is included so the numbers can be re-checked.

---

## 0. Summary

The goal is a cleaner, faster, more readable gkmSVM that keeps **bit-identical DNA results and
the published user-facing contract**, and that opens the door to per-position alphabets
(two-block, then general **B**).

Three issues were raised, and the review confirms all three plus a set of latent memory bugs that
should be fixed on the way:

| # | Issue | Status found in the code |
|---|---|---|
| 1 | Reliance on unique sequence names | Worse than "an error message": duplicate names inside one file **silently merge two sequences into one kernel row** *(measured: a 60-record positive file with one repeated name produces `npos 59`)*. Names are truncated to 100 characters for matching, and names ≥100 chars **overflow a `new char[100]` buffer**. The kernel matrix itself carries **no identifiers at all** — row order is the only link to the sequences. |
| 2 | Text data files | The kernel matrix is `%e`-formatted text. *(measured, n=5000: 155 MB text vs 47.7 MB binary float32; R load 16.5 s via `read.table` vs 0.25 s via `readBin` — **3.2× smaller, 65× faster**.)* A `-b`/`OutputBinary` flag is already parsed but does nothing. |
| 3 | DNA-only optimisation | Not merely "optimised for b=4": with the shipped `MAX_ALPHABET_SIZE 4`, any alphabet larger than 4 **segfaults** *(measured: `-A` with a 5-letter and a 20-letter alphabet both exit 139)*. `readAlphabetFile` prints an error and `return`s **without aborting**, so execution continues writing past the end of every trie node. From R this kills the whole session. |

Three further findings shape the plan:

* **The GitHub tree and CRAN have diverged.** CRAN ships **0.83.0 (2023-08-20, maintainer Mike Beer)**;
  this repo is **0.80 (2018)**. CRAN has the packaging/`snprintf` hardening, `useDynLib(..., .registration=TRUE)`
  and Bioconductor deps demoted to `Suggests`; **it does not have** this repo's multithreading
  (`std::thread`, `std::atomic`, `-T`), the duplicate-ID checks, or `normalizePath`. The useful CRAN
  changes get folded into this tree as a low-priority parallel track (Phase 0b) — it does not gate
  the refactor, but it should not be forgotten either.
* **The standalone CLI cannot be built from this repo.** `src/` has no `main()` and no Makefile — it
  produces only the R shared object. The published `gkmsvm_kernel` / `gkmsvm_classify` binaries come
  from a separate `gkmsvm-2.0.tar.gz`, i.e. a third lineage. Since those binaries stay supported,
  Phase 0 adds `main()` shims and a Makefile so one tree builds the package *and* the CLI — which
  also gives the golden tests something to drive.
* **`CCalcWmML::calcKernel` reads out of bounds on every single run.** `wm` is `new double[K+1]`
  but the loop indexes `wm[i]` for `i ≤ m ≤ L` (`src/CalcWmML.cpp:30` vs `:145-150`). With the
  defaults L=10, K=6 that is a 4-double overread, ASAN-confirmed on a plain DNA run. It is normally
  masked because `dCombinations(L-m, K-i)` returns 0 for `i>K` — but `0 × NaN = NaN`, so it is a
  latent silent-corruption bug, not only a formality. *(measured: bounding the loop with
  `i<=m && i<=K` makes the DNA output **byte-identical** and the plain DNA run **ASAN-clean** — so it
  is a zero-risk one-line fix. With it applied, the next finding to surface is the 100-character name
  overflow at `mainGkmKernel.cpp:696`, which ASAN otherwise never reaches.)*

The plan is **nine phases (plus a low-priority CRAN track)**, ordered so that every phase is independently reviewable, keeps the
test suite green, and can be released on its own. Phases 0–2 change no behaviour.

---

## 1. What must not break (the compatibility contract)

Anything in this list is covered by a golden test in Phase 0 and may only change with an explicit
version bump and a documented migration.

**R API** — `gkmsvm_kernel`, `gkmsvm_train`, `gkmsvm_trainCV`, `gkmsvm_classify`, `gkmsvm_delta`,
`genNullSeqs`: names, argument order, defaults, and the positional call forms used in the tutorial
(`gkmsvm_kernel(posfn, negfn, kernelfn)` etc.).

**File formats**
* Kernel matrix: tab-separated **lower triangle including the diagonal**, `%e` off-diagonal, the
  literal string `1.0` on the diagonal, one row per sequence, rows = positives in file order then
  negatives in file order, no header, no names. Read back by
  `read.table(fill=TRUE, col.names=paste("V",1:nseq))`.
* `{prefix}_svalpha.out`: `name<TAB>%11.6e alpha`, negative-class alphas negated.
* `{prefix}_svseq.fa`: SV sequences in kernel-row order, named by their FASTA names.
* classify output: `name<TAB>%f score`, one line per test sequence, test-file order.
* `gkmsvm_delta` output: space-separated `name score`, subtraction **by row position**.
* `genNullSeqs` bed/fasta outputs and their `chr_start_end_{pos,neg}_i` naming.

**CLI** — the standalone `gkmsvm_kernel` / `gkmsvm_classify` flag set
(`-l -k -d -m -n -t -a -b -M -L -A -T -R -p` + positional files), which mirrors the published
gkmsvm-2.0 binaries.

**Numerics** — for DNA, every kernel value must match the current output to the last printed digit.

---

## 2. Target architecture

Today `src/` is 63 files / ~13 k lines in one flat directory, of which a large majority is
unreachable from the two R entry points. The proposed layout:

```
src/
  core/                 # pure computation: no file I/O, no globals, no R headers
    alphabet.{h,cpp}          Alphabet: symbols, b, code table, optional complement
    weights.{h,cpp}           GkmWeights (was CCalcWmML): wm, kernel, cTr, h, c, wildcard, mismatch
    lmer_trie.h               template<int B> LmerTrie  +  LmerTrieDyn
    passes.{h,cpp}            iDL pass construction (was CiDLPasses + CbinMMtree/CbinMMtable)
    profile.{h,cpp}           mismatch-profile / kernel accumulation, thread pool
    kernel_matrix.{h,cpp}     the gkmsvm_kernel computation
    scorer.{h,cpp}            the gkmsvm_classify computation
  io/
    fasta.{h,cpp}             streaming reader, no statics, bounded, std::string names
    sequence_set.{h,cpp}      SeqRecord/SequenceSet: ids, names, labels (see §4)
    kernel_file.{h,cpp}       text + binary v1, auto-detecting reader
    model_file.{h,cpp}        SV model: legacy pair + new single-file format
  cli/
    main_kernel.cpp           getopt front-end, byte-compatible with today
    main_classify.cpp
    main_train.cpp            C-SVC training on the kernel matrix (Phase 4b)
  libsvm/                     vendored LIBSVM 3.33 (BSD-3): svm.cpp, svm.h, COPYRIGHT
  r/
    RcppExports.cpp
    api.cpp                   R -> core via a params struct (no argv round-trip)
  legacy/                     quarantined, compiled only with -DGKM_LEGACY (see §3.2)
```

Three cross-cutting rules for the new code:

1. **No mutable globals.** `gMMProfile`, `gMAXMM`, `gLM1`, `gDFSlistT[1000]`, `gDFSMMlist[1000]`,
   `gDFSMMtree[1000]`, `gGTreeLeaves*`, `globalConverter`, `globtmpstr` all become members of a
   `KernelContext` / `ScoreContext` passed down the recursion. This is what makes the code
   re-entrant, thread-safe, and safe to call repeatedly from one R session. (Today the alphabet
   loaded by `-A` *persists into the next call* in the same R session — a real bug.)
2. **RAII everywhere.** `std::vector` / `std::unique_ptr` / a bump allocator for trie nodes;
   no raw `new`/`delete`, no fixed-size arrays sized by `#define`.
3. **One output sink.** A `Reporter` interface with an `Rprintf` implementation for the package and
   a `printf` implementation for the CLI. Today `RPACKAGE` is commented out
   (`src/gkmCommonLib.cpp:208`), so the *package* writes to stdout with `printf` and calls `rand()` —
   both are R CMD check violations.

---

## 3. Phase plan

Each phase = one PR. "Gate" = what must hold before merging.

### Phase 0 — Baseline and safety net (no code changes to the algorithms)

1. **Build the standalone CLI binaries from this repo.** Today `src/` contains **no `main()` and no
   Makefile** — it only builds the R shared object, so `gkmsvm_kernel` / `gkmsvm_classify` cannot be
   produced here at all; the published binaries come from the separate `gkmsvm-2.0.tar.gz` on
   beerlab.org, i.e. a *third* lineage that will keep drifting. Add `src/cli/main_kernel.cpp` and
   `src/cli/main_classify.cpp` (thin `main()` shims over the existing `mainGkmKernel` /
   `mainSVMclassify`) plus a `Makefile`, so one source tree produces both the R package and the
   binaries (three, once Phase 4b adds `gkmsvm_train`). `dev/baseline.sh` already prototypes
   exactly this shim.
   This lands **first** because the golden tests below need a way to drive the code from a shell,
   and because it makes the CLI a permanently supported artefact rather than a fork.
2. **Golden-test corpus.** Fixture FASTAs (DNA small/medium, duplicate names, names >100 chars,
   empty records, lowercase, `N`s, CRLF, b=2 alphabet) × a parameter grid
   (L ∈ {6,8,10,12}, K ∈ {4,6,8}, d ∈ {-1,2,3,4}, t ∈ {0,1,2,3,4}, alg ∈ {1,2}, ±RC, ±pseudocount).
   Freeze the current outputs as reference; compare exactly for text and to 1e-12 for values.
   `testthat` for the R layer, the CLI binaries for the core.
3. **Independent numerical oracle.** Cross-check `c[m]`, `h[m]`, `kernel[m]`, `kernelTruncated[m]`
   and whole small kernels against the exact-arithmetic pure-Python implementation in
   `gkmsvm3/generalb/generalb_gkm.py` (uses `fractions.Fraction`, already unit-tested). This is a
   genuinely independent oracle, and it is also the bridge to Phase 7.
4. **CI**: GitHub Actions matrix (macOS/Linux, gcc/clang, R release/devel) running R CMD check,
   the golden tests, a CLI build + smoke test, an **ASAN+UBSAN** job, and a **benchmark gate**
   (DNA kernel wall-clock, fail on >2 % regression).
5. **Baseline measurements** committed (`dev/baseline.sh`, this document's *(measured)* numbers).

*Gate:* golden tests reproduce today's output byte-for-byte; both CLI binaries build and run in CI.

### Phase 0b — CRAN 0.83.0 reconciliation *(parallel track, low priority, non-blocking)*

CRAN ships **0.83.0 (2023-08-20, maintainer Mike Beer)** while this repo is **0.80 (2018)**; the two
have diverged (see §0). The useful CRAN changes should be folded in here, but this does **not** gate
the refactor and can happen whenever convenient:

* `snprintf` hardening throughout, `useDynLib(gkmSVM, .registration = TRUE)`, Bioconductor
  dependencies demoted from `Imports` to `Suggests`, the duplicated `GenomeInfoDb` entry removed,
  and the `.Rd` fixes — all straightforwardly good and worth taking.
* Where the two conflict, **this tree wins**: CRAN 0.83.0 has no multithreading (`int ***gMMProfile`,
  no `std::thread`, no `-T`), no duplicate-ID handling and no `normalizePath`. Those are this repo's
  improvements and must survive the merge.
* Best done right after Phase 1, so the CRAN hardening lands on top of the latent-bug fixes rather
  than being re-litigated during them. Worth a note to Mike Beer about who publishes the next CRAN
  release once the refactor produces one.

### Phase 1 — Build hygiene and latent-bug fixes (behaviour-preserving)

Confirmed defects to fix, each with a regression test. None of these change correct-input results
(except where noted), so they can land before any restructuring:

| Bug | Location | Fix |
|---|---|---|
| `wm[i]` read past end for `i>K` (ASAN-confirmed, **every run**) | `CalcWmML.cpp:145-150` vs `:30` | bound the inner loop by `min(m,K)` — *verified byte-identical output* |
| `h[L-K+1 .. L]` never initialised, read on all `useTgkm==2` paths | `CalcWmML.cpp:97-104`, `LList.cpp:503,575` | zero-fill `h` |
| `seqBaseId[-1]` written on every zero-length read (i.e. at every EOF) | `Sequence.cpp:248` | early-return on `length==0` (as `readBasic`/`readString` already do) |
| Sequence longer than `maxseqlen` overflows 4 heap buffers | `Sequence.cpp:222-227` | bounds check + clear error |
| Names ≥100 chars overflow `new char[100]` (ASAN-confirmed once the row above is fixed) | `mainGkmKernel.cpp:695,747`; `SequenceNames.cpp:109` | `std::string` (removed entirely in Phase 3) |
| More sequences than `-n maxnumseq` overruns `seqsB/seqname/norm/...` | all four read loops | bounds check + error (removed in Phase 2 by using `std::vector`) |
| `b > MAX_ALPHABET_SIZE` prints an error then **continues** → segfault | `Converter.cpp:160-163` | return a status; abort the call with an R-level error |
| classify passes literal `4` instead of `b` to wildcard/mismatch weights | `mainSVMclassify.cpp:503,509` | pass `b` (**changes results** for `-t 4` with b≠4 — currently wrong) |
| `b==4` `calcnorm` sums all `m` while the score truncates at `maxnmm` | `mainSVMclassify.cpp:744` vs `:723-726` | truncate both consistently (**changes results** when `-d` is small; the b≠4 path is already consistent) |
| Out-of-bounds write when a classify batch has zero unique L-mers | `mainSVMclassify.cpp:663,670` | skip empty batches |
| `svmClassifySimple` never `fclose`s its output (data can be lost in a long-lived R session) | `mainSVMclassify.cpp:409` | RAII file handle |
| Leaks: SV names per batch, `CiDLPasses::passTrees`, `seqsTS` root, `sgi`, whole `svmClassifySimple` | several | RAII |
| `delete` vs `delete[]` mismatch | `LTreeS.cpp:2136`, `SequenceNames.cpp:41` | `delete[]` / containers |
| `RPACKAGE` commented out → package `printf`s to stdout and calls `rand()` | `gkmCommonLib.cpp:208` | the `Reporter` sink + `R_unif_index` |
| Signed-`char` indexing of 256-entry tables | `Sequence.cpp:222,245` | `unsigned char` |
| Alphabet from `-A` persists into the next call in the same R session | `global.cpp:21` | context object (Phase 2); interim: reset per call |
| **Repeated calls abort the R session**: ~38 consecutive `gkmsvm_kernel()` calls in one R process die with SIGABRT (heap corruption from the rows above; each call alone succeeds) *(measured in Phase 0, macOS/R 4.5.3)* | consequence of the reader bugs above | fixed by the rows above; regression test = `tests/testthat` "many calls in one session" (skipped until then) |
| 120-character names **abort R** (exit 134) although the CLI survives the same input by luck *(measured in Phase 0)* | `mainGkmKernel.cpp:695` | same as the ≥100-char row |

*Gate:* golden tests unchanged except the two rows marked **changes results**, which get their own
before/after documentation.

### Phase 2 — Core extraction: kill the globals, one path per job

* Introduce `KernelContext`/`ScoreContext`; delete the `g*` globals and the `[1000]` scratch arrays.
* Move the Phase 0 CLI shims into `cli/` proper and give them the `Options` struct below.
* Replace the **argv round-trip** (`R list → sprintf into 50×5000-char buffers → getopt`) with a
  direct `struct Options` call. The CLI keeps `getopt` and fills the same struct. This removes the
  silent path-length overflow at 5000 chars and makes parameter validation reportable to R
  (today an invalid `K>L` just prints usage to stdout and returns 0 — R sees success).
* Replace `CSequenceNames`' `MAXNSeqs=1000000` fixed arrays (**~16 MB per object**) with vectors.
* Rewrite the FASTA reader without `static` look-ahead state (`Sequence.cpp:183-186`), so two files
  can be read concurrently and the reader is thread-safe.
* **Quarantine dead code.** Unreachable from both R entry points, with evidence:
  `CountKLmers`, `CountKLmersGeneral`, `CountKLmersH`,
  `MLEstimKLmers`, `MLEstimKLmersLog`, `KLmer`, `EstimLogRatio`, `SequenceData` (empty stubs),
  `LKTree`, `CLTreeMemorize`, `GTree`/`GTreeLeafData` (never constructed), `GTree2`/`GTreeLeafData2`
  (guarded by the literal `int UseGTree = 0;`), `CLTreeS::DFST/DFSTn` (dead **and** containing an
  abandoned depth heuristic at `LTreeS.cpp:1580-1587` that would give wrong profiles), plus ~70 % of
  `LTreeS.cpp` that is commented-out history. Move to `legacy/` in this phase (git keeps the history),
  delete in Phase 6. Expected: **~13 k → ~5 k lines**.
* Fold `FAST_TRACK` / `USE_GLOBAL` (never defined) and the `-b` no-op flag out of the live path.
* `mainSVMtrain.cpp` is **kept and rewired** rather than quarantined: Phase 4b replaces its solver
  (`SVMtrain.{h,cpp}` is deleted) while keeping its CLI and output format.

*Gate:* golden tests bit-identical; benchmark within 2 %; ASAN/UBSAN clean.

### Phase 3 — Sequence identity (**issue 1**)

**Design.** Identity becomes explicit and positional; names become metadata.

```cpp
struct SeqRecord {
    int         index;   // 0..n-1  == kernel row
    std::string id;      // assigned at read time: "pos1".."posN", "neg1".."negM", "seq1".. for test
    std::string name;    // original FASTA name, may repeat, may be empty, any length
    int8_t      label;   // +1 positive, -1 negative, 0 test
    int64_t     nlmers;
};
class SequenceSet { std::vector<SeqRecord> records; ... };
```

* Ids are generated on read (`pos1…`, `neg1…`) and are unique **by construction**. Nothing downstream
  ever compares FASTA names to establish identity.
* **`mergeByName` becomes an explicit option, default OFF.** Today duplicate names inside one file
  silently merge sequences (`find_str`, `mainGkmKernel.cpp:683-692,735-745`). That is a legitimate
  feature (multi-exon regions) but a dangerous default. Default OFF = "one FASTA record = one kernel
  row"; `mergeByName=TRUE` restores today's behaviour. *This is the one intentional behaviour change
  in the plan and needs a NEWS entry.* (Users reaching gkmSVM through R cannot be relying on it
  today: the R wrapper `stop()`s on duplicates before ever calling C++.)
* **Row identity is written out.** Binary kernels embed the id/name/label table (§4). For text
  kernels — whose format must not change — a sidecar `<outfile>.index` is written:
  `row<TAB>id<TAB>name<TAB>label<TAB>length<TAB>nlmers`. `gkmsvm_train`/`trainCV` prefer it and fall
  back to today's `read.fasta`-based counting when it is absent.
* **Classify stops matching by name.** The new model file (§4) stores SVs by id/index, so lookup is
  O(1) instead of the current O(records × Nseqs) linear scan with 100-char truncated comparisons.
  For legacy `svalpha.out`/`svseq.fa` pairs the reader keeps name matching but uses a hash map and
  **fails loudly** on a duplicate or missing name instead of silently dropping SVs (today: "the
  sequences for only %d out of %d ... were found", then it scores with a partial model).
* **R side:** remove the three `stop("Error: duplicated sequence ID")` blocks
  (`gkmsvm_kernel.R:40-64`, `gkmsvm_train.R:25-39`, `gkmsvm_trainCV.R:118-132`) — they are no longer
  needed and they cost a full `seqinr::read.fasta` of both files on every call. Emit an informational
  message instead when duplicates are seen.

*Gate:* new tests — duplicate names across pos/neg, duplicate names within a file (both
`mergeByName` settings), 200-character names, empty names, non-ASCII names, names with tabs.
Golden DNA results unchanged for unique-name inputs.

### Phase 4 — Binary formats (**issue 2**)

**`.gkmk` v1 kernel file**

```
offset  field
0       magic       "GKMK" + uint8 version(=1)
5       uint8  dtype        0=float32, 1=float64
6       uint8  layout       0=lower triangle incl. diagonal, 1=full, 2=lower excl. diagonal
7       uint8  flags        bit0 = names table present
8       int32  n, npos
16      int32  L, K, maxnmm, useTgkm, b, addRC, usePseudocnt      (provenance)
        char[] alphabet (length-prefixed)
        names table: n × (id, name, label, nlmers)  (length-prefixed UTF-8)
        payload: n(n+1)/2 values, row-major lower triangle
        uint32 crc32 of the payload
```

* **float32 by default** — the current text format prints `%e`, i.e. ~7 significant digits, so
  float32 (24-bit mantissa, ~7.2 digits) *loses nothing relative to today's file*; float64 available
  via `dtype`. *(measured, n=5000: 155 MB → 47.7 MB, and R load 16.5 s → 0.25 s.)*
* Reader auto-detects text vs binary by magic, so `gkmsvm_train(kernelfn, ...)` works unchanged with
  either.
* New R helpers `read_gkm_kernel(file)` / `write_gkm_kernel(x, file)`.
* `gkmsvm_kernel(..., format = c("text","binary"))`. **Default stays `"text"`** for this release
  (implementing the long-dead `-b` flag along the way); the default flips to `"binary"` in the
  release that bumps the minor version, with `format="text"` kept forever.
* The same treatment for the SV model: a single `.gkmmodel` file (header + ids + alphas + sequences)
  replacing the `{prefix}_svalpha.out` + `{prefix}_svseq.fa` pair, with the legacy pair still written
  and still readable.
* Optional gzip (zlib is already available to R packages) for the text path.

*Gate:* round-trip tests (text→binary→text is identity to `%e` precision); a golden binary file
committed as a fixture so future format changes are caught; cross-endian read test.

### Phase 4b — SVM training in C++ (`gkmsvm_train`)

**Why replace the current trainer.** `CSVMtrain::train` (`src/SVMtrain.cpp:41-101`) is not an SVM
solver. It is the iterative heuristic of Jaakkola, Diekhans & Haussler (2000): a fixed number of
coordinate sweeps (`niter20 = 5`, i.e. 100·N updates) with the multipliers clipped to `[0,1]`,
started from `myrandom()` values. It therefore has

* **no `C` parameter** — the box is hard-coded to `[0,1]`, so the regularisation cannot be changed
  (the CLI's only option is `-n niter20`);
* **no bias/intercept**;
* **no convergence criterion** — it stops after a fixed iteration count whatever the KKT state;
* **no portability of results** — the initialisation calls `rand()`, which is never seeded, so a
  given build is deterministic but the answer depends on the platform's libc `rand()`;
* **a dense model** — *(measured on the 60+60 fixture: all 120 training sequences come out with
  |alpha| > 1e-10, i.e. every sequence is written to `svseq.fa` as a "support vector")*.

That is materially weaker than what R users get from `kernlab::ksvm(type="C-svc", C=...)`.

**What already works.** `mainSVMtrain` is otherwise a complete CLI: it reads the kernel matrix and
both FASTA files and writes `{prefix}_svalpha.out` + `{prefix}_svseq.fa` in **exactly the format the
R/kernlab path produces** (name TAB alpha, negative-class alphas negated, SVs thresholded at
`|alpha| > 1e-10`). So this is a **solver swap inside an existing, format-compatible tool**, not a
new pipeline. It is dead today only because nothing can call it (no `main()`; fixed in Phase 0).

**Library: LIBSVM, vendored into `src/libsvm/`.** Vendor the **current release 3.37**
(2025-12-29); the spike below used 3.33, whose API and file sizes are identical.

| requirement | libsvm 3.33 |
|---|---|
| license | **BSD-3-Clause** — compatible with **both** GPL-2 and GPL-3, so `GPL (>= 2)` is preserved. Obligations: retain `COPYRIGHT`, state that LIBSVM is used |
| footprint | `svm.cpp` 3 312 lines + `svm.h` 105 lines + `COPYRIGHT`; **no dependencies**, no build system needed |
| precomputed kernel | `kernel_type = PRECOMPUTED` — exactly our case, since gkmSVM computes its own kernel |
| output capture | `svm_set_print_string_function()` — routes solver output to `Rprintf`, fixing one of the current R CMD check violations |
| cross-validation | `svm_cross_validation()` — gives `gkmsvm_trainCV` a C++ equivalent for free |
| class weights | `nr_weight`/`weight_label`/`weight` — useful for imbalanced positive/negative sets |
| validation | `svm_check_parameter()` — real error messages instead of silent misbehaviour |
| SV recovery | `svm_get_sv_indices()` returns 1-based indices into the training set — how we map SVs back to sequences |

**Four libsvm gotchas to design around** (none is a blocker, all bite if ignored):

* **`svm_save_model` is useless under `PRECOMPUTED`** — it writes only `0:<index>` per support
  vector, so the model degenerates to (alpha, training index) with no way to score a new sequence.
  gkmSVM must therefore keep writing its own model (the `.gkmmodel` of Phase 4, carrying the SV
  *sequences*), and use `svm_get_sv_indices()` to map back. This is expected, not a surprise.
* **`sizeof(svm_node) == 16`**, so materialising a dense n×n problem for `-t 4` costs ~16·n² bytes —
  **~1.6 GB at n = 10 000, ~6.4 GB at n = 20 000**, i.e. *twice* the raw `double` matrix. The
  precomputed route is therefore the more memory-hungry of the two, which strengthens the case for
  4b-2 below.
* **`rand()`** is used in `svm_cross_validation`'s shuffling and in the probability estimates —
  must be replaced with R's RNG (`unif_rand` under `GetRNGstate`/`PutRNGstate`) for CRAN.
* **`svm_set_print_string_function` is not thread-safe** — one process-wide static pointer. Fine for
  us (training is single-threaded) but it must not be reset per call from concurrent threads.

*Spike (verified, not hypothetical):* a ~60-line driver reading the lower-triangle kernel file that
`gkmsvm_kernel` writes today, feeding it to libsvm as a `PRECOMPUTED` kernel, trains C-SVC on the
motif fixture and returns alphas correctly bounded by ±C and a proper `rho`. `svm_cross_validation`
on the same input gives **5-fold accuracy 0.725 at C = 0.1 and 0.775 at C = 1 and C = 10** — i.e.
the regularisation parameter demonstrably matters on this data, and it is exactly the knob the
current solver cannot express. The integration is a small, well-understood piece of work.

Alternatives considered and rejected:

* **ThunderSVM** — its `PRECOMPUTED` enum is a **stub that silently computes RBF instead**
  (`case SvmParam::PRECOMPUTED://precomputed uses rbf as default`), so it cannot do our job at all;
  dormant since 2019; CMake+Eigen+CUDA. Separately, Apache-2.0 is compatible with GPL-3 but **not
  GPL-2**, and `GPL (>= 2)` does *not* rescue that — it offers a disjunction, so a licensee may
  elect GPL-2.0-only, which we could not grant for Apache-derived code. Taking it would force a
  move to `GPL (>= 3)`, a licence change CRAN requires be flagged.
* **LIBLINEAR** (BSD-3) — no kernels at all; `parameter` has no `kernel_type`. By design.
* **mlpack** (BSD-3) — has **no kernel SVM**, only `linear_svm`.
* **dlib** (BSL-1.0) — viable (precomputed via a custom kernel functor, header-only) but far heavier
  than 3.4k lines for no gain here.
* **Shark** — LGPL-3+, unmaintained, needs compiled Boost.
* **SVMlight** — disqualified outright: "free for scientific use … must not be further distributed
  without prior permission" is both a field-of-use and a no-redistribution restriction, GPL-
  incompatible and incompatible with CRAN's perpetual-distribution requirement.

**CRAN precedent.** `e1071` vendors libsvm (3.23) but does **not** expose the precomputed kernel.
The packages that vendor *and* expose `-t 4` are **WeightSVM** (CRAN, GPL-2 | GPL-3, libsvm 3.23)
and **kebabs** (Bioconductor, libsvm 3.20) — the latter is a *sequence-kernel* package doing
essentially our use case. `kernlab` vendors a 2002 libsvm/BSVM fork and exposes precomputed via
`as.kernelMatrix`. gkmSVM 0.83.0 itself vendors no SVM code; it `Imports: kernlab`.

**Integration in two stages.**

* **4b-1 — precomputed matrix (drop-in).** Read the kernel (text or the Phase-4 binary format),
  hand it to libsvm as `PRECOMPUTED`, write both the legacy two-file output and the new single-file
  model. Memory stays O(n²) as today, and libsvm's `svm_node` representation *doubles* it:
  *(computed)* the `double**` matrix alone is 0.2 GB at n = 5 000 and **3.0 GB at n = 20 000**,
  while the `svm_node` problem handed to `-t 4` adds ~16·n² bytes on top (1.6 GB at n = 10 000,
  6.4 GB at n = 20 000). Feeding libsvm row-by-row from a memory-mapped binary kernel avoids
  holding both at once; beyond ~20 000 sequences only 4b-2 is viable.
* **4b-2 — kernel on the fly (later, with Phase 6).** Give libsvm a kernel callback backed by the
  L-mer trie plus its LRU cache, so the n² matrix never has to exist. This is the same move as
  dropping the mismatch profile from the hot path in Phase 6, and together they remove the memory
  wall that currently caps the method at a few thousand sequences.
  **Prior art to follow: LS-GKM** (Dongwon Lee), the successor to gkmsvm-2.0, does exactly this — it
  embeds a modified libsvm 3.18, *removes* `PRECOMPUTED` entirely, replaces `svm_node` with a
  sequence type, and computes kernel rows on demand through a new `kernel_function_batch` driven by
  a k-mer tree over the training set, called from `SVC_Q::get_Q` while keeping libsvm's `Cache`.
  It also already ships a `gkmmatrix` tool that dumps the Gram matrix. Read it before designing
  4b-2. **Licence note:** LS-GKM is GPL-3.0-or-later, so its *code* cannot be copied into a
  `GPL (>= 2)` package without moving to `GPL (>= 3)` — borrow the design, not the source, unless we
  accept that change.

**Compatibility decisions.**

* **Bias.** libsvm returns `rho`; today's pipeline has no intercept and `gkmsvm_classify` applies
  none. Carry `rho` in the new model format and apply it there; keep the legacy two-file output
  bias-free. Ranking metrics (AUC, and the `gkmsvm_delta` difference) are unaffected by a constant
  offset, so this is safe.
* **kernlab and libsvm are the *same solver*, so agreement should be tight.** `ksvm(type="C-svc")`
  dispatches to `.Call(smo_optim, …)`, which is libsvm's `Solver` class; a source diff shows only two
  substantive differences (a dropped `info()` warning and one extra disjunct in
  `select_working_set`'s stopping test). Second-order working-set selection, shrinking,
  `reconstruct_gradient`, `calculate_rho` and `Qfloat=float` are identical. Measured head-to-head on
  the same Gram matrix with `shrinking=FALSE` on both sides, the coefficient vectors agree to
  **7.2e-07 at `tol=1e-8` and 1.5e-09 at `tol=1e-10`** — the gap scales with the stopping slack,
  which is the signature of one solver run twice, not two solvers. So expect **agreement to the
  stopping tolerance**, not merely "similar".
  *Caveat:* where the dual optimum is non-unique (many bounded SVs) individual alphas can differ a
  lot while the model does not — in one linear-kernel case `max|Δcoef| = 0.70` yet decision values
  correlated 0.9999992. **Validate on decision values and AUC, not on alphas.**
* **Sign/label convention differs.** libsvm makes the first-seen label `+1`; kernlab makes the higher
  sorted factor level `+1`. kernlab's reported `b` is `rho` (its docs call it "the negative
  intercept"). The migration must pin the convention explicitly or the whole model flips sign.
* **`shrinking` needs care — and there may be a kernlab bug here.** Both libraries share a
  byte-identical `Kernel::swap_index` that permutes the `x` pointers only. libsvm's `PRECOMPUTED`
  kernel is permutation-safe by construction, because row identity travels *inside the data*
  (`kernel_precomputed` reads `x[i][(int)(x[j][0].value)]`). kernlab's `kernelMatrix` path instead
  indexes positionally (`*(K + m*i + j)`) into an R matrix that `swap_index` never permutes, so with
  shrinking enabled the solver appears to read permuted indices into an unpermuted `K`. Measured on
  the same data, `ksvm(K, shrinking=TRUE)` gave a different objective, `b`, nSV and training error
  than the feature-space fit, while `shrinking=FALSE` reproduced it to 6.4e-15.
  **This is source reading plus one measurement, not a confirmed upstream bug — verify it
  independently before relying on it, and if it holds, report it to the kernlab maintainers.**
  Two consequences either way: (i) any "reproduce the old numbers" test must keep
  `shrinking=FALSE` on the kernlab side; (ii) gkmSVM already sets `shrinking=FALSE` as its own
  default at both call sites (`gkmsvm_train.R:2,48`, `gkmsvm_trainCV.R:177`) — against kernlab's
  default of `TRUE` — so existing results are sound, but every training run has been paying the full
  no-shrinking cost. libsvm's `PRECOMPUTED` path can safely turn shrinking back **on**, which is a
  free speed-up kernlab cannot currently offer.
* Because the solvers agree, `gkmsvm_train` still gains `backend = c("kernlab", "libsvm")`
  defaulting to `"kernlab"` for one release — but the bar for flipping the default is low, and once
  it flips the heavy `kernlab` dependency can leave `Imports` entirely.
* **CLI.** `gkmsvm_train` keeps its positional arguments and gains `-c C` (and `-w` class weights,
  `-v nfold` for CV). `-n niter20` is accepted and ignored with a deprecation notice.

*Gate:* on the tutorial CTCF dataset the libsvm backend reproduces the kernlab AUC within noise;
the R backend and the CLI produce the same model from the same kernel; model-file round-trip tests;
`CSVMtrain` deleted rather than quarantined.

### Phase 5 — Alphabet generalisation, single B (**issue 3, part 1**)

The constraint is explicit: **DNA performance must not regress**. The design achieves that by making
the alphabet size a *compile-time constant on the fast path* and a runtime value only on the slow path.

```cpp
template<int B> class LmerTrie { std::array<Child, B> child; ... };   // B=2,4 instantiated
class LmerTrieDyn { std::vector<Child> child; ... };                  // any b

// single runtime dispatch at the top of the computation:
switch (alphabet.size()) {
  case 2:  return runKernel<LmerTrie<2>>(ctx);
  case 4:  return runKernel<LmerTrie<4>>(ctx);   // <- today's code, same constants, same codegen
  default: return runKernel<LmerTrieDyn>(ctx);
}
```

* DNA takes the `B=4` instantiation: identical loop bounds and identical node layout to today, so the
  benchmark gate (≤2 %) is expected to pass trivially. Today's `MAX_ALPHABET_SIZE=24` workaround for
  protein is exactly what this replaces — it would make **every** node carry 24 pointers even for DNA.
* `Alphabet` validates at runtime and returns an error instead of crashing: the current
  `readAlphabetFile` behaviour (print, `return`, continue, segfault) is the direct cause of the
  measured `b=5`/`b=20` crashes.
* Reverse complement becomes a property of the alphabet, not the assumption `complement(i)=b-1-i`
  (`Converter.cpp:65`). DNA/RNA get built-in maps; an alphabet file may declare pairs; alphabets
  without a complement reject `addRC=TRUE` with a clear message rather than silently disabling it.
* The XOR/`CLList` algorithm (`alg=1`) stays **b=4 only** and is documented as such — it already is,
  via the forced `alg=2` for b≠4, and its 2-bit packing (`<<2`, `&3`, `HamDist` tables) is not worth
  generalising.
* `gkmsvm_kernel(alphabetFN=)` gains a first-class R interface: `alphabet="dna"|"rna"|"protein"|<file>`.

*Gate:* protein (b=20) and b=5 runs **complete and are numerically verified against
`generalb_gkm.py`**; DNA benchmark within 2 %; DNA golden tests bit-identical.

### Phase 6 — Performance

Measured starting point *(n=500 sequences × 500 bp, L=10 K=6 d=3)*:

| threads | wall | CPU | efficiency |
|---|---|---|---|
| 1 | 6.55 s | 6.53 s | — |
| 4 | 2.04 s | 7.62 s | 80 % |
| 20 | 0.97 s | 16.62 s | **34 %** |

* **Atomics are not the bottleneck.** The comment at `global.h:25` claims removing `MULTI_THREAD_SAFE`
  is ~2× faster; *measured*, single-threaded with and without atomics: 6.54 s vs 6.75 s — no
  difference. The scaling loss comes from (a) the per-pass tree clone, (b) three heap allocations per
  internal node visit in `DFSTiDL` (`LTreeS.cpp:1662-1666`), and (c) static `j % nThreads` pass
  assignment with unequal pass costs. Fix with an arena allocator for DFS lists and a work-stealing
  queue over passes. Then re-evaluate whether atomics can be dropped in favour of per-thread tiles.
* **Drop the mismatch profile from the hot path.** `gMMProfile[n][d+1][n]` is the memory wall:
  *(computed)* 0.47 GB at n=5000 and **7.45 GB at n=20 000**. At a leaf the mismatch count `m` is
  known, so `K[i][j] += c[m]` can be accumulated directly, removing the `(d+1)` factor and the
  second pass entirely; keep the full profile only when `usePseudocnt` or when a caller asks for it.
  With lower-triangle-only storage: 7.45 GB → **0.8 GB** at n=20 000.
* **Tile the computation** so peak memory is bounded independently of n (the classify path already
  batches; the kernel path does not).
* Only then consider micro-optimisation (SoA node layout, prefetch, `__builtin_popcount`).

*Gate:* golden tests bit-identical (the accumulation reorder must be validated for
floating-point associativity — use a fixed summation order); benchmark improves; memory ceiling
documented.

### Phase 7 — Two-block, then general B (**issue 3, part 2**)

Only after Phases 3–6 are released and stable.

* Generalise the coefficient tables from a single alphabet size to a vector
  **B = (b₁,…,b_ℓ)**, reusing the theory and the *already-tested* exact implementations in the
  `gkmsvm3` repo (`generalb/generalb_gkm.py`, `twoblock/`): `H(u,w) = (1/Πbᵢ) Σ_j e_j(w)` computed in
  O(ℓk) per pair, block-wise mismatch classes, the truncated filter as a Kronecker contraction.
  Note that today's `CCalcWmML::calcc` is precisely the single-block specialisation of the six-fold
  sum in `twoblock/README`, and `dCombinations` already implements the generalised-binomial
  convention that the general case needs — so the two code bases meet cleanly.
* The trie generalises to a **per-position** alphabet size: `child` count depends on depth. Two-block
  first (`B = (b₁,)·L₁ + (b₂,)·L₂`), which covers DNA+methylation and RNA+structure, then general B.
* Validation: every new coefficient table is checked against `generalb_gkm.py` at exact
  (`Fraction`) precision, and end-to-end kernels are checked against the pure-Python
  `kernel_matrix()` on small inputs. The `gkmsvm3/experiments/` harness then provides real-data
  regression tests.

*Gate:* general-B kernels match the exact Python implementation; DNA path untouched and still
bit-identical.

---

## 4. Effort and sequencing

| Phase | Content | Rough size | Risk | Depends on |
|---|---|---|---|---|
| 0 | CLI build, golden tests, CI, ASAN, numerical oracle | ~1 week | low | — |
| 0b | Fold in the useful CRAN 0.83.0 changes *(low priority, non-blocking)* | ~2 days | low | 1 |
| 1 | Latent-bug fixes | ~3 days | low | 0 |
| 2 | De-globalise, RAII, quarantine dead code | ~1–2 weeks | medium (large diff, zero behaviour change) | 0,1 |
| 3 | Sequence identity | ~1 week | medium (one deliberate default change) | 2 |
| 4 | Binary formats | ~1 week | low (additive) | 3 |
| 4b | C++ SVM training via vendored libsvm | ~1 week | low–medium | 0, 4 |
| 5 | Runtime alphabet, templated trie | ~1–2 weeks | medium | 2 |
| 6 | Performance | ~1–2 weeks | medium (numerics must be pinned) | 2,5 |
| 7 | Two-block → general B | project-sized | high | 5,6 |

Phases 3, 4 and 5 are independent of each other once 2 has landed and can be done in parallel or
reordered by priority. Phase 0b is independent of everything after Phase 1.

---

## 5. Decisions taken (2026-08-29)

All six open questions are resolved; the phases above already reflect them.

1. **CRAN reconciliation — do it, but it is not high priority.** Fold the useful and reasonable
   0.83.0 changes into this tree (Phase 0b), where this tree's own improvements — multithreading,
   duplicate-ID handling, `normalizePath` — take precedence on conflict. It does **not** gate the
   refactor; best scheduled right after Phase 1.
2. **`mergeByName` default OFF** — one FASTA record = one kernel row. Today's silent merging of
   same-named records becomes opt-in via `mergeByName=TRUE`, with a NEWS entry. (No R user can be
   relying on it today: the wrapper `stop()`s on duplicates before C++ is reached.)
3. **Binary output is opt-in first, default later.** Implement the long-dead `-b` flag and add
   `format = c("text","binary")` with `"text"` as the default for this release; the auto-detecting
   reader lands at the same time, so flipping the default at the next minor bump is a one-line
   change and old kernels keep working forever.
4. **Fix both result-changing bugs** — the classify path passing a literal `4` instead of `b` to the
   wildcard/mismatch weights, and the `b==4` norm summing all mismatch levels while the score
   truncates at `maxnmm`. Both are defects; document the before/after in NEWS.
5. **Move to C++17** — `Makevars` with `CXX_STD = CXX17` and `SystemRequirements: C++17` in
   DESCRIPTION (CRAN accepts it unconditionally). `if constexpr`, `std::string_view` and structured
   bindings make the templated-trie code in Phase 5 substantially cleaner.
6. **Keep shipping the standalone CLI** — `gkmsvm_kernel`, `gkmsvm_classify` and (per decision 7)
   `gkmsvm_train` stay supported artefacts so the C++ version can be run without R. Because they
   cannot be built from this repo at all today (no `main()`, no Makefile), adding them is the
   *first* item of Phase 0, where they also serve as the driver for the golden tests.

7. **Train in C++ too** — `gkmsvm_train` becomes a third supported binary, with the hand-rolled
   `CSVMtrain` heuristic replaced by **vendored LIBSVM** (BSD-3, no dependencies). See Phase 4b.

## 6. Reproducing the measurements

```bash
dev/baseline.sh            # builds a standalone driver from src/, runs the measurements above
```

It builds `src/*.cpp` with a small `main()` shim, then reports: DNA kernel timing vs thread count,
atomic vs non-atomic, duplicate-name merging, the `b>4` crash, name-length overflow (under ASAN),
and text-vs-binary kernel size and R load time.

## 7. Execution log

* **Phase 0 — done (PR: `phase0/baseline-safety-net`).** `Makefile` + `src/cli/main_{kernel,classify,train}.cpp`
  build the three binaries from this tree (`make`, `make ASAN=1`). Golden corpus: `tests/golden/`
  (22 deterministic fixtures from `make_fixtures.py`, 89 cases in `cases.tsv`, 86 outputs frozen from
  `222cc50` byte-for-byte, 3 known crashes tracked as `xfail-phase1`/`xfail-phase5` tags that fail the
  suite the moment they stop crashing). Oracle: `tests/oracle/oracle_check.py` against a vendored copy of
  `gkmsvm3/generalb/generalb_gkm.py` in `dev/oracle/` (the upstream repo is private, so CI cannot clone
  it) — 206 coefficient values and 4 068 kernel entries agree to the printed precision, for b ∈ {2,4},
  t ∈ {0,1,2}, ±RC, both algorithms. R layer: `tests/testthat/` replays the corpus through
  `gkmsvm_kernel`/`gkmsvm_classify` (one `Rscript` per case until Phase 1, see the two new rows in the
  Phase 1 table), plus format tests for `gkmsvm_train`/`gkmsvm_delta`; `dev/scratch_install.sh`
  installs the package with the Bioconductor `Imports` trimmed so this runs without BSgenome.
  CI: `.github/workflows/ci.yml` — CLI+golden+oracle on ubuntu gcc/clang and macOS clang, ASAN+UBSAN
  (advisory until Phase 1), the R scratch install + testthat, and a benchmark gate (`dev/bench.sh`,
  best-of-5 single-thread, PR vs. base on the same runner, 2 %). A full `R CMD check` job is deferred to
  Phase 0b when the Bioconductor deps become `Suggests`. Baseline on this machine: 6.33 s.
* **Phase 1 — done (PR: `phase1/latent-bug-fixes`, stacked on Phase 0).** Every row of the Phase 1 table
  is fixed (see `NEWS.md` for the user-facing list and the before/after numbers of the two
  result-changing fixes). Golden: all kernel outputs byte-identical; exactly the seven classify cases
  with `alg=2, d=3` changed (the approved norm truncation), verified independently by a new classify
  section of the oracle (exact scores for d ∈ {3, L}, t ∈ {0,1,2}, ±RC, DNA and b=2). `k_emptyrec` runs
  and is frozen; the b=5 alphabet and the new `-m`/`-n`/missing-alphabet cases are `expect-error`
  golden cases (clean exit 1 + ERROR message, no output file, and an R `stop()`). Whole corpus
  ASAN+UBSAN-clean (`detect_leaks=0`; the sanitizer CI job is now required). C++17 + `src/Makevars`
  (`-DRPACKAGE`: `Rprintf`, `R_unif_index`). Worker threads no longer print (Rprintf is main-thread
  only). The R golden tests run in-process again, plus new regression tests (60 calls in one session,
  alphabet does not leak between calls, alg=1 classify output complete).
* **Phase 2a — done (PR: `phase2/core-extraction`, stacked on Phase 1).** Dead code quarantined into
  `src/legacy/` (13 translation units, the `UseGTree=0` branch, `CLTreeS::addToGTree`, the GTree
  globals; live source 13.2k → 10.1k lines); the never-defined `FAST_TRACK`/`USE_GLOBAL` blocks removed
  (comment-aware — a naive unifdef unbalanced one of `LTreeS.cpp`'s comment blocks and silently changed
  results, which the golden corpus caught); `src/gkmOptions.h` with `OptsGkmKernel`/`OptsSVMClassify`,
  `gkmKernelRun()`/`svmClassifyRun()` shared by the CLI (getopt) and the R bindings (direct struct fill —
  no more `sprintf` into 50×5000-char buffers); parameter validation is now an error (exit 1 / R
  `stop()`, six new `expect-error` cases) instead of "print usage, return 0"; `CSequenceNames` on
  `std::vector` (was two 1 M-entry fixed arrays, 16 MB per object); the FASTA reader's look-ahead state
  is per object (was function-static). Gates: golden 108/108 byte-identical, ASAN/UBSAN clean, oracle OK,
  testthat 10/10, bench 6.41 s (+1.2 % vs Phase 0).
* **Phase 2b — TODO (not started).** The remaining Phase 2 items: `KernelContext`/`ScoreContext` replacing
  `gMMProfile`, `gMMProfile0`, `gMAXMM`, `gLM1`, `gDFSlistT/gDFSMMlist/gDFSMMtree/gDFSlist[1000]`,
  `globalConverter`, `globtmpstr` (≈150 use sites, ~110 of them inside the `CLTreeS`/`CLTreef` recursion in
  `LTreeS.cpp`/`LTreef.cpp` — needs a context parameter threaded through `DFST*`/`cloneReorder`); deleting
  the dead `CLTreeS::DFST/DFSTn` and the commented-out history in `LTreeS.cpp`; the `-b` flag is kept as a
  parsed-but-inert option for Phase 4 to implement. Do 2b before Phase 6 (the thread-pool work needs the
  context object); Phases 3 and 4 do not depend on it.
* **Phase 3 — done (PR: `phase3/sequence-identity`, stacked on Phase 2a).** `src/SequenceSet.h`
  (`SeqRecord`, `seqRecordId`, `writeIndexSidecar`); `mergeByName` option (default off, `-N`, R
  argument, documented) — `-N` reproduces the previous merged output byte-for-byte (`k_dupnames_merge`
  == the old `k_dupnames`); `<outfile>.index` written by both kernel paths and frozen for eight golden
  cases (the runner compares it whenever one exists); `CSequenceNames::nextSeq` rewritten on an
  `unordered_map` with loud failures (duplicate in the alpha file, missing from or repeated in the
  FASTA — three `expect-error` cases); R: `mergeByName=`, duplicate-name `stop()`s replaced by a
  `message()`, `gkmsvm_train`/`trainCV` prefer the sidecar and fall back to unique-name counting for
  old kernels, ids replace names in the legacy model pair when names are not unique. Name edge cases
  (empty, tab, non-ASCII, 200 chars, cross-file duplicates) frozen with their sidecars. Not done: the
  model file that stores SVs by id is Phase 4's `.gkmmodel`; the full `SequenceSet` class (readers
  returning records) is folded into Phase 2b's reader rewrite.
