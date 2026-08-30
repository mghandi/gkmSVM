# Phase 7 — gkm-SVM over heterogeneous alphabets (two-block, then general B)

Status: **plan approved for execution 2026-08-29**; sub-phases 7a → 7e below, one branch + PR each,
stacked on `master` (0.90.0). Companion to `dev/REFACTORING_PLAN.md` §3 "Phase 7"; this file is the
detailed design and the record of what was decided and measured.

Theory and exact reference implementations come from the (private) `gkmsvm3` repository:
`generalb/generalb_gkm.py` (any B; vendored verbatim as `dev/oracle/generalb_gkm.py`),
`twoblock/twoblock_gkm.py` (two alphabet sizes; an *independent* second reference, to be vendored
as `dev/oracle/twoblock_gkm.py`) and `experiments/01..07` (the application ladder).

---

## 1. What is being built

Today an input record is one string over one alphabet of size b, and a kernel between two
sequences is `Σ_m N_m · c[m]` where `N_m` is the number of L-mer pairs with `m` mismatches
(the *mismatch profile*, filled by the trie DFS) and `c[m]` a coefficient table
(`CCalcWmML`: `c` = full filter, `cTr` = truncated filter, `h` = gapped k-mer counts).

Phase 7 makes an input record **T aligned tracks** (e.g. DNA + methylation flag, residue + secondary
structure, several SAX-symbolised sensor channels), each over its own alphabet α_t. A window of
length L over T tracks is an ℓ-mer with ℓ = T·L whose positions carry alphabet sizes

    B = (|α_1|,)·L + (|α_2|,)·L + … + (|α_T|,)·L          (track-major concatenation)

Positions with the same alphabet size form a **block** (blocks are defined by *size*, not by track:
two DNA tracks form one block of 2L positions). With r blocks of lengths ℓ_1..ℓ_r, a pair of ℓ-mers
falls into a **mismatch class** `m = (m_1,…,m_r)`, `0 ≤ m_j ≤ ℓ_j`; there are `C = Π(ℓ_j+1)` classes.
The kernel becomes

    K(S_1,S_2) = Σ_m N(m) · c(m)

with the r-block profile `N(m)` and, per kernel type (`-t`):

| `-t` | name | c(m) | today's r = 1 table |
|---|---|---|---|
| 0 | filter | `H(m)` — Theorem 2 of Mohammad-Noori, Ghareghani & Ghandi: `H(u,w) = (1/Πb_i) Σ_{j≤k} e_j(w)`, `w_i = b_i−1` at matches, `−1` at mismatches | `c` |
| 1 | truncated filter | `⟨g, (E_1(m_1)⊗…⊗E_r(m_r)) g⟩` with `g` = H with the dominance rule (`H(m) ≤ 0` zeroes every `a ≥ m` componentwise) | `cTr` |
| 2 | gkm counts | `C(ℓ−Σm_j, k)` | `h` |
| 3, 4 | wildcard / mismatch (Leslie–Kuang) | **not generalised** (error for T > 1) | — |

For r = 1 every formula reduces to today's (`CCalcWmML::calcc` *is* the one-block six-fold sum; the
dominance rule *is* the 2014 "first non-positive coefficient" truncation; `dCombinations` already uses
the generalised-binomial convention `C(−1,m) = (−1)^m` that the diagonal needs). This is the
consistency test that anchors the new tables to the old ones.

### Why the trie needs almost no change

The DFS (`CLTreeS::DFSTiDL/DFSTnIDL`, `CLTreef::DFST/DFSTn`) carries one integer per live path,
`curMismatchCnt`, incremented by 1 at a mismatch and used as the profile row at the leaf. The total
mismatch bound (`-d`) is enforced separately, by the pass tree (`CbinMMtree`, child1 = "one more
mismatch allowed"). Therefore replacing

    curMismatchCnt + 1     →     curMismatchCnt + step[depth]

with `step[depth] = stride_{block(position at depth)}` and `stride_j = Π_{j'<j}(ℓ_{j'}+1)` turns the
counter into the **class index** `Σ m_j·stride_j`, and the rest of the DFS is unchanged. For r = 1
every `step` is 1, so the DNA path is *the same code with the same values* (golden tests must remain
byte-identical; the benchmark must not move). The trie node type (`daughter[MAX_ALPHABET_SIZE]`,
instantiated for 4 and 32) is reused: a per-position alphabet only needs `max_t |α_t|` children, so
DNA + binary flags runs on the **b4 trie**.

---

## 2. Design decisions (D1–D14)

* **D1 — input format.** Multi-track FASTA (`.mfa`, the format of `generalb_gkm.read_mfa` and of
  `twoblock`'s `.2fa`): a header line followed by exactly T non-empty lines of equal length, one per
  track. A symbol outside its track's alphabet makes every window covering it **skipped** (alignment
  must be preserved; the single-track reader *drops* such characters instead — that behaviour is
  untouched for T = 1). `N` is therefore "skip" for DNA tracks, as in gkmsvm3. Records with a
  duplicate name follow the Phase 3 rules (one record = one row; `-N` merges).
* **D2 — alphabet specification.** `-A spec[,spec…]`, one spec per track: a keyword (`dna`, `rna`,
  `protein`), a path to an alphabet file (today's format), or `=SYMBOLS` for a literal alphabet
  (`=01`, `=NUM`, `=abc`). T = number of specs; a single spec is exactly today's `-A`. In R:
  `alphabets = c("dna", "01")` — an entry that is a keyword or an existing file is passed through,
  anything else is sent as a literal (`=…`).
* **D3 — position order and blocks.** ℓ-mer = track-major concatenation of the T windows, blocks by
  alphabet size in order of first appearance (identical to `generalb_gkm.Blocks(B)`), class index
  `Σ m_j·stride_j`. The iDL pass permutations act on the concatenated ℓ-mer, so the per-depth step
  vector is computed per pass as `step[d] = stride[block(passOrder[j][d])]`.
* **D4 — profile rows.** `KernelContext`/`ScoreContext` get `nclasses`; row `idx` is allocated only if
  its class is *reachable* (`Σm_j ≤ d`), unreachable rows are `NULL` and never touched. For r = 1 this
  is exactly today's `maxmm+1` rows. `calcinnerprod`/`svmScoreunorm` sum over reachable rows.
  Tiling (`-r`, `tileMemoryMB`) counts reachable rows.
* **D5 — `-d` semantics.** Bound on the **total** number of mismatches `Σm_j` (the quantity the pass
  trees bound). Automatic (`-1`): filter → ℓ, gkm → ℓ−k, truncated → largest `|m|` with `c(m) ≠ 0`,
  which for r = 1 equals the current `2·(kernelTruncatedLength−1)`.
* **D6 — reverse complement.** As in gkmsvm3: track 1 complemented **and** reversed, every other
  track reversed (a methylation track must mark both bases of a CpG, `meth2s`, for this to be
  consistent). Requires track 1 to have complements (`hasComplement`), otherwise the existing "use
  -R" error. `-R` disables.
* **D7 — kernel types.** `-t 0/1/2` generalised; `-t 3/4` refused for T > 1 with a clear message.
  Pseudocounts (`-p`) refused for T > 1 (their `b^L` closed form is single-block).
* **D8 — coefficient tables.** New module `src/GeneralB.{h,cpp}` (blocks, strides, `H_table`,
  `gkm_table`, `truncate` (dominance rule), `block_enum`, `contract`, `truncated_table`, auto-`d`),
  in `double` (`long double` accumulation where cheap). **T = 1 keeps `CCalcWmML`** as the source of
  truth (DNA byte-identical by construction); T ≥ 2 always uses `GeneralB`, even when r = 1
  (two DNA tracks). A unit test asserts `GeneralB` == `CCalcWmML` for r = 1 to 1e-12 relative, and
  the oracle asserts `GeneralB` == `Fraction` tables to 1e-9.
* **D9 — trie instantiation** is selected by `bmax = max_t |α_t|` (b4 when every track has ≤ 4
  symbols, b32 otherwise; > 32 is an error as today).
* **D10 — passes.** The greedy iDL design enumerates all mismatch patterns with `≤ d` ones
  (`CbinMMtable`, `Σ_{i≤d} C(ℓ,i)` rows × ℓ ints). That is fine up to ~2^20 patterns (ℓ = 20, d = 20:
  80 MB) but not for ℓ = 30–60 with a large `d` (gkmsvm3's experiments use `kind='filter'`, i.e.
  `d = ℓ`, and small k for joint words). Above a threshold (measured in 7d, provisionally 2^20
  patterns) switch to **prefix-split passes**: fix the match/mismatch pattern of the first q positions
  (one pass per pattern with `≤ d` ones, e.g. q = 6 → ≤ 64 passes), and represent the remaining
  levels by an *implicit* bound tree — a DAG of `(depth, mismatches used)` nodes with `child1` present
  iff `used < d` — that plugs into the DFS through the same `child0/child1` pointers. Every pattern is
  enumerated exactly once, the passes are independent (parallel with the existing scheduler), and no
  table is materialised. DNA defaults (ℓ = 10, d = 3: 176 patterns) never take this path.
  Optional refinement, measured before keeping: per-position match probability `p_i = 1/b_i` in
  `CiDLPasses::calcCost` instead of the single `p = 1/b`.
* **D11 — model files.** `.gkmmodel` gets a header line `#alphabets dna,=01`; SV records are
  multi-track (header + T lines). `gkmsvm_classify` takes the alphabets from the model when `-A` is
  absent, and refuses a conflicting `-A`. The legacy `svseq.fa` + `svalpha.out` pair is also written
  in `.mfa` form for T > 1. Text kernels are unchanged; `.gkmk` stores the spec list in its existing
  `alphabet` string field and `b = bmax` (format version stays 1; the R reader already returns the
  string).
* **D12 — memory.** Profile = reachable classes × n(n+1)/2 × 4 bytes (e.g. n = 4 000, ℓ = 20, all 121
  classes: 3.9 GB → 4 tiles at the default 1 GB budget; each tile re-runs the passes). Accumulating the
  kernel directly at the leaf ("kernel on the fly", plan item 4b-2) would make this n²·8 bytes
  regardless of C; it is the follow-up if 7e's measurements show tiling dominating.
* **D13 — R interface.** `alphabets=` on `gkmsvm_kernel`, `gkmsvm_classify`, `gkmsvm_train`
  (which must copy `.mfa` records for the SVs: new exported helpers `read_mfa()` / `write_mfa()`),
  `gkmsvm_trainCV`. `gkmsvm_delta` stays DNA-only (documented). `genNullSeqs` unchanged.
* **D14 — experiments harness.** `gkmsvm3/experiments/common/methods.py` gets an optional C++
  backend (environment variable `GKMSVM_BIN` = directory with `gkmsvm_kernel`): it writes the samples
  as `.mfa`, runs the binary, reads the normalised kernel. The pure-Python path stays the reference;
  each rerun reports both. Changes to gkmsvm3 go to a branch (`cpp-backend`) there, not to `main`.

Out of scope for Phase 7: a user-specified per-position B that is not derived from tracks (e.g. a
different alphabet at each position of one track); `gkmsvm_delta` for multi-track; the XOR/hash
algorithm (`-a 1`, DNA-only by design); wildcard/mismatch kernels for T > 1.

---

## 3. Sub-phases

Each sub-phase: branch off the previous one (stacked PRs), the standard gates (golden corpus
byte-identical for every existing case, oracle OK, ASAN/UBSAN clean, `testthat`, `R CMD check
--as-cran`, CI green, benchmark within noise), numbers in the PR body from commands actually run.

### 7a — Foundation (no behaviour change)

1. `src/GeneralB.{h,cpp}`: `AlphabetVector` (B, blocks, strides, `nclasses`, `classOf(pos)`),
   `H_table`, `gkm_table`, `truncate`, `block_enum`, `contract`, `truncated_table`, `autoMaxmm`,
   reachable-class enumeration. Pure functions, no I/O.
2. Test binary `build/gkm_selftest` (not shipped): `GeneralB` vs `CCalcWmML` for
   b ∈ {2,4,5,20}, (L,K) grid, all three tables; exits non-zero on mismatch. Wired into `make test`.
3. `tests/oracle/oracle_check.py`: new section comparing `gkmsvm_kernel --print-coeffs` (new hidden
   flag: print the table for `-A spec,…` and exit) with `generalb_gkm.kernel_coefficients` for
   B ∈ {(4,)·6+(2,)·6, (4,)·4+(2,)·4+(3,)·4, (20,)·4+(3,)·4, (2,)·8, (4,)·10}, kinds filter/truncated/gkm.
4. DFS plumbing: `const int *step` argument through `DFSTiDL/DFSTnIDL` (kernel) and `DFST/DFSTn`
   (classify); `KernelContext.nclasses`, `ScoreContext.nclasses`; `calcinnerprod`/`svmScoreunorm`
   over reachable rows; tiling accounting by reachable rows. For T = 1: `step ≡ 1`, `nclasses = maxmm+1`.
   Benchmark before/after (`make bench`, 3 runs each); if the indirection costs > 2 %, template the
   DFS on a `Uniform` bool instead.
5. Gate: 151/151 golden byte-identical, oracle (old sections + new tables) OK, ASAN clean, bench flat.

### 7b — Two blocks: kernel (DNA + one annotation track)

1. `src/MultiTrack.{h,cpp}`: `TrackAlphabets` (parse `-A` list per D2; per-track `CConverter`),
   `readMfaRecord` (D1), window enumeration into a flat `int[ℓ]` per position (skipping windows with
   out-of-alphabet symbols), reverse complement per D6. **Restricted to T = 2 in this sub-phase**
   (T > 2 → "not yet supported"), so that `twoblock_gkm.py` is an exact independent oracle.
2. `gkmKernelSuffixTree`: when T > 1, build the trie from the flat windows (`addSeq` per window; the
   `.index` sidecar records `nlmers` = windows kept), tables from `GeneralB`, `d` per D5, kernel
   types per D7, output as today (text / `.gkmk`, D11). Dispatch on `bmax` (D9).
3. Fixtures: `twoblock/examples/{pos,neg}.2fa` (vendored), synthetic `.mfa` with `N`s, with/without
   RC, `-t 0/1/2`, `-d` bounded and unbounded. Oracle: entry-wise vs `generalb_gkm.kernel_matrix`
   **and** `twoblock_gkm.kernel_matrix` (rtol 2e-6, as for DNA). Golden cases frozen from the C++
   only after the oracle agrees.
4. Experiments (first pass, via D14): `twoblock/examples/run_example.py` (LOO centroid AUC table:
   gkm/filter/truncated 1.000, DNA-only 0.502, methylation-only 0.713), gkmsvm3 experiment 01
   (studies A and B) and 02 (study B) with the C++ backend — AUCs must match the Python run to the
   float tolerance of the kernel, runtimes recorded side by side.

### 7c — Two blocks: classify, train, R

1. `svmClassifySuffixTree` for T = 2: SV trie from the `.gkmmodel`/`svseq` multi-track records,
   test records from `.mfa`, `ScoreContext` over classes, `calcnorm` via the general self-profile
   (build the sequence's own windows as a trie, run the class-aware DFS against itself) — used for
   every T > 1 case; the b = 4 legacy/`legacyNorm` paths untouched for T = 1.
2. `.gkmmodel` header `#alphabets` (D11), `gkmsvm_train` (C++ and R) copies SV records from the
   `.mfa` inputs by row index (the sidecar), `read_mfa`/`write_mfa` in R.
3. R: `alphabets=` on kernel/classify/train/trainCV; Rd pages; `testthat` (`test-multitrack.R`:
   kernel == oracle fixture, train → classify round trip, `gkmsvm_train(backend="libsvm")` on a
   multi-track kernel).
4. Oracle: classify scores vs exact `Σ_j α_j K(x_i,x_j)` for the two-block SV fixture.

### 7d — General B

1. Lift the T = 2 restriction; `three_tracks.mfa` and a residue + 2-state + 3-state fixture
   (experiment 04's shape); blocks with r = 3 and with two tracks of equal size (r = 1, ℓ = 2L).
2. Prefix-split passes with the implicit bound DAG (D10) when the pattern count exceeds the
   threshold; unit test that the union of the pass patterns is exactly the `≤ d` set (small ℓ), and
   a kernel equality test between the two pass designs on a fixture where both are feasible.
3. Per-position `p_i` in the greedy cost (D10, optional): keep only if the pass count or wall clock
   improves on a two-track benchmark.
4. Oracle for r = 3 and for `.mfa` with T = 4; golden cases; ASAN.
5. Benchmark table in the PR: DNA (unchanged), DNA+binary (ℓ = 20, `-t 0/2`, d auto), three tracks
   (ℓ = 30), n = 250+250 and 1 000+1 000, T = 1/4/8 threads; Python reference timings from
   `generalb/EXAMPLES.md` and `experiments/README.md` for comparison.

### 7e — Experiments at scale and documentation

1. gkmsvm3 reruns with the C++ backend: 01, 02 (exact agreement), 03 `BasicMotions`/`RacketSports`
   (per-channel SAX alphabets, r up to 3, product alphabet n/a — the positive result), 06 CTCF +
   methylation at n = 300 (agreement with the Python numbers, 0.640–0.659 GM12878) and then at the
   full prepared set (`sites_GM.tsv`: 400 + 400) and, data permitting, ≥ 2 000 + 2 000 — the scale the
   Python harness could not reach; 07 (r = 3) at full size. Record runtimes in gkmsvm3's
   `results/RESULTS.md` through its `Runner` (never typed by hand).
2. `tutorials/gkmsvm-multitrack-tutorial.md`: DNA + methylation end to end (kernel → train →
   classify → per-window scores), with the encoding advice from experiment 02 (`meth2s`).
3. `NEWS.md`, `README.md`, `CLAUDE.md` §2, `dev/REFACTORING_PLAN.md` §7 execution log.

Rough effort: 7a 1 day, 7b 1–2 days, 7c 1 day, 7d 1–2 days, 7e 1 day of compute + writing.

---

## 4. Validation matrix

| level | reference | inputs | tolerance |
|---|---|---|---|
| tables, r = 1 | `CCalcWmML` (C++) | b ∈ {2,4,5,20}, (L,K) grid | 1e-12 rel |
| tables, any r | `generalb_gkm.kernel_coefficients` (exact) | 5 alphabet vectors × 3 kinds | 1e-9 rel |
| kernels, r = 2 | `twoblock_gkm.kernel_matrix` **and** `generalb_gkm.kernel_matrix` | `pos/neg.2fa`, synthetic with `N`, ± RC, `-t 0/1/2`, `-d` auto and bounded | 2e-6 rel, 1e-7 abs |
| kernels, r = 3, T = 4 | `generalb_gkm.kernel_matrix` | `three_tracks.mfa`, 4-track synthetic | same |
| classify | exact `Σ α_j K(x_i,x_j)/(‖x_i‖‖x_j‖)` | two-block SV fixture | 2e-6 rel |
| DNA regression | golden corpus (151 cases) | unchanged | byte-identical |
| memory/UB | ASAN + UBSAN | every fixture | clean |
| applications | gkmsvm3 experiments 01, 02, 03, 06, 07 | Python vs C++ backend | AUC equal to 1e-3 where both run; runtimes reported |

## 5. Risks and how they are handled

* **DNA regression through the DFS change** — mitigated by `step ≡ 1` for T = 1, golden corpus, and
  the benchmark gate (fallback: template on `Uniform`).
* **Class count blow-up** (`C = 2^ℓ` when all sizes differ) — reachable-row allocation (D4), tiling,
  and a hard error with the class count in the message when even one row set exceeds the budget;
  documented guidance: few distinct alphabet *sizes*, small k.
* **Pass table blow-up for long ℓ-mers** — D10 prefix-split passes; threshold measured, not guessed.
* **Profile counter overflow** — `aint` is 32-bit; a pair of sequences contributes at most
  (windows × strands)² per class; error out if `nlmers_max² ≥ 2^31`.
* **Double-precision tables** — every table is checked against exact rationals; `Πb_i` up to
  `4^20·2^20 ≈ 1e18` is still exact in `long double` accumulation, and the e_j sums are integers < 2^53
  for ℓ ≤ 60, k ≤ 8 (asserted at table construction).
* **Ambiguous alphabet spec** (a literal that is also a file name) — D2's `=` prefix on the C++ side
  is unambiguous; only the R convenience layer guesses, and it says so in the Rd page.
