# Making the general-B kernel fast at n ≥ 10 000 — measured diagnosis and plan

*2026-08-30. Follows the full-scale experiment reruns in gkmsvm3 (`experiments/README.md`,
"C++ full-scale reruns"): at 4 000–12 000 sequences the joint kernels took 20 min – 2 h each and the
runs had to be cut. Everything below is grounded in a profile of a representative kernel
(experiment 06 joint, 4 000 records × 40 bp, `-l 10 -k 6 -d -1 -t 0 -A dna,=01`) taken with
`sample` on the M3 Ultra, and in the source (`kernel_tree_impl.h`, `LTreeS_impl.h`).*

## 1. Where the time actually goes (measured)

1. **One pass runs alone for most of the wall clock.** The prefix-split design makes 64 passes;
   28 threads start, but the pass costs are wildly unequal, and after ~1 min the profile shows
   exactly **one worker thread** in `DFSTiDL` while the other 27 cores idle (`sample`: main thread
   in `join`, one thread at 100 %). The "single-threaded n² stage" observed in the overnight run is
   this tail pass, not a serial phase after the passes.
2. **~80 % of the tail is internal-node traversal, not leaf counting.** Top-of-stack: `DFSTiDL`
   (list building at internal nodes) 4 828 samples vs `DFSTnIDL` (leaf-pair profile increments)
   1 240. With `-d -1` the mismatch-pattern tree never prunes — every compatible subtree is carried
   in the DFS list at every depth — so the traversal itself is the quadratic cost.
3. **Tiling multiplies the traversal.** The per-class profile needs `nclasses × 4 B × n` per row
   (121 classes for the 06 joint kernel), tiled into 1 GB: at n = 7 824 that is ~28 tiles, and
   **each tile re-runs all passes**; the leaf work partitions cleanly across tiles (seqID pruning)
   but the dominant internal-node traversal is largely repeated per tile. This is why kernel time
   grows visibly faster than n² past ~4 000 sequences.
4. `MULTI_THREAD_SAFE` makes every profile counter a `std::atomic<int>` (~2× on the leaf part —
   secondary, given 1–3).
5. Hardware headroom: Apple M3 Ultra — 28 CPU cores (20 P + 8 E; the kernel caps itself at
   `-T = 2·l = 20`), **60-core GPU (Metal 4, ~28 TFLOP/s fp32)**, 256 GB unified memory, AMX units
   reachable through Accelerate (`cblas_ssyrk`, ~1–2 TFLOP/s effective). The current code uses at
   most one CPU core during the dominant phase.

## 2. Plan A — fix the tree path (days of work, est. 50–300× at n ≥ 8 000)

**A1. Kernel-mode accumulation (biggest single win).** The harness only ever needs the kernel, not
the per-class profile. At the leaf, the class index is known, so accumulate
`K[i][j] += c[class]` into one fp64 (or fp32 + fp64-reduce) triangle instead of
`++profile[i][class][j]` ints. Effects:
  * memory per row drops by the factor `nclasses` (121× for 06) → **no tiles at all** up to
    n ≈ 50 000 → the ~28× repeated traversal disappears;
  * the atomics question disappears if each worker owns a private triangle (n = 16 000 →
    0.5 GB fp32 per thread, 10 GB for 20 threads — fine in 256 GB) reduced once at the end.
  Keep the profile path for `-d ≥ 0` diagnostics and the truncated filter (whose coefficients are
  not linear in the class counts until the dominance mask is applied — it stays on the old path).
**A2. Split the heavy passes.** Recursive prefix-split: estimate each pass's cost (pattern-tree
  mass × trie node counts at the split depth), and split any pass above ~1/(4·threads) of the total
  by extending its prefix. Removes the one-thread tail; with A1's private accumulators there is no
  contention. Target: ≥ 15× parallel efficiency on 20 P-cores.
**A3. Raise the `-T` default** to `hardware_concurrency()` once A2 exists.

Estimate for the measured cases: H1 all-sites joint kernel (7 824 seqs) 7 339 s today →
`/~28` (no tile repetition) `/~15–20` (parallel passes) ≈ **15–30 s**. Even if the tile factor is
half that, it lands in minutes. A1 and A2 are independent, verifiable against the current path
bit-for-bit (A1 in fp64 vs profile×c) and against `tests/`/`dev/oracle`.

## 3. Plan B — exact pattern-Gram (GEMM) reformulation → AMX / Metal GPU (1–2 weeks, exact, scales to n ≥ 30 000)

**Verified identity** (exact-arithmetic check in gkmsvm3, 2026-08-30, for the 05/06/07 block
configurations): the filter coefficient `H(m)` is a **polynomial of total degree ≤ k in the block
mismatch counts** m = (m₁..m_r), i.e. `H(m) = Σ_{|t|≤k} β_t Π_b C(m_b, t_b)` exactly; the gkm
coefficient `C(ℓ−|m|, k)` likewise (degree k in |m|). Consequently, for `-t 0` (filter) and `-t 2`
(gkm) with `-d -1`:

    K = Σ_{patterns p, |p| ≤ k} γ_p · M_p,        M_p[i][j] = Σ_w X_p[w,i] · X_p[w,j]

where p ranges over position subsets of size ≤ k, X_p[w,i] = number of windows of sequence i whose
projection onto p equals the |p|-mer w (key space Π_{pos∈p} b_pos ≤ 4^k), and γ_p follows from β by
inclusion–exclusion (γ depends only on p's block weights). **No pair enumeration, no mismatch
bound, exact.** Each M_p is a rank-`b^|p|` symmetric update: `cblas_ssyrk` on a dense
(b^|p| × n) count matrix, batched over patterns; K accumulates in fp64.

Cost = Σ_p (rows_p · n²/2) MACs. For the experiment configurations:
  * 06 at the probe optimum k = 3–4 (ℓ = 20): ~6 200 patterns, ≈ 7·10⁵ · n²/2 MACs → n = 13 400:
    ~65 TFLOP → **seconds on the 60-core GPU, ~1 min on AMX**. Today: hours.
  * 05 (ℓ = 12, k = 4) at n = 30 000 (the full CLIP sets, currently untouchable): ~800 patterns →
    minutes on AMX/GPU, and the n² memory is one fp32 triangle (1.8 GB).
  * The k = 6, ℓ = 20 point is the expensive corner (~4·10⁴ patterns, ~1.7 PFLOP at n = 13 400 →
    ~1 h on the GPU): still bounded, but Plan A's tree remains competitive there. Rule of thumb:
    pattern-Gram wins for k ≤ 4–5 (exactly where the (L,k) probes put every joint model).

Implementation route: CPU first (`Accelerate` syrk, OpenMP-style batching over patterns), then the
same batches as Metal compute (MPS matrix multiplication or a simd-group-matmul kernel; unified
memory means no copies). Side benefit: X_p **is** an explicit sparse feature map — it gives an
LS-GKM-style linear-SVM training path and cheap deltaSVM-style scoring at any n.

## 4. Plan C — the `-d` study: **done (2026-08-30), and it changes the priorities**

Run first at the owner's request (gkmsvm3 `experiments/*/dstudy.py`, `results_dstudy/`; harness knob
`GKMSVM_D`). Facts established:

* A capped kernel **equals** the exact `-d -1` kernel with the coefficients of the classes |m| > d
  zeroed — verified against the exact Python path with zeroed coefficients (max diff 1.8e-8), so the
  already-measured exact AUCs served as baselines with no recomputation.
* The exact AUC is reproduced at **d = 3–4** for ℓ = 12 joint models (experiment 05: 30–130× kernel
  speedup; the 12 000-sequence joint ctx5 kernel: 65 min → 31 s) and ℓ = 20 with `meth2s`
  (experiment 06: 0.7498 vs 0.7484, 5 s vs 1 328 s ≈ 270×), and at **d = 6** for `meth3` and for the
  ℓ = 24 3-block model (experiment 07: 0.8044 vs 0.8042, 235 s vs 7 620 s ≈ 32×). DNA-only and
  product-alphabet kernels are d-insensitive. The needed d grows with the joint word length.
* Caveats: capping changes the normalisation, so the normalised kernel does **not** converge
  monotonically in d — always confirm d against d + 2; and one SVM instability was seen on a single
  capped kernel (05 product alphabet at exactly d = 3), so check AUC too, not just kernel diffs.

**Revised priorities.** With `GKMSVM_D` set per configuration every kernel in the current
experiments runs in seconds-to-minutes, so the capped tree path is the workhorse for the
applications paper. Plans A and B remain worth doing, in this order, for (i) the exact `-d -1`
reference at scale, (ii) `meth3` / long-word configurations where the converged d (≈ 6) still costs
minutes-to-hours at n ≥ 10⁴, and (iii) n ≥ 30 000 (full CLIP sets), where A1's untiled kernel-mode
and B's GEMM shape are the difference between feasible and not.

## 5. The capped-d regime (d = 3–6) — measured profile and re-weighted priorities (2026-08-30)

Since capped d reproduces the exact AUC at these sizes, **d = 3–6 is the regime users will actually
run**, and the optimisation targets there differ from the `-d -1` picture of section 1. Measured on
the M3 Ultra (idle machine), representative configurations:

| config | `-T 1` | default (28 thr) | speedup | CPU used | **`-P 2` (prefix-split)** |
|---|---|---|---|---|---|
| 06 joint meth2s d = 4, n = 13 398 (ℓ = 20) | 458 s | 95 s | 4.8× | 500 s | **34 s** (18 eff. cores) |
| 07 joint 3-block d = 6, n = 4 000 (ℓ = 24) | 800 s | 240 s | 3.3× | 1 571 s (!) | **61 s** (17.8 eff. cores) |

* Mid-run profile of the 06 case (greedy design): leaf counting 31 %, **`cloneReorder` 31 %**,
  traversal 27 %, idle 5 %. The greedy iDL design clones + reorders the whole trie *per pass per
  tile* — 11 tiles × 40 passes = 440 clones of a 329 394-l-mer trie. Invisible at `-d -1`
  (traversal dominates there), a third of the work at capped d.
* The 07 case burns 2× the CPU multithreaded vs `-T 1` (1 571 s vs 798 s): `MULTI_THREAD_SAFE`
  atomic profile counters under contention.
* **Prefix-split (`-P 2`) at capped d: identical kernel (verified to 0.0 on both cases), 2.8–3.9×
  faster wall.** No per-pass clone (it traverses the shared trie), 64 passes balance far better
  (18 effective cores vs ~5). Its total CPU is ~20–35 % above greedy, so greedy remains the right
  choice for 1–2 threads; the `-P 0` auto rule should become threads-aware (use prefix-split
  whenever the thread budget is ≳ 4, regardless of the pattern-table size). The gkmsvm3 harness now
  passes `-P 2` whenever `GKMSVM_D ≥ 0`.

Thread-efficiency sweep over d (same 06 configuration, n = 13 398, ℓ = 20, 28 threads; "eff." =
CPU-seconds / wall-seconds):

| d | greedy (`-P 1`) wall / CPU (eff.) | prefix-split (`-P 2`) wall / CPU (eff.) |
|---|---|---|
| 2 | 8.6 s / 66 s (7.7) | 7.3 s / **26 s** (3.6) |
| 3 | 14.4 s / 134 s (9.4) | 11.6 s / 140 s (12.0) |
| 4 | 95 s / 500 s (5.3) | **34 s** / 618 s (18.2) |
| −1 | 1 thread for most of the wall (section 1) | — |

The imbalance exists at every d, mildest at d = 3–4 with `-P 2`; at d = 2 the greedy design's CPU
is dominated by its fixed per-pass cost (the trie clone, independent of d — 440 clones here), which
is why prefix-split wins on *total CPU* there despite worse balance. Tiling multiplies the tail:
each of the 11 tiles re-runs the pass set behind a join barrier, so the straggler penalty is paid
per tile.

History note: the imbalance is *not* a refactoring regression. v0.80 (2018) assigned passes
statically (`thread j0 takes passes j ≡ j0 mod nThreads`, no stealing) and even rounded nThreads
toward a divisor of M (`nThreads = M/(int)(M/nThreads)`); Phase 6 (2026) replaced this with a
dynamic atomic work queue — bit-identical output, strictly better schedule. What is new since
v0.80 are the regimes that expose the remaining tail: multi-track words of ℓ = 20–32, `-d -1`,
n ≥ 10⁴, and the prefix-split design itself. Divisibility of the thread count into M is irrelevant
under dynamic pulling; the floor is the heaviest single pass, which only splitting it can lower.

**Priorities for the capped-d workhorse path, in order:**

1. *(shipped, zero code)* `-P 2` when d is capped and threads are available: 2.8–3.9×. Flip the
   `-P 0` heuristic in `kernel_tree_impl.h` accordingly.
2. **A1 kernel-mode accumulation** (section 2) pays off just as much here: even with `-P 2` the
   6–11 tiles repeat the traversal per tile (the 618 s CPU of the 06 `-P 2` run still contains
   ~11× duplicated traversal), and private per-thread accumulators remove the 2× atomic-contention
   tax. Estimate: 06 case 34 s → ~8–12 s; grows with n (54 tiles at n = 30 000).
3. **A2a — makespan-aware greedy assignment (cheapest structural fix).** `newGreedyIDLPasses`
   assigns every mismatch pattern to the pass where *its own* estimated DFS cost (`calcCost`, node
   pairs surviving per depth) is minimal — the objective is total work, with no term for the
   *maximum* per-pass load, which is exactly why one pass ends up dominating. Keeping the same cost
   model, assign each pattern to the least-loaded pass among those within (1 + ε) of its cheapest
   (ε ≈ 0.2, patterns visited in decreasing min-cost order, LPT-style): ~10 lines in
   `CiDLPasses.cpp`, no change to correctness (any assignment is exact), directly lowers the
   straggler floor. The same idea applies to the prefix-split design by splitting heavy prefixes one
   level deeper (A2 proper).
4. **A2b — cost-model fix for general B (measured 2026-08-30, gkmsvm3 experiment 06 data,
   830 676 joint 20-mers).** The pass-design cost `C_pi(delta)` uses a single mean match
   probability (`p = mean(1/b_i)`) and depth weights `w_i` from the identity-order trie. Both are
   fine for constant B: pooled sliding l-mers make the reordered-trie node counts
   permutation-invariant to within 3 % (measured over rotations and scrambled orders), and the
   empirical per-position match probabilities are 1/3.84 vs the 1/4 assumed. For **general B** both
   break: the methylation track's *effective* alphabet is 1.02 (not 2), and true node counts differ
   between orders by up to **1 700×** at mid depths (meth-first order: 182 nodes at depth 10 vs
   312 499 for DNA-first) — the model is blind to which positions are cheap. Re-running the greedy
   assignment (6 196 patterns, the real 40-order dictionary) with per-position empirical match
   probabilities and per-order w lowers the predicted total DFS cost 6.8× and the predicted
   heaviest-pass load 12× (34 % → 19 % of total; evaluated under the improved cost model, so treat
   the exact factors as indicative). Implementation cost: per-position p_j = one pass over the
   l-mers (free); per-order w by *subsampled exact counting* (hash ~50 k reordered prefixes per
   order, ~ms each in C++) — the closed-form independence estimate `min(N, prod 1/p_j)` is up to
   148× off on the correlated methylation track and should only be a fallback. Pairs naturally
   with A2a: better per-order costs concentrate patterns onto genuinely cheap passes, which makes
   the load-aware tie-breaking more, not less, important.
5. **A2 pass splitting** then lifts the remaining tail (prefix-split already gets 18/28 cores; the
   last passes still straggle).
6. Plan B (pattern-Gram GEMM) does **not** produce capped kernels (the capped coefficient is no
   longer a degree-k polynomial) — it stays the route to *exact* kernels at k ≤ 4–5, which, where it
   applies, supersedes capping altogether.

## 6. Plan G — a GPU gkmSVM: reducing the mismatch profile to matrix products

The question is whether the core — the mismatch profile N_m[i][j] (number of window pairs of
sequences i, j whose per-block mismatch counts are m) — can be written as matrix multiplication.
Three exact reductions exist; together they cover every kernel we run. Notation: W windows in
total, U distinct l-mers (06 all-sites: W = 830 676, U = 329 394; 07 at 4 000 windows: U = 113 709),
n sequences, cnt_u ∈ ℕⁿ the per-sequence multiplicity of l-mer u.

**G1 — exact filter/gkm kernel (any d, k ≤ 5–6) = batched syrk over pattern-count matrices.**
Because H(m) and C(ℓ−|m|,k) are polynomials of degree ≤ k in the block mismatch counts (verified
exactly, section 3), K = Σ_{|p|≤k} γ_p · X_pᵀX_p with X_p the (b^|p| × n) window-projection count
matrix of position subset p. Dense, tensor-core-shaped; streams pattern batches; fp16/bf16 inputs
(counts ≤ 32 767) with fp32/fp64 accumulation. Cost Σ_p b^|p| · n²/2 MACs: 06 at k = 4 → 65 TFLOP
(seconds on H100/A100, ~1 min on the M3 Ultra GPU); 06 at k = 6 → 1.7 PFLOP (≈ 2–5 s on H100 dense
bf16, ~1 h on M3 Ultra) vs ≈ 2 h CPU exact today. GaKCo (Singh et al. 2017) is the CPU counting
precedent for the gkm case. Cannot produce *capped* kernels (a capped coefficient is not a
polynomial) — but where it applies it makes capping unnecessary.

**G2 — capped-d profile = sparse Gram of "deletion-index" matrices + binomial inversion.**
For patterns q with t ≤ d wildcard positions, M_q = X_qᵀX_q counts window pairs agreeing on all
non-wildcard positions; a pair at block distance s is counted in Π_b C(ℓ_b−s_b, t_b−s_b) such
patterns and a pair at distance > d in none, so A_t := Σ_{|q̄|=t} M_q = Σ_{s≤t} [Π_b C(ℓ_b−s_b,t_b−s_b)] N_s
(block-resolved), triangular and integer-invertible: the whole capped profile from d + 1 sparse
Grams. GPU shape: pack the projected windows into ≤ 128-bit keys, radix-sort / hash-join per
pattern, segmented reduction, scatter of postings outer products. Work = patterns × W key
insertions (06 d = 4: 6 196 × 830 676 = 5·10⁹ keys — seconds on a GPU sort) **plus the covering
multiplicity on near pairs**, which is the catch: identical/near-identical windows (the CTCF motif
core in every 06 window) are re-counted in thousands of patterns. Best when near-duplicate windows
are rare (typical genome-wide negative sets), poor on motif-centred or repetitive sets.

**G3 — capped-d profile = one Gram over distinct l-mers with a fused classify epilogue (the
flagship).** Encode each distinct l-mer one-hot per block: O_b ∈ {0,1}^{(ℓ_b·b_b) × U}. Then
(O_bᵀO_b)[u][v] = number of matching positions of u and v in block b, so the per-block mismatch
counts of *every* l-mer pair are the entries of r small Gram matrices — literally matrix
multiplication (int8/one-hot: 06 has Σ ℓ_b b_b = 60 rows). Tiled U × U (never materialised), the
epilogue computes the class m per pair, drops pairs with |m| > d, and emits (u, v, c(m)) — or, in
the bit-parallel variant, XOR + per-field popcount on packed codes (≈ 10 integer ops per pair, no
tensor cores needed, the natural Metal version). The surviving near pairs are then scattered as
c(m)·cnt_u ⊗ cnt_v into the n × n kernel (kernel-mode) or into per-class n × n counters (profile
mode); most cnt vectors have one entry, heavy l-mers (motif cores, poly-runs) are batched as small
dense outer products. Each pair is classified exactly once — no covering multiplicity, no
d-dependence in the dominant stage, no data-dependent pathologies, trivially multi-GPU.
Cost U²/2 pair classifications: 06 all-sites 5.4·10¹⁰ → ≈ 5–20 s on the M3 Ultra GPU (custom
Metal), ≈ 1–3 s on an H100 (int8 GEMM ~3·10¹² MACs); 07 (U = 113 709) → < 1 s H100, ~3 s M3 Ultra vs
61 s CPU at d = 6; scatter is output-sensitive (near pairs only). Above U ≈ 10⁶ (full CLIP sets)
add a multi-index-hashing prefilter (split positions into d + 1 parts; a pair with ≤ d mismatches
is exact on some part) to cut the quadratic — the standard GPU all-pairs-Hamming design.

What G3 does *not* do: the full, uncapped profile with all classes — its scatter is Σ over all
l-mer pairs of postings products = W², inherent to the object (that is why the exact kernel goes
through G1's polynomial identity instead). Nobody needs that object.

**Coverage.** Exact kernel, k ≤ 6 → G1. Capped d, any k, any B → G3 (G2 as the sort-based
alternative on near-duplicate-free data). Both give bit-exact integers/verified coefficients and
are checked against `-a 0` and the oracle corpus like every other path. Expected end state on a
CUDA card: exact 06 kernel (today ≈ 2 h CPU) in seconds; capped kernels (today 30–60 s CPU, often
2–4 min under load) in ≈ 1 s; the classify/score path (deltaSVM-style) reuses G3's pair
classification against the support-vector l-mers.

**Roadmap (Plan G).**
1. Prototype G3 kernel-mode on Metal (M3 Ultra) and G1 on Accelerate/AMX — one week each — with the
   oracle corpus as the acceptance test; keep the tile/epilogue design portable (Metal ↔ CUDA are
   both plain SIMT kernels here).
2. CUDA backend (`gkmsvm_kernel -a 5`, R `useGPU=TRUE`) with cuBLASLt int8/bf16 for G1 and a
   custom popcount/one-hot tile kernel for G3; multi-index prefilter for U > 10⁶.
3. Only then revisit the CPU tree path: with G1/G3 available, the CPU items A1/A2/A2a/A2b become
   the fallback path for machines without a GPU, still worth doing but no longer on the critical
   path.

## 7. What not to do

* Do **not** port the DFS/trie to the GPU — irregular pointer chasing, no coalescing; the GPU win
  comes through the Gram formulations of Plan G (G1 for exact kernels, G3 for capped d).
* Micro-optimising `addmmprof`/branchless leaf code buys ≤ 20 % (leaf work is the minor term).
* More threads on the current pass design do nothing — the tail pass is one unit of work.

## 8. Suggested order

1. Done: C's d-study (section 4) and `-P 2` at capped d (section 5).
2. Threads-aware `-P 0` auto rule (one heuristic in `kernel_tree_impl.h`).
3. A1 kernel-mode (no tiles, private accumulators) — the biggest remaining lever in *both* regimes;
   verify vs oracle, measure at capped d and `-d -1`.
4. A2 pass splitting — measure parallel efficiency at n = 8 000–30 000.
5. B on CPU/AMX for k ≤ 4, then the Metal batch; wire into `gkmsvm_kernel -a 3` and the harness
   (exact kernels without capping where it applies).
