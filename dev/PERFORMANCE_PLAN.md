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

## 4. Plan C — an empirical `-d` study (half a day, possibly makes the corner cases moot)

The 2014/2016 tools defaulted to `-d 3`; only the gkmsvm3 harness insists on `-d -1` to match the
exact Python reference. With `-d` bounded the pattern tree prunes hard and the existing tree path
is fast. One systematic comparison (kernel Δ and AUC Δ at d = 4, 6, 8 vs −1 on 06 and 05 data)
would tell whether the applications paper can simply run d-capped kernels; if AUC is unchanged at
d ≈ 6, most large-n pain disappears with zero new code. Report it either way — it is an
approximation knob, not a bug fix, and the exact paths (A/B) are still worth having.

## 5. What not to do

* Do **not** port the DFS/trie to the GPU — irregular pointer chasing, no coalescing; the GPU win
  comes only through the Plan-B GEMM shape.
* Micro-optimising `addmmprof`/branchless leaf code buys ≤ 20 % (leaf work is the minor term).
* More threads on the current pass design do nothing — the tail pass is one unit of work.

## 6. Suggested order

1. A1 kernel-mode (no tiles, private accumulators) — verify vs oracle, measure.
2. A2 pass splitting — measure parallel efficiency at n = 8 000–16 000.
3. C's d-study on real data (uses the now-fast tree path).
4. B on CPU/AMX for k ≤ 4, then the Metal batch; wire into `gkmsvm_kernel -a 3` and the harness.
