# metal_g3 — Route 2 ("G3") prototype of the gkm / general-B kernel on Apple Metal

Prototype of section 6 of `dev/PERFORMANCE_PLAN.md` and section 4 of
`gkmsvm3/theory/gkm_kernels_constantB.pdf`: the mismatch profile is computed by classifying every
pair of **distinct** l-mers on the GPU and scattering `cnt_u ⊗ cnt_v` into per-class n×n integer
counters; the kernel `K = Σ_m c(m) N_m` is then formed on the CPU in double precision, normalised
and written as `.gkmk` for direct comparison with `gkmsvm_kernel`.

```
make
./gkm_metal -l L -k K -d D [-t 0|2] -A dna,=01 [-R] [-o out.gkmk] [-M GB] [-g groups] [-L light] [-C chunk] in.mfa
```
Same conventions as `gkmsvm_kernel` (multi-track `.mfa`, `-A` specs `dna|rna|protein|=SYMBOLS`,
reverse complement on by default: track 1 complemented + reversed, other tracks reversed; windows
containing a symbol outside the alphabet are skipped; `-d -1` = all classes). Kinds: `-t 0` filter,
`-t 2` gkm (the truncated filter is not polynomial and is not implemented here).

## Design

* **Encoding.** Each distinct l-mer is packed into a 64-bit code, `f = ⌈log2 max b⌉` bits per
  position (prototype limit `ℓ·f ≤ 64`, ≤ 8 blocks). Per-block mismatch counts of a pair:
  `x = a ^ b`, OR the `f` bit-planes of every field into the field's low bit, mask per block,
  `popcount` — about ten integer operations per pair, no tensor cores.
* **Classification.** 2-D thread grid: thread `(r, c)` handles row `u = rowLo + r` and the 512
  columns `v ∈ [512c, min(512(c+1), u))` of the strict lower triangle (plus the self pair when
  `c = 0`). Rows are processed in chunks (`-C`, default 4096 rows) so that the deferred-pair list
  stays bounded. The class index (mixed radix over blocks) is mapped through a table to a compact
  index or "unreachable" (`|m| > d`).
* **Scatter.** Pairs whose posting-list product is ≤ `-L` (default 16) are scattered inline with
  32-bit atomics into `N[class][tri(i,j)]` (lower triangle, `i ≥ j`; the two ordered window-pair
  terms of a symmetric cell both come out of the same double loop, the diagonal cell gets factor 2).
  Larger products are appended to an overflow list and processed after each chunk by a second
  kernel with one thread per (entry, posting of the longer list).
* **Kernel.** `K = Σ_c coef[c] N_c` over the reachable classes in double; `H(m)` from the
  elementary symmetric polynomials, gkm from binomials; normalised; `.gkmk` float32 with CRC-32.
* **Class groups (2026-08-31).** The counter buffer holds `Rg ≤ R` classes at a time (`-M GB`, default
  60 % of `maxBufferLength`; or `-g` groups): the pair space is classified once per group with a
  class table that maps the other classes to "skip", each group's counts are folded into `K` on the
  CPU (classes in increasing compact order, so the double sum is the same as with one group), and
  the buffer is reused. Classification is cheap, so this costs little: exact 06 kernel at 13 398
  sites 36.8 s with one 43 GB group, 37.3 s with 6 groups of 8 GB, 51 s with 25 groups of 2 GB —
  all bit-identical. Memory is no longer the limit on n.

## Results (Apple M3 Ultra, 60-core GPU; CPU = `gkmsvm_kernel -P 2`, 28 threads; all bit-exact)

| input | ℓ, blocks | d | Metal total | CPU | speed-up |
|---|---|---|---|---|---|
| exp. 06, n = 4 000 (U = 140 k) | 20, (4¹⁰\|2¹⁰) | 4 | 0.73 s | 2.5 s | 3.5× |
| same | | 6 | 0.98 s | 15.6 s | 16× |
| same | | **−1 (exact)** | **2.9 s** | 427 s | **150×** |
| exp. 07, n = 4 000 (U = 114 k) | 24, (4⁸\|2⁸\|3⁸) | 6 | 3.0 s | 57.7 s | 19× |
| exp. 06, n = 13 398 (U = 329 k) | 20 | 4 | 4.0 s | 33.6 s | 8× |
| same | | 6 | 10.1 s | 395 s | 39× |
| same | | **−1 (exact)** | **36.8 s** (43 GB profile) | not finished on CPU (est. ≥ 6 h) | — |

Additional checks on 2026-08-31 (every one 0 differing entries vs `gkmsvm_kernel`, whose kernel-mode
path was itself checked against the profile path): class groups `-g 4` / `-g 121` on the exact 4 000-site
06 kernel (3.05 s / 18.9 s vs 3.18 s), `-g 3` on 07 d = 6, `-M 8` / `-M 2` on the exact 13 398-site
kernel (37.3 s / 51.3 s). Through the gkmsvm3 harness (`GKMSVM_GPU=1`) the prototype reproduced every
CPU result of experiments 05–07 (largest |ΔAUC| 1e-4, LIBSVM noise) and ran the configurations that
had been cut on the CPU: 06 at all matched sites for six cell types (60 analyses, 44 min), 06
transfers (81, 7 min), 05 all six CLIP sets at 12 000 sequences (106, 68 min), 07 all stages (18,
3 min).

Classification runs at 1–3·10¹⁰ pairs/s when few pairs are near; the run time is dominated by the
scatter of near pairs (inline atomics and the heavy-pair kernel), which is why cost grows with `d`
and is highest for the exact case, where every pair contributes. The comparison files were produced
with `gkmsvm_kernel ... -P 2 -b` and compared with `cpp_backend.read_gkmk` (max |Δ| = 0.0 in every
run: the profile is integer-exact on both sides and the same double arithmetic forms the kernel).

## Limitations / next steps

* `ℓ·f ≤ 64` bits (extend to 128-bit codes for ℓ = 32 with 5-symbol tracks), ≤ 8 blocks,
  `-t 1` (truncated filter) not supported, single GPU.
* The profile buffer, `reachable classes × n(n+1)/2 × 4 B`, is bounded by `-M` through class
  groups (above). A kernel-mode variant that accumulates `c(m)·cnt` directly into one triangle is
  not available on Metal: the shading language has no 64-bit `fetch_add` and no double, and 32-bit
  digit-split accumulators would need 3–4 atomics per pair where the scatter already dominates.
* Threadgroup-shared code tiles were tried on 2026-08-31 and removed: the classify stage is not
  memory-bound (0.35 vs 0.34 s at 4 000 sites, 2.5 vs 2.3 s exact), the scatter is what costs, and
  the first version was wrong for partial threadgroups — a reminder that every variant must be
  compared entry for entry before use. Simdgroup compaction of the deferred pairs would help the
  classify stage only, for the same reason. The exact-case scatter (all pairs) can be reformulated
  as Route 1 (pattern-Gram GEMM) when k ≤ 5.
* Multi-index prefilter for U ≫ 10⁶; a CUDA port is a straight translation (the kernel is plain SIMT).
