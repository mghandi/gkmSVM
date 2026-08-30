# Prompt for the gkmsvm3 session: use the gkmSVM C++ kernel for the experiments

Paste the text below into the Claude Code session that works in `../gkmsvm3`.

---

**Context: a fast C++ implementation of the general-B kernel now exists in the gkmSVM repo — use
it for the experiments.**

The gkmSVM package at `/Users/mghandi/CDS/gkmSVM` (https://github.com/mghandi/gkmSVM, `master`,
version 0.90.0) implements the gapped k-mer kernel over heterogeneous alphabets ("Phase 7": design
in `dev/PHASE7_PLAN.md`, user guide in `tutorials/gkmsvm-multitrack-tutorial.md`). It computes
exactly what `generalb/generalb_gkm.py` computes: its test suite checks 33 096 multi-track kernel
entries and 276 classify scores against `generalb_gkm.py` and `twoblock_gkm.py` (vendored in
`dev/oracle/`) to 2e-6. It is 20× or more faster than the pure-Python path and multithreaded.

**Building and the binaries.** `make` produces
`build/gkmsvm_kernel`, `build/gkmsvm_train`, `build/gkmsvm_classify` (`make test` runs the golden
corpus, `make oracle` the exact-arithmetic checks). The R package (`R CMD INSTALL .`) exposes the
same through `gkmsvm_kernel(..., alphabets = c("dna","01"))`, `gkmsvm_train`, `gkmsvm_classify`,
`gkmsvm_trainCV`, `read_mfa`/`write_mfa`.

**Input format** = our multi-track FASTA (`.mfa` / `.2fa`): a header line, then one line per
track, all of equal length. Alphabets are given per track with `-A spec[,spec...]` where a spec
is `dna`, `rna`, `protein`, an alphabet file (one symbol per line), or `=SYMBOLS` for a literal
alphabet (`=01`, `=NUM`, `=abcdef`). Windows are track-major words of T·L positions; positions of
equal alphabet size form one block (same convention as `generalb_gkm.Blocks`). A symbol outside
its alphabet (e.g. `N`) skips the windows covering it. Reverse complement = track 1 complemented +
reversed, other tracks reversed (`-R` disables it; it is off automatically if track 1 has no
complements).

**Kernel:** `gkmsvm_kernel -l L -k k -d D -t T -A dna,=01 [-R] [-T threads] pos.mfa neg.mfa out.txt`

* `-t 0` = `kind='filter'` (H), `-t 1` = `'truncated'` (dominance rule), `-t 2` = `'gkm'`
  (C(ℓ−|m|,k)).
* `-d` = bound on the *total* mismatches (`-1` = all classes, i.e. exactly our Python kernels).
* Output: normalised lower-triangle text kernel (rows = pos records then neg records;
  `out.txt.index` sidecar with names/row ids), `-b` for the binary `.gkmk`. Options must come
  *before* the file names (BSD getopt).
* `k` may be up to T·L; with one track K > L is also accepted (as in `generalb`).
* Long words: above 750 000 mismatch patterns it switches automatically to "prefix-split" passes
  (`-P 2` forces them, `-P 1` the classical iDL table); results are identical.
* Limits: ≤ 32 symbols per alphabet; symbols ` >#;,=` cannot be used literally (remap them);
  memory for the class profile is Π(block length+1) × n² × 4 bytes, tiled automatically within
  1 GB (`-r rows`, `tileMemoryMB`).

**Training / scoring:** `gkmsvm_train -q -c C -A dna,=01 kernel.txt pos.mfa neg.mfa prefix`
(LIBSVM C-SVC on the precomputed kernel; writes `prefix.gkmmodel` with `#rho` and `#alphabets`,
plus the legacy `prefix_svalpha.out`/`prefix_svseq.fa`);
`gkmsvm_classify -l L -k k -d D -t T test.mfa prefix.gkmmodel prefix.gkmmodel scores.txt`
(score = Σ α_j⟨x,x_j⟩/(|x||x_j|) − rho; the model carries its alphabets).

**The harness already has a backend for it — on the gkmsvm3 branch `cpp-backend` (not merged
into `main` yet; merge or rebase it first):** `experiments/common/cpp_backend.py`, and
`methods.py`/`classifiers.py` switch on environment variables:

* `GKMSVM_BIN=/Users/mghandi/CDS/gkmSVM/build` → every `kernel_for_method` kernel comes from
  `gkmsvm_kernel` (unsafe symbols are remapped automatically; alphabets of > 32 symbols fall back
  to Python with a note).
* `GKMSVM_SVM=libsvm` → `kernel_svm()` trains through `gkmsvm_train` (exact C-SVC) instead of the
  simplified SMO, which is what dominates the runtime above a few hundred sequences.
* `EXP_RESULTS_DIR=results_cpp` → writes next to the Python `results/`;
  `common/compare_results.py A.json B.json` compares two runs.

Already rerun this way (numbers in `experiments/README.md`, section "C++ backend"): 01, 02, 03,
06, 07 reproduce the Python results (largest |ΔAUC| 0.03, on a flat single-track kernel) at 4–23×
the speed, and 06 at 1 000 + 1 000 sites took 15 min (`N_PER=1000 python3 prepare.py;
GKMSVM_SVM=libsvm python3 run.py --n 1000`, results in `06_ctcf_methylation/results_cpp_n1000/`).

**What to do with it:** run the experiments that were subsampled for lack of speed at full scale —
06 with all matched sites and more cell types, 05 and 07 without subsampling, the k/L probes on
more than 250 samples — always with `GKMSVM_BIN` and `GKMSVM_SVM=libsvm` set, through `Runner` so
runtimes are recorded, and keep the pure-Python path as the reference (spot-check a small
configuration in both modes with `compare_results.py` before trusting a new setup). Record in each
README's "Decisions" section that the C++ backend was used.
