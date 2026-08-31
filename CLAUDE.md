# CLAUDE.md — handover for working in the gkmSVM repository

Read this first, then `dev/REFACTORING_PLAN.md` (the plan is the source of truth; this file is the
operating manual for executing it).

## 1. What this repo is

`gkmSVM`: the R package wrapping the gkm-SVM (gapped k-mer SVM) C++ implementation of Ghandi et al.
(*PLoS Comput Biol* 2014 / *Bioinformatics* 2016). **Public repository that people download and use**
(https://github.com/mghandi/gkmSVM, 17 stars), and the source of the CRAN package. Owner: Mahmoud
Ghandi (mghandi). CRAN maintainer since 0.83.0: Mike Beer <mbeer@jhu.edu>.

Related repo on this machine: `../gkmsvm3` — the generalisation of gkm-SVM to per-position alphabets
(different alphabet sizes per position). Its `generalb/generalb_gkm.py` is an **exact-arithmetic
pure-Python implementation** (`fractions.Fraction`, unit-tested) that serves as the independent
numerical oracle for this refactor, and its theory is the target of Phase 7.

## 2. Current state (2026-08-29)

* `master` = **0.80**, plus the refactor scaffolding merged from PR #5 (`ebd50e9`):
  `dev/REFACTORING_PLAN.md`, `dev/baseline.sh`, `CLAUDE.md`, `.claude/settings.json`,
  `.Rbuildignore`, `.gitignore`. `src/` and `R/` are untouched at `222cc50`.
* **Branch stack (none merged): `phase0/baseline-safety-net` (PR #7) ← `phase1/latent-bug-fixes`
  (PR #8) ← `phase2/core-extraction` (PR #9, Phase 2a) ← `phase3/sequence-identity` (PR #10) ← `phase4/binary-formats` (PR #11) ← `phase4b/libsvm-train` (PR #12) ← `phase5/alphabet-generalisation` (PR #13) ← `phase6/performance` (PR #14) ← `phase0b/cran-reconciliation` (PR #15) ← `phase2b/kernel-context` (PR #16) ← `phase6b/tiled-profile` (PR #17) ← `phase2c/alphabet-per-call` (PR #18) ← `phase2d/reporter` (PR #19) ← `phase8/legacy-norm`.** CI (ubuntu gcc/clang,
  macOS, ASAN+UBSAN, R, benchmark gate) is green on every branch from #8 up after the Linux
  portability fixes were forward-merged through the stack; only #7 is red on Linux, by design
  (see the plan §7, "Linux CI reconciliation"). If you change something on a lower branch, merge it
  forward through the chain (`git merge -X ours <parent>` from each child) — an independently applied
  duplicate makes the child PR "CONFLICTING" and GitHub then runs no CI for it at all. See §7 of the plan for what each contains and what
  is left (Phase 2b = the context object), and `NEWS.md` for user-visible changes. Merge in order and
  retarget each PR to `master` after its base merges. Key commands: `make`, `make test` (golden), `make oracle`,
  `make bench`, `dev/scratch_install.sh` then `Rscript -e 'testthat::test_dir("tests/testthat",
  package="gkmSVM", load_package="installed")'`. Later phases branch from `master` but need Phase 0's
  files — until it is merged, branch from the Phase 0 branch and say so in the PR.
* The stale branch `refactor/plan` can be ignored (its content is on `master`).
* **The plan is approved.** All seven decisions are recorded in §5 of the plan. Execute the phases
  in order starting with Phase 0.
* CRAN ships **0.83.0 (2023)**, a diverged lineage: it has packaging hardening (`snprintf`,
  `.registration=TRUE`, Bioc deps in `Suggests`) but **not** this tree's multithreading, duplicate-ID
  checks or `normalizePath`. Folding the good CRAN changes in is **Phase 0b: low priority,
  non-blocking**, best done after Phase 1. A copy of the CRAN tarball is not committed; fetch from
  `https://cran.r-project.org/src/contrib/gkmSVM_0.83.0.tar.gz` when you get to it.

## 3. Where to start: Phase 0

Phase 0 is fully specified in the plan. Its first item unblocks everything else:

1. **`main()` shims + Makefile** so `gkmsvm_kernel`, `gkmsvm_classify` (and later `gkmsvm_train`)
   build from this tree. There is **no `main()` and no Makefile** today — the repo only builds the R
   shared object, and the published CLI binaries come from a separate `gkmsvm-2.0.tar.gz`.
   `dev/baseline.sh` already contains a working shim; lift it into `src/cli/`.
2. Golden-test corpus + `testthat`, 3. numerical oracle vs `../gkmsvm3/generalb/generalb_gkm.py`,
4. CI (R CMD check, golden tests, ASAN/UBSAN, benchmark gate), 5. baseline numbers (already in
`dev/baseline.sh`).

## 4. Verified facts you can rely on (do not re-derive)

Measured on this machine; `dev/baseline.sh` reproduces all of them.

* `CCalcWmML::calcKernel` (`src/CalcWmML.cpp:150`) reads `wm[i]` past the end of `new double[K+1]`
  on **every run**. Fix = bound the loop with `i<=m && i<=K`; **verified byte-identical output** and
  it makes a plain DNA run ASAN-clean.
* Duplicate FASTA names inside one file **silently merge** two records into one kernel row
  (60 records → `npos 59`).
* Names ≥100 chars overflow `new char[100]` (`mainGkmKernel.cpp:696`) — ASAN-confirmed *after* the
  CalcWmML fix (before it, ASAN aborts earlier and never reaches this).
* Alphabet size > `MAX_ALPHABET_SIZE` (compiled as 4) → **SIGSEGV** (`-A` with b=5 and b=20 both
  exit 139), because `readAlphabetFile` prints an error and `return`s without aborting.
* Kernel matrix at n=5000: 155 MB text vs 47.7 MB binary float32; R load 16.5 s vs 0.25 s (65×).
* Threading: 250+250 seqs × 500 bp, L=10 K=6 d=3 → T=1 6.55 s, T=4 2.04 s, T=20 0.97 s
  (**34 % efficiency**). Atomics are **not** the cause: 6.54 s vs 6.75 s single-threaded with
  `MULTI_THREAD_SAFE` on/off, contradicting the comment in `global.h:25`.
* `gMMProfile` is n²·(d+1)·4 bytes → 0.47 GB at n=5000, **7.45 GB at n=20 000**.
* libsvm 3.33: BSD-3, `svm.cpp` 3312 + `svm.h` 105 lines, no deps, has `PRECOMPUTED`,
  `svm_set_print_string_function`, `svm_cross_validation`. A spike training C-SVC on a real gkmSVM
  kernel file works (5-fold accuracy 0.725 at C=0.1 vs 0.775 at C=1).
* kernlab's `ksvm(type="C-svc")` **is** libsvm's `Solver`; coefficients agree to 7.2e-07 at
  `tol=1e-8`. Validate on decision values/AUC, not alphas (dual degeneracy).
* Before Phase 1: ~38 consecutive `gkmsvm_kernel()` calls in one R session aborted R (heap
  corruption), 120-char names aborted R, an empty FASTA record crashed with SIGBUS. All fixed in
  Phase 1 and covered by golden/testthat cases; the whole corpus is ASAN+UBSAN clean since then.
* **Linux differs from macOS in three places that the golden corpus caught on CI** (all fixed from
  Phase 2a up): `min`/`max` macros in `global.h` vs libstdc++ `<vector>`; `snprintf` of a buffer
  onto itself empties it under glibc; the sign of exact zeros differs between gcc and clang. Run the
  corpus on Linux (CI) before trusting a change that only passed here.

## 5. Environment and tooling on this Mac

* R 4.5.3 at `/usr/local/bin/R`. Installed and ready: `Rcpp`, `kernlab`, `seqinr`, `ROCR`,
  `testthat`. **Not** installed: the Bioconductor stack (`BSgenome`, `GenomicRanges`, …), so
  `genNullSeqs` cannot run. Since Phase 0b those are `Suggests`, so `R CMD INSTALL .` and
  `R CMD check --as-cran` (with `_R_CHECK_FORCE_SUGGESTS_=false`) work directly;
  `dev/scratch_install.sh` still works too.
* `clang++` (Apple) and `g++` available; `-fsanitize=address` works. `cmake` is **not** installed.
* `gh` authenticated as `mghandi` with **ADMIN** on the repo — you can open PRs.
* Scratch space: use the session scratchpad dir, never the repo, for build artefacts.

## 6. Conventions

* **Never commit to `master` and never push to it.** Every phase = its own branch + PR, branched
  from `master`. A `pre-push` hook in `.git/hooks/` on this clone refuses pushes to `master`
  outright; if you hit it, that is working as intended — open a PR.
* Commit messages in imperative mood, ending with the `Co-Authored-By: Claude` and `Claude-Session`
  trailers used in the existing history.
* **Do not `git add -A`** — it has already swept in `.DS_Store` and a pandoc-rendered
  `dev/REFACTORING_PLAN.html` once. Add named paths.
* Keep `dev/` out of the built package (`.Rbuildignore` already does this).
* Every claim in a PR description that is a number must come from a command you actually ran.
* The CI benchmark job is advisory (shared runners vary ±10 % run to run); the 2 % gate is
  `dev/bench.sh` best-of-5 on a quiet local machine, quoted in the PR.

## 7. Unattended (overnight) runs — decision defaults

When running autonomously, **do not stop to ask questions**. Decide with these defaults, write the
decision into the PR description or the plan, and keep going.

* **Order:** Phase 0 → 1 → 2 → 3 → 4 → 4b → 5 → 6. (Phase 0b whenever convenient after 1; Phase 7 is
  project-sized — do not start it unattended.) Finish a phase (code, tests green, PR opened, pushed)
  before starting the next. If a phase is blocked, say so in its PR, mark it in the plan, and move on.
* **Never merge a PR** and never push to `master`. Opening PRs is fine and expected.
* **Behaviour changes:** only the ones §5 of the plan already approved (`mergeByName` default off;
  the two result-changing bug fixes; binary opt-in; C++17; libsvm backend opt-in). Anything else
  that would change DNA results: don't — record it as a question in the PR instead.
* **Golden tests are the gate.** If a refactor changes DNA output, it is wrong unless the plan says
  that specific change was approved. Do not "update the expected values" to make a test pass.
* **Scope control:** a phase that is ballooning should be split — open the PR with what is done and
  working rather than carrying a huge uncommitted diff.
* **Time budget:** roughly 1–2 h per phase unattended. Prefer finishing Phase 0 and Phase 1 properly
  over starting Phase 2 in a hurry — the safety net is what makes everything after it cheap.
* **If tests need data:** generate fixtures (as `dev/baseline.sh` does); do not download genomes.
* **At the end:** update the status section of this file and post a summary comment on the newest PR.

## 8. Performance work and the GPU prototype (2026-08-30)

* `dev/PERFORMANCE_PLAN.md` is the plan for large-n kernels: measured diagnosis of the tree path
  (pass imbalance, tiling, atomics, greedy trie clones), the `-d` study (capped d reproduces the
  exact AUC at d = 4–6 for the gkmsvm3 experiments), and Plan G — the kernel and the capped
  mismatch profile as Gram matrices. The mathematics is written up in
  `../gkmsvm3/theory/gkm_kernels_as_gram_matrices.pdf` (general alphabets) and
  `gkm_kernels_constantB.pdf` (single alphabet, every matrix printed).
* `dev/metal_g3/` is a working Metal prototype of Route 2 (`make`; `gkm_metal` takes the
  `gkmsvm_kernel` flags and writes `.gkmk`): bit-exact against `gkmsvm_kernel` on every
  configuration tried, 3.5–150× faster than 28 CPU threads (numbers in its README). It is on branch
  `phase7/e-experiments-docs`, not integrated into `gkmsvm_kernel` or the R package.
* Harness knobs on the gkmsvm3 side: `GKMSVM_BIN`, `GKMSVM_SVM=libsvm`, `GKMSVM_D` (cap; the
  harness then passes `-P 2`, which is 2.8–3.9× faster than the greedy design on many cores and
  bit-identical), planned `GKMSVM_GPU=1`.
* Any change to the C++ kernel path must keep `make test` (181 golden cases) and `make oracle`
  green; any new kernel path is first compared with `gkmsvm_kernel -b` output on a small input.
