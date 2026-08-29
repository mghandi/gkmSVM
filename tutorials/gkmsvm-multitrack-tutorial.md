# gkmSVM with several tracks per sequence (DNA + methylation, …)

Since version 0.90 the gapped k-mer kernel is defined for *multi-track* sequences: every record
carries, besides its DNA, one or more aligned annotation tracks over their own alphabets — a
methylation flag per base, a chromatin state, a predicted structure, several sensor channels —
and a word of length L covers all tracks at once. The theory is that of
Mohammad-Noori, Ghareghani & Ghandi, *Generalized Gapped k-mer Filters for Robust Frequency
Estimation* (positions with different alphabet sizes; see [`dev/PHASE7_PLAN.md`](../dev/PHASE7_PLAN.md)
for how it is implemented). This page shows the workflow on a small synthetic DNA + methylation
example; `tutorials/run_multitrack_tutorial.R` runs all of it unattended in well under a minute.

## Input format: multi-track FASTA

A record is a header line followed by **one line per track**, all of the same length:

```
>pos1
ATGACAGGCCGGCACGTGCGAGAAAACATGACTCAATGCC      track 1: DNA
1010010010011111110010000001000000011011      track 2: methylation flag (0/1) per base
```

The alphabets are given per track, in the same order, with `alphabets = c("dna", "01")`: a
keyword (`dna`, `rna`, `protein`), the path of an alphabet file (one symbol per line), or simply
the symbols themselves (`"01"`, `"NUM"`, `"abc"`). A symbol outside its track's alphabet (an `N`
in the DNA) makes every window covering it skipped. `read_mfa()` / `write_mfa()` read and write
the format from R (a list with one character vector per record, one string per track).

Practical notes for methylation: mark **both** bases of a CpG (`meth2s` in the gkmsvm3
experiments), so that the reverse complement of a window — DNA complemented and reversed, every
other track reversed — is consistent; a 2-state flag gave better estimates and classification
than a 3-state encoding with an explicit non-CpG symbol.

## What the kernel computes

A window of length L over T tracks is a word of T·L positions. Positions with the same alphabet
size form a *block* (all DNA positions; all 0/1 positions), and two windows are compared by their
numbers of mismatches **per block**, so the model can gap a base without gapping its flag and
vice versa. `K` is the number of informative positions of the whole word (up to T·L), `maxnmm`
bounds the total number of mismatches (`-1` = all), and `useTgkm` chooses the full filter (0),
the truncated filter (1) or gapped k-mer counts (2) as before. Everything else — kernel matrix,
`gkmsvm_trainCV`, `gkmsvm_train`, `gkmsvm_classify`, text or binary kernels — works as for DNA.

## A worked example

The synthetic set (`run_multitrack_tutorial.R`, 60 + 60 sequences of 40 bp) contains two motifs
in every sequence, `TGACTCA` and `CACGTG`; in the positives the first is unmethylated and the
second methylated, in the negatives it is the other way round. Neither track alone can tell the
classes apart; the joint (base, methylation) pattern can.

```r
library(gkmSVM)
# kernel matrix from the two tracks (L = 6, K = 4, all mismatch levels, full filter)
gkmsvm_kernel('pos.mfa', 'neg.mfa', 'kernel_joint.txt', L = 6, K = 4, maxnmm = -1, useTgkm = 0,
              alphabets = c('dna', '01'))
# cross-validated performance, then a model
gkmsvm_trainCV('kernel_joint.txt', 'pos.mfa', 'neg.mfa', svmfnprfx = 'model_joint', nCV = 5,
               alphabets = c('dna', '01'))
gkmsvm_train('kernel_joint.txt', 'pos.mfa', 'neg.mfa', 'model_joint', alphabets = c('dna', '01'))
# score new multi-track sequences: the model file records its alphabets
gkmsvm_classify('test.mfa', 'model_joint', 'scores.txt', L = 6, K = 4, maxnmm = -1, useTgkm = 0)
```

`model_joint.gkmmodel` is a multi-track FASTA of the support vectors with their weights, plus
`#rho` and `#alphabets dna,=01` header lines; `gkmsvm_classify` takes the alphabets from it.

For comparison, the DNA track alone is the 2014 gkm-SVM (`gkmsvm_kernel('pos_dna.fa', …)` with
no `alphabets`), and the methylation track alone is a single-track kernel over a 2-letter
alphabet (`alphabets = '01'`, `addRC = FALSE`).

## Test run (2026-08-29, macOS, version 0.90.0)

`Rscript tutorials/run_multitrack_tutorial.R` (5-fold CV, C = 1):

| model | AUROC |
|---|---|
| joint (DNA + methylation) | 1.000 (AUPRC 1.000) |
| DNA only | 0.530 |
| methylation only | 0.680 |

The same construction in the gkmsvm3 reference implementation (`twoblock/examples`, 20 + 20
sequences, leave-one-out centroid classifier) gives 1.000 for the joint kernels, 0.502 for DNA
only and 0.713 for methylation only; the C++ kernels reproduce those numbers exactly (max
|ΔK| = 5·10⁻⁸).

## Command line

```
gkmsvm_kernel   -l 6 -k 4 -d -1 -t 0 -A dna,=01 pos.mfa neg.mfa kernel_joint.txt
gkmsvm_train    -A dna,=01 kernel_joint.txt pos.mfa neg.mfa model_joint
gkmsvm_classify -l 6 -k 4 -d -1 -t 0 test.mfa model_joint.gkmmodel model_joint.gkmmodel scores.txt
```

(`=SYMBOLS` is the literal-alphabet form of `-A`; a comma-separated list gives one alphabet per
track.) Long words — several tracks, or all mismatch levels of two tracks with L = 10 — switch
automatically to the *prefix-split* pass design, which needs no table of mismatch patterns
(`-P` / `passDesign` selects a design explicitly; the results are identical).

## Limits

* Filter types 3 and 4 (wildcard and mismatch kernels), pseudocounts and `gkmsvm_delta` are
  single-alphabet only.
* The number of mismatch classes is Π(block length + 1); with many *distinct* alphabet sizes it
  grows quickly, so keep the number of alphabet sizes small (several tracks may share a size).
* A multi-track file given to a run *without* its `alphabets` is read as single-track FASTA (the
  annotation lines contain no DNA symbols and are dropped) — there is no way to detect it.
