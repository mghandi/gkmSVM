# Vendored verbatim from https://github.com/mghandi/gkmsvm3 (private) generalb/generalb_gkm.py @ 471623c
# Same author as this package. Used only by tests/oracle/oracle_check.py as an exact-arithmetic reference;
# not part of the R package (dev/ is in .Rbuildignore). Re-copy when the upstream file changes.
#!/usr/bin/env python3
"""
generalb_gkm.py -- gapped k-mer filters, mismatch profiles and kernels for a GENERAL alphabet vector
B = (b_1, ..., b_l): position i of an l-mer uses an alphabet of size b_i (any b_i >= 2, any pattern).

Key facts used (see README.md / EXAMPLES.md for derivations and worked examples)
--------------------------------------------------------------------------------
1. The gkm-filter H = A^+ A has the closed form (Theorem 2 of Mohammad-Noori, Ghareghani, Ghandi)

       H(u,w) = (1 / prod_i b_i) * sum_{G subset [l], |G| >= l-k} (-1)^{|S \\ G|} prod_{i in R \\ G} (b_i - 1)

   (R = matching positions, S = mismatching positions).  Writing w_i = b_i - 1 for i in R and w_i = -1 for
   i in S, the sum over G equals the sum over the complements T = [l] \\ G with |T| <= k of prod_{i in T} w_i,
   i.e. the sum of the elementary symmetric polynomials e_0 + e_1 + ... + e_k of (w_1, ..., w_l):

       H(u,w) = (1 / prod_i b_i) * sum_{j=0}^{k} e_j(w_1, ..., w_l).            (*)

   e_0..e_k are the coefficients of prod_i (1 + w_i z) truncated at z^k: O(l*k) per pair, no 2^l sum.

2. H(u,w) depends only on the numbers of mismatches inside each *block* of equal alphabet size.  Positions
   are grouped into r blocks by their b_i value; a mismatch class is a vector m = (m_1, ..., m_r) with
   0 <= m_j <= l_j (l_j = block length).  The number of classes is prod_j (l_j + 1).

3. Kernels: <f_1, f_2> = sum_m N(m) c(m), N(m) = number of l-mer pairs (one from each sequence) in class m.
     kind 'gkm'       : c(m) = C(l - |m|, k)                      (shared gapped k-mers, |m| = sum_j m_j)
     kind 'filter'    : c(m) = H(m)                               (H symmetric idempotent => <Hx1,Hx2> = x1' H x2)
     kind 'truncated' : c(m) = sum_u g(d(u,u1)) g(d(u,u2))  for the truncated filter g,
                        = < g , (E_1(m_1) (x) ... (x) E_r(m_r)) g >   (Kronecker/tensor contraction)
   where E_j(m_j)[a, b] = number of words of block j with a mismatches to the first and b to the second
   of two block-words that differ in m_j positions (per-block enumeration, PLoS 2014 eq. 14).

All coefficient tables are exact (fractions.Fraction).

Command line
------------
  python3 generalb_gkm.py hpair  --B 3,2,2 -k 2 -u 000 -w 010          # H(u,w) with the e_j steps shown
  python3 generalb_gkm.py coeffs --B 4,2,3,2,4 -k 3 [--kind filter|gkm|truncated] [--exact]
  python3 generalb_gkm.py kernel --L 4 -k 3 --alphabets ACGT,01,abc [--kind ...] [--rc] seqs.mfa -o K.txt

Multi-track FASTA (".mfa"): header line followed by one line per track, all of equal length:
    >seq1
    ACGTACGT
    01100101
    abcabcab
A window of length L over T tracks is an l-mer with l = T*L and B = (|alpha_1|,)*L + ... + (|alpha_T|,)*L.
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from itertools import product
from math import comb, sqrt
from typing import Dict, Iterable, List, Sequence, Tuple

Word = Tuple[int, ...]                 # an l-mer as a flat tuple of integers, symbol i in range(b_i)
Cls = Tuple[int, ...]                  # mismatch class (m_1, ..., m_r)
Table = Dict[Cls, Fraction]


# ============================================================================ blocks
class Blocks:
    """Decomposition of B into blocks of equal alphabet size (in order of first appearance)."""

    def __init__(self, B: Sequence[int]):
        self.B = tuple(int(b) for b in B)
        if any(b < 2 for b in self.B):
            raise ValueError("all alphabet sizes must be >= 2")
        self.sizes: List[int] = []                     # x_j
        self.positions: List[List[int]] = []           # positions of block j
        for i, b in enumerate(self.B):
            if b in self.sizes:
                self.positions[self.sizes.index(b)].append(i)
            else:
                self.sizes.append(b)
                self.positions.append([i])
        self.lengths = [len(p) for p in self.positions]          # l_j
        self.r = len(self.sizes)
        self.l = len(self.B)
        self.dims = [lj + 1 for lj in self.lengths]
        self.n_classes = 1
        for d in self.dims:
            self.n_classes *= d

    def classes(self) -> List[Cls]:
        return [m for m in product(*[range(d) for d in self.dims])]

    def mismatch_class(self, u: Word, w: Word) -> Cls:
        return tuple(sum(u[i] != w[i] for i in pos) for pos in self.positions)

    def representative(self, m: Cls) -> Tuple[Word, Word]:
        """u = 0...0 and w with symbol 1 at the first m_j positions of block j."""
        u = (0,) * self.l
        w = [0] * self.l
        for pos, mj in zip(self.positions, m):
            for i in pos[:mj]:
                w[i] = 1
        return u, tuple(w)

    def count_in_class(self, m: Cls) -> int:
        """Number of words w in class m relative to a fixed u:  prod_j C(l_j, m_j) (x_j - 1)^{m_j}."""
        n = 1
        for lj, xj, mj in zip(self.lengths, self.sizes, m):
            n *= comb(lj, mj) * (xj - 1) ** mj
        return n

    def describe(self) -> str:
        return "; ".join(f"block {j+1}: x={x}, positions {pos} (l_{j+1}={len(pos)})"
                         for j, (x, pos) in enumerate(zip(self.sizes, self.positions)))


# ============================================================================ H
def esp_prefix(ws: Sequence[int], k: int) -> List[int]:
    """[e_0, e_1, ..., e_k] of the numbers ws via the truncated product prod (1 + w z).  O(len(ws) * k)."""
    e = [1] + [0] * k
    for x in ws:
        for j in range(min(k, len(ws)), 0, -1):
            e[j] += x * e[j - 1]
    return e


def H_pair(B: Sequence[int], k: int, u: Word, w: Word) -> Fraction:
    """H(u,w) for arbitrary B by formula (*): (1/prod b_i) * sum_{j<=k} e_j(w_1..w_l)."""
    ws = [(b - 1) if a == c else -1 for b, a, c in zip(B, u, w)]
    den = 1
    for b in B:
        den *= b
    return Fraction(sum(esp_prefix(ws, k)), den)


def H_pair_subsets(B: Sequence[int], k: int, u: Word, w: Word) -> Fraction:
    """Theorem 2 evaluated literally as a sum over subsets G (2^l terms) -- reference implementation."""
    from itertools import combinations
    l = len(B)
    R = [i for i in range(l) if u[i] == w[i]]
    S = [i for i in range(l) if u[i] != w[i]]
    tot, den = 0, 1
    for b in B:
        den *= b
    for r in range(l - k, l + 1):
        for G in combinations(range(l), r):
            sgn = (-1) ** len([i for i in S if i not in G])
            p = 1
            for i in R:
                if i not in G:
                    p *= B[i] - 1
            tot += sgn * p
    return Fraction(tot, den)


def H_table(B: Sequence[int], k: int, blocks: Blocks | None = None) -> Table:
    """H on every mismatch class m (one O(l k) evaluation per class)."""
    bl = blocks or Blocks(B)
    return {m: H_pair(bl.B, k, *bl.representative(m)) for m in bl.classes()}


def gkm_table(B: Sequence[int], k: int, blocks: Blocks | None = None) -> Table:
    bl = blocks or Blocks(B)
    return {m: Fraction(comb(bl.l - sum(m), k) if bl.l - sum(m) >= k else 0) for m in bl.classes()}


def truncate_table(H: Table) -> Table:
    """Dominance rule: whenever H(m) <= 0, every class a with a_j >= m_j for all j is set to zero."""
    killers = [m for m, v in H.items() if v <= 0]
    return {a: (Fraction(0) if any(all(aj >= mj for aj, mj in zip(a, m)) for m in killers) else v)
            for a, v in H.items()}


# ============================================================================ truncated-filter coefficients
def block_enum(L: int, x: int, m: int, a: int, b: int) -> int:
    """E_{L,x}(m; a, b): number of words of a block (length L, alphabet x) with a mismatches to v1 and b to v2,
    where v1, v2 differ in m positions.
       sum_t C(L-m, t)(x-1)^t C(m, a-t) C(a-t, r)(x-2)^r,  r = a+b-2t-m,
       max(0, a-m) <= t <= min(a, L-m),  0 <= r <= a-t   (0^0 = 1 for x = 2)."""
    tot = 0
    for t in range(max(0, a - m), min(a, L - m) + 1):
        r = a + b - 2 * t - m
        if r < 0 or r > a - t:
            continue
        tot += comb(L - m, t) * (x - 1) ** t * comb(m, a - t) * comb(a - t, r) * (x - 2) ** r
    return tot


def block_matrices(blocks: Blocks) -> List[List[List[List[int]]]]:
    """E[j][m][a][b] for every block j and every m."""
    return [[[[block_enum(L, x, m, a, b) for b in range(L + 1)] for a in range(L + 1)] for m in range(L + 1)]
            for L, x in zip(blocks.lengths, blocks.sizes)]


def _strides(dims: Sequence[int]) -> List[int]:
    s, acc = [], 1
    for d in reversed(dims):
        s.append(acc)
        acc *= d
    return list(reversed(s))


def contract(blocks: Blocks, E, m: Cls, g_flat: List[Fraction]) -> List[Fraction]:
    """(E_1(m_1) (x) ... (x) E_r(m_r)) g, computed one mode at a time.  Cost O(n_classes * sum_j (l_j+1))."""
    dims, strides = blocks.dims, _strides(blocks.dims)
    G = list(g_flat)
    for j in range(blocks.r):
        Ej, dj, sj = E[j][m[j]], dims[j], strides[j]
        new = [Fraction(0)] * len(G)
        for idx in range(len(G)):
            aj = (idx // sj) % dj
            base = idx - aj * sj
            new[idx] = sum((Ej[aj][b] * G[base + b * sj] for b in range(dj) if Ej[aj][b]), Fraction(0))
        G = new
    return G


def filter_kernel_table(B: Sequence[int], g: Table, blocks: Blocks | None = None) -> Table:
    """c(m) = sum_u g(d(u,u1)) g(d(u,u2)) = < g, (E_1(m_1) (x) ... (x) E_r(m_r)) g > for every class m."""
    bl = blocks or Blocks(B)
    E = block_matrices(bl)
    cls = bl.classes()
    g_flat = [g[m] for m in cls]
    out: Table = {}
    for m in cls:
        Tg = contract(bl, E, m, g_flat)
        out[m] = sum((a * b for a, b in zip(g_flat, Tg)), Fraction(0))
    return out


def kernel_coefficients(B: Sequence[int], k: int, kind: str = "filter") -> Table:
    bl = Blocks(B)
    if kind == "gkm":
        return gkm_table(B, k, bl)
    H = H_table(B, k, bl)
    if kind == "filter":
        return H
    if kind == "truncated":
        return filter_kernel_table(B, truncate_table(H), bl)
    raise ValueError(f"unknown kernel kind {kind!r}")


# ============================================================================ sequences
DEFAULT_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class MultiTrackSeq:
    """T aligned tracks over alphabets alphabets[t]; a window of length L is the l-mer obtained by
    concatenating the T track windows (l = T*L, B = (|alpha_1|,)*L + ... + (|alpha_T|,)*L)."""

    def __init__(self, tracks: Sequence[str], alphabets: Sequence[str], name: str = ""):
        if len(tracks) != len(alphabets):
            raise ValueError(f"{name}: {len(tracks)} tracks but {len(alphabets)} alphabets")
        if len({len(t) for t in tracks}) != 1:
            raise ValueError(f"{name}: track lengths differ")
        self.name, self.tracks, self.alphabets = name, [t for t in tracks], list(alphabets)
        self.sizes = [len(a) for a in alphabets]
        self._codes = [[a.index(c) if c in a else -1 for c in t] for t, a in zip(tracks, alphabets)]

    def __len__(self):
        return len(self.tracks[0])

    def B(self, L: int) -> Tuple[int, ...]:
        return tuple(x for x in self.sizes for _ in range(L))

    def lmers(self, L: int, revcomp: bool = False, complement=DEFAULT_COMPLEMENT) -> List[Word]:
        out: List[Word] = []
        for i in range(len(self) - L + 1):
            wins = [c[i:i + L] for c in self._codes]
            if any(min(wn) < 0 for wn in wins):
                continue
            out.append(tuple(s for wn in wins for s in wn))
            if revcomp:      # track 1 complemented + reversed, other tracks reversed
                a0 = self.alphabets[0]
                rc0 = [a0.index(complement[a0[s]]) for s in reversed(wins[0])]
                rest = [list(reversed(wn)) for wn in wins[1:]]
                out.append(tuple(s for wn in [rc0] + rest for s in wn))
        return out


def read_mfa(path: str, alphabets: Sequence[str]) -> List[MultiTrackSeq]:
    seqs, name, lines = [], None, []
    def flush():
        if name is not None:
            if len(lines) != len(alphabets):
                raise ValueError(f"record {name}: expected {len(alphabets)} track lines, got {len(lines)}")
            seqs.append(MultiTrackSeq(lines, alphabets, name))
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name, lines = (line[1:].split() or [""])[0], []
            else:
                lines.append(line)
    flush()
    return seqs


# ============================================================================ mismatch profile & kernels
def _masks(words: Iterable[Word], blocks: Blocks) -> Dict[Tuple[int, ...], int]:
    """One one-hot bitmask per block; mismatches in block j = l_j - popcount(mask_a & mask_b)."""
    counts: Dict[Tuple[int, ...], int] = {}
    for w in words:
        key = []
        for x, pos in zip(blocks.sizes, blocks.positions):
            mk = 0
            for t, i in enumerate(pos):
                mk |= 1 << (t * x + w[i])
            key.append(mk)
        key = tuple(key)
        counts[key] = counts.get(key, 0) + 1
    return counts


def mismatch_profile(words_a: Sequence[Word], words_b: Sequence[Word], blocks: Blocks) -> Dict[Cls, int]:
    """N(m): number of pairs (u in a, w in b) in mismatch class m (with multiplicity)."""
    A, Bm = _masks(words_a, blocks), _masks(words_b, blocks)
    N: Dict[Cls, int] = {}
    for ka, ca in A.items():
        for kb, cb in Bm.items():
            m = tuple(lj - (x & y).bit_count() for lj, x, y in zip(blocks.lengths, ka, kb))
            N[m] = N.get(m, 0) + ca * cb
    return N


def inner_product(words_a, words_b, coeffs: Table, blocks: Blocks) -> Fraction:
    N = mismatch_profile(words_a, words_b, blocks)
    return sum((Fraction(n) * coeffs[m] for m, n in N.items()), Fraction(0))


def kernel_matrix(seqs: Sequence[MultiTrackSeq], L: int, k: int, kind: str = "filter", revcomp: bool = False,
                  normalize: bool = True) -> List[List[float]]:
    if not seqs:
        return []
    B = seqs[0].B(L)
    bl = Blocks(B)
    coeffs = kernel_coefficients(B, k, kind)
    words = [s.lmers(L, revcomp) for s in seqs]
    n = len(seqs)
    G = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            G[i][j] = G[j][i] = inner_product(words[i], words[j], coeffs, bl)
    if not normalize:
        return [[float(v) for v in row] for row in G]
    return [[(float(G[i][j]) / sqrt(float(G[i][i] * G[j][j])) if G[i][i] * G[j][j] > 0 else 0.0)
             for j in range(n)] for i in range(n)]


# ============================================================================ brute-force helpers (tests / examples)
def all_words(B: Sequence[int]) -> List[Word]:
    return [w for w in product(*[range(b) for b in B])]


def gapped_kmer_counts(words: Sequence[Word], k: int) -> Dict[tuple, int]:
    from itertools import combinations
    feats: Dict[tuple, int] = {}
    for w in words:
        for keep in combinations(range(len(w)), k):
            key = tuple(w[i] if i in keep else None for i in range(len(w)))
            feats[key] = feats.get(key, 0) + 1
    return feats


def filtered_counts(words: Sequence[Word], B: Sequence[int], g: Table, blocks: Blocks) -> Dict[Word, Fraction]:
    return {u: sum((g[blocks.mismatch_class(u, w)] for w in words), Fraction(0)) for u in all_words(B)}


# ============================================================================ CLI
def _parse_B(s: str) -> Tuple[int, ...]:
    return tuple(int(t) for t in s.split(","))


def _parse_word(s: str) -> Word:
    return tuple(int(c) for c in s)


def explain_H_pair(B, k, u, w) -> str:
    ws = [(b - 1) if a == c else -1 for b, a, c in zip(B, u, w)]
    lines = [f"B = {tuple(B)}, k = {k}, u = {''.join(map(str,u))}, w = {''.join(map(str,w))}",
             f"match/mismatch weights w_i = b_i-1 (match) or -1 (mismatch): {ws}"]
    e = [1] + [0] * k
    for i, x in enumerate(ws):
        for j in range(min(k, len(ws)), 0, -1):
            e[j] += x * e[j - 1]
        lines.append(f"  after position {i+1} (w={x:>2d}): e_0..e_{k} = {e}")
    den = 1
    for b in B:
        den *= b
    lines.append(f"H(u,w) = (e_0 + ... + e_{k}) / prod(b_i) = {sum(e)} / {den} = {Fraction(sum(e), den)}")
    lines.append(f"check (Theorem 2 subset sum): {H_pair_subsets(B, k, u, w)}")
    return "\n".join(lines)


def format_table(tab: Table, blocks: Blocks, exact: bool) -> str:
    rows = [f"# blocks: {blocks.describe()}", "# class m = (" + ", ".join(f"m_{j+1}" for j in range(blocks.r)) + ")   value"]
    for m in blocks.classes():
        v = tab[m]
        rows.append(f"{str(m):<{4*blocks.r+4}s} {str(v) if exact else f'{float(v):.6g}'}")
    return "\n".join(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("hpair"); h.add_argument("--B", required=True); h.add_argument("-k", type=int, required=True)
    h.add_argument("-u", required=True); h.add_argument("-w", required=True)
    c = sub.add_parser("coeffs"); c.add_argument("--B", required=True); c.add_argument("-k", type=int, required=True)
    c.add_argument("--kind", choices=["filter", "gkm", "truncated"], default="filter"); c.add_argument("--exact", action="store_true")
    kk = sub.add_parser("kernel"); kk.add_argument("files", nargs="+"); kk.add_argument("--L", type=int, required=True)
    kk.add_argument("-k", type=int, required=True); kk.add_argument("--alphabets", required=True, help="comma separated, e.g. ACGT,01")
    kk.add_argument("--kind", choices=["filter", "gkm", "truncated"], default="filter"); kk.add_argument("--rc", action="store_true")
    kk.add_argument("--no-normalize", action="store_true"); kk.add_argument("-o", "--out", default="-")
    a = ap.parse_args(argv)
    if a.cmd == "hpair":
        print(explain_H_pair(_parse_B(a.B), a.k, _parse_word(a.u), _parse_word(a.w)))
    elif a.cmd == "coeffs":
        B = _parse_B(a.B)
        print(f"# {a.kind} coefficients, B={B}, k={a.k}")
        print(format_table(kernel_coefficients(B, a.k, a.kind), Blocks(B), a.exact))
    elif a.cmd == "kernel":
        alphas = a.alphabets.split(",")
        seqs = [s for f in a.files for s in read_mfa(f, alphas)]
        K = kernel_matrix(seqs, a.L, a.k, a.kind, a.rc, normalize=not a.no_normalize)
        out = sys.stdout if a.out == "-" else open(a.out, "w")
        names = [s.name for s in seqs]
        print("\t".join([""] + names), file=out)
        for name, row in zip(names, K):
            print("\t".join([name] + [f"{v:.6f}" for v in row]), file=out)
        if out is not sys.stdout:
            out.close()
            print(f"wrote {a.out} ({len(seqs)} sequences, L={a.L}, k={a.k}, kind={a.kind}, B={seqs[0].B(a.L)})", file=sys.stderr)


if __name__ == "__main__":
    main()
