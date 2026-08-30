# Vendored verbatim from https://github.com/mghandi/gkmsvm3 (private) twoblock/twoblock_gkm.py @ 58adde2
# Same author as this package. Independent (of generalb_gkm.py) reference for the two-block case; used only
# by tests/oracle/oracle_check.py. Not part of the R package (dev/ is in .Rbuildignore).
#!/usr/bin/env python3
"""
twoblock_gkm.py -- gapped k-mer mismatch profiles, filters and kernels for the TWO-BLOCK alphabet case.

Model
-----
An l-mer has l = l1 + l2 positions: the first l1 positions ("block 1") use an alphabet of size x1, the
remaining l2 positions ("block 2") an alphabet of size x2.  Typical use: a DNA window (x1 = 4) together
with an aligned binary annotation window such as methylation status (x2 = 2), so l1 = l2 = L.

For two l-mers u, w let (m1, m2) be the numbers of mismatching positions in block 1 / block 2.
Everything needed for kernels depends only on (m1, m2):

  * gapped k-mer kernel      <f1,f2> = sum_{m1,m2} N(m1,m2)(S1,S2) * h(m1,m2),   h = C(l-m1-m2, k)
  * gkm-filter (robust count estimate) kernel:
                             <Hx1,Hx2> = x1' H H x2 = x1' H x2   (H symmetric idempotent)
                             => c(m1,m2) = H(m1,m2),  H from Theorem 3 of dif-b-app-7
  * truncated filter kernel (if H(m1,m2) <= 0, all entries with >= m1 AND >= m2 mismatches are zeroed):
                             c(m1,m2) = sum over all l-mers u of g(mm(u,u1)) g(mm(u,u2)),
                             computed with the 6-fold sum of Cm1m2.docx (per-block enumeration)

where N(m1,m2)(S1,S2) is the two-block mismatch profile: the number of pairs (u in S1, w in S2) of
l-mers with exactly m1 block-1 mismatches and m2 block-2 mismatches.

Exact rational arithmetic is used for all coefficient tables; kernel values are returned as floats
unless exact=True.

Command line
------------
  python3 twoblock_gkm.py coeffs  --L1 6 --L2 6 --x1 4 --x2 2 -k 4 [--kind filter|gkm|truncated]
  python3 twoblock_gkm.py profile --L 6 -a ACGTACGT:01100101 -b ACGTTCGA:01100001
  python3 twoblock_gkm.py kernel  --L 6 -k 4 --kind filter [--rc] [--no-normalize] pos.2fa neg.2fa -o kernel.txt

Two-track FASTA (".2fa") format: each record is a header line, a track-1 line (e.g. DNA) and a track-2
line (e.g. 0/1 methylation flags) of equal length:
    >seq1
    ACGTACGTAC
    0110010011
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from itertools import product
from math import comb, factorial, sqrt
from typing import Dict, Iterable, List, Sequence, Tuple

Lmer = Tuple[Tuple[int, ...], Tuple[int, ...]]        # (block-1 word, block-2 word) as integer tuples
Table = Dict[Tuple[int, int], Fraction]               # (m1, m2) -> coefficient


# ============================================================================ coefficient tables
def gbinom(r: int, m: int) -> int:
    """Generalised binomial C(r,m) = r(r-1)...(r-m+1)/m!  (integer r, may be negative); 0 for m < 0.
    Theorem 3 needs C(-1, m) = (-1)^m on the diagonal p = l."""
    if m < 0:
        return 0
    num = 1
    for j in range(m):
        num *= (r - j)
    return num // factorial(m)


def f_coef(l: int, k: int, p: int, n: int) -> int:
    """f(l,k,p,n) = C(l-p-1, l-n) + (-1)^k C(l-p-1, l-k-n-1)   (dif-b-app-7, Theorem 3)."""
    return gbinom(l - p - 1, l - n) + (-1) ** k * gbinom(l - p - 1, l - k - n - 1)


def H_entry(L1: int, L2: int, x1: int, x2: int, k: int, m1: int, m2: int) -> Fraction:
    """H(u,w) for l-mers with m1 / m2 mismatches in block 1 / block 2 (Theorem 3 of dif-b-app-7).
    p_i = L_i - m_i matching positions, p = p1 + p2:
       H = 1/(x1^L1 x2^L2) sum_{n=0}^{p} sum_{c1+c2=p-n} f(l,k,p,n) C(p1,c1) C(p2,c2) (1-x1)^c1 (1-x2)^c2
    """
    l, p1, p2 = L1 + L2, L1 - m1, L2 - m2
    p = p1 + p2
    tot = 0
    for n in range(p + 1):
        f = f_coef(l, k, p, n)
        if f == 0:
            continue
        for c1 in range(p1 + 1):
            c2 = p - n - c1
            if 0 <= c2 <= p2:
                tot += f * comb(p1, c1) * comb(p2, c2) * (1 - x1) ** c1 * (1 - x2) ** c2
    return Fraction(tot, x1 ** L1 * x2 ** L2)


def H_table(L1: int, L2: int, x1: int, x2: int, k: int) -> Table:
    """The full gkm-filter H indexed by (m1, m2)."""
    return {(m1, m2): H_entry(L1, L2, x1, x2, k, m1, m2) for m1 in range(L1 + 1) for m2 in range(L2 + 1)}


def truncate_table(H: Table) -> Table:
    """Truncated filter (dominance rule): whenever H(m1,m2) <= 0, every entry (a, b) with a >= m1 AND b >= m2
    is set to zero -- pairs with at least as many mismatches in *both* blocks are expected to contribute less
    to the similarity than (m1,m2), so they are dropped as well.  All surviving entries are > 0, hence the
    count estimates are non-negative.  If no entry is <= 0, H is returned unchanged."""
    killers = [(m1, m2) for (m1, m2), v in H.items() if v <= 0]
    return {(a, b): (Fraction(0) if any(a >= m1 and b >= m2 for (m1, m2) in killers) else v)
            for (a, b), v in H.items()}


def gkm_table(L1: int, L2: int, k: int) -> Table:
    """Gapped k-mer kernel weights h(m1,m2) = C(l - m1 - m2, k): number of gapped k-mers shared by two
    l-mers with m1+m2 mismatches (the k informative positions must avoid all mismatches)."""
    l = L1 + L2
    return {(m1, m2): Fraction(comb(l - m1 - m2, k) if l - m1 - m2 >= k else 0)
            for m1 in range(L1 + 1) for m2 in range(L2 + 1)}


def _block_enum(L: int, x: int, m: int, a: int, b: int) -> int:
    """Number of block words v (length L, alphabet x) with exactly a mismatches to v1 and b mismatches
    to v2, where v1, v2 have m mismatches between them.  (Ghandi et al. 2014, eq. 14, one block.)
    t   = mismatches of v placed on the L-m positions where v1 = v2 (differs from both: x-1 choices)
    a-t = mismatches placed on the m positions where v1 != v2; of those, r differ from both (x-2 choices),
          the other a-t-r equal v2.   =>  b = t + (m - (a-t)) + r  =>  r = a + b - 2t - m.
    """
    tot = 0
    for t in range(0, min(a, L - m) + 1):
        s = a - t                 # mismatches on the m "v1 != v2" positions
        if s > m:
            continue
        r = a + b - 2 * t - m
        if r < 0 or r > s:
            continue
        tot += comb(L - m, t) * (x - 1) ** t * comb(m, s) * comb(s, r) * (x - 2) ** r
    return tot


def filter_kernel_table(L1: int, L2: int, x1: int, x2: int, g: Table) -> Table:
    """Kernel coefficients c(m1,m2) = sum_{u} g(mm(u,u1)) g(mm(u,u2)) for an arbitrary (e.g. truncated)
    filter g, via the six-fold sum of Cm1m2.docx:
       c(m1,m2) = sum_{a1,b1} sum_{a2,b2} E1(m1;a1,b1) E2(m2;a2,b2) g(a1,a2) g(b1,b2)
    with E_j the per-block enumeration above (summed over t).  For the full filter g = H this returns H."""
    E1 = {(m, a, b): _block_enum(L1, x1, m, a, b) for m in range(L1 + 1) for a in range(L1 + 1) for b in range(L1 + 1)}
    E2 = {(m, a, b): _block_enum(L2, x2, m, a, b) for m in range(L2 + 1) for a in range(L2 + 1) for b in range(L2 + 1)}
    out: Table = {}
    for m1 in range(L1 + 1):
        for m2 in range(L2 + 1):
            tot = Fraction(0)
            for a1 in range(L1 + 1):
                for b1 in range(L1 + 1):
                    e1 = E1[(m1, a1, b1)]
                    if e1 == 0:
                        continue
                    for a2 in range(L2 + 1):
                        for b2 in range(L2 + 1):
                            e2 = E2[(m2, a2, b2)]
                            if e2 == 0:
                                continue
                            tot += e1 * e2 * g[(a1, a2)] * g[(b1, b2)]
            out[(m1, m2)] = tot
    return out


def kernel_coefficients(L1: int, L2: int, x1: int, x2: int, k: int, kind: str = "filter") -> Table:
    """c(m1,m2) for kind in {'gkm', 'filter', 'truncated'}."""
    if kind == "gkm":
        return gkm_table(L1, L2, k)
    H = H_table(L1, L2, x1, x2, k)
    if kind == "filter":
        return H                         # <Hx1, Hx2> = x1' H x2 because H = H' = H^2
    if kind == "truncated":
        return filter_kernel_table(L1, L2, x1, x2, truncate_table(H))
    raise ValueError(f"unknown kernel kind {kind!r}")


# ============================================================================ sequences and l-mers
DEFAULT_ALPHABET1 = "ACGT"
DEFAULT_ALPHABET2 = "01"
DEFAULT_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class TwoTrackSeq:
    """An aligned pair of tracks, e.g. DNA + methylation flags, over alphabets alpha1 / alpha2."""

    def __init__(self, track1: str, track2: str, name: str = "", alpha1: str = DEFAULT_ALPHABET1,
                 alpha2: str = DEFAULT_ALPHABET2):
        if len(track1) != len(track2):
            raise ValueError(f"{name}: track lengths differ ({len(track1)} vs {len(track2)})")
        self.name, self.track1, self.track2 = name, track1.upper(), track2
        self.alpha1, self.alpha2 = alpha1, alpha2
        self.x1, self.x2 = len(alpha1), len(alpha2)
        self._t1 = [alpha1.index(c) if c in alpha1 else -1 for c in self.track1]
        self._t2 = [alpha2.index(c) if c in alpha2 else -1 for c in self.track2]

    def __len__(self):
        return len(self.track1)

    def lmers(self, L: int, revcomp: bool = False, complement: Dict[str, str] = DEFAULT_COMPLEMENT) -> List[Lmer]:
        """All windows of length L (l1 = l2 = L); windows containing symbols outside the alphabets
        (e.g. 'N') are skipped.  With revcomp=True the reverse complement of each window is added too
        (track 1 complemented and reversed, track 2 reversed)."""
        out: List[Lmer] = []
        for i in range(len(self) - L + 1):
            b1, b2 = self._t1[i:i + L], self._t2[i:i + L]
            if min(b1) < 0 or min(b2) < 0:
                continue
            out.append((tuple(b1), tuple(b2)))
            if revcomp:
                rc1 = tuple(self.alpha1.index(complement[self.alpha1[c]]) for c in reversed(b1))
                out.append((rc1, tuple(reversed(b2))))
        return out


def read_2fa(path: str, alpha1: str = DEFAULT_ALPHABET1, alpha2: str = DEFAULT_ALPHABET2) -> List[TwoTrackSeq]:
    """Read a two-track FASTA file (header line, track-1 line, track-2 line per record)."""
    seqs, name, lines = [], None, []
    def flush():
        if name is not None:
            if len(lines) != 2:
                raise ValueError(f"record {name}: expected 2 sequence lines, got {len(lines)}")
            seqs.append(TwoTrackSeq(lines[0], lines[1], name, alpha1, alpha2))
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name, lines = line[1:].split()[0] if line[1:].split() else "", []
            else:
                lines.append(line)
    flush()
    return seqs


# ============================================================================ mismatch profiles
def _masks(lmers: Iterable[Lmer], x1: int, x2: int) -> Dict[Tuple[int, int], int]:
    """Encode each l-mer as a pair of one-hot bitmasks (bit i*x + s for symbol s at position i) and
    count multiplicities.  Mismatches between two encoded words = L - popcount(mask_a & mask_b)."""
    counts: Dict[Tuple[int, int], int] = {}
    for (b1, b2) in lmers:
        m1 = 0
        for i, s in enumerate(b1):
            m1 |= 1 << (i * x1 + s)
        m2 = 0
        for i, s in enumerate(b2):
            m2 |= 1 << (i * x2 + s)
        counts[(m1, m2)] = counts.get((m1, m2), 0) + 1
    return counts


def mismatch_profile(lmers_a: Sequence[Lmer], lmers_b: Sequence[Lmer], L1: int, L2: int, x1: int, x2: int
                     ) -> Dict[Tuple[int, int], int]:
    """Two-block mismatch profile N(m1,m2): number of pairs (u in a, w in b) with m1 block-1 and m2 block-2
    mismatches.  Pairs are counted with multiplicity (every window pair counts)."""
    A, B = _masks(lmers_a, x1, x2), _masks(lmers_b, x1, x2)
    N: Dict[Tuple[int, int], int] = {}
    for (a1, a2), ca in A.items():
        for (b1, b2), cb in B.items():
            m1 = L1 - (a1 & b1).bit_count()
            m2 = L2 - (a2 & b2).bit_count()
            N[(m1, m2)] = N.get((m1, m2), 0) + ca * cb
    return N


def inner_product(lmers_a: Sequence[Lmer], lmers_b: Sequence[Lmer], coeffs: Table, L1: int, L2: int,
                  x1: int, x2: int) -> Fraction:
    """<f_a, f_b> = sum_{m1,m2} N(m1,m2) c(m1,m2)   (exact)."""
    N = mismatch_profile(lmers_a, lmers_b, L1, L2, x1, x2)
    return sum((Fraction(n) * coeffs[key] for key, n in N.items()), Fraction(0))


def kernel_matrix(seqs: Sequence[TwoTrackSeq], L: int, k: int, kind: str = "filter", revcomp: bool = False,
                  normalize: bool = True, exact: bool = False) -> List[List[float]]:
    """Kernel matrix K[i][j] between all sequences; K = <fi,fj> / sqrt(<fi,fi><fj,fj>) if normalize."""
    if not seqs:
        return []
    x1, x2 = seqs[0].x1, seqs[0].x2
    coeffs = kernel_coefficients(L, L, x1, x2, k, kind)
    lm = [s.lmers(L, revcomp) for s in seqs]
    n = len(seqs)
    G = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            G[i][j] = G[j][i] = inner_product(lm[i], lm[j], coeffs, L, L, x1, x2)
    if normalize:
        K = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                d = G[i][i] * G[j][j]
                if exact:
                    K[i][j] = G[i][j] / _exact_sqrt(d) if d > 0 else Fraction(0)
                else:
                    K[i][j] = float(G[i][j]) / sqrt(float(d)) if d > 0 else 0.0
        return K
    return G if exact else [[float(v) for v in row] for row in G]


def _exact_sqrt(q: Fraction) -> Fraction:
    from math import isqrt
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    if n * n != q.numerator or d * d != q.denominator:
        raise ValueError("kernel normaliser is not a perfect square; use exact=False")
    return Fraction(n, d)


# ============================================================================ brute-force references (for tests)
def all_lmers(L1: int, L2: int, x1: int, x2: int) -> List[Lmer]:
    return [(b1, b2) for b1 in product(range(x1), repeat=L1) for b2 in product(range(x2), repeat=L2)]


def mismatches(u: Lmer, w: Lmer) -> Tuple[int, int]:
    return (sum(a != b for a, b in zip(u[0], w[0])), sum(a != b for a, b in zip(u[1], w[1])))


def gapped_kmer_counts(lmers: Sequence[Lmer], L1: int, L2: int, k: int) -> Dict[tuple, int]:
    """Explicit gapped k-mer feature vector (for brute-force checks)."""
    from itertools import combinations
    l = L1 + L2
    feats: Dict[tuple, int] = {}
    for (b1, b2) in lmers:
        word = b1 + b2
        for keep in combinations(range(l), k):
            key = tuple(word[i] if i in keep else None for i in range(l))
            feats[key] = feats.get(key, 0) + 1
    return feats


def filtered_counts(lmers: Sequence[Lmer], L1: int, L2: int, x1: int, x2: int, g: Table) -> Dict[Lmer, Fraction]:
    """Explicit robust count estimates  xhat(u) = sum_{w in lmers} g(mm(u,w))  over ALL l-mers u."""
    return {u: sum((g[mismatches(u, w)] for w in lmers), Fraction(0)) for u in all_lmers(L1, L2, x1, x2)}


# ============================================================================ CLI
def _fmt_table(tab: Table, L1: int, L2: int, exact: bool) -> str:
    rows = ["m1\\m2 " + " ".join(f"{m2:>12d}" for m2 in range(L2 + 1))]
    for m1 in range(L1 + 1):
        cells = [str(tab[(m1, m2)]) if exact else f"{float(tab[(m1, m2)]):.6g}" for m2 in range(L2 + 1)]
        rows.append(f"{m1:>5d} " + " ".join(f"{c:>12s}" for c in cells))
    return "\n".join(rows)


def _parse_pair(s: str) -> TwoTrackSeq:
    t1, t2 = s.split(":")
    return TwoTrackSeq(t1, t2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("coeffs", help="print the coefficient table c(m1,m2)")
    c.add_argument("--L1", type=int, required=True); c.add_argument("--L2", type=int, required=True)
    c.add_argument("--x1", type=int, default=4); c.add_argument("--x2", type=int, default=2)
    c.add_argument("-k", type=int, required=True)
    c.add_argument("--kind", choices=["filter", "gkm", "truncated"], default="filter")
    c.add_argument("--exact", action="store_true", help="print exact fractions")

    p = sub.add_parser("profile", help="mismatch profile N(m1,m2) between two track pairs")
    p.add_argument("--L", type=int, required=True)
    p.add_argument("-a", required=True, help="TRACK1:TRACK2, e.g. ACGTAC:011001")
    p.add_argument("-b", required=True)
    p.add_argument("--rc", action="store_true")

    kk = sub.add_parser("kernel", help="kernel matrix for the sequences in one or more .2fa files")
    kk.add_argument("files", nargs="+")
    kk.add_argument("--L", type=int, required=True); kk.add_argument("-k", type=int, required=True)
    kk.add_argument("--kind", choices=["filter", "gkm", "truncated"], default="filter")
    kk.add_argument("--rc", action="store_true", help="add reverse complements of every window")
    kk.add_argument("--no-normalize", action="store_true")
    kk.add_argument("-o", "--out", default="-")

    a = ap.parse_args(argv)
    if a.cmd == "coeffs":
        tab = kernel_coefficients(a.L1, a.L2, a.x1, a.x2, a.k, a.kind)
        print(f"# {a.kind} coefficients, L1={a.L1} L2={a.L2} x1={a.x1} x2={a.x2} k={a.k}")
        print(_fmt_table(tab, a.L1, a.L2, a.exact))
    elif a.cmd == "profile":
        sa, sb = _parse_pair(a.a), _parse_pair(a.b)
        N = mismatch_profile(sa.lmers(a.L, a.rc), sb.lmers(a.L, a.rc), a.L, a.L, sa.x1, sa.x2)
        print("m1 m2 N")
        for (m1, m2), n in sorted(N.items()):
            print(m1, m2, n)
    elif a.cmd == "kernel":
        seqs = [s for f in a.files for s in read_2fa(f)]
        K = kernel_matrix(seqs, a.L, a.k, a.kind, a.rc, normalize=not a.no_normalize)
        out = sys.stdout if a.out == "-" else open(a.out, "w")
        names = [s.name for s in seqs]
        print("\t".join([""] + names), file=out)
        for name, row in zip(names, K):
            print("\t".join([name] + [f"{v:.6f}" for v in row]), file=out)
        if out is not sys.stdout:
            out.close()
            print(f"wrote {a.out} ({len(seqs)} sequences, L={a.L}, k={a.k}, kind={a.kind})", file=sys.stderr)


if __name__ == "__main__":
    main()
