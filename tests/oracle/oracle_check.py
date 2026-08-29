#!/usr/bin/env python3
"""Independent numerical oracle for gkmSVM.

Cross-checks the C++ implementation against dev/oracle/generalb_gkm.py, an exact-arithmetic
(fractions.Fraction) pure-Python implementation of the gapped k-mer filters and kernels that shares
no code with src/. Two levels:

  1. coefficient tables  c[m] (t=0, full filter), cTr[m] (t=1, truncated filter), h[m] (t=2, gkm
     counts) as printed by gkmsvm_kernel, for a grid of (L, K) and alphabet sizes b in {2, 4, 5, 20};
  2. whole kernel matrices on the small oracle fixtures (6+6 x 60 bp), with and without reverse
     complement, both algorithms, t in {0, 1, 2}, compared entry by entry to the normalised exact
     kernel; the same for b in {2, 5, 20} without reverse complement (Phase 5 gate: b=5 and protein
     complete and match the exact implementation);
  2b. (Phase 7) the general-B coefficient tables for several alphabet vectors B, every mismatch
     class and the three kinds, printed by `gkmsvm_selftest coeffs`, against the exact tables;
  3. gkmsvm_classify scores for the fixed SV model (sv_svseq.fa / sv_svalpha.out) against the exact
     score  sum_j alpha_j <x_i,x_j>_d / (|x_i|_d |x_j|_d), where <.,.>_d keeps mismatch levels
     m <= d only (the -d option; norm and score truncated alike, Phase 1 decision 4), for d in {3, L}.

Tolerance: gkmsvm_kernel prints %e (7 significant digits), so agreement is required to a relative
1e-6 for tables and 2e-6 for kernel entries (absolute 1e-7 near zero).

  python3 tests/oracle/oracle_check.py --bin build
"""
import argparse, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FIX = os.path.join(ROOT, "tests", "golden", "fixtures")
sys.path.insert(0, os.path.join(ROOT, "dev", "oracle"))
import generalb_gkm as gb  # noqa: E402

KINDS = {0: "filter", 1: "truncated", 2: "gkm"}


def run_kernel(bindir, L, K, t, pos, neg, out, rc=True, alg=2, alphabet=None):
    argv = [os.path.join(bindir, "gkmsvm_kernel"), "-l", str(L), "-k", str(K), "-d", "-1", "-t", str(t), "-a", str(alg)]
    if not rc:
        argv.append("-R")
    if alphabet:
        argv += ["-A", alphabet]
    argv += [pos, neg, out]
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(argv)} failed: {p.stderr}")
    coef = {int(m): float(v) for m, v in re.findall(r"c\[(\d+)\] = (\S+)", p.stdout)}
    return coef, p.stdout


def rel_close(a, b, rtol, atol):
    return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))


def check_tables(bindir, tmp, failures):
    n = 0
    grid = [(b, L, K) for b in (4, 2) for (L, K) in ((4, 2), (6, 4), (8, 5), (10, 6), (12, 8), (7, 3))]
    grid += [(20, L, K) for (L, K) in ((4, 2), (5, 3), (6, 4))] + [(5, L, K) for (L, K) in ((6, 4), (8, 5))]
    inputs = {4: (FIX + "/oracle_pos.fa", FIX + "/oracle_neg.fa", None),
              2: (FIX + "/b2_pos.fa", FIX + "/b2_neg.fa", FIX + "/alphabet_b2.txt"),
              20: (FIX + "/protein_pos.fa", FIX + "/protein_neg.fa", FIX + "/alphabet_protein.txt"),
              5: (FIX + "/b5_pos.fa", FIX + "/b5_neg.fa", FIX + "/alphabet_b5.txt")}
    for b, L, K in grid:
        pos, neg, alpha = inputs[b]
        for t, kind in KINDS.items():
            coef, _ = run_kernel(bindir, L, K, t, pos, neg, os.path.join(tmp, "k.txt"), alphabet=alpha)
            table = {m[0]: float(v) for m, v in gb.kernel_coefficients((b,) * L, K, kind).items()}
            for m, v in coef.items():
                n += 1
                if not rel_close(v, table[m], 1e-6, 1e-12):
                    failures.append(f"table b={b} L={L} K={K} t={t}: c[{m}] C++ {v!r} vs exact {table[m]!r}")
            # every non-zero exact coefficient within the printed range must be present in C++ output
            for m, v in table.items():
                if v != 0 and m not in coef:
                    failures.append(f"table b={b} L={L} K={K} t={t}: exact c[{m}]={v} but C++ did not use m={m}")
    return n


def check_generalb_tables(bindir, tmp, failures):
    """Phase 7: the general-B coefficient tables (src/GeneralB.cpp, printed by `gkmsvm_selftest coeffs B k`)
    against the exact tables for every mismatch class, several alphabet vectors and all three kinds."""
    n = 0
    grid = [((4,) * 6 + (2,) * 6, 4), ((4,) * 4 + (2,) * 4 + (3,) * 4, 5), ((20,) * 4 + (3,) * 4, 3),
            ((2,) * 8, 5), ((4,) * 10, 6), ((4, 2, 3, 2, 4), 3), ((4,) * 3 + (2,) * 3 + (4,) * 3, 4),
            ((5,) * 5 + (2,) * 5, 3)]
    exe = os.path.join(bindir, "gkmsvm_selftest")
    for B, k in grid:
        p = subprocess.run([exe, "coeffs", ",".join(map(str, B)), str(k)], capture_output=True, text=True)
        if p.returncode != 0:
            failures.append(f"gkmsvm_selftest coeffs {B} {k} failed: {p.stderr}")
            continue
        rows = {}
        autod = {}
        for line in p.stdout.splitlines():
            if line.startswith("m="):
                f = dict(x.split("=", 1) for x in line.split("\t"))
                rows[tuple(int(v) for v in f["m"].split(","))] = f
            elif line.startswith("autod"):
                autod = dict((x.split("=")[0], int(x.split("=")[1])) for x in line.split()[1:])
        for kind in ("filter", "truncated", "gkm"):
            table = gb.kernel_coefficients(B, k, kind)
            scale = max(abs(float(v)) for v in table.values())
            for m, v in table.items():
                n += 1
                got = float(rows[m][kind])
                if abs(got - float(v)) > 1e-9 * scale:
                    failures.append(f"generalB table B={B} k={k} {kind}: c{m} C++ {got!r} vs exact {v} ({float(v)!r})")
        # dominance-rule truncated impulse response and the automatic -d bounds
        g = gb.truncate_table(gb.H_table(B, k))
        scale = max(abs(float(v)) for v in g.values())
        for m, v in g.items():
            if abs(float(rows[m]["g"]) - float(v)) > 1e-9 * scale:
                failures.append(f"generalB g B={B} k={k}: g{m} C++ {rows[m]['g']} vs exact {v}")
        ctr = gb.kernel_coefficients(B, k, "truncated")
        exp_tr = max((sum(m) for m, v in ctr.items() if v != 0), default=0)
        l = len(B)
        if autod != {"filter": l, "truncated": exp_tr, "gkm": l - k}:
            failures.append(f"generalB auto -d B={B} k={k}: C++ {autod} vs expected filter={l} truncated={exp_tr} gkm={l - k}")
    return n


def check_kernels(bindir, tmp, failures):
    n = 0
    configs = [(4, 6, 4), (4, 8, 5), (4, 10, 6)]
    for b, L, K in configs:
        if b == 4:
            pos, neg, alpha, alphabets = FIX + "/oracle_pos.fa", FIX + "/oracle_neg.fa", None, ["ACGT"]
        seqs = gb.read_mfa(pos, alphabets) + gb.read_mfa(neg, alphabets)
        for rc in (True, False):
            for t, kind in KINDS.items():
                exact = gb.kernel_matrix(seqs, L, K, kind, revcomp=rc)
                for alg in (1, 2):
                    out = os.path.join(tmp, "k.txt")
                    run_kernel(bindir, L, K, t, pos, neg, out, rc=rc, alg=alg)
                    rows = [ln.split() for ln in open(out).read().splitlines() if ln.strip()]
                    if len(rows) != len(seqs):
                        failures.append(f"kernel L={L} K={K} t={t} rc={rc} alg={alg}: {len(rows)} rows, expected {len(seqs)}")
                        continue
                    for i, row in enumerate(rows):
                        if len(row) != i + 1:
                            failures.append(f"kernel L={L} K={K} t={t} rc={rc} alg={alg}: row {i} has {len(row)} entries")
                            break
                        for j, s in enumerate(row):
                            n += 1
                            v = float(s)
                            if not rel_close(v, exact[i][j], 2e-6, 1e-7):
                                failures.append(f"kernel L={L} K={K} t={t} rc={rc} alg={alg}: K[{i}][{j}] C++ {v!r} vs exact {exact[i][j]!r}")
    # non-DNA alphabets, tree algorithm only (the XOR algorithm is b=4 only), no reverse complement:
    # b=2, b=5 (the plan's "crashes today" case) and protein b=20 (Phase 5 gate)
    others = [(2, "AB", "b2", "alphabet_b2.txt", ((6, 4), (8, 5))),
              (5, "ACGTN", "b5", "alphabet_b5.txt", ((6, 4), (8, 5))),
              (20, "ACDEFGHIKLMNPQRSTVWY", "protein", "alphabet_protein.txt", ((4, 2), (5, 3), (6, 4)))]
    for b, letters, stem, afile, LKs in others:
        pos, neg = FIX + f"/{stem}_pos.fa", FIX + f"/{stem}_neg.fa"
        seqs = gb.read_mfa(pos, [letters]) + gb.read_mfa(neg, [letters])
        for L, K in LKs:
            for t, kind in KINDS.items():
                exact = gb.kernel_matrix(seqs, L, K, kind, revcomp=False)
                out = os.path.join(tmp, "k.txt")
                run_kernel(bindir, L, K, t, pos, neg, out, rc=False, alg=2, alphabet=FIX + "/" + afile)
                rows = [ln.split() for ln in open(out).read().splitlines() if ln.strip()]
                if len(rows) != len(seqs):
                    failures.append(f"kernel b={b} L={L} K={K} t={t}: {len(rows)} rows, expected {len(seqs)}")
                    continue
                for i, row in enumerate(rows):
                    for j, s in enumerate(row):
                        n += 1
                        if not rel_close(float(s), exact[i][j], 2e-6, 1e-7):
                            failures.append(f"kernel b={b} L={L} K={K} t={t}: K[{i}][{j}] C++ {s} vs exact {exact[i][j]!r}")
    return n



def exact_multitrack_kernel(seqs, L, K, kind, rc, d):
    """Normalised kernel with the mismatch bound d (classes with |m| > d dropped), exact arithmetic."""
    from fractions import Fraction
    from math import sqrt
    B = seqs[0].B(L)
    bl = gb.Blocks(B)
    coeffs = gb.kernel_coefficients(B, K, kind)
    if d is not None:
        coeffs = {m: (v if sum(m) <= d else Fraction(0)) for m, v in coeffs.items()}
    words = [s.lmers(L, rc) for s in seqs]
    n = len(seqs)
    G = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            G[i][j] = G[j][i] = gb.inner_product(words[i], words[j], coeffs, bl)
    # the diagonal is written as 1.0 by gkmsvm_kernel even for a record without any window (norm 0)
    return [[(1.0 if i == j else float(G[i][j]) / sqrt(float(G[i][i] * G[j][j])) if i == j or G[i][i] * G[j][j] > 0 else 0.0)
             for j in range(n)] for i in range(n)]


def check_multitrack_kernels(bindir, tmp, failures):
    """Phase 7: two-track kernels (DNA + a second track) against generalb_gkm (any B) and, for the
    two-block case, against the independent twoblock_gkm implementation; with and without reverse
    complement, t in {0,1,2}, automatic and bounded -d; a fixture with N's, a 3-state second track,
    a record shorter than the window and a record whose every window contains an N."""
    import twoblock_gkm as tb
    n = 0
    cases = [("meth", ["ACGT", "01"], "dna,=01", 6, 4, [None, 4]),
             ("mt3", ["ACGT", "NUM"], "dna,=NUM", 5, 3, [None, 3])]
    for stem, alphabets, spec, L, K, ds in cases:
        pos, neg = FIX + f"/{stem}_pos.mfa", FIX + f"/{stem}_neg.mfa"
        seqs = gb.read_mfa(pos, alphabets) + gb.read_mfa(neg, alphabets)
        for rc in (True, False):
            for t, kind in KINDS.items():
                for d in ds:
                    exact = exact_multitrack_kernel(seqs, L, K, kind, rc, d)
                    out = os.path.join(tmp, "k.txt")
                    argv = [os.path.join(bindir, "gkmsvm_kernel"), "-l", str(L), "-k", str(K), "-d", str(-1 if d is None else d), "-t", str(t), "-A", spec]
                    if not rc:
                        argv.append("-R")
                    argv += [pos, neg, out]
                    p = subprocess.run(argv, capture_output=True, text=True)
                    if p.returncode != 0:
                        failures.append(f"{' '.join(argv)} failed: {p.stdout[-300:]}")
                        continue
                    rows = [ln.split() for ln in open(out).read().splitlines() if ln.strip()]
                    if len(rows) != len(seqs):
                        failures.append(f"multitrack {stem} t={t} rc={rc} d={d}: {len(rows)} rows, expected {len(seqs)}")
                        continue
                    for i, row in enumerate(rows):
                        for j, sv in enumerate(row):
                            n += 1
                            if not rel_close(float(sv), exact[i][j], 2e-6, 1e-7):
                                failures.append(f"multitrack {stem} t={t} rc={rc} d={d}: K[{i}][{j}] C++ {sv} vs exact {exact[i][j]!r}")
                    # the independent two-block implementation (automatic d only: it has no bound)
                    if stem == "meth" and d is None:
                        tseqs = tb.read_2fa(pos) + tb.read_2fa(neg)
                        tk = tb.kernel_matrix(tseqs, L, K, kind, revcomp=rc)
                        for i, row in enumerate(rows):
                            for j, sv in enumerate(row):
                                n += 1
                                if not rel_close(float(sv), tk[i][j], 2e-6, 1e-7):
                                    failures.append(f"multitrack {stem} t={t} rc={rc} (twoblock): K[{i}][{j}] C++ {sv} vs twoblock {tk[i][j]!r}")
    return n

def exact_scores(test, svs, alphas, L, K, kind, d, rc, B, rho=0.0):
    bl = gb.Blocks(B)
    coeffs = {m: (v if sum(m) <= d else gb.Fraction(0)) for m, v in gb.kernel_coefficients(B, K, kind).items()}
    wt = [s.lmers(L, rc) for s in test]
    ws = [s.lmers(L, rc) for s in svs]
    nsv = [gb.inner_product(w, w, coeffs, bl) for w in ws]
    out = []
    for w in wt:
        nt = gb.inner_product(w, w, coeffs, bl)
        s = sum((gb.Fraction(a) * gb.inner_product(w, wj, coeffs, bl) / gb.Fraction(float((nt * nj) ** 0.5)) if nt * nj > 0 else 0
                 for wj, nj, a in zip(ws, nsv, alphas)), gb.Fraction(0))
        out.append(float(s) - rho)
    return out


def check_classify_multitrack(bindir, tmp, failures):
    """Phase 7: two-track classify scores (the .gkmmodel with rho and #alphabets, and the legacy pair)
    against the exact scores; t in {0,1,2}, +-RC, automatic and bounded -d; the test set has a record
    with an N and a record shorter than the window (score -rho)."""
    n = 0
    alphabets = ["ACGT", "01"]
    test = gb.read_mfa(FIX + "/test_meth.mfa", alphabets)
    svs = gb.read_mfa(FIX + "/sv_meth_svseq.mfa", alphabets)
    alphas = [float(ln.split()[1]) for ln in open(FIX + "/sv_meth_svalpha.out") if ln.strip()]
    rho = [float(ln.split()[1]) for ln in open(FIX + "/sv_meth.gkmmodel") if ln.startswith("#rho")][0]
    L, K = 6, 4
    B = (4,) * L + (2,) * L
    for model in ("pair", "gkmmodel"):
        for rc in (True, False):
            for t, kind in KINDS.items():
                for d in (4, 2 * L):
                    exact = exact_scores(test, svs, alphas, L, K, kind, d, rc, B, rho if model == "gkmmodel" else 0.0)
                    out = os.path.join(tmp, "s.txt")
                    argv = [os.path.join(bindir, "gkmsvm_classify"), "-l", str(L), "-k", str(K), "-d", str(d if d < 2 * L else -1), "-t", str(t)]
                    if not rc:
                        argv.append("-R")
                    if model == "pair":
                        argv += ["-A", "dna,=01", FIX + "/test_meth.mfa", FIX + "/sv_meth_svseq.mfa", FIX + "/sv_meth_svalpha.out", out]
                    else:
                        argv += [FIX + "/test_meth.mfa", FIX + "/sv_meth.gkmmodel", FIX + "/sv_meth.gkmmodel", out]
                    p = subprocess.run(argv, capture_output=True, text=True)
                    if p.returncode != 0:
                        failures.append(f"classify multitrack {model} t={t} d={d} rc={rc}: exit {p.returncode}: {p.stdout[-300:]}")
                        continue
                    got = [float(ln.split()[1]) for ln in open(out) if ln.strip()]
                    if len(got) != len(exact):
                        failures.append(f"classify multitrack {model} t={t} d={d} rc={rc}: {len(got)} scores, expected {len(exact)}")
                        continue
                    for i, (g, e) in enumerate(zip(got, exact)):
                        n += 1
                        if abs(g - e) > 2e-6:
                            failures.append(f"classify multitrack {model} t={t} d={d} rc={rc}: seq {i} C++ {g} vs exact {e:.7f}")
    return n


def check_classify(bindir, tmp, failures):
    n = 0
    setups = [("test_seqs.fa", "sv_svseq.fa", "sv_svalpha.out", ["ACGT"], None, [(10, 6), (8, 5)]),
              ("test_b2.fa", "sv_b2_svseq.fa", "sv_b2_svalpha.out", ["AB"], FIX + "/alphabet_b2.txt", [(8, 5)])]
    for tf, svf, af, alphabets, alpha_file, LKs in setups:
        test = gb.read_mfa(FIX + "/" + tf, alphabets)
        svs = gb.read_mfa(FIX + "/" + svf, alphabets)
        alphas = [float(ln.split()[1]) for ln in open(FIX + "/" + af) if ln.strip()]
        for L, K in LKs:
            for rc in ((True, False) if alpha_file is None else (False,)):
                for t, kind in KINDS.items():
                    for d in (3, L):
                        for alg in (1, 2):
                            exact = exact_scores(test, svs, alphas, L, K, kind, d, rc, (len(alphabets[0]),) * L)
                            out = os.path.join(tmp, "s.txt")
                            argv = [os.path.join(bindir, "gkmsvm_classify"), "-l", str(L), "-k", str(K), "-d", str(d), "-t", str(t), "-a", str(alg)]
                            if not rc:
                                argv.append("-R")
                            if alpha_file:
                                argv += ["-A", alpha_file]
                            argv += [FIX + "/" + tf, FIX + "/" + svf, FIX + "/" + af, out]
                            p = subprocess.run(argv, capture_output=True, text=True)
                            if p.returncode != 0:
                                failures.append(f"classify {tf} L={L} K={K} t={t} d={d} rc={rc} alg={alg}: exit {p.returncode}")
                                continue
                            got = [float(ln.split()[1]) for ln in open(out) if ln.strip()]
                            if alg == 1 and d < L:
                                continue  # the XOR/hash-table path (-a 1) sums all mismatch levels regardless of -d; documented
                            for i, (g, e) in enumerate(zip(got, exact)):
                                n += 1
                                if abs(g - e) > 2e-6:  # scores are printed with %f (6 decimals)
                                    failures.append(f"classify {tf} L={L} K={K} t={t} d={d} rc={rc} alg={alg}: seq {i} C++ {g} vs exact {e:.7f}")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default="build")
    a = ap.parse_args()
    bindir = os.path.abspath(a.bin)
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        nt = check_tables(bindir, tmp, failures)
        print(f"coefficient tables: {nt} values compared")
        ng = check_generalb_tables(bindir, tmp, failures)
        print(f"general-B tables:   {ng} values compared")
        nk = check_kernels(bindir, tmp, failures)
        print(f"kernel entries:     {nk} values compared")
        nm = check_multitrack_kernels(bindir, tmp, failures)
        print(f"multi-track kernel: {nm} values compared")
        nc = check_classify(bindir, tmp, failures)
        print(f"classify scores:    {nc} values compared")
        ncm = check_classify_multitrack(bindir, tmp, failures)
        print(f"multi-track scores: {ncm} values compared")
    for f in failures[:50]:
        print("FAIL", f)
    if len(failures) > 50:
        print(f"... {len(failures) - 50} more")
    print("oracle: OK" if not failures else f"oracle: {len(failures)} mismatches")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
