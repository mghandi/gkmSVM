#!/usr/bin/env python3
"""Independent numerical oracle for gkmSVM.

Cross-checks the C++ implementation against dev/oracle/generalb_gkm.py, an exact-arithmetic
(fractions.Fraction) pure-Python implementation of the gapped k-mer filters and kernels that shares
no code with src/. Two levels:

  1. coefficient tables  c[m] (t=0, full filter), cTr[m] (t=1, truncated filter), h[m] (t=2, gkm
     counts) as printed by gkmsvm_kernel, for a grid of (L, K) and alphabet sizes b in {2, 4};
  2. whole kernel matrices on the small oracle fixtures (6+6 x 60 bp), with and without reverse
     complement, both algorithms, t in {0, 1, 2}, compared entry by entry to the normalised exact
     kernel;
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
    for b, L, K in grid:
        pos, neg, alpha = (FIX + "/oracle_pos.fa", FIX + "/oracle_neg.fa", None) if b == 4 else \
                          (FIX + "/b2_pos.fa", FIX + "/b2_neg.fa", FIX + "/alphabet_b2.txt")
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
    # two-letter alphabet, tree algorithm only (the XOR algorithm is b=4 only), no reverse complement
    pos, neg = FIX + "/b2_pos.fa", FIX + "/b2_neg.fa"
    seqs = gb.read_mfa(pos, ["AB"]) + gb.read_mfa(neg, ["AB"])
    for L, K in ((6, 4), (8, 5)):
        for t, kind in KINDS.items():
            exact = gb.kernel_matrix(seqs, L, K, kind, revcomp=False)
            out = os.path.join(tmp, "k.txt")
            run_kernel(bindir, L, K, t, pos, neg, out, rc=False, alg=2, alphabet=FIX + "/alphabet_b2.txt")
            rows = [ln.split() for ln in open(out).read().splitlines() if ln.strip()]
            for i, row in enumerate(rows):
                for j, s in enumerate(row):
                    n += 1
                    if not rel_close(float(s), exact[i][j], 2e-6, 1e-7):
                        failures.append(f"kernel b=2 L={L} K={K} t={t}: K[{i}][{j}] C++ {s} vs exact {exact[i][j]!r}")
    return n


def exact_scores(test, svs, alphas, L, K, kind, d, rc, B):
    bl = gb.Blocks(B)
    coeffs = {m: (v if m[0] <= d else gb.Fraction(0)) for m, v in gb.kernel_coefficients(B, K, kind).items()}
    wt = [s.lmers(L, rc) for s in test]
    ws = [s.lmers(L, rc) for s in svs]
    nsv = [gb.inner_product(w, w, coeffs, bl) for w in ws]
    out = []
    for w in wt:
        nt = gb.inner_product(w, w, coeffs, bl)
        s = sum((gb.Fraction(a) * gb.inner_product(w, wj, coeffs, bl) / gb.Fraction(float((nt * nj) ** 0.5)) if nt * nj > 0 else 0
                 for wj, nj, a in zip(ws, nsv, alphas)), gb.Fraction(0))
        out.append(float(s))
    return out


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
        nk = check_kernels(bindir, tmp, failures)
        print(f"kernel entries:     {nk} values compared")
        nc = check_classify(bindir, tmp, failures)
        print(f"classify scores:    {nc} values compared")
    for f in failures[:50]:
        print("FAIL", f)
    if len(failures) > 50:
        print(f"... {len(failures) - 50} more")
    print("oracle: OK" if not failures else f"oracle: {len(failures)} mismatches")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
