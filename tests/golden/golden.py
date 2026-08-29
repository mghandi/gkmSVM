#!/usr/bin/env python3
"""Golden tests for the gkmSVM command-line binaries (and, through tests/testthat, the R wrappers).

  python3 tests/golden/golden.py check  --bin build          # run every case, compare byte-for-byte
  python3 tests/golden/golden.py freeze --bin build          # (re)write tests/golden/expected/*.out
  python3 tests/golden/golden.py list [--tag T] [--filter RE]

Cases live in cases.tsv (one per line, explicit columns; see build_argv() for the CLI mapping and
tests/testthat/helper-golden.R for the R mapping). Expected outputs are frozen from master @ 222cc50
(gkmSVM 0.80) and must NOT be re-frozen to make a failing test pass unless dev/REFACTORING_PLAN.md
approved that exact behaviour change; in that case re-freeze only the affected cases with --filter
and explain the diff in the PR.

Tags: `numeric` compares numbers to a relative 1e-5 and everything else exactly (SVM solver output);
`xfail-<phase>` marks a case that crashes today (passes while it crashes, fails as XPASS once it
does not); `expect-error` marks invalid input that must be rejected with a clean non-zero exit and an
ERROR message (no crash, no output file).

Comparison is exact on the bytes of the output file. --tol X additionally reports, for a failing
case, the largest relative numeric difference so a formatting-only change can be told apart from a
numeric one (it never turns a failure into a pass).
"""
import argparse, csv, os, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
EXP = os.path.join(HERE, "expected")
ACT = os.path.join(HERE, "actual")


def load_cases(path=os.path.join(HERE, "cases.tsv")):
    with open(path, newline="") as f:
        cases = list(csv.DictReader(f, delimiter="\t"))
    # case names become file names: they must also be unique on a case-insensitive file system
    # (macOS): two cases differing only by case once overwrote each other's frozen output
    seen = {}
    for c in cases:
        k = c["name"].lower()
        if k in seen:
            sys.exit(f"cases.tsv: names {seen[k]!r} and {c['name']!r} differ only by letter case")
        seen[k] = c["name"]
    return cases


def build_argv(c, bindir, outfile):
    """CLI argv for one case. Only options with a value are passed, so that defaults are exercised."""
    if c["tool"] == "train":
        # positional: kernel pos neg prefix; the compared output is <prefix>_svalpha.out (see run_case)
        argv = [os.path.join(bindir, "gkmsvm_train"), "-q"]
        for col, flag in (("C", "-c"), ("w", "-w"), ("v", "-v")):
            if c.get(col, ""):
                argv += [flag, c[col]]
        if c.get("S", "") == "1":
            argv.append("-S")
        argv += [os.path.join(EXP, c["kernel"]), os.path.join(FIX, c["pos"]), os.path.join(FIX, c["neg"]), outfile]
        return argv
    exe = os.path.join(bindir, "gkmsvm_kernel" if c["tool"] == "kernel" else "gkmsvm_classify")
    argv = [exe]
    for col, flag in (("L", "-l"), ("K", "-k"), ("d", "-d"), ("t", "-t"), ("alg", "-a"),
                      ("M", "-M"), ("lam", "-L"), ("T", "-T"), ("batch", "-b"), ("m", "-m"), ("n", "-n"), ("r", "-r")):
        if c.get(col, "") != "":
            argv += [flag, c[col]]
    if c["A"]:
        spec = c["A"]
        if spec not in ("dna", "rna", "protein") and "," not in spec and not spec.startswith("="):
            spec = os.path.join(FIX, spec)  # an alphabet file; keywords and multi-track specs (Phase 7) pass through
        argv += ["-A", spec]
    if c["RC"] == "0":
        argv.append("-R")
    if c["p"] == "1":
        argv.append("-p")
    if c.get("legacy", "") == "1" and c["tool"] == "classify":
        argv.append("-y")
    if c.get("N", "") == "1":
        argv.append("-N")
    if c.get("bin", "") == "1" and c["tool"] == "kernel":
        argv.append("-b")
    if c["tool"] == "kernel":
        argv += [os.path.join(FIX, c["pos"]), os.path.join(FIX, c["neg"]), outfile]
    else:
        argv += [os.path.join(FIX, c["seq"]), os.path.join(FIX, c["svseq"]), os.path.join(FIX, c["svalpha"]), outfile]
    return argv


def run_case(c, bindir):
    d = os.path.join(ACT, c["name"])
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "out.txt")
    argv = build_argv(c, bindir, out if c["tool"] != "train" else os.path.join(d, "model"))
    t0 = time.time()
    with open(os.path.join(d, "stdout.txt"), "wb") as so, open(os.path.join(d, "stderr.txt"), "wb") as se:
        rc = subprocess.call(argv, stdout=so, stderr=se)
    if c["tool"] == "train" and os.path.exists(os.path.join(d, "model_svalpha.out")):
        os.rename(os.path.join(d, "model_svalpha.out"), out)
    return dict(name=c["name"], rc=rc, out=out, argv=argv, secs=time.time() - t0, dir=d)


NUM = re.compile(rb"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def numeric_equal(a, b, rtol):
    """Same non-numeric text and all numbers within rtol (for solver output whose last digits vary by compiler)."""
    if NUM.sub(b"#", a) != NUM.sub(b"#", b):
        return False
    return max_rel_diff(a, b) <= rtol


def max_rel_diff(a, b):
    xa, xb = NUM.findall(a), NUM.findall(b)
    if len(xa) != len(xb):
        return float("inf")
    worst = 0.0
    for u, v in zip(xa, xb):
        u, v = float(u), float(v)
        worst = max(worst, abs(u - v) / max(abs(u), abs(v), 1e-300))
    return worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["check", "freeze", "list"])
    ap.add_argument("--bin", default="build", help="directory holding gkmsvm_kernel / gkmsvm_classify")
    ap.add_argument("--tag", action="append", default=[], help="only cases carrying this tag (repeatable, AND)")
    ap.add_argument("--filter", help="regex on the case name")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--tol", type=float, default=None, help="report max relative numeric diff on failure")
    ap.add_argument("--index", action="store_true", help="freeze: also freeze the <out>.index sidecar")
    a = ap.parse_args()

    cases = load_cases()
    if a.tag:
        cases = [c for c in cases if all(t in c["tags"].split() for t in a.tag)]
    if a.filter:
        cases = [c for c in cases if re.search(a.filter, c["name"])]
    if a.cmd == "list":
        for c in cases:
            print(c["name"], "\t", c["tags"], "\t", " ".join(build_argv(c, a.bin, "<out>")[1:]))
        return
    bindir = os.path.abspath(a.bin)
    for exe in ("gkmsvm_kernel", "gkmsvm_classify"):
        if not os.access(os.path.join(bindir, exe), os.X_OK):
            sys.exit(f"missing binary {os.path.join(bindir, exe)} (run make first)")
    os.makedirs(EXP, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        results = list(ex.map(lambda c: run_case(c, bindir), cases))

    failed = []
    tags = {c["name"]: c["tags"].split() for c in cases}
    for r in results:
        exp = os.path.join(EXP, r["name"] + ".out")
        xfail = [t for t in tags[r["name"]] if t.startswith("xfail")]
        crashed = r["rc"] != 0 or not os.path.exists(r["out"])
        if "expect-error" in tags[r["name"]]:
            # must fail cleanly: a non-zero exit status (not a signal), an ERROR message, no output file
            msgs = open(os.path.join(r["dir"], "stdout.txt"), "rb").read() + open(os.path.join(r["dir"], "stderr.txt"), "rb").read()
            if 0 < r["rc"] < 128 and b"ERROR" in msgs and not os.path.exists(r["out"]):
                print(f"xerror {r['name']:32s} exit {r['rc']}: {msgs.decode(errors='replace').strip().splitlines()[-1][:80]}")
            else:
                failed.append((r, f"expected a clean error exit, got exit {r['rc']}, output file {'present' if os.path.exists(r['out']) else 'absent'}"))
            continue
        if xfail:
            # a known crash (see the tag for the phase that fixes it): passes while it still crashes,
            # and fails loudly once it stops crashing so the case gets frozen and the tag removed
            if crashed:
                print(f"xfail  {r['name']:32s} exit {r['rc']} ({xfail[0]})")
            else:
                failed.append((r, f"XPASS: now runs (exit 0) - freeze it and drop the {xfail[0]} tag"))
            continue
        if crashed:
            tail = b""
            for fn in ("stdout.txt", "stderr.txt"):
                t = open(os.path.join(r["dir"], fn), "rb").read()
                if t.strip():
                    tail += b"\n       [" + fn.encode() + b"] " + b" | ".join(t.strip().splitlines()[-3:])[:300]
            failed.append((r, f"exit code {r['rc']}" + tail.decode(errors="replace")))
            continue
        actual = open(r["out"], "rb").read()
        if a.cmd == "freeze":
            with open(exp, "wb") as f:
                f.write(actual)
            if a.index and os.path.exists(r["out"] + ".index"):
                with open(exp + ".index", "wb") as f:
                    f.write(open(r["out"] + ".index", "rb").read())
            print(f"froze  {r['name']:32s} {len(actual):8d} bytes  {r['secs']:.2f}s")
            continue
        if not os.path.exists(exp):
            failed.append((r, "no expected output (run freeze)"))
            continue
        expected = open(exp, "rb").read()
        if actual != expected and "numeric" in tags[r["name"]] and numeric_equal(actual, expected, 1e-5):
            print(f"ok~    {r['name']:32s} {r['secs']:.2f}s (numeric, rel 1e-5)")
            continue
        if actual != expected:
            msg = "output differs"
            if a.tol is not None:
                msg += f" (max relative numeric diff {max_rel_diff(actual, expected):.3e})"
            failed.append((r, msg))
            continue
        # the .index sidecar (row identity) is compared too whenever one was frozen for this case
        if os.path.exists(exp + ".index") and open(r["out"] + ".index", "rb").read() != open(exp + ".index", "rb").read():
            failed.append((r, "the .index sidecar differs"))
            continue
        print(f"ok     {r['name']:32s} {r['secs']:.2f}s")
    if a.cmd == "freeze":
        print(f"froze {len(results)} cases into {EXP}")
        return
    for r, msg in failed:
        print(f"FAIL   {r['name']:32s} {msg}\n       {' '.join(r['argv'])}")
    print(f"\n{len(results) - len(failed)} passed, {len(failed)} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
