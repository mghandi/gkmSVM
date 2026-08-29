#!/usr/bin/env bash
# DNA kernel wall-clock benchmark (the Phase 6 gate in dev/REFACTORING_PLAN.md).
#
#   dev/bench.sh [bindir]                 -> prints "bench_seconds <best-of-N>" for the reference workload
#   dev/bench.sh --compare BASE NEW [pct] -> exits 1 if NEW is slower than BASE by more than pct % (default 2)
#
# Workload: 250+250 sequences x 500 bp, -l 10 -k 6 -d 3 (the plan's reference), single-threaded so the
# number measures the algorithm rather than the scheduler. Best-of-N (BENCH_RUNS, default 5) to
# suppress noise. Fixtures are generated deterministically into a scratch dir; nothing is downloaded.
set -euo pipefail
cd "$(dirname "$0")/.."
if [ "${1:-}" = "--compare" ]; then
  base=$2; new=$3; pct=${4:-2}
  awk -v b="$base" -v n="$new" -v p="$pct" 'BEGIN{
    r=(n-b)/b*100; printf "baseline %.3fs  new %.3fs  change %+.2f%% (gate %s%%)\n", b, n, r, p;
    exit (r > p) ? 1 : 0 }'
  exit $?
fi
BIN=${1:-build}
RUNS=${BENCH_RUNS:-5}
WORK=${BENCH_WORK:-$(mktemp -d)}
python3 - "$WORK" <<'PY'
import random, sys, os
out = sys.argv[1]; rng = random.Random(0)
def write(fn, n, length, prefix, motif=None):
    with open(os.path.join(out, fn), "w") as f:
        for i in range(n):
            s = [rng.choice("ACGT") for _ in range(length)]
            if motif:
                p = rng.randrange(0, length - len(motif)); s[p:p+len(motif)] = list(motif)
            f.write(f">{prefix}{i}\n{''.join(s)}\n")
write("bench_pos.fa", 250, 500, "bp", motif="TGACTCA")
write("bench_neg.fa", 250, 500, "bn")
PY
best=""
for i in $(seq "$RUNS"); do
  s=$(python3 -c 'import time;print(time.time())')
  "$BIN/gkmsvm_kernel" -l 10 -k 6 -d 3 -T 1 "$WORK/bench_pos.fa" "$WORK/bench_neg.fa" "$WORK/bench_kernel.txt" >/dev/null 2>&1
  e=$(python3 -c "import time;print(time.time()-$s)")
  best=$(awk -v a="$best" -v b="$e" 'BEGIN{ if (a=="" || b<a) print b; else print a }')
done
echo "bench_seconds $best"
