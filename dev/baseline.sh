#!/usr/bin/env bash
# Baseline measurements for dev/REFACTORING_PLAN.md.
#
# Builds a standalone driver from src/ (a tiny main() that dispatches to mainGkmKernel /
# mainSVMclassify), generates synthetic FASTA fixtures, and reports the numbers quoted in the plan:
#   * DNA kernel wall-clock vs thread count
#   * atomic (MULTI_THREAD_SAFE) vs non-atomic, single-threaded
#   * duplicate sequence names silently merging two records into one kernel row
#   * alphabet size > MAX_ALPHABET_SIZE crashing instead of erroring
#   * heap overflow for FASTA names >= 100 characters (needs ASAN, see --asan)
#   * kernel matrix size and R load time, text vs binary
#
# Nothing here is part of the R package; it only reads src/ and writes to a scratch directory.
# Usage: dev/baseline.sh [outdir]      (default: a temp dir; set --asan for the sanitizer run)
set -euo pipefail
cd "$(dirname "$0")/.."
SRC=$PWD/src
OUT=${1:-$(mktemp -d)}; OUT=${OUT%--asan}
ASAN=""; for a in "$@"; do [ "$a" = "--asan" ] && ASAN="-fsanitize=address -O1 -g"; done
mkdir -p "$OUT"/{obj,obj_noatomic}
echo "== work dir: $OUT"
CXX=${CXX:-clang++}
CXXFLAGS=${CXXFLAGS:-"-std=c++11 -O2"}
[ -n "$ASAN" ] && CXXFLAGS="-std=c++11 $ASAN"

cat > "$OUT/driver.cpp" <<'DRV'
#include <cstring>
#include <cstdio>
int mainGkmKernel(int argc, char** argv);
int mainSVMclassify(int argc, char** argv);
int main(int argc, char** argv){
  if(argc<2){ printf("usage: gkm kernel|classify <args...>\n"); return 1; }
  if(!strcmp(argv[1],"kernel"))   return mainGkmKernel(argc-1, argv+1);
  if(!strcmp(argv[1],"classify")) return mainSVMclassify(argc-1, argv+1);
  printf("unknown mode %s\n", argv[1]); return 1;
}
DRV

build() {  # build <objdir> <exe> [extra include dir first on the path]
  local objdir=$1 exe=$2 inc=${3:-$SRC}
  for f in "$SRC"/*.cpp; do
    case $(basename "$f") in RcppExports.cpp|gkmsvm_*.cpp) continue;; esac
    $CXX $CXXFLAGS -I"$inc" -I"$SRC" -c "$f" -o "$objdir/$(basename "${f%.cpp}").o" 2>/dev/null \
      || { echo "FAILED to compile $f"; return 1; }
  done
  $CXX $CXXFLAGS "$OUT/driver.cpp" "$objdir"/*.o -o "$exe"
}

echo "== building (this compiles src/ with a main() shim; ~30 translation units)"
build "$OUT/obj" "$OUT/gkm"
# second build with MULTI_THREAD_SAFE disabled, to measure the cost of the atomic counters
mkdir -p "$OUT/inc_noatomic"
sed 's|^#define MULTI_THREAD_SAFE|//#define MULTI_THREAD_SAFE|' "$SRC/global.h" > "$OUT/inc_noatomic/global.h"
build "$OUT/obj_noatomic" "$OUT/gkm_noatomic" "$OUT/inc_noatomic"

echo "== generating fixtures"
python3 - "$OUT" <<'PY'
import random, sys, os
out = sys.argv[1]; rng = random.Random(0)
def seq(n, alpha="ACGT"): return "".join(rng.choice(alpha) for _ in range(n))
def write(fn, n, length, names, motif=None, alpha="ACGT"):
    with open(os.path.join(out, fn), "w") as f:
        for i in range(n):
            s = list(seq(length, alpha))
            if motif:
                p = rng.randrange(0, length - len(motif)); s[p:p+len(motif)] = list(motif)
            f.write(f">{names(i)}\n{''.join(s)}\n")
write("pos.fa",     60, 300, lambda i: f"pos{i}", motif="TGACTCA")
write("neg.fa",     60, 300, lambda i: f"neg{i}")
write("big_pos.fa",250, 500, lambda i: f"bp{i}")
write("big_neg.fa",250, 500, lambda i: f"bn{i}")
# one repeated name inside the positive file: records 0 and 1 share "posX"
write("pos_dup.fa", 60, 300, lambda i: "posX" if i < 2 else f"pos{i}")
# names longer than the 100-char buffer in mainGkmKernel
write("pos_longname.fa", 20, 300, lambda i: "N"*120 + f"_{i}")
write("prot_pos.fa", 30, 200, lambda i: f"pp{i}", alpha="ACDEFGHIKLMNPQRSTVWY")
write("prot_neg.fa", 30, 200, lambda i: f"pn{i}", alpha="ACDEFGHIKLMNPQRSTVWY")
open(os.path.join(out,"alphabet_protein.txt"),"w").write("\n".join("ACDEFGHIKLMNPQRSTVWY")+"\n")
open(os.path.join(out,"alphabet_b5.txt"),"w").write("A\nC\nG\nT\nN\n")
PY

cd "$OUT"

asan_run() {  # asan_run <posfile>: run once under ASAN, print only the sanitizer verdict
  # the child aborts on the first finding; run it in a nested shell whose job report we discard
  bash -c "ASAN_OPTIONS=detect_leaks=0 '$OUT/gkm' kernel -l 10 -k 6 -d 3 -T 1 '$1' neg.fa /dev/null > '$OUT/asan.log' 2>&1" 2>/dev/null || true
  { grep -m1 "ERROR: AddressSanitizer" "$OUT/asan.log"; grep -m1 "    #0 " "$OUT/asan.log"; } \
      2>/dev/null | sed 's/^/      /'
  grep -q "AddressSanitizer" "$OUT/asan.log" || echo "      (clean)"
}

if [ -n "$ASAN" ]; then
  # Timing sections are meaningless (and abort) under the sanitizer: the CalcWmML overread below
  # is fatal with ASAN, which is the point of this mode.
  echo
  echo "== AddressSanitizer report (timing sections skipped)"
  echo "   NOTE: ASAN aborts on the FIRST finding. CalcWmML.cpp:150 fires on every run and masks"
  echo "         the others; fix it first (Phase 1) to surface the name-buffer overflow below."
  echo "   -- plain DNA run, default -l 10 -k 6:"
  asan_run pos.fa
  echo "   -- FASTA names >= 100 characters:"
  asan_run pos_longname.fa
  echo
  echo "== done (sanitizer mode). artefacts in $OUT"
  exit 0
fi

echo
echo "== 1. DNA kernel, 250+250 seqs x 500 bp, -l 10 -k 6 -d 3: wall clock vs threads"
for T in 1 4 20; do
  /usr/bin/time -p ./gkm kernel -l 10 -k 6 -d 3 -T $T big_pos.fa big_neg.fa /dev/null >/dev/null 2>t.$T
  printf "   -T %-3s real %ss  user %ss\n" $T "$(awk '/real/{print $2}' t.$T)" "$(awk '/user/{print $2}' t.$T)"
done

echo
echo "== 2. atomic vs non-atomic counters, single thread"
/usr/bin/time -p ./gkm          kernel -l 10 -k 6 -d 3 -T 1 big_pos.fa big_neg.fa ka.txt >/dev/null 2>ta
/usr/bin/time -p ./gkm_noatomic kernel -l 10 -k 6 -d 3 -T 1 big_pos.fa big_neg.fa kn.txt >/dev/null 2>tn
printf "   MULTI_THREAD_SAFE on : real %ss\n" "$(awk '/real/{print $2}' ta)"
printf "   MULTI_THREAD_SAFE off: real %ss\n" "$(awk '/real/{print $2}' tn)"
cmp -s ka.txt kn.txt && echo "   (results identical single-threaded)"

echo
echo "== 3. duplicate sequence names inside one file (records 0 and 1 share a name)"
echo "   60 records in pos_dup.fa; rows actually produced:"
./gkm kernel -l 10 -k 6 -d 3 pos_dup.fa neg.fa k_dup.txt 2>&1 | grep -E "npos|nneg" | tr -d '\n'; echo
echo "   -> two records were merged into one kernel row, silently"

echo
echo "== 4. alphabet larger than MAX_ALPHABET_SIZE (compiled as 4)"
for A in alphabet_b5.txt alphabet_protein.txt; do
  # run in a child shell so its crash report does not reach this script's stderr
  set +e; bash -c "'$OUT/gkm' kernel -l 6 -k 4 -d 2 -A '$A' prot_pos.fa prot_neg.fa 'k_$A.txt' >/dev/null 2>&1" 2>/dev/null; rc=$?; set -e
  printf "   -A %-22s exit=%s %s\n" "$A" "$rc" "$([ $rc -eq 139 ] && echo '(SIGSEGV)')"
done

echo
echo "== 5. (re-run with --asan for the AddressSanitizer findings)"

echo
echo "== 6. kernel matrix on disk and in R: text vs binary (n = 5000, synthetic values)"
python3 - <<'PY'
import random, struct
n = 5000; rng = random.Random(0)
with open("synth_kernel.txt", "w") as f:
    for i in range(n):
        f.write("\t".join(["%e" % rng.random() for _ in range(i)] + ["1.0"]) + "\t\n")
with open("synth_kernel.bin", "wb") as f:
    f.write(struct.pack("<ii", n, 0))
    for i in range(n):
        f.write(struct.pack("<%df" % (i+1), *([rng.random() for _ in range(i)] + [1.0])))
PY
ls -l synth_kernel.txt synth_kernel.bin | awk '{printf "   %-20s %8.1f MB\n", $9, $5/1048576}'
command -v Rscript >/dev/null && Rscript -e '
n <- 5000
t1 <- system.time(m <- data.matrix(utils::read.table("synth_kernel.txt", fill=TRUE, col.names=paste("V",1:n))))
t2 <- system.time({f <- file("synth_kernel.bin","rb"); readBin(f,"integer",2); readBin(f,"numeric",n*(n+1)/2,size=4); close(f)})
cat(sprintf("   read.table (today)  %6.2f s\n   readBin  (proposed) %6.2f s\n   speedup             %6.1fx\n",
            t1[["elapsed"]], t2[["elapsed"]], t1[["elapsed"]]/max(t2[["elapsed"]],1e-3)))' 2>/dev/null \
  || echo "   (Rscript not found; skipping the R load-time comparison)"

echo
echo "== done. artefacts in $OUT"
