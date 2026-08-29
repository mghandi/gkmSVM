#!/usr/bin/env bash
# Install the package into a scratch library with the Bioconductor Imports trimmed away, so the
# C++ entry points and the R wrappers can be tested on a machine (or CI runner) without the ~1 GB
# BSgenome stack. genNullSeqs is not usable from such an install; everything else is.
#
#   dev/scratch_install.sh [libdir]      (default: $SCRATCH_LIB or a temp dir; prints the lib path)
#
# Nothing in the repository is modified: the package is copied to a temp dir first.
set -euo pipefail
cd "$(dirname "$0")/.."
LIB=${1:-${SCRATCH_LIB:-$(mktemp -d)/lib}}
mkdir -p "$LIB"
TMP=$(mktemp -d)
mkdir -p "$TMP/gkmSVM"
cp -R DESCRIPTION NAMESPACE R src man "$TMP/gkmSVM/"
rm -rf "$TMP/gkmSVM/src/cli" "$TMP/gkmSVM/src/"*.o "$TMP/gkmSVM/src/"*.so
# trim Imports to what the tested code paths use
python3 - "$TMP/gkmSVM/DESCRIPTION" <<'PY'
import re, sys
p = sys.argv[1]; s = open(p).read()
s = re.sub(r"^Imports:.*$", "Imports: Rcpp, kernlab, seqinr, utils, ROCR", s, flags=re.M)
open(p, "w").write(s)
PY
sed -i.bak '/importMethodsFrom(GenomicRanges/d' "$TMP/gkmSVM/NAMESPACE" && rm "$TMP/gkmSVM/NAMESPACE.bak"
R CMD INSTALL --no-docs --no-html --library="$LIB" "$TMP/gkmSVM" >"$TMP/install.log" 2>&1 || { cat "$TMP/install.log"; exit 1; }
echo "$LIB"
