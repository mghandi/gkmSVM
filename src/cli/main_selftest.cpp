/* gkmsvm_selftest : consistency checks that need no input files (not shipped with the R package).
 *
 *   gkmsvm_selftest                 run the checks, exit 1 on the first failure
 *   gkmsvm_selftest coeffs B k      print the general-B coefficient tables (t=0,1,2) for the alphabet
 *                                   vector B ("4,4,4,2,2,2") -- compared with dev/oracle/generalb_gkm.py
 *                                   by tests/oracle/oracle_check.py
 *
 * Check 1: for a single alphabet size (r = 1) the general-B tables (GeneralB.h) must reproduce
 * CCalcWmML's c, cTr and h, which the DNA path keeps using. Tolerance 1e-10 relative to the largest
 * entry of each table (c and cTr are sums with cancellation; CCalcWmML computes c through the
 * impulse response and the six-fold sum, GeneralB through H directly).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "global.h"
#include "CalcWmML.h"
#include "GeneralB.h"

static int nfail = 0;
static void report(const char *what, int L, int K, int b, int m, double a, double g, double scale)
{
	printf("FAIL %s L=%d K=%d b=%d m=%d: CCalcWmML %.17g  GeneralB %.17g  (scale %.3g)\n", what, L, K, b, m, a, g, scale);
	nfail++;
}

static void checkUniform(int L, int K, int b)
{
	CCalcWmML w(L, K, b);
	AlphabetVector av(std::vector<int>(L, b));
	GeneralBTables t(av, K);
	if (av.r != 1 || av.nclasses != L + 1) { printf("FAIL blocks L=%d b=%d: r=%d nclasses=%d\n", L, b, av.r, av.nclasses); nfail++; return; }
	double sc = 0, scTr = 0, sh = 0;
	for (int m = 0; m <= L; m++) { sc = fmax(sc, fabs(w.c[m])); scTr = fmax(scTr, fabs(w.cTr[m])); sh = fmax(sh, fabs(w.h[m])); }
	for (int m = 0; m <= L; m++) {
		if (fabs(w.c[m] - t.H[m]) > 1e-10 * sc) report("c", L, K, b, m, w.c[m], t.H[m], sc);
		if (fabs(w.cTr[m] - t.cTr[m]) > 1e-10 * scTr) report("cTr", L, K, b, m, w.cTr[m], t.cTr[m], scTr);
		if (fabs(w.h[m] - t.h[m]) > 1e-10 * sh) report("h", L, K, b, m, w.h[m], t.h[m], sh);
		if (fabs(w.kernelTruncated[m] - t.g[m]) > 1e-10 * sc) report("g", L, K, b, m, w.kernelTruncated[m], t.g[m], sc);
	}
	int autoTr = 2 * (w.kernelTruncatedLength - 1); if (autoTr > L) autoTr = L;
	// CCalcWmML truncates at the first impulse-response value < 1e-50, GeneralB at H <= 0 exactly (integer
	// arithmetic). When H(m) is exactly 0, CCalcWmML may keep a rounding-noise entry (~1e-17) and its rule
	// then over-estimates the support by 2; the coefficients it keeps beyond GeneralB's bound must be ~0.
	if (t.autoMaxmm(1) > autoTr) { printf("FAIL auto -d (t=1) L=%d K=%d b=%d: CCalcWmML rule %d, GeneralB %d\n", L, K, b, autoTr, t.autoMaxmm(1)); nfail++; }
	for (int m = t.autoMaxmm(1) + 1; m <= autoTr; m++) if (fabs(w.cTr[m]) > 1e-10 * scTr) { printf("FAIL auto -d (t=1) L=%d K=%d b=%d: cTr[%d] = %g beyond GeneralB bound %d\n", L, K, b, m, w.cTr[m], t.autoMaxmm(1)); nfail++; }
	if (t.autoMaxmm(2) != L - K) { printf("FAIL auto -d (t=2) L=%d K=%d b=%d: %d\n", L, K, b, t.autoMaxmm(2)); nfail++; }
}

static void checkIndexing()
{
	// B = (4,2,3,2,4): blocks x=4 {0,4}, x=2 {1,3}, x=3 {2}; 3*3*2 = 18 classes (generalb EXAMPLES.md 2.1)
	int Barr[] = {4, 2, 3, 2, 4};
	AlphabetVector av(std::vector<int>(Barr, Barr + 5));
	if (av.r != 3 || av.nclasses != 18 || av.lengths[0] != 2 || av.lengths[1] != 2 || av.lengths[2] != 1) { printf("FAIL blocks of (4,2,3,2,4): %s\n", av.describe().c_str()); nfail++; }
	for (int idx = 0; idx < av.nclasses; idx++) if (av.classIndex(av.classOf(idx)) != idx) { printf("FAIL class index round trip %d\n", idx); nfail++; }
	int expectStep[] = {1, 3, 9, 3, 1};
	for (int i = 0; i < 5; i++) if (av.step[i] != expectStep[i]) { printf("FAIL step[%d] = %d (expected %d)\n", i, av.step[i], expectStep[i]); nfail++; }
	if ((int)av.reachable(1).size() != 4) { printf("FAIL reachable(1) of (4,2,3,2,4): %d classes (expected 4)\n", (int)av.reachable(1).size()); nfail++; }
	// H idempotence check through the contraction: with g = H the truncated table must equal H when no class is zeroed.
	// (4,4,4,2,2,2), k=4: generalb EXAMPLES.md shows H > 0 on every class? Not necessarily; instead check the row sums of E:
	for (int L = 1; L <= 6; L++) for (int x = 2; x <= 5; x++) for (int m = 0; m <= L; m++) {
		double s = 0; for (int a = 0; a <= L; a++) for (int b = 0; b <= L; b++) s += gkmBlockEnum(L, x, m, a, b);
		if (fabs(s - pow((double)x, L)) > 1e-9 * pow((double)x, L)) { printf("FAIL sum_ab E_{%d,%d}(%d) = %g != %g\n", L, x, m, s, pow((double)x, L)); nfail++; }
	}
}

static int printCoeffs(const char *Bspec, int k)
{
	std::vector<int> B;
	const char *p = Bspec;
	while (*p) { char *e; long v = strtol(p, &e, 10); if (e == p || v < 2) { fprintf(stderr, "bad B: %s\n", Bspec); return 2; } B.push_back((int)v); p = (*e == ',') ? e + 1 : e; }
	AlphabetVector av(B);
	GeneralBTables t(av, k);
	printf("B=%s k=%d %s\n", Bspec, k, av.describe().c_str());
	for (int idx = 0; idx < av.nclasses; idx++)
		printf("m=%s\t|m|=%d\tfilter=%.17g\ttruncated=%.17g\tgkm=%.17g\tg=%.17g\n", av.classLabel(idx).c_str(), av.totalMismatches(idx), t.H[idx], t.cTr[idx], t.h[idx], t.g[idx]);
	printf("autod filter=%d truncated=%d gkm=%d\n", t.autoMaxmm(0), t.autoMaxmm(1), t.autoMaxmm(2));
	return 0;
}

int main(int argc, char **argv)
{
	if (argc == 4 && strcmp(argv[1], "coeffs") == 0) return printCoeffs(argv[2], atoi(argv[3]));
	if (argc != 1) { fprintf(stderr, "usage: %s [coeffs B k]\n", argv[0]); return 2; }
	int bs[] = {2, 4, 5, 20};
	int LK[][2] = {{4, 2}, {6, 4}, {7, 3}, {8, 5}, {10, 6}, {12, 8}, {15, 7}, {5, 5}, {6, 1}};
	int n = 0;
	for (size_t ib = 0; ib < sizeof bs / sizeof *bs; ib++)
		for (size_t i = 0; i < sizeof LK / sizeof *LK; i++) { checkUniform(LK[i][0], LK[i][1], bs[ib]); n++; }
	checkIndexing();
	if (nfail) { printf("gkmsvm_selftest: %d failure(s)\n", nfail); return 1; }
	printf("gkmsvm_selftest: OK (%d single-alphabet table comparisons, block indexing, E row sums)\n", n);
	return 0;
}
