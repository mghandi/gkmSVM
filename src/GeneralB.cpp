/* GeneralB.cpp : see GeneralB.h.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */
#include "GeneralB.h"
#include "global.h"
#include <math.h>
#include <stdio.h>

AlphabetVector::AlphabetVector(const std::vector<int> &B_) : B(B_)
{
	l = (int)B.size();
	blockOf.assign(l, -1);
	for (int i = 0; i < l; i++) {
		int j = -1;
		for (size_t s = 0; s < sizes.size(); s++) if (sizes[s] == B[i]) { j = (int)s; break; }
		if (j < 0) { j = (int)sizes.size(); sizes.push_back(B[i]); positions.push_back(std::vector<int>()); }
		positions[j].push_back(i);
		blockOf[i] = j;
	}
	r = (int)sizes.size();
	lengths.resize(r); dims.resize(r); strides.resize(r);
	nclasses = 1;
	for (int j = 0; j < r; j++) {
		lengths[j] = (int)positions[j].size();
		dims[j] = lengths[j] + 1;
		strides[j] = nclasses;
		nclasses *= dims[j];
	}
	step.resize(l);
	for (int i = 0; i < l; i++) step[i] = strides[blockOf[i]];
}

int AlphabetVector::classIndex(const std::vector<int> &m) const
{
	int idx = 0;
	for (int j = 0; j < r; j++) idx += m[j] * strides[j];
	return idx;
}

std::vector<int> AlphabetVector::classOf(int idx) const
{
	std::vector<int> m(r);
	for (int j = 0; j < r; j++) m[j] = (idx / strides[j]) % dims[j];
	return m;
}

int AlphabetVector::totalMismatches(int idx) const
{
	int t = 0;
	for (int j = 0; j < r; j++) t += (idx / strides[j]) % dims[j];
	return t;
}

std::string AlphabetVector::classLabel(int idx) const
{
	std::vector<int> m = classOf(idx);
	std::string s;
	char buf[16];
	for (int j = 0; j < r; j++) { snprintf(buf, sizeof buf, "%s%d", j ? "," : "", m[j]); s += buf; }
	return s;
}

std::string AlphabetVector::describe() const
{
	std::string s;
	char buf[64];
	for (int j = 0; j < r; j++) {
		snprintf(buf, sizeof buf, "%sblock %d: alphabet size %d, %d position%s", j ? "; " : "", j + 1, sizes[j], lengths[j], lengths[j] == 1 ? "" : "s");
		s += buf;
	}
	snprintf(buf, sizeof buf, " (%d mismatch classes)", nclasses);
	s += buf;
	return s;
}

std::vector<int> AlphabetVector::reachable(int maxmm) const
{
	std::vector<int> out;
	for (int idx = 0; idx < nclasses; idx++) if (totalMismatches(idx) <= maxmm) out.push_back(idx);
	return out;
}

// ---------------------------------------------------------------------------------------------- H
#ifdef __SIZEOF_INT128__
__extension__ typedef __int128 gkmWide;   // e_j sums are integers: exact up to ~1.7e38 (__extension__: no -pedantic warning)
#else
typedef long double gkmWide;
#endif

// H(m) = (e_0 + ... + e_k)(w) / prod b_i for the representative pair of class m: mismatches at the
// first m_j positions of block j (any representative gives the same value).
static double hOfClass(const AlphabetVector &av, int k, const std::vector<int> &m)
{
	std::vector<gkmWide> e(k + 1, 0);
	e[0] = 1;
	std::vector<int> used(av.r, 0);
	for (int i = 0; i < av.l; i++) {
		int j = av.blockOf[i];
		int w = (used[j] < m[j]) ? -1 : av.B[i] - 1;   // mismatch : match
		used[j]++;
		for (int q = k; q >= 1; q--) e[q] += (gkmWide)w * e[q - 1];
	}
	gkmWide num = 0;
	for (int q = 0; q <= k; q++) num += e[q];
	long double den = 1.0L;
	for (int i = 0; i < av.l; i++) den *= av.B[i];
	return (double)((long double)num / den);
}

double gkmBlockEnum(int L, int x, int m, int a, int b)
{
	double tot = 0;
	int tlo = (a - m > 0) ? a - m : 0;
	int thi = (a < L - m) ? a : L - m;
	for (int t = tlo; t <= thi; t++) {
		int rr = a + b - 2 * t - m;
		if (rr < 0 || rr > a - t) continue;
		// (x-2)^rr with 0^0 = 1 (binary blocks): pow(0,0) is 1 in C
		tot += dCombinations(L - m, t) * pow((double)(x - 1), (double)t) * dCombinations(m, a - t) * dCombinations(a - t, rr) * pow((double)(x - 2), (double)rr);
	}
	return tot;
}

GeneralBTables::GeneralBTables(const AlphabetVector &av_, int k_) : av(av_), k(k_)
{
	const int C = av.nclasses;
	H.assign(C, 0.0); g.assign(C, 0.0); cTr.assign(C, 0.0); h.assign(C, 0.0);
	std::vector<std::vector<int> > cls(C);
	for (int idx = 0; idx < C; idx++) {
		cls[idx] = av.classOf(idx);
		H[idx] = hOfClass(av, k, cls[idx]);
		int rem = av.l - av.totalMismatches(idx);
		h[idx] = (rem >= k) ? dCombinations(rem, k) : 0.0;
	}
	// dominance rule: H(m) <= 0 zeroes every class a >= m componentwise
	g = H;
	for (int idx = 0; idx < C; idx++) {
		if (H[idx] > 0) continue;
		for (int a = 0; a < C; a++) {
			bool dom = true;
			for (int j = 0; j < av.r && dom; j++) if (cls[a][j] < cls[idx][j]) dom = false;
			if (dom) g[a] = 0.0;
		}
	}
	// E[j][m][a][b]
	std::vector<std::vector<std::vector<std::vector<double> > > > E(av.r);
	for (int j = 0; j < av.r; j++) {
		int Lj = av.lengths[j], xj = av.sizes[j];
		E[j].resize(Lj + 1);
		for (int m = 0; m <= Lj; m++) {
			E[j][m].assign(Lj + 1, std::vector<double>(Lj + 1, 0.0));
			for (int a = 0; a <= Lj; a++) for (int b = 0; b <= Lj; b++) E[j][m][a][b] = gkmBlockEnum(Lj, xj, m, a, b);
		}
	}
	// cTr(m) = < g, (E_1(m_1) (x) ... (x) E_r(m_r)) g >, one mode at a time
	std::vector<double> G(C), Gn(C);
	maxTotalTruncated = 0;
	for (int idx = 0; idx < C; idx++) {
		G = g;
		for (int j = 0; j < av.r; j++) {
			const std::vector<std::vector<double> > &Ej = E[j][cls[idx][j]];
			int dj = av.dims[j], sj = av.strides[j];
			for (int q = 0; q < C; q++) {
				int aj = (q / sj) % dj;
				int base = q - aj * sj;
				double s = 0;
				for (int b = 0; b < dj; b++) s += Ej[aj][b] * G[base + b * sj];
				Gn[q] = s;
			}
			G.swap(Gn);
		}
		double s = 0;
		for (int q = 0; q < C; q++) s += g[q] * G[q];
		cTr[idx] = s;
		if (s != 0.0) { int t = av.totalMismatches(idx); if (t > maxTotalTruncated) maxTotalTruncated = t; }
	}
}

const double *GeneralBTables::table(int useTgkm) const
{
	if (useTgkm == 0) return H.data();
	if (useTgkm == 1) return cTr.data();
	if (useTgkm == 2) return h.data();
	return NULL;
}

int GeneralBTables::autoMaxmm(int useTgkm) const
{
	if (useTgkm == 1) return maxTotalTruncated;
	if (useTgkm == 2) return av.l - k;
	return av.l;
}
