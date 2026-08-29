/* GeneralB.h : gapped k-mer filter / kernel coefficient tables for an ARBITRARY alphabet vector
 * B = (b_1, ..., b_l) -- position i of an l-mer uses an alphabet of size b_i (Phase 7 of
 * dev/REFACTORING_PLAN.md; design in dev/PHASE7_PLAN.md).
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Theory (Mohammad-Noori, Ghareghani, Ghandi, "Generalized Gapped k-mer Filters for Robust
 * Frequency Estimation"; exact reference implementation dev/oracle/generalb_gkm.py):
 *
 *   H(u,w) = (1 / prod_i b_i) * sum_{j=0}^{k} e_j(w_1..w_l),   w_i = b_i - 1 (match) or -1 (mismatch),
 *
 * with e_j the elementary symmetric polynomials. H depends only on the number of mismatches inside
 * each BLOCK (positions of equal alphabet size), so a pair of l-mers has a mismatch class
 * m = (m_1..m_r), 0 <= m_j <= l_j, and there are prod_j (l_j + 1) classes. The kernel between two
 * sequences is sum_m N(m) c(m) with the class profile N(m) and
 *   t=0 (filter)           c(m) = H(m)
 *   t=1 (truncated filter) c(m) = < g, (E_1(m_1) (x) ... (x) E_r(m_r)) g >, g = H with the dominance rule
 *   t=2 (gkm counts)       c(m) = C(l - |m|, k)
 * For r = 1 these are exactly CCalcWmML's c, cTr and h (tested by src/cli/main_selftest.cpp).
 *
 * Classes are addressed by a single integer index sum_j m_j * stride_j, stride_0 = 1,
 * stride_j = prod_{j'<j} (l_j' + 1): the trie DFS adds step[position] = stride_{block(position)} at a
 * mismatch instead of 1, and the profile row is the class index (for r = 1 every step is 1 and the
 * index is the mismatch count, i.e. today's behaviour).
 */
#pragma once
#include <vector>
#include <string>

struct AlphabetVector {
	std::vector<int> B;                       // b_i per position; l = B.size()
	std::vector<int> sizes;                   // x_j: alphabet size of block j (order of first appearance)
	std::vector<std::vector<int> > positions; // positions of block j
	std::vector<int> lengths;                 // l_j
	std::vector<int> dims;                    // l_j + 1
	std::vector<int> strides;                 // stride_j
	std::vector<int> blockOf;                 // block of each position
	std::vector<int> step;                    // strides[blockOf[pos]] for each position
	int r = 0, l = 0, nclasses = 0;

	AlphabetVector() {}
	explicit AlphabetVector(const std::vector<int> &B);
	int classIndex(const std::vector<int> &m) const;
	std::vector<int> classOf(int idx) const;      // m_1..m_r
	int totalMismatches(int idx) const;           // |m| = sum_j m_j
	bool uniform() const { return r == 1; }
	std::string classLabel(int idx) const;        // "m1,m2,..."
	std::string describe() const;                 // blocks, sizes, class count
	// The classes with |m| <= maxmm (the rows of the profile that the DFS can reach when the total
	// number of mismatches is bounded by maxmm); in increasing index order.
	std::vector<int> reachable(int maxmm) const;
};

struct GeneralBTables {
	const AlphabetVector &av;
	int k;
	std::vector<double> H;    // t=0: c(m) = H(m)              (indexed by class index)
	std::vector<double> g;    // truncated impulse response: H with the dominance rule
	std::vector<double> cTr;  // t=1: <g, (E_1(m_1) (x) ... (x) E_r(m_r)) g>
	std::vector<double> h;    // t=2: C(l - |m|, k)
	int maxTotalTruncated;    // largest |m| with cTr(m) != 0 (the automatic -d for t=1)

	GeneralBTables(const AlphabetVector &av, int k);
	// The table for a filter type (0, 1, 2); NULL for the types that are not generalised (3, 4).
	const double *table(int useTgkm) const;
	// Automatic mismatch bound (-d -1) for a filter type: l (t=0), maxTotalTruncated (t=1), l-k (t=2).
	int autoMaxmm(int useTgkm) const;
};

// E_{L,x}(m; a, b): number of words of a block of length L over x symbols with a mismatches to the
// first and b mismatches to the second of two block-words that differ in m positions
// (PLoS Comput Biol 2014 eq. 14; twoblock/README.md).
double gkmBlockEnum(int L, int x, int m, int a, int b);
