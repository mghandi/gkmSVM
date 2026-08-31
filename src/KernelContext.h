/* KernelContext.h : per-call state of the kernel and classify computations (Phase 2b).
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * These replace the process-wide globals gMMProfile, gMMProfile0, gMAXMM and gLM1 (and the
 * per-level DFS scratch arrays gDFSlistT/gDFSMMlist), so that two computations can never see each
 * other's state and the trie recursion is re-entrant. Only size-independent types live here; the
 * per-alphabet-size scratch (which holds CLTreeS pointers) is ScoreScratch in LTreef.h.
 */
#pragma once
#include "global.h"
#include <atomic>

struct KernelContext {
	aint ***mmProfile;   // mmProfile[i][m][j], j <= i (lower triangle): number of L-mer pairs of sequences i, j with m mismatches; NULL outside the current tile
	int maxmm;           // maximum number of mismatches recorded (was gMAXMM)
	int LM1;             // L - 1 (was gLM1)
	int rowLo, rowHi;    // the tile: only rows i in [rowLo, rowHi] are recorded (Phase 6: bounded memory)
	int nclasses;        // rows per sequence: mismatch classes (Phase 7). Single alphabet: maxmm+1, row m = m mismatches; general B: prod (l_j+1), row = class index sum m_j stride_j, unreachable rows NULL
	const int *stepPos;  // stepPos[position] added to the class index at a mismatch (all 1 for a single alphabet); the DFS receives it permuted into pass order
	// Kernel-mode accumulation (Plan A1, dev/PERFORMANCE_PLAN.md): when kacc is set the leaf adds the integer-scaled
	// coefficient kcoef[class] straight into kacc[i(i+1)/2 + j] (j <= i) instead of counting the pair in mmProfile;
	// one int64 lower triangle instead of nclasses of them, hence no tiles. kaccAtomic: several threads share this
	// accumulator (relaxed fetch_add); otherwise it is private to one thread (plain load/add/store).
	std::atomic<long long> *kacc;
	const long long *kcoef;
	bool kaccAtomic;
	KernelContext() : mmProfile(NULL), maxmm(0), LM1(0), rowLo(0), rowHi(0x7fffffff), nclasses(0), stepPos(NULL), kacc(NULL), kcoef(NULL), kaccAtomic(false) {}
};

struct ScoreContext {
	myFlt **mmProfile0;  // mmProfile0[m][i]: weighted count of support-vector L-mers at m mismatches from sequence i
	int maxmm;           // bound on the total number of mismatches
	int LM1;
	int nclasses;        // rows of mmProfile0 (see KernelContext)
	const int *step;     // per trie depth (identity order in classify): class-index increment at a mismatch
	ScoreContext() : mmProfile0(NULL), maxmm(0), LM1(0), nclasses(0), step(NULL) {}
};
