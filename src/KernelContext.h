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

struct KernelContext {
	aint ***mmProfile;   // mmProfile[i][m][j], j <= i (lower triangle): number of L-mer pairs of sequences i, j with m mismatches
	int maxmm;           // maximum number of mismatches recorded (was gMAXMM)
	int LM1;             // L - 1 (was gLM1)
	KernelContext() : mmProfile(NULL), maxmm(0), LM1(0) {}
};

struct ScoreContext {
	myFlt **mmProfile0;  // mmProfile0[m][i]: weighted count of support-vector L-mers at m mismatches from sequence i
	int maxmm;
	int LM1;
	ScoreContext() : mmProfile0(NULL), maxmm(0), LM1(0) {}
};
