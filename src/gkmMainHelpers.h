/* gkmMainHelpers.h : small helpers shared by the front ends and their per-alphabet tree drivers. */
#pragma once
#include <cstdio>
#include <string>
#include "global.h"
#include "gkmOptions.h"
#include "KernelFile.h"

inline int gkmCannotOpen(const char *fn)
{
	gkmMsg("\n ERROR: cannot open file %s\n", fn);
	return 1;
}

inline int gkmTooManySequences(int limit)
{
	gkmMsg("\n ERROR: more than %d sequences (set -n / maxnumseq to at least the number of sequences).\n", limit);
	return 1;
}

inline void gkmFillGkmkHeader(GkmkWriter &w, const OptsGkmKernel &opt, int maxnmm, int n, int npos, int b, const std::string &alphabet)
{
	w.hdr.n = n; w.hdr.npos = npos;
	w.hdr.L = opt.L; w.hdr.K = opt.K; w.hdr.maxnmm = maxnmm; w.hdr.useTgkm = opt.useTgkm;
	w.hdr.b = b; w.hdr.addRC = opt.addRC ? 1 : 0; w.hdr.usePseudocnt = opt.usePseudocnt ? 1 : 0;
	w.hdr.alphabet = alphabet; // one alphabet: its symbols; multi-track (Phase 7): the -A spec list, b = the largest track alphabet
}
inline void gkmFillGkmkHeader(GkmkWriter &w, const OptsGkmKernel &opt, int maxnmm, int n, int npos, const CConverter &conv)
{
	gkmFillGkmkHeader(w, opt, maxnmm, n, npos, conv.b, std::string(conv.alphabet, conv.alphabet + conv.b));
}
