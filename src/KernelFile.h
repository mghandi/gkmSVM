/* KernelFile.h : the .gkmk v1 binary kernel format (Phase 4 of dev/REFACTORING_PLAN.md).
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * All integers little-endian. Layout (offsets in bytes):
 *   0   "GKMK"                       magic
 *   4   uint8  version   (1)
 *   5   uint8  dtype     0 = float32, 1 = float64
 *   6   uint8  layout    0 = lower triangle including the diagonal, row-major (the text format's layout)
 *   7   uint8  flags     bit0 = names table present
 *   8   int32  n
 *   12  int32  npos
 *   16  int32  L, K, maxnmm, useTgkm, b, addRC, usePseudocnt          (7 x int32: provenance)
 *   44  uint32 alphabet length, then that many bytes
 *       names table (if flags&1), n entries: uint32 idlen, id, uint32 namelen, name, int8 label,
 *                                            int64 length, int64 nlmers
 *       payload: n(n+1)/2 values of dtype
 *       uint32 crc32 of the payload bytes
 * The R reader is R/gkm_kernel_io.R (read_gkm_kernel / write_gkm_kernel).
 */
#pragma once
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <cstring>
#include "SequenceSet.h"
#include "crc32.h"

struct GkmkHeader {
	int n = 0, npos = 0;
	int L = 0, K = 0, maxnmm = 0, useTgkm = 0, b = 4, addRC = 1, usePseudocnt = 0;
	std::string alphabet;
};

namespace gkmk_detail {
inline void put8(std::vector<unsigned char> &o, uint8_t v) { o.push_back(v); }
inline void put32(std::vector<unsigned char> &o, uint32_t v) { for (int i = 0; i < 4; i++) o.push_back((v >> (8 * i)) & 0xFF); }
inline void put64(std::vector<unsigned char> &o, uint64_t v) { for (int i = 0; i < 8; i++) o.push_back((v >> (8 * i)) & 0xFF); }
inline void putstr(std::vector<unsigned char> &o, const std::string &s) { put32(o, (uint32_t)s.size()); o.insert(o.end(), s.begin(), s.end()); }
}

// Collects the lower triangle row by row (the same order the text writer prints it) and writes the
// file at the end. `values` holds n(n+1)/2 doubles; dtype 0 stores them as float32.
class GkmkWriter {
public:
	GkmkHeader hdr;
	std::vector<double> values;
	int dtype = 0;
	void add(double v) { values.push_back(v); }
	int write(const std::string &fn, const std::vector<SeqRecord> &recs) const {
		using namespace gkmk_detail;
		std::vector<unsigned char> h;
		h.insert(h.end(), {'G', 'K', 'M', 'K'});
		put8(h, 1); put8(h, (uint8_t)dtype); put8(h, 0); put8(h, recs.empty() ? 0 : 1);
		put32(h, (uint32_t)hdr.n); put32(h, (uint32_t)hdr.npos);
		for (int v : {hdr.L, hdr.K, hdr.maxnmm, hdr.useTgkm, hdr.b, hdr.addRC, hdr.usePseudocnt}) put32(h, (uint32_t)v);
		putstr(h, hdr.alphabet);
		for (const SeqRecord &r : recs) {
			putstr(h, r.id); putstr(h, r.name); put8(h, (uint8_t)(int8_t)r.label);
			put64(h, (uint64_t)r.length); put64(h, (uint64_t)r.nlmers);
		}
		std::vector<unsigned char> payload;
		payload.reserve(values.size() * (dtype ? 8 : 4));
		for (double d : values) {
			if (dtype == 0) { float f = (float)d; uint32_t u; memcpy(&u, &f, 4); put32(payload, u); }
			else { uint64_t u; memcpy(&u, &d, 8); put64(payload, u); }
		}
		FILE *f = fopen(fn.c_str(), "wb");
		if (f == NULL) return 1;
		uint32_t crc = gkm_crc32(payload.data(), payload.size());
		std::vector<unsigned char> tail; put32(tail, crc);
		bool ok = fwrite(h.data(), 1, h.size(), f) == h.size()
		       && fwrite(payload.data(), 1, payload.size(), f) == payload.size()
		       && fwrite(tail.data(), 1, tail.size(), f) == tail.size();
		fclose(f);
		return ok ? 0 : 1;
	}
};

namespace gkmk_detail {
inline bool get8(FILE *f, uint8_t &v) { int c = fgetc(f); if (c == EOF) return false; v = (uint8_t)c; return true; }
inline bool get32(FILE *f, uint32_t &v) { unsigned char b[4]; if (fread(b, 1, 4, f) != 4) return false; v = (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24); return true; }
inline bool get64(FILE *f, uint64_t &v) { uint32_t lo, hi; if (!get32(f, lo) || !get32(f, hi)) return false; v = ((uint64_t)hi << 32) | lo; return true; }
inline bool getstr(FILE *f, std::string &s) { uint32_t n; if (!get32(f, n)) return false; s.resize(n); return n == 0 || fread(&s[0], 1, n, f) == n; }
}

// Reads a .gkmk file whose 4 magic bytes have already been consumed (see readKernel in
// mainSVMtrain.cpp). Fills hdr, records and values (row-major lower triangle as doubles).
class GkmkReader {
public:
	GkmkHeader hdr;
	std::vector<SeqRecord> records;
	std::vector<double> values;
	std::string error;
	int read(FILE *f) {
		using namespace gkmk_detail;
		uint8_t version, dtype, layout, flags;
		if (!get8(f, version) || !get8(f, dtype) || !get8(f, layout) || !get8(f, flags)) return fail("truncated header");
		if (version != 1) return fail("unsupported .gkmk version");
		if (layout != 0) return fail("unsupported layout");
		uint32_t u;
		if (!get32(f, u)) return fail("truncated header"); hdr.n = (int)u;
		if (!get32(f, u)) return fail("truncated header"); hdr.npos = (int)u;
		int *prov[7] = {&hdr.L, &hdr.K, &hdr.maxnmm, &hdr.useTgkm, &hdr.b, &hdr.addRC, &hdr.usePseudocnt};
		for (int i = 0; i < 7; i++) { if (!get32(f, u)) return fail("truncated header"); *prov[i] = (int)u; }
		if (!getstr(f, hdr.alphabet)) return fail("truncated header");
		if (flags & 1) {
			records.resize(hdr.n);
			for (int i = 0; i < hdr.n; i++) {
				SeqRecord &r = records[i]; r.index = i; uint8_t lab; uint64_t v;
				if (!getstr(f, r.id) || !getstr(f, r.name) || !get8(f, lab) || !get64(f, v)) return fail("truncated names table");
				r.label = (int8_t)lab; r.length = (long)v;
				if (!get64(f, v)) return fail("truncated names table"); r.nlmers = (long)v;
			}
		}
		size_t m = (size_t)hdr.n * (hdr.n + 1) / 2, size = dtype == 0 ? 4 : 8;
		std::vector<unsigned char> payload(m * size);
		if (m && fread(payload.data(), 1, payload.size(), f) != payload.size()) return fail("truncated payload");
		uint32_t crc; if (!get32(f, crc)) return fail("missing checksum");
		if (crc != gkm_crc32(payload.data(), payload.size())) return fail("payload checksum mismatch (file corrupted or truncated)");
		values.resize(m);
		for (size_t k = 0; k < m; k++) {
			if (dtype == 0) { uint32_t b = (uint32_t)payload[4*k] | ((uint32_t)payload[4*k+1] << 8) | ((uint32_t)payload[4*k+2] << 16) | ((uint32_t)payload[4*k+3] << 24); float x; memcpy(&x, &b, 4); values[k] = x; }
			else { uint64_t b = 0; for (int i = 7; i >= 0; i--) b = (b << 8) | payload[8*k+i]; double x; memcpy(&x, &b, 8); values[k] = x; }
		}
		return 0;
	}
private:
	int fail(const char *msg) { error = msg; return 1; }
};
