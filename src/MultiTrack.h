/* MultiTrack.h : multi-track input for gkm-SVM over heterogeneous alphabets (Phase 7,
 * dev/PHASE7_PLAN.md D1-D3, D6).
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * A record of a multi-track FASTA file (".mfa") is a header line followed by exactly T lines of
 * equal length, one per track:
 *     >seq1
 *     ACGTACGT          track 1 (e.g. DNA)
 *     01100101          track 2 (e.g. methylation flags)
 * A window of length L over the T tracks is an l-mer of length l = T*L, the track-major
 * concatenation of the T track windows, whose position i uses the alphabet of track i / L. A symbol
 * outside its track's alphabet makes every window covering it skipped (the tracks must stay aligned,
 * so such symbols are not dropped as the single-track reader does).
 *
 * The alphabets are given per track as "-A spec[,spec...]": a keyword (dna, rna, protein), the path
 * of an alphabet file (one symbol per line, optional complement), or "=SYMBOLS" for a literal
 * alphabet ("=01", "=NUM"). One spec = the single-track behaviour of earlier versions.
 */
#pragma once
#include <string>
#include <vector>
#include <cstdio>
#include "Converter.h"

class TrackAlphabets {
public:
	std::vector<CConverter *> conv;   // one converter per track (owned)
	std::vector<std::string> specs;   // the specs as given
	TrackAlphabets() {}
	~TrackAlphabets();
	int T() const { return (int)conv.size(); }
	int bmax() const;                 // largest alphabet size over the tracks (selects the trie instantiation)
	// Parse the spec list; 0 on success, 1 on error (a message was printed). maxAlphabet = GKM_MAX_ALPHABET.
	int parse(const std::string &list, int maxAlphabet);
	std::string canonical() const;    // "dna,=01": what the model / kernel files record
	std::vector<int> alphabetVector(int L) const; // B = (b_1,)*L + ... + (b_T,)*L
private:
	TrackAlphabets(const TrackAlphabets &);
	TrackAlphabets &operator=(const TrackAlphabets &);
};

struct MultiTrackRecord {
	std::string name;
	std::vector<std::string> tracks;
	int length() const { return tracks.empty() ? 0 : (int)tracks[0].size(); }
};

// Reads one record: 1 = a record was read, 0 = end of file, -1 = format error (message printed).
// `pending` carries the look-ahead header between calls (start with an empty string).
int readMfaRecord(FILE *f, int T, MultiTrackRecord &rec, std::string &pending, int maxseqlen);

// Encodes every window of length L as a flat l-mer (l = T*L, track-major), skipping windows with a
// symbol outside its track's alphabet; with addRC also the reverse complement of each window (track 1
// complemented and reversed, every other track reversed). Appends to `out`; returns the number of
// l-mers appended (each occupies T*L ints).
int encodeWindows(const MultiTrackRecord &rec, const TrackAlphabets &ta, int L, bool addRC, std::vector<int> &out);
