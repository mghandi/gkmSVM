/* MultiTrack.cpp : see MultiTrack.h.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */
#include "MultiTrack.h"
#include "global.h"
#include <string.h>
#include <ctype.h>

TrackAlphabets::~TrackAlphabets()
{
	for (size_t t = 0; t < conv.size(); t++) delete conv[t];
}

int TrackAlphabets::bmax() const
{
	int b = 0;
	for (size_t t = 0; t < conv.size(); t++) if (conv[t]->b > b) b = conv[t]->b;
	return b;
}

int TrackAlphabets::parse(const std::string &list, int maxAlphabet)
{
	for (size_t t = 0; t < conv.size(); t++) delete conv[t];
	conv.clear(); specs.clear();
	size_t start = 0;
	while (start <= list.size()) {
		size_t comma = list.find(',', start);
		std::string spec = list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
		start = (comma == std::string::npos) ? list.size() + 1 : comma + 1;
		if (spec.empty()) {
			if (list.empty() && conv.empty()) { conv.push_back(new CConverter()); specs.push_back("dna"); break; } // no -A: the built-in DNA alphabet, silently
			gkmMsg("\n ERROR: empty alphabet specification in '%s'.\n", list.c_str()); return 1;
		}
		CConverter *cv = new CConverter();
		conv.push_back(cv); specs.push_back(spec);
		if (spec[0] == '=') {
			std::string sym = spec.substr(1);
			if ((int)sym.size() < 2) { gkmMsg("\n ERROR: alphabet '%s' must have at least two symbols.\n", spec.c_str()); return 1; }
			if ((int)sym.size() > maxAlphabet) { gkmMsg("\n ERROR: alphabet '%s' has more than %d symbols.\n", spec.c_str(), maxAlphabet); return 1; }
			for (size_t i = 0; i < sym.size(); i++) {
				if (isspace((unsigned char)sym[i]) || sym[i] == '>' || sym[i] == '#' || sym[i] == ';') { gkmMsg("\n ERROR: symbol '%c' is not allowed in alphabet '%s'.\n", sym[i], spec.c_str()); return 1; }
				for (size_t j = 0; j < i; j++) if (toupper((unsigned char)sym[i]) == toupper((unsigned char)sym[j])) { gkmMsg("\n ERROR: symbol '%c' appears twice in alphabet '%s'.\n", sym[i], spec.c_str()); return 1; }
			}
			if (cv->setAlphabet(sym.c_str(), (int)sym.size(), NULL) != 0) return 1;
			gkmMsg("Alphabet: %s (%d symbols)\n", spec.c_str(), cv->b);
		} else {
			if (cv->readAlphabetFile(spec.c_str(), maxAlphabet) != 0) return 1;
		}
		if (list.empty()) break;
	}
	return 0;
}

std::string TrackAlphabets::canonical() const
{
	std::string s;
	for (size_t t = 0; t < specs.size(); t++) { if (t) s += ","; s += specs[t]; }
	return s;
}

std::vector<int> TrackAlphabets::alphabetVector(int L) const
{
	std::vector<int> B;
	for (size_t t = 0; t < conv.size(); t++) for (int i = 0; i < L; i++) B.push_back(conv[t]->b);
	return B;
}

// Reads a whole line (any length) without the newline; false at end of file with nothing read.
static bool readLine(FILE *f, std::string &line)
{
	line.clear();
	char buf[4096];
	bool any = false;
	while (fgets(buf, sizeof buf, f) != NULL) {
		any = true;
		size_t n = strlen(buf);
		if (n && buf[n - 1] == '\n') { line.append(buf, n - 1); break; }
		line.append(buf, n);
	}
	while (!line.empty() && (line[line.size() - 1] == '\r' || line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t')) line.erase(line.size() - 1);
	return any;
}

int readMfaRecord(FILE *f, int T, MultiTrackRecord &rec, std::string &pending, int maxseqlen)
{
	rec.name.clear(); rec.tracks.clear();
	std::string line;
	if (pending.empty()) {
		// find the next header
		for (;;) {
			if (!readLine(f, line)) return 0;
			if (line.empty() || line[0] == '#' || line[0] == ';') continue;
			if (line[0] == '>') break;
			gkmMsg("\n ERROR: expected a '>' header line in the multi-track FASTA file, found '%.40s'.\n", line.c_str());
			return -1;
		}
		pending = line;
	}
	{
		size_t i = 1; while (i < pending.size() && isspace((unsigned char)pending[i])) i++;
		size_t j = i; while (j < pending.size() && !isspace((unsigned char)pending[j])) j++;
		rec.name = pending.substr(i, j - i);
		if (rec.name.empty()) rec.name = "NA";
	}
	pending.clear();
	while ((int)rec.tracks.size() < T) {
		if (!readLine(f, line)) break;
		if (line.empty() || line[0] == '#' || line[0] == ';') continue;
		if (line[0] == '>') { pending = line; break; }
		rec.tracks.push_back(line);
	}
	if ((int)rec.tracks.size() != T) {
		gkmMsg("\n ERROR: record '%s' has %d track line(s); %d alphabets were given (-A), so every record needs exactly %d lines after its header.\n", rec.name.c_str(), (int)rec.tracks.size(), T, T);
		return -1;
	}
	for (int t = 1; t < T; t++) if (rec.tracks[t].size() != rec.tracks[0].size()) {
		gkmMsg("\n ERROR: record '%s': track %d has length %d but track 1 has length %d (tracks must be aligned).\n", rec.name.c_str(), t + 1, (int)rec.tracks[t].size(), (int)rec.tracks[0].size());
		return -1;
	}
	if (rec.length() > maxseqlen) {
		gkmMsg("\n ERROR: sequence '%s' is longer than the maximum sequence length %d (set -m / maxseqlen to at least its length).\n", rec.name.c_str(), maxseqlen);
		return -1;
	}
	return 1;
}

int encodeWindows(const MultiTrackRecord &rec, const TrackAlphabets &ta, int L, bool addRC, std::vector<int> &out)
{
	const int T = ta.T(), n = rec.length();
	if (n < L) return 0;
	// per-track codes, -1 outside the alphabet; bad[p] = number of invalid positions in [0, p)
	std::vector<std::vector<int> > codes(T, std::vector<int>(n));
	std::vector<int> bad(n + 1, 0);
	for (int p = 0; p < n; p++) {
		int invalid = 0;
		for (int t = 0; t < T; t++) {
			unsigned char ch = (unsigned char)rec.tracks[t][p];
			if (ta.conv[t]->isInAlphabet[ch]) codes[t][p] = ta.conv[t]->cidx[ch];
			else { codes[t][p] = -1; invalid = 1; }
		}
		bad[p + 1] = bad[p] + invalid;
	}
	int added = 0;
	const CConverter &c0 = *ta.conv[0];
	for (int p = 0; p + L <= n; p++) {
		if (bad[p + L] - bad[p] != 0) continue;
		for (int t = 0; t < T; t++) for (int i = 0; i < L; i++) out.push_back(codes[t][p + i]);
		added++;
		if (addRC) {
			for (int i = 0; i < L; i++) out.push_back(c0.bidcompl[codes[0][p + L - 1 - i]]);
			for (int t = 1; t < T; t++) for (int i = 0; i < L; i++) out.push_back(codes[t][p + L - 1 - i]);
			added++;
		}
	}
	return added;
}
