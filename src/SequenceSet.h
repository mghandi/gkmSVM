/* SequenceSet.h : explicit, positional sequence identity (Phase 3 of dev/REFACTORING_PLAN.md).
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * A kernel row is identified by its position; the FASTA name is metadata that may repeat, be empty
 * or be arbitrarily long. Ids ("pos1".., "neg1"..) are assigned at read time and are unique by
 * construction. Because the text kernel format carries no identifiers, gkmsvm_kernel writes the
 * table next to it as <outfile>.index (tab separated, one header line):
 *     row  id  name  label  length  nlmers
 */
#pragma once
#include <string>
#include <vector>
#include <cstdio>

struct SeqRecord {
	int index;          // 0..n-1 == kernel row
	std::string id;     // pos1.., neg1..
	std::string name;   // original FASTA name (first whitespace-delimited token of the header)
	int label;          // +1 positive, -1 negative
	long length;        // total residues (sum over merged records when mergeByName is on)
	long nlmers;        // L-mers contributed to the kernel (both strands when addRC)
};

inline std::string seqRecordId(int index, int npos) {
	char buf[32];
	if (index < npos) snprintf(buf, sizeof buf, "pos%d", index + 1);
	else snprintf(buf, sizeof buf, "neg%d", index - npos + 1);
	return buf;
}

// Writes <outfile>.index; returns 0 on success.
inline int writeIndexSidecar(const std::string &outfile, const std::vector<SeqRecord> &recs) {
	std::string fn = outfile + ".index";
	FILE *f = fopen(fn.c_str(), "w");
	if (f == NULL) return 1;
	fprintf(f, "row\tid\tname\tlabel\tlength\tnlmers\n");
	for (size_t i = 0; i < recs.size(); i++) {
		const SeqRecord &r = recs[i];
		fprintf(f, "%d\t%s\t%s\t%d\t%ld\t%ld\n", r.index, r.id.c_str(), r.name.c_str(), r.label, r.length, r.nlmers);
	}
	fclose(f);
	return 0;
}
