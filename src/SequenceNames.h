/* SequenceNames.h
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "global.h"
#include "Sequence.h"

#include <vector>
#include <string>
#include <unordered_map>

class CSequenceNames
{
public:
	CSequenceNames(void);
	~CSequenceNames(void);
	int Nseqs; 
	std::vector<char *> seqNames; 
	std::vector<double> weight; 

	int readSeqNames(const char *seqNamesFN); 
	int readSeqNamesandWeights(const char *seqNamesFN); 

	void openSeqFile( const char *seqFN, int maxSeqLength, const CConverter *conv); 
	CSequence *nextSeq(); 

	CSequence *curSeq; 
private: 
	FILE *seqf; 
	int nSeqsRead;  // number of times a sequence was returned by nextSeq()
	int nextSeqtoRead;  // (unused since Phase 3; kept for layout stability)
	std::unordered_map<std::string, int> nameIndex; // name -> position in seqNames/weight
	std::vector<char> used;                         // which alpha entries have been matched
public:
	int error;      // set when a name in the alpha file is duplicated, missing from the FASTA, or repeated in the FASTA
	double rho;     // bias from a .gkmmodel file (0 for the legacy alpha file): score = sum_j alpha_j K(x, x_j) - rho
	int isModelFile; // 1 when the alpha file was a .gkmmodel (then the same file also holds the sequences)
};

