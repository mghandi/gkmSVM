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

#define MAXNSeqs 1000000 
#define MAXSeqnameLENGTH 100 

class CSequenceNames
{
public:
	CSequenceNames(void);
	~CSequenceNames(void);
	int Nseqs; 
	char *seqNames[MAXNSeqs]; 

	double weight[MAXNSeqs]; 

	int readSeqNames(char *seqNamesFN); 
	int readSeqNamesandWeights(char *seqNamesFN); 

	void openSeqFile( char *seqFN, int maxSeqLength); 
	CSequence *nextSeq(); 

	CSequence *curSeq; 
private: 
	FILE *seqf; 
	int nSeqsRead;  // number of times a sequence was returned by nextSeq()
	int nextSeqtoRead;  // index to the next potential sequence 
};

