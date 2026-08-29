/* SequenceNames.cpp
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

#include "SequenceNames.h"
#include <string.h>


CSequenceNames::CSequenceNames(void)
{
	Nseqs=0; 
		
	seqf = NULL; 

	nSeqsRead=0;
	nextSeqtoRead=0;
	curSeq=NULL; 
	error=0;
	rho=0;
	isModelFile=0;
}


CSequenceNames::~CSequenceNames(void)
{
	int i; 
	if (Nseqs!=0)
	{
		for(i=0;i<Nseqs;i++)
		{
			delete []seqNames[i]; 
		}
		seqNames.clear(); weight.clear(); used.clear(); nameIndex.clear();
		Nseqs = 0; 		
	}

	if (seqf!=NULL)
	{
		fclose(seqf); 
		seqf = NULL; 
	}

	if (this->curSeq!=NULL) delete curSeq; 

}

int CSequenceNames::readSeqNames(const char *seqNamesFN)
{
	int i; 
	if (Nseqs!=0)
	{
		for(i=0;i<Nseqs;i++)
		{
			delete []seqNames[i]; 
		}
		seqNames.clear(); weight.clear(); used.clear(); nameIndex.clear();
		Nseqs = 0; 		
	}
	
	char stmp[10000]; 
	
	FILE *f = fopen(seqNamesFN, "r") ; 
	if (f == NULL) { gkmMsg("\n ERROR: cannot open %s\n", seqNamesFN); return 0; }
	while (!feof(f))
	{
		if (fgets(stmp, 10000-5, f))
		{
			if (stmp[0]!=0) {
				char *name = new char[strlen(stmp)+1];  // a %s token can never be longer than the line
				if (sscanf(stmp, "%s", name)!=1) { delete []name; continue; }
				seqNames.push_back(name); weight.push_back(0.0);
				Nseqs++; 
			}
		}
	}
	fclose(f); 
	return Nseqs; 
}


int CSequenceNames::readSeqNamesandWeights(const char *seqNamesFN)
{
	int i; 
	if (Nseqs!=0)
	{
		for(i=0;i<Nseqs;i++)
		{
			delete []seqNames[i]; 
		}
		seqNames.clear(); weight.clear(); used.clear(); nameIndex.clear();
		Nseqs = 0; 		
	}
	
	char stmp[10000]; 
	
	FILE *f = fopen(seqNamesFN, "r") ; 
	if (f == NULL) { gkmMsg("\n ERROR: cannot open %s\n", seqNamesFN); return 0; }
	// A .gkmmodel file (Phase 4) is a FASTA file whose headers are ">id<TAB>alpha" and whose
	// comment lines "#key value" carry the bias: "#gkmmodel 1" must be the first line.
	while (!feof(f))
	{
		if (fgets(stmp, 10000-5, f)) 
		{
			if (Nseqs==0 && !isModelFile && strncmp(stmp, "#gkmmodel", 9)==0) { isModelFile = 1; continue; }
			if (isModelFile) {
				if (stmp[0]=='#') { double v; if (sscanf(stmp, "#rho %lf", &v)==1) rho = v; continue; }
				if (stmp[0]!='>') continue; // sequence lines
				memmove(stmp, stmp+1, strlen(stmp)); // drop the '>' and parse "id alpha" below
			}
			if (stmp[0]!=0) {
				char *name = new char[strlen(stmp)+1];  // a %s token can never be longer than the line (was a fixed 100 bytes)
				double w = 0;
				if (sscanf(stmp, "%s%lf", name, &w)!=2) { delete []name; continue; } // blank/malformed line
				if (nameIndex.count(name)) {
					gkmMsg("\n ERROR: sequence name '%s' appears more than once in %s\n", name, seqNamesFN);
					delete []name; error = 1; break;
				}
				nameIndex[name] = Nseqs;
				seqNames.push_back(name); weight.push_back(w); used.push_back(0);
				Nseqs++; 
			}
		}

	}
	fclose(f); 
	return Nseqs; 
}

void CSequenceNames::openSeqFile( const char *seqFN,  int maxSeqLength, const CConverter *conv)
{
	this->seqf = fopen(seqFN, "r"); 
	if (this->seqf == NULL) { gkmMsg("\n ERROR: cannot open %s\n", seqFN); }

	if (this->curSeq!=NULL) delete curSeq; 

	this->curSeq = new CSequence(maxSeqLength, conv); 
}

// Returns the next FASTA record (in file order) whose name is listed in the alpha file, with its
// weight attached; NULL at the end of the file or on error (check `error`). Every alpha entry must
// be matched exactly once: a listed name missing from the FASTA, or present twice, is an error
// instead of silently scoring with a partial model (pre-Phase-3 behaviour).
CSequence *CSequenceNames::nextSeq()
{
	while (seqf!=NULL && !feof(seqf))
	{
		if (curSeq->readFsa(this->seqf) < 0) { fclose(seqf); seqf = NULL; error = 1; return NULL; }
		if (curSeq->getLength() <= 0) continue;
		std::unordered_map<std::string, int>::iterator it = nameIndex.find(curSeq->getName());
		if (it == nameIndex.end()) continue; // not a support vector
		int k = it->second;
		if (used[k]) {
			gkmMsg("\n ERROR: support vector '%s' appears more than once in the sequence file\n", curSeq->getName());
			fclose(seqf); seqf = NULL; error = 1; return NULL;
		}
		used[k] = 1;
		curSeq->setWeight(weight[k]); 
		curSeq->setNameLink(seqNames[k]); 
		nSeqsRead++; 
		return curSeq; 
	}

	if (seqf!=NULL) fclose(seqf); 
	seqf = NULL; 
	if (nSeqsRead < Nseqs) {
		int firstMissing = 0; while (firstMissing < Nseqs && used[firstMissing]) firstMissing++;
		gkmMsg("\n ERROR: the sequences for only %d out of %d support vectors were found in the sequence file (first missing: '%s')\n", nSeqsRead, Nseqs, seqNames[firstMissing]);
		error = 1;
	}
	return NULL; 
}
