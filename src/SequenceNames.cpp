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
		seqNames.clear(); weight.clear();
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
		seqNames.clear(); weight.clear();
		Nseqs = 0; 		
	}
	
	char stmp[10000]; 
	
	FILE *f = fopen(seqNamesFN, "r") ; 
	if (f == NULL) { sprintf(globtmpstr,"\n ERROR: cannot open %s\n", seqNamesFN); Printf(globtmpstr); return 0; }
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
		seqNames.clear(); weight.clear();
		Nseqs = 0; 		
	}
	
	char stmp[10000]; 
	
	FILE *f = fopen(seqNamesFN, "r") ; 
	if (f == NULL) { sprintf(globtmpstr,"\n ERROR: cannot open %s\n", seqNamesFN); Printf(globtmpstr); return 0; }
	while (!feof(f))
	{
		if (fgets(stmp, 10000-5, f)) 
		{
			if (stmp[0]!=0) {
				char *name = new char[strlen(stmp)+1];  // a %s token can never be longer than the line (was a fixed 100 bytes)
				double w = 0;
				if (sscanf(stmp, "%s%lf", name, &w)!=2) { delete []name; continue; } // blank/malformed line
				seqNames.push_back(name); weight.push_back(w);
				Nseqs++; 
			}
		}

	}
	fclose(f); 
	return Nseqs; 
}

void CSequenceNames::openSeqFile( const char *seqFN,  int maxSeqLength)
{
	this->seqf = fopen(seqFN, "r"); 
	if (this->seqf == NULL) { sprintf(globtmpstr,"\n ERROR: cannot open %s\n", seqFN); Printf(globtmpstr); }

	if (this->curSeq!=NULL) delete curSeq; 

	this->curSeq = new CSequence(maxSeqLength); 
}

CSequence *CSequenceNames::nextSeq()
{
	while (seqf!=NULL && !feof(seqf))
	{
		if (this->nextSeqtoRead ==0)
		{
			if (curSeq->readFsa(this->seqf) < 0) { fclose(seqf); seqf = NULL; return NULL; }  // read a new sequence; -1 = too long
		}

		while (nextSeqtoRead<this->Nseqs)
		{
			if (stringcompare(this->seqNames[nextSeqtoRead], curSeq->getName(), MAX_LINE_WIDTH)) // full-length comparison (was truncated at 100 chars)			//if (strcmp(this->seqNames[nextSeqtoRead], curSeq->getName())==0) 
			{
				curSeq->setWeight(weight[nextSeqtoRead]); 
				curSeq->setNameLink(seqNames[nextSeqtoRead]); 

				nextSeqtoRead++; 
				nSeqsRead++; 

				if (nSeqsRead==Nseqs)
				{
					fclose(seqf);
					seqf = NULL; 
				}

				return curSeq; 
			}
			nextSeqtoRead++; 
		}
		nextSeqtoRead = 0; 
	}

	if (seqf!=NULL) fclose(seqf); 
	seqf = NULL; 
	return NULL; 
}
