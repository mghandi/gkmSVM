/* Sequence.cpp: implementation of the CSequence class.
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

// I have changed this class to partially support generalized alphabet, originally it was designed only for DNA --
// WARNING: some of the functionality (such as reverse complement, dinucl, etc) will not work for general alphabet (globalConverter does not support them)
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"

#include "Sequence.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSequence::CSequence(int maxLength,  CSequence *sCopyFrom )
{
	static int serialnumber = 0; 
	
	this->maxLength = maxLength;
	seqName = new char[MAX_LINE_WIDTH]; 
	seqLabel = new char[MAX_LINE_WIDTH]; 
	seq = new char[maxLength]; 
	subseq = new char[maxLength]; 
	seqBaseId = new baseId[maxLength];
	length = 0; 
	dinucl = new dinuclId[maxLength]; 
	reverseComplement = NULL; 
	sline = new char[MAX_LINE_WIDTH+3];
	nextName = new char[MAX_LINE_WIDTH];
	hasNextName = 0;

//	seqName = new char(maxLength); 
	sprintf(seqName,"seq_%d",serialnumber++);	 
	sprintf(seq,"%s","");
	sprintf(seqLabel,"%s","");
	//bendingData = NULL; 
	//this->occupancyData = NULL;
	weight = 0; 
	NameLink = NULL; 
	readError = 0; 

	if (sCopyFrom!=NULL)
	{
		// does not copy reverseCompl, if needed later, will recalculate that
		CSequence *s = sCopyFrom;
		length = s->getLength(); 
		sprintf(seqName,"%s", s->getName());
		sprintf(seqLabel,"%s", s->getLabel());
		int i; 
		char *sseq = s->getSeq(); 
		char *ssubseq = s->getSubseq(); 
		dinuclId *sdinucl = s->getDinucl(); 
		baseId *sseqBaseId = s->getSeqBaseId();
		weight = s->weight; 
		NameLink = s->NameLink; 
		for(i=0;i<length;i++)
		{
			seq[i] = sseq[i]; 
			subseq[i] = ssubseq[i]; 
			seqBaseId[i] = sseqBaseId[i]; 
			dinucl[i] = sdinucl[i]; 
		}
	}

}

CSequence::~CSequence()
{

	if (this->reverseComplement!=NULL)
	{
		delete reverseComplement; 
		this->reverseComplement = NULL; 
	}

	freeMem(seq); 
	freeMem(subseq); 
	freeMem(seqName); 
	freeMem(seqLabel); 
	freeMem(dinucl); 
	freeMem(seqBaseId); 
	freeMem(sline); 
	freeMem(nextName); 

	length = 0; 
}


char *CSequence::getSeq()
{
	return this->seq;
}

char *CSequence::getSubseq()
{
	return this->subseq;
}

double CSequence::getWeight()
{
	return this->weight;
}

void CSequence::setWeight(double w)
{
	this->weight=w;
}

char * CSequence::getNameLink()
{
	return this->NameLink;
}

void CSequence::setNameLink(char *sp)
{
	this->NameLink=sp;
}


char *CSequence::getSubseq(int p1,int p2)
{
	if(p1<0) {p1=0; }
	if(p2>length-1) {p2 = length-1; }
	for(int i=p1;i<=p2;i++)
	{
		subseq[i-p1]=seq[i]; 
	}
	subseq[p2-p1+1]=0; 
	return this->subseq;
}

int *CSequence::getSubseqBaseId(int p1,int p2, int *obid)
{
	if(p1<0) {p1=0; }
	if(p2>(length-1)) {p2 = length-1; }
	for(int i=p1;i<=p2;i++)
	{
		obid[i-p1]=this->seqBaseId[i]; 
	}
	return obid;
}


char *CSequence::getName()
{
	return this->seqName;
}

char *CSequence::getLabel()
{
	return this->seqLabel;
}

int CSequence::getLength() 
{
	return this->length;
}

dinuclId *CSequence::getDinucl()
{
	return this->dinucl;
}

baseId *CSequence::getSeqBaseId()
{
	return this->seqBaseId;
}


int CSequence::readFsa(FILE *f, int SkipAlphabetCheck)  // reads one sequence from already opened file f;
{
	length = 0; 
	readError = 0; 
	
	if (f==NULL) 
	{
		return 0; 
	}
	int overflow = 0; 

	if (fgets(sline, MAX_LINE_WIDTH, f) == NULL) sline[0] = 0; 
	if (sline[0]=='>') //first line
	{
		if (sscanf(sline+1, "%s", nextName) != 1) nextName[0] = 0; // ">" with no name
		if (fgets(sline, MAX_LINE_WIDTH, f) == NULL) sline[0] = 0;
		hasNextName = 1;  
	}

	sprintf(seqName,"%s",(hasNextName?nextName:"NA")); 
	hasNextName = 0; 
	
	while(true)
	{
		if (feof(f) || (sline[0]=='>'))
		{
			break;
		}
		int i=0; 

		if (sline[0]!=';' && sline[0]!='#') // ignore comment lines starting with ";" or "#" (the .gkmmodel header)
		{
			while (sline[i]!=0)
			{
/*				if ((ucase(sline[i])=='A')||
					(ucase(sline[i])=='C')||
					(ucase(sline[i])=='G')||
					(ucase(sline[i])=='T'))
					*/

				if ((globalConverter.isInAlphabet[(unsigned char)sline[i]])||SkipAlphabetCheck)
				{
					if (length >= maxLength-1) { overflow = 1; break; } // leave room for the terminator
					this->seq[this->length] = sline[i]; 
					length++;
				}
				i++;
			}
		}
		if (fgets(sline, MAX_LINE_WIDTH, f) == NULL) sline[0] = 0; 
	}
	if (sline[0]=='>') //first line
	{
		if (sscanf(sline+1, "%s", nextName) != 1) nextName[0] = 0; 
		hasNextName =1; 
	}
	
	seq[length] = 0; 

	if (overflow)
	{
		sprintf(globtmpstr,"\n ERROR: sequence '%s' is longer than the maximum sequence length %d (set -m / maxseqlen to at least its length).\n", seqName, maxLength-3); Printf(globtmpstr);
		readError = 1; 
		length = 0; 
		return -1; 
	}
	if (length==0) // empty record, or EOF: nothing to convert (seqBaseId[-1] used to be written here)
	{
		return 0; 
	}

	int i; 
	int *cidx = globalConverter.cidx;
	for(i=0;i<length-1;i++)
	{
		this->seqBaseId[i] = cidx[(unsigned char)seq[i]];
		this->dinucl[i] = globalConverter.dnidx(seq+i);  // for i==length-1 it is meaningless
	}
	this->seqBaseId[length-1] = cidx[(unsigned char)seq[length-1]];

	if (TALK)
	{
		//printf("\n%s\t%s\t%d\n%s",seqName,seqLabel, length,seq);
	}
			
	return length; 
}

void CSequence::writeFsa(FILE *f)  
{
	if (f==NULL)
	{
		Printf("\n cannot write to file (file not open)");
		return;
	}
	fprintf(f,">%s\t%d\t%s",this->seqName, length, this->seqLabel);
	int i=0; 
	while (i<length)
	{
		if (i%60==0)
		{
			fprintf(f,"\n");
		}
		fprintf(f,"%c",seq[i]); 
		i++; 
	}
	fprintf(f,"\n");
}


int CSequence::readBasic(FILE *f)  // reads one sequence from already opened file f; format: name\tlabel\tseq 
{

	length = 0; 
	
	if (f==NULL) 
	{
		return 0; 
	}

	if (fgets(sline, MAX_LINE_WIDTH, f) == NULL) sline[0] = 0; 

	sscanf(sline,"%s%s%s", this->seqName, this->seqLabel, this->seq);
	
	this->length = strlength(seq); 

	if (length==0)
	{
		return length; 
	}
	
	int i; 
	int *cidx = globalConverter.cidx;
	for(i=0;i<length-1;i++)
	{
		this->seqBaseId[i] = cidx[(unsigned char)seq[i]]; 
		this->dinucl[i] = globalConverter.dnidx(seq+i);  // for i==length-1 it is meaningless
	}
	this->seqBaseId[length-1] = cidx[(unsigned char)seq[length-1]];

	if (TALK)
	{
		//printf("\n%s\t%s\t%d\n%s",seqName,seqLabel, length,seq);
	}
				
	return length; 
}

int CSequence::readString(char *s)  // reads sequence from a string (coverts chat * to sequence)
{

	length = 0; 
	
	sscanf(s,"%s", this->seq);
	
	this->length = strlength(seq); 

	if (length==0)
	{
		return length; 
	}
	
	int i; 
	int *cidx = globalConverter.cidx;
	for(i=0;i<length-1;i++)
	{
		this->seqBaseId[i] = cidx[(unsigned char)seq[i]]; 
		this->dinucl[i] = globalConverter.dnidx(seq+i);  // for i==length-1 it is meaningless
	}
	this->seqBaseId[length-1] = cidx[(unsigned char)seq[length-1]];

	return length; 
}

void CSequence::writeBasic(FILE *f)  
{
	if (f==NULL)
	{
		Printf("\n cannot write to file (file not open)");
		return;
	}
	fprintf(f,"%s\t%s\t%s\n",this->seqName,this->seqLabel,this->seq);
}

void CSequence::mutateOneBase(int pos, baseId nwbid)
{
	if (pos>=length)
	{
		sprintf(globtmpstr,"\n error : cannot mutate pos %d while length is %d",pos, length); Printf(globtmpstr);
		return; 
	}

	this->seq[pos]= globalConverter.icidx[nwbid]; 
	this->seqBaseId[pos] = nwbid; 
	if (pos>0)
	{
		this->dinucl[pos-1] = globalConverter.dnidx(&seq[pos-1]);
	}
	if (pos<length-1)
	{
		this->dinucl[pos] = globalConverter.dnidx(&seq[pos]);
	}

}


CSequence *CSequence::getReverseComplement()
{
	if (this->reverseComplement==NULL)
	{
		this->reverseComplement = new CSequence(this->maxLength, this); 
	}
	else
	{
		this->reverseComplement->length = length; 
		// copy the name/label to the reverse-complement object (the original code copied seqName onto
		// itself: undefined for overlapping buffers, and reported by ASAN's sprintf interceptor on Linux)
		sprintf(reverseComplement->seqName,"%s", seqName);
		sprintf(reverseComplement->seqLabel,"%s",seqLabel); 
	}
	int i;


	char *rseq = reverseComplement->getSeq();  
	
	dinuclId *rdinucl = reverseComplement->getDinucl();  	
	baseId *rseqBaseId = reverseComplement->getSeqBaseId();  	
	for(i=0;i<length;i++)
	{
		rseq[i] = globalConverter.bcompl[(unsigned char)seq[length-1-i]];
	}

	rseq[length] =0; 





	for(i=0;i<length-1;i++)
	{
		rseqBaseId[i] = globalConverter.cidx[(unsigned char)rseq[i]]; 
		rdinucl[i] = globalConverter.dnidx(rseq+i);  // for i==length-1 it is meaningless
	}
	rseqBaseId[length-1] = globalConverter.cidx[(unsigned char)rseq[length-1]];

	
	return reverseComplement; 
}



int countKLmerHitsNDCONVUPPERC(char *KLmerseq, int L, char *s, int size)
{
	int N = size-L+1; 
	int i,j;
	
	for(i=0;i<L;i++) KLmerseq[i]=toupper(KLmerseq[i]);
	for(i=0;i<size;i++) s[i]=toupper(s[i]);


	int cnt = 0; 
	for (i=0;i<N;i++)
	{
		for(j=0;j<L;j++)
		{
			if ((KLmerseq[j]!='.')&&(s[j]!=KLmerseq[j]))
			{
				break; 
			}
		}
		if (j==L) cnt++;
		s++;
	}
	return cnt; 
}
