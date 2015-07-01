/* LTreef.cpp
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

#include "LTreef.h"
#include "LTreeS.h"
#include "Sequence.h"

CLTreef::CLTreef(void)
{
//	daughter[0].f=daughter[1].f=daughter[2].f=daughter[3].f=0;
//	daughter[0].p=daughter[1].p=daughter[2].p=daughter[3].p=NULL;
	for(int i=0;i<MAX_ALPHABET_SIZE;i++){
		daughter[i].f=0;
		daughter[i].p=NULL;
	}
	nonEmptyDaughterCnt=0;

}


CLTreef::~CLTreef(void)
{
}

void CLTreef::addSeq(int *bid, int n, int cnt)
{
	if (n==1)
	{
		if (this->daughter[*bid].i==0){
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		this->daughter[*bid].i += cnt;
	}
	else
	{
		if (this->daughter[*bid].p == NULL)
		{
			this->daughter[*bid].p = new CLTreef(); 
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		daughter[*bid].p->addSeq(bid+1, n-1, cnt); 
	}
}

void CLTreef::addSeq(int *bid, int n, myFlt w)
{
	if (n==1)
	{
		//check the given bid has already been added.  debugged by Dongwon Lee 1/11/15
		int found = 0;
		for (int i=0;i<this->nonEmptyDaughterCnt;i++) 
		{
			if (this->nonEmptyDaughterIdxs[i] == *bid) 
			{
				found = 1;
				break;
			}
		}
		if (found == 0){
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		this->daughter[*bid].f +=w; 
	}
	else
	{
		if (this->daughter[*bid].p == NULL)
		{
			this->daughter[*bid].p = new CLTreef(); 
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		daughter[*bid].p->addSeq(bid+1, n-1, w); 
	}
}


int CLTreef::icount(int *bid, int n) //returns the number of times the sequence is in found in the tree
{
	if (n==1)
	{
		return daughter[*bid].i; 
	}
	else
	{
		if (this->daughter[*bid].p == NULL)
		{
			return 0; 
		}
		return daughter[*bid].p->icount(bid+1, n-1); 
	}
}


myFlt CLTreef::fcount(int *bid, int n) //returns the number of times the sequence is in found in the tree
{
	if (n==1)
	{
		return daughter[*bid].f; 
	}
	else
	{
		if (this->daughter[*bid].p == NULL)
		{
			return 0; 
		}
		return daughter[*bid].p->fcount(bid+1, n-1); 
	}
}

int CLTreef::leavesCount(int withMultiplicity, int n, int alphabetSize)  //returns the number of sequences in the tree.  //call with n=L from outside
{
	int nleaves = 0; 
	for (int i=0;i<alphabetSize;i++)
	{
		if (daughter[i].p!=NULL)
		{
			if (n==1)
			{
				if (withMultiplicity)
				{
					int frq= daughter[i].i; 
					nleaves +=frq; 
				}
				else 
				{
					nleaves++; 					
				}
			}
			else
			{
				nleaves+=daughter[i].p->leavesCount(withMultiplicity, n-1, alphabetSize);
			}
		}
	}
	return nleaves; 
}



void CLTreef::deleteTree(int n, int alphabetSize)
{
	if (n>1)
	{
		for (int i=0;i<alphabetSize;i++)
		{
			if (daughter[i].p!=NULL)
			{
				daughter[i].p->deleteTree(n-1,alphabetSize);
				delete daughter[i].p;
			}
		}
	}
}


int CLTreef::addSequence(int *bid, int n, int L)  //adds all the L-subseqs 
{
	n = n-L+1;
	for(int i=0;i<n;i++)
	{
		addSeq(bid,L,1);
		bid++;
	}
	return n; 
}
int CLTreef::addSequence(int *bid, int n, int L, myFlt w) //adds all the L-subseqs w times 
{
	n = n-L+1;
	for(int i=0;i<n;i++)
	{
		addSeq(bid,L,w);
		bid++;
	}
	return n; 
}




void  CLTreef::addSequences(char *FsaFileName, int L, int maxSequenceLength, int addrevcompl, int numberOfCVPartitions, int selectPartitionNumber) // adds all the sequence from a FSA file 
{
	if (numberOfCVPartitions==0)
	{
		numberOfCVPartitions = 1; 
	}
	selectPartitionNumber = selectPartitionNumber % numberOfCVPartitions; 

	FILE *fi;

	fi = fopen(FsaFileName,"r"); 
		
	CSequence *s = new CSequence(maxSequenceLength+3);

	int counter = 0; 

	while (!feof(fi))
	{
		s->readFsa(fi); 
		if(s->getLength()>0)
		{
			counter++; 
			if (counter % numberOfCVPartitions==selectPartitionNumber)
			{ 
				this->addSequence(s->getSeqBaseId(), s->getLength(), L); 

				if (addrevcompl)
				{
					this->addSequence(s->getReverseComplement()->getSeqBaseId(), s->getLength(), L); 
				}
			}
		}
	}

	fclose(fi); 

	delete s; 
}




void CLTreef::imismatchCount(int *bid, int n, int *cnt) // fills the mismatch count array (how many sequences with m mismatches to a given sequence exist)
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=daughter[i].i;
		cnt++; 
		i=(i+1)&3;  //(i+1)%4
		*cnt+=daughter[i].i;
		i=(i+1)&3; 
		*cnt+=daughter[i].i;
		i=(i+1)&3; 
		*cnt+=daughter[i].i;
	}
	else
	{
		n--; 
		bid++; 
		
		if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt);

		cnt++; 
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt);

	}
}


void CLTreef::fmismatchCount(int *bid, int n, myFlt *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=daughter[i].f;
		cnt++; 
		i=(i+1)&3;  //(i+1)%4
		*cnt+=daughter[i].f;
		i=(i+1)&3; 
		*cnt+=daughter[i].f;
		i=(i+1)&3; 
		*cnt+=daughter[i].f;
	}
	else
	{
		n--; 
		bid++; 
		
		if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt);

		cnt++; 
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt);

	}
}


void CLTreef::imismatchCount(int *bid, int n, int *cnt, int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=daughter[i].i;
		if (maxmm!=0)
		{
			cnt++; 
			i=(i+1)&3;  //(i+1)%4
			*cnt+=daughter[i].i;
			i=(i+1)&3; 
			*cnt+=daughter[i].i;
			i=(i+1)&3; 
			*cnt+=daughter[i].i;
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--; 
			bid++; 
		
			if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt,maxmm);

			maxmm--;
			cnt++; 
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->imismatchCount(bid,n,cnt,maxmm);
		}
		else
		{   // fast track for maxmm==0
			CLTreef *cur = this;  
			while (--n)
			{
				if ((cur = cur->daughter[*bid++].p)==NULL) return;  
			}
			*cnt+=cur->daughter[*bid].i; 
		}

	}
}

void CLTreef::fmismatchCount(int *bid, int n, myFlt *cnt ,int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=daughter[i].f;
		if (maxmm!=0)
		{
			cnt++; 
			i=(i+1)&3;  //(i+1)%4
			*cnt+=daughter[i].f;
			i=(i+1)&3; 
			*cnt+=daughter[i].f;
			i=(i+1)&3; 
			*cnt+=daughter[i].f;
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--; 
			bid++; 
		
			if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt,maxmm);

			maxmm--;
			cnt++; 
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i].p!=NULL) daughter[i].p->fmismatchCount(bid,n,cnt,maxmm);
		}
		else
		{   // fast track for maxmm==0
			CLTreef *cur = this;  
			while (--n)
			{
				if ((cur = cur->daughter[*bid++].p)==NULL) return;  
			}
			*cnt+=cur->daughter[*bid].f; 
		}

	}
}



void CLTreef::fimismatchCount(class CLTreef *iTree, int n, myFlt *cnt ,int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given set of sequences exist) 
{
	fintptr_t df0 = daughter[0];
	fintptr_t df1 = daughter[1];
	fintptr_t df2 = daughter[2];
	fintptr_t df3 = daughter[3];

	fintptr_t *idaughter = iTree->daughter;

	fintptr_t di0 = idaughter[0];
	fintptr_t di1 = idaughter[1];
	fintptr_t di2 = idaughter[2];
	fintptr_t di3 = idaughter[3];

	if (n==1)
	{
		*cnt+=df0.f * di0.i+
			  df1.f * di1.i+
			  df2.f * di2.i+
			  df3.f * di3.i;

		if (maxmm!=0)
		{
			cnt++;
#ifdef FTHENI
			if (df0.p!=NULL)
				*cnt+=  df0.f*(di1.i+di2.i+di3.i);
			if (df1.p!=NULL)
				*cnt+=  df1.f*(di0.i+di2.i+di3.i);
			if (df2.p!=NULL)
				*cnt+=  df2.f*(di1.i+di0.i+di3.i);
			if (df3.p!=NULL)
				*cnt+=  df3.f*(di1.i+di2.i+di0.i);
#endif
#ifndef FTHENI
			if (di0.p!=NULL)
				*cnt+=  di0.i*(df1.f+df2.f+df3.f);
			if (di1.p!=NULL)
				*cnt+=  di1.i*(df0.f+df2.f+df3.f);
			if (di2.p!=NULL)
				*cnt+=  di2.i*(df1.f+df0.f+df3.f);
			if (di3.p!=NULL)
				*cnt+=  di3.i*(df1.f+df2.f+df0.f);
#endif
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--;
			myFlt *cnt1 = cnt+1; 
			int maxmm1 = maxmm-1;

#ifdef FTHENI
			if (df0.p!=NULL)
			{
				if (di0.p!=NULL) df0.p->fimismatchCount(di0.p,n,cnt,maxmm);
				if (di1.p!=NULL) df0.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df0.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df0.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df1.p!=NULL)
			{
				if (di0.p!=NULL) df1.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df1.p->fimismatchCount(di1.p,n,cnt,maxmm);
				if (di2.p!=NULL) df1.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df1.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df2.p!=NULL)
			{
				if (di0.p!=NULL) df2.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df2.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df2.p->fimismatchCount(di2.p,n,cnt,maxmm);
				if (di3.p!=NULL) df2.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df3.p!=NULL)
			{
				if (di0.p!=NULL) df3.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df3.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df3.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df3.p->fimismatchCount(di3.p,n,cnt,maxmm);
			}
#endif
#ifndef FTHENI
			if (di0.p!=NULL)
			{
				if (df0.p!=NULL) df0.p->fimismatchCount(di0.p,n,cnt,maxmm);
				if (df1.p!=NULL) df1.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
				if (df2.p!=NULL) df2.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
				if (df3.p!=NULL) df3.p->fimismatchCount(di0.p,n,cnt1,maxmm1);
			}
			if (di1.p!=NULL)
			{
				if (df0.p!=NULL) df0.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
				if (df1.p!=NULL) df1.p->fimismatchCount(di1.p,n,cnt,maxmm);
				if (df2.p!=NULL) df2.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
				if (df3.p!=NULL) df3.p->fimismatchCount(di1.p,n,cnt1,maxmm1);
			}
			if (di2.p!=NULL)
			{
				if (df0.p!=NULL) df0.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
				if (df1.p!=NULL) df1.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
				if (df2.p!=NULL) df2.p->fimismatchCount(di2.p,n,cnt,maxmm);
				if (df3.p!=NULL) df3.p->fimismatchCount(di2.p,n,cnt1,maxmm1);
			}
			if (di3.p!=NULL)
			{
				if (df0.p!=NULL) df0.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
				if (df1.p!=NULL) df1.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
				if (df2.p!=NULL) df2.p->fimismatchCount(di3.p,n,cnt1,maxmm1);
				if (df3.p!=NULL) df3.p->fimismatchCount(di3.p,n,cnt,maxmm);
			}
#endif

		}
		else
		{
			// maybe we can implement fast track to calculate it faster in a non-recursive form maybe

			n--;

#ifdef FTHENI
			if (df0.p!=NULL)
			{
				if (di0.p!=NULL) df0.p->fimismatchCount(di0.p,n,cnt,maxmm);
			}
			if (df1.p!=NULL)
			{
				if (di1.p!=NULL) df1.p->fimismatchCount(di1.p,n,cnt,maxmm);
			}
			if (df2.p!=NULL)
			{
				if (di2.p!=NULL) df2.p->fimismatchCount(di2.p,n,cnt,maxmm);
			}
			if (df3.p!=NULL)
			{
				if (di3.p!=NULL) df3.p->fimismatchCount(di3.p,n,cnt,maxmm);
			}
#endif
#ifndef FTHENI
			if (di0.p!=NULL)
			{
				if (df0.p!=NULL) df0.p->fimismatchCount(di0.p,n,cnt,maxmm);
			}
			if (di1.p!=NULL)
			{
				if (df1.p!=NULL) df1.p->fimismatchCount(di1.p,n,cnt,maxmm);
			}
			if (di2.p!=NULL)
			{
				if (df2.p!=NULL) df2.p->fimismatchCount(di2.p,n,cnt,maxmm);
			}
			if (di3.p!=NULL)
			{
				if (df3.p!=NULL) df3.p->fimismatchCount(di3.p,n,cnt,maxmm);
			}
#endif
		}
	}
}


void CLTreef::iimismatchCount(class CLTreef *iTree, int n, int *cnt ,int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given set of sequences exist) 
{
	fintptr_t df0 = daughter[0];
	fintptr_t df1 = daughter[1];
	fintptr_t df2 = daughter[2];
	fintptr_t df3 = daughter[3];

	fintptr_t *idaughter = iTree->daughter;

	fintptr_t di0 = idaughter[0];
	fintptr_t di1 = idaughter[1];
	fintptr_t di2 = idaughter[2];
	fintptr_t di3 = idaughter[3];

	if (n==1)
	{
		*cnt+=df0.i * di0.i+
			  df1.i * di1.i+
			  df2.i * di2.i+
			  df3.i * di3.i;

		if (maxmm!=0)
		{
			cnt++;

			if (df0.p!=NULL)
				*cnt+=  df0.i*(di1.i+di2.i+di3.i);
			if (df1.p!=NULL)
				*cnt+=  df1.i*(di0.i+di2.i+di3.i);
			if (df2.p!=NULL)
				*cnt+=  df2.i*(di1.i+di0.i+di3.i);
			if (df3.p!=NULL)
				*cnt+=  df3.i*(di1.i+di2.i+di0.i);
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--;
			int *cnt1 = cnt+1; 
			int maxmm1 = maxmm-1;

			if (df0.p!=NULL)
			{
				if (di0.p!=NULL) df0.p->iimismatchCount(di0.p,n,cnt,maxmm);
				if (di1.p!=NULL) df0.p->iimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df0.p->iimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df0.p->iimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df1.p!=NULL)
			{
				if (di0.p!=NULL) df1.p->iimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df1.p->iimismatchCount(di1.p,n,cnt,maxmm);
				if (di2.p!=NULL) df1.p->iimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df1.p->iimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df2.p!=NULL)
			{
				if (di0.p!=NULL) df2.p->iimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df2.p->iimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df2.p->iimismatchCount(di2.p,n,cnt,maxmm);
				if (di3.p!=NULL) df2.p->iimismatchCount(di3.p,n,cnt1,maxmm1);
			}
			if (df3.p!=NULL)
			{
				if (di0.p!=NULL) df3.p->iimismatchCount(di0.p,n,cnt1,maxmm1);
				if (di1.p!=NULL) df3.p->iimismatchCount(di1.p,n,cnt1,maxmm1);
				if (di2.p!=NULL) df3.p->iimismatchCount(di2.p,n,cnt1,maxmm1);
				if (di3.p!=NULL) df3.p->iimismatchCount(di3.p,n,cnt,maxmm);
			}

		}
		else
		{
			// maybe we can implement fast track to calculate it faster in a non-recursive form maybe

			n--;


			if (df0.p!=NULL)
			{
				if (di0.p!=NULL) df0.p->iimismatchCount(di0.p,n,cnt,0);
			}
			if (df1.p!=NULL)
			{
				if (di1.p!=NULL) df1.p->iimismatchCount(di1.p,n,cnt,0);
			}
			if (df2.p!=NULL)
			{
				if (di2.p!=NULL) df2.p->iimismatchCount(di2.p,n,cnt,0);
			}
			if (df3.p!=NULL)
			{
				if (di3.p!=NULL) df3.p->iimismatchCount(di3.p,n,cnt,0);
			}

		}
	}
}


void CLTreef::iimismatchCountGeneral(class CLTreef *iTree, int n, int *cnt ,int maxmm, int alphabetSize) // fills the mismatch count array (howmany sequences with m mismatches to a given set of sequences exist)
{
	fintptr_t *idaughter = iTree->daughter;
	if (n==1)
	{
		double s,si,s0; s=si=s0=0;
		for (int i=0;i<alphabetSize;i++){
			s+=daughter[i].i;
			si+=idaughter[i].i;
			s0+=daughter[i].i*idaughter[i].i;
		}
		*cnt+=s0;
		if (maxmm!=0)
		{
			cnt++;
			*cnt+=(s*si-s0);
		}
	}else{
		if(maxmm!=0){
			n--;
			int *cnt1 = cnt+1;
			int maxmm1 = maxmm-1;

			for (int i=0;i<alphabetSize;i++){
				if (daughter[i].p!=NULL){
					for (int j=0;j<alphabetSize;j++){
						if (idaughter[j].p!=NULL){
							if (i==j){
								daughter[i].p->iimismatchCountGeneral(idaughter[j].p,n,cnt,maxmm,alphabetSize);
							}else{
								daughter[i].p->iimismatchCountGeneral(idaughter[j].p,n,cnt1,maxmm1,alphabetSize);
							}
						}
					}
				}
			}
		}else{
			n--;
			for (int i=0;i<alphabetSize;i++){
				if ((daughter[i].p!=NULL)&&(idaughter[i].p!=NULL)){
					daughter[i].p->iimismatchCountGeneral(idaughter[i].p,n,cnt,maxmm,alphabetSize);
				}
			}
		}
	}
}



double CLTreef::calciScore(int *bid, int L, double *kernel, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->imismatchCount(bid, L, tmpcnt); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}


double CLTreef::calcfScore(int *bid, int L, double *kernel, myFlt *tmpcnt)//calculates the score. tmpcnt is myFlt[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->fmismatchCount(bid, L, tmpcnt); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}
double CLTreef::calcfScore(class CLTreef *iTree, int L, double *kernel,int maxmm, myFlt *tmpcnt)//calculates the score. tmpcnt is myFlt[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->fimismatchCount(iTree, L, tmpcnt, maxmm); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}
double CLTreef::calciScore(class CLTreef *iTree, int L, double *kernel,int maxmm, int *tmpcnt)
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->iimismatchCount(iTree, L, tmpcnt, maxmm); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}




void CLTreef::addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid, int alphabetSize) // used by LList. adds sequences to lis from tree
{
	for (int i=0;i<alphabetSize;i++)
	{
		if (daughter[i].p!=NULL)
		{
			tmpbid[n] =i; 
			if (n==Lm1)
			{
				int frq= daughter[i].i; 
				if ((frq==1)==single)
				{
					list->addSeq(tmpbid, frq); 
				}
			}
			else
			{
				daughter[i].p->addToList(list, n+1, Lm1, single, tmpbid, alphabetSize);
			}
		}
	}
}



/*

void CLTreef::DFSTn(CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize) // without nonEmptyDaughterCnt
{
	int i,j,k;

	for(int bid=0;bid<alphabetSize;bid++)
	{
		if(daughter[bid].p==NULL) continue;
		myFlt nodeival = daughter[bid].f; //LTreeSnodeData *nodei=daughter[bid].node;
//		if (nodei->n==1)
//		{
//		int curnodeid = nodei->seqIDs.i;
		//int **mmprofile=gMMProfile[curnodeid];
		myFlt **mmprofile=gMMProfile0;
		for(int fbid=0;fbid<alphabetSize;fbid++)
		{

			if (bid==fbid)
			{
				for(int i=0;i<listlen;i++)
				{
					if(matchingLmers[i][fbid].node!=NULL)
					{
						if (matchingLmers[i][fbid].node->n==1)
						{
							mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.i]+=nodeival;
						}
						else
						{
							for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
							{
								//if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break;
								mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.p[j]]+=nodeival;
							}
						}
					}
				}
			}
			else
			{
				for(int i=0;i<listlen;i++)
				{
					if(matchingLmers[i][fbid].node!=NULL)
					{
						if (curMismatchCnt[i]<gMAXMM)
						{
							if (matchingLmers[i][fbid].node->n==1)
							{
								mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.i]+=nodeival;
							}
							else
							{
								for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
								{
									//if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break;
									mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.p[j]]+=nodeival;
								}
							}
						}
					}
				}

			}
		}
	}
//	}

}

void CLTreef::DFST( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize) // without nonEmptyDaughterCnt
{
	//LTreeSnodeData **matchingLmers = gDFSlist[pos];
	//int *curMismatchCnt= gDFSMMlist[pos];

	if(pos==gLM1) //LM1 is L-1
	{
		DFSTn(matchingLmers, listlen, curMismatchCnt, alphabetSize); // process the node.
	}
	else
	{
		CLTreeSptr **newlist = gDFSlistT[pos+1];
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0;
		CLTreeSptr **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		for(int bid=0;bid<alphabetSize;bid++)
		{
			if(daughter[bid].p==NULL) continue;
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;
			//int daughter_maxSeqID = daughter[bid].t->maxSeqID;
			for(int fbid=0;fbid<alphabetSize;fbid++) //  foreign bid
			{
				if (bid==fbid)
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].t!=NULL)
						{
							//if (matchingLmers[i][fbid].t->minSeqID >daughter_maxSeqID) continue;
							*newlistnewlistlen=matchingLmers[i][fbid].t->daughter;
							newlistnewlistlen++;
							*newMismatchCntnewlistlen=curMismatchCnt[i];
							newMismatchCntnewlistlen++;
							newlistlen++;
						}
					}
				}
				else
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].t!=NULL)
						{
							if (curMismatchCnt[i]<gMAXMM)
							{
								//if (matchingLmers[i][fbid].t->minSeqID >daughter_maxSeqID) continue;

								*newlistnewlistlen=matchingLmers[i][fbid].t->daughter;
								newlistnewlistlen++;
								*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
								newMismatchCntnewlistlen++;
								newlistlen++;
							}
						}
					}
				}
			}

			if (newlistlen!=0)
			{
				daughter[bid].p->DFST(newlist,newlistlen,newMismatchCnt,pos+1, alphabetSize);
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
*/


void CLTreef::DFSTn(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize) // with nonEmptyDaughterCnt
{
	int bid;
	for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
	{
		bid = this->nonEmptyDaughterIdxs[ibid];
//		if(daughter[bid].p==NULL) continue;
		myFlt nodeival = daughter[bid].f; //LTreeSnodeData *nodei=daughter[bid].node;
//		if (nodei->n==1)
//		{
//		int curnodeid = nodei->seqIDs.i; 
		//int **mmprofile=gMMProfile[curnodeid];
		myFlt **mmprofile=gMMProfile0;

		for(int i=0;i<listlen;i++)
		{
			CLTreeS *imatchingLmer =matchingLmers[i];
			int fbid;
			for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
			{
				fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];

				if (bid==fbid)
				{

					LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
					if (nodej->n==1)
					{
						mmprofile[curMismatchCnt[i]][nodej->seqIDs.i]+=nodeival;
					}
					else
					{
						for(int j=0;j<nodej->n;j++)
						{
//							if (nodej->seqIDs.p[j]>curnodeid) break;
							mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]+=nodeival;
						}
					}
				}
				else
				{
					if (curMismatchCnt[i]<gMAXMM){
						LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
						if (nodej->n==1)
						{
							mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.i]+=nodeival;
						}
						else
						{
							for(int j=0;j<nodej->n;j++)
							{
	//							if (nodej->seqIDs.p[j]>curnodeid) break;
								mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.p[j]]+=nodeival;
							}
						}

					}

				}
			}
		}

	}

}

void CLTreef::DFST( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize) // with nonEmptyDaughterCnt
{
	//LTreeSnodeData **matchingLmers = gDFSlist[pos]; 
	//int *curMismatchCnt= gDFSMMlist[pos];
		
	if(pos==gLM1) //LM1 is L-1
	{
		DFSTn(matchingLmers, listlen, curMismatchCnt, alphabetSize); // process the node.
	}
	else
	{
		CLTreeS **newlist = gDFSlistT[pos+1];
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		CLTreeS **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
//		for(int bid=0;bid<alphabetSize;bid++)
//		{
		int bid;
		for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
		{
			bid = this->nonEmptyDaughterIdxs[ibid];

//			if(daughter[bid].p==NULL) continue;
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;
			//int daughter_maxSeqID = daughter[bid].t->maxSeqID; 

			for(int i=0;i<listlen;i++)
			{
				//CLTreeSptr *imatchingLmer =matchingLmers[i];
				CLTreeS *imatchingLmer =matchingLmers[i];
				int fbid;
				for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
				{
					fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];

					if (bid==fbid){
						CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
						//if (newnode->minSeqID >daughter_maxSeqID) continue;

						*newlistnewlistlen=newnode;
						newlistnewlistlen++;
						*newMismatchCntnewlistlen=curMismatchCnt[i];
						newMismatchCntnewlistlen++;
						newlistlen++;

					} else {
						if (curMismatchCnt[i]<gMAXMM)
						{
							CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
							//if (newnode->minSeqID >daughter_maxSeqID) continue;

							*newlistnewlistlen=newnode;
							newlistnewlistlen++;
							*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
							newMismatchCntnewlistlen++;
							newlistlen++;
						}

					}

				}
			}

			if (newlistlen!=0)
			{
				daughter[bid].p->DFST(newlist,newlistlen,newMismatchCnt,pos+1,alphabetSize);
			}

		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
