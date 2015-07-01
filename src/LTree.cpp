/* LTree.cpp
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

#include "LTree.h"
#include "Sequence.h"

CLTree::CLTree(void)
{
//	daughter[0]=daughter[1]=daughter[2]=daughter[3]=NULL;
	for(int i=0;i<MAX_ALPHABET_SIZE;i++){
		daughter[i]=NULL;
	}

}

CLTree::~CLTree(void)
{
}

void CLTree::addSeq(int *bid, int n, int cnt)
{
	if (n==1)
	{
		this->daughter[*bid] = (CLTree *)((intptr_t)(daughter[*bid])+cnt); 
	}
	else
	{
		if (this->daughter[*bid] == NULL)
		{
			this->daughter[*bid] = new CLTree(); 
		}
		daughter[*bid]->addSeq(bid+1, n-1, cnt); 
	}
}

int CLTree::count(int *bid, int n) //returns the number of times the sequence is in found in the tree
{
	if (n==1)
	{
		return (int)(intptr_t)(daughter[*bid]);
	}
	else
	{
		if (this->daughter[*bid] == NULL)
		{
			return 0; 
		}
		return daughter[*bid]->count(bid+1, n-1); 
	}
}

int CLTree::leavesCount(int withMultiplicity, int n)  //returns the number of sequences in the tree.  //call with n=L from outside
{
	int nleaves = 0; 
	for (int i=0;i<MAX_ALPHABET_SIZE;i++)
	{
		if (daughter[i]!=NULL)
		{
			if (n==1)
			{
				if (withMultiplicity)
				{
					int frq= (int)(intptr_t)(daughter[i]);
					nleaves +=frq; 
				}
				else 
				{
					nleaves++; 					
				}
			}
			else
			{
				nleaves+=daughter[i]->leavesCount(withMultiplicity, n-1); 
			}
		}
	}
	return nleaves; 
}

void CLTree::deleteTree(int n)
{
	if (n>1)
	{
		for (int i=0;i<MAX_ALPHABET_SIZE;i++)
		{
			if (daughter[i]!=NULL)
			{
				daughter[i]->deleteTree(n-1);
				delete daughter[i];
			}
		}
	}
}

int CLTree::addSequence(int *bid, int n, int L)  //adds all the L-subseqs 
{
	n = n-L+1;
	for(int i=0;i<n;i++)
	{
		addSeq(bid,L,1);
		bid++;
	}
	return n; 
}

void  CLTree::addSequences(char *FsaFileName, int L, int maxSequenceLength, int addrevcompl, int numberOfCVPartitions, int selectPartitionNumber) // adds all the sequence from a FSA file 
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

void CLTree::mismatchCount(int *bid, int n, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=(intptr_t)(daughter[i]);
		cnt++; 
		i=(i+1)&3;  //(i+1)%4
		*cnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*cnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*cnt+=(intptr_t)(daughter[i]);
	}
	else
	{
		n--; 
		bid++; 
		
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);

		cnt++; 
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);

	}
}

void CLTree::mismatchCount(int *bid, int n, int *cnt, int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=(intptr_t)(daughter[i]);
		if (maxmm!=0)
		{
			cnt++; 
			i=(i+1)&3;  //(i+1)%4
			*cnt+=(intptr_t)(daughter[i]);
			i=(i+1)&3; 
			*cnt+=(intptr_t)(daughter[i]);
			i=(i+1)&3; 
			*cnt+=(intptr_t)(daughter[i]);
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--; 
			bid++; 
		
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);

			maxmm--;
			cnt++; 
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
		}
		else
		{   // fast track for maxmm==0
			CLTree *cur = this;  
			while (--n)
			{
				if ((cur = cur->daughter[*bid++])==NULL) return;  
			}
			*cnt+=(intptr_t)(cur->daughter[*bid]); 
		}
	}
}

double CLTree::calcScore(int *bid, int L, double *kernel, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->mismatchCount(bid, L, tmpcnt); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}

double CLTree::calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->mismatchCount(bid, L, tmpcnt,maxmm); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}

double CLTree::calcScore(int *bid,int *bidrc, int L, int slen, double *kernel, int maxmm, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	slen = slen-L+1; 
	for(i=0;i<slen;i++)
	{
		this->mismatchCount(bid, L, tmpcnt,maxmm); 
		bid++; 
	}
	if (bidrc!=NULL)
	{
		for(i=0;i<slen;i++)
		{
			this->mismatchCount(bidrc, L, tmpcnt,maxmm); 
			bidrc++; 
		}
	}
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}


void CLTree::addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid) // used by LList. adds sequences to lis from tree 
{
	for (int i=0;i<MAX_ALPHABET_SIZE;i++)
	{
		if (daughter[i]!=NULL)
		{
			tmpbid[n] =i; 
			if (n==Lm1)
			{
				int frq= (int)(intptr_t)(daughter[i]);
				if ((frq==1)==single)
				{
					list->addSeq(tmpbid, frq); 
				}
			}
			else
			{
				daughter[i]->addToList(list, n+1, Lm1, single, tmpbid); 
			}
		}
	}
}


CLTreeMemorize::CLTreeMemorize(int UseLookupTable, int LookupTableMaxSize, class CLTree *tree)
{
	this->tree = tree; 
	this->UseLookupTable = UseLookupTable; // = (L<=15); // use a look up table for calcScore to speed up. (I used 15 tentatively)
	this->LookupTableMaxSize = LookupTableMaxSize; // = 1<<26;
}

CLTreeMemorize::~CLTreeMemorize(void)
{
	this->lookuptable.clear();
}

double CLTreeMemorize::calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int key=0; 
	if (UseLookupTable)
	{
		key = convert2int(bid,L);
		Mymap::iterator it; 
		it = lookuptable.find(key); 
		if (it!=lookuptable.end())
		{
			return it->second;
		}
	}

	double res = this->tree->calcScore(bid, L,kernel, maxmm, tmpcnt); 
	
	if (UseLookupTable)
	{
		if (lookuptable.size()<LookupTableMaxSize)
		{
			lookuptable.insert(Mymap::value_type(key, res));
		}
	}
	return res; 
}

