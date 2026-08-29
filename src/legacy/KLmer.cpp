/* KLmer.cpp
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

#include "KLmer.h"
#include "global.h" 

CKLmer::CKLmer(int L, int K)
{
	this->K = K; 
	this->L = L; 
	seq = new char[L+1]; 
	
	seq[L]=0;  
	iseq = new int[L]; 


}

void CKLmer::readKLmer(char *s)
{
	for (int i=0;i<L;i++)
	{
		seq[i] = toupper(s[i]); 
		iseq[i] = 1<<::globalConverter.cidx[s[i]];
		if (seq[i]=='.') {iseq[i]=15; }
	}
}

int CKLmer::countHits(char *s, int size)
{
	int N = size-L+1; 
	int i,j;
	int cnt = 0; 
	for (i=0;i<N;i++)
	{
		for(j=0;j<L;j++)
		{
			if ((seq[j]!='.')&&(s[j]!=seq[j]))
			{
				break; 
			}
		}
		if (j==L) cnt++;
		s++;
	}
	return cnt; 
}

int CKLmer::commonKMerCnt(CKLmer *klmerj)
{
	int ndot=0; 
	int *jseq = klmerj->iseq; 
	for(int i=0;i<L;i++)
	{
		if (((iseq[i])&(jseq[i]))==0)
		{
			return 0; 
		}
		if (((iseq[i])&(jseq[i]))==15)
		{
			ndot++; 
		}
	}
	return 1<<(2*ndot); 
}

static const double fbgi[]={0,log(5.0/4),-log(5.0/4),0,-log(5.0/4), 0,0, 0,log(5.0/4),0,0,0,0,0,0,0};  // A,T: +log(5/4) and C,G: -log(5/4)
double CKLmer::calcfbg()
{
	double z = 0; 
	for(int i=0;i<L; i++)
	{
		z+=fbgi[iseq[i]]; 
	}
	return z; 
}

CKLmer::~CKLmer(void)
{
	delete []seq; 
	delete []iseq; 
}
