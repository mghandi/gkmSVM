/* LList.cpp
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

#include "LList.h"


int cmmcnt(int v); // return number of mismatches between x and y. (v = x^y)

CLList::CLList(int L, int maxnLmers, int *extHdist)
{

	if (USEMMLOOKUPTABLE)
	{
		L12THR = L12THRM;
	}
	else
	{
		L12THR = 15;
	}


	twidth= (L+L12THR-1)/L12THR; 
	L12THR = ((L+twidth-1)/twidth);

	int i,j; 

	this->L = L; 
	nmultiples= nsingles= 0; 

	if (twidth<2) 
	{
		twidth=2; //min
	}

	freq = new int [maxnLmers]; 

	table = new int*[twidth]; 
	for(i =0;i<twidth;i++) 
	{
		table[i] = new int [maxnLmers]; 
		for(j=0;j<maxnLmers;j++)
		{
			table[i][j] = freq[j] = 0; 
		}
	}
	
	tmpinttw = new int[twidth]; 


	/// hamdist 
	if (extHdist)
	{
		this->HamDist = this->extHamdist = extHdist; 
	}
	else 
	{

		if (USEMMLOOKUPTABLE)
		{
			HamDist = new int[(1<<(2*L12THR+1))]; 
			int Lm = (2*L12THR+1)/2 +2; 
			if (L<Lm) Lm = L; 
			for(i=0;i<(1<<(2*L12THR+1));i++) 
			{
				HamDist[i] = 0; 
				int mask=3; 
				for(int k=0;k<Lm;k++) 
				{
					if ((mask&i)!=0) HamDist[i]++; 				
					mask <<=2; 
				}
			}
		}
		else
		{
			HamDist =NULL;
		}
	}
	///

	UseLookupTable = (L<=15); // use a look up table for calcScore to speed up. (I used 15 tentatively)
	LookupTableMaxSize = 1<<26;


}


CLList::~CLList(void)
{

	int i; 
	delete []freq; 

	for(i =0;i<twidth;i++) 
	{
		delete []table[i];
	}
	
	delete []table; 
	delete []tmpinttw; 
	if (HamDist!=extHamdist) 
	{
		delete []HamDist; 
	}
	
}

void CLList::mismatchCount2(int *bid, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; for width=2
{
	int B1,B2;
	B1 = convert(bid,0); 
	B2 = convert(bid,1); 

	int *L1, *L2; 
	L1 =table[0]; 
	L2 =table[1]; 
	int *Frq = this->freq; 
		
	int *hamdist=HamDist; 

	int i=this->nsingles; 
	Frq+=i; 
	while(i!=0)
	{	
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)]++; 
		L1++; L2++; 
		i--; 
	}

	i=this->nmultiples; 
	while(i!=0)
	{
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)]+=*Frq; 
		L1++; L2++; Frq++; 
		i--; 
	}
}

void CLList::mismatchCount3(int *bid, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; ; for width=3
{
	int B1,B2,B3; 
	B1 = convert(bid,0); 
	B2 = convert(bid,1); 
	B3 = convert(bid,2); 

	int *L1,*L2,*L3; 
	L1 =table[0]; 
	L2 =table[1]; 
	L3 =table[2]; 
	int *Frq = this->freq; 
		
	int *hamdist=HamDist; 

	int i=nsingles; 
	Frq+=i; 
	while(i!=0)
	{
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)+HAMDIST((*L3)^B3)]++; 
		L1++; L2++; L3++; 
		i--; 
	}

	i=nmultiples; 
	while(i!=0)
	{
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)+HAMDIST((*L3)^B3)]+=*Frq; 
		L1++; L2++; L3++; Frq++; 
		i--; 
	}
}

void CLList::mismatchCount4(int *bid, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; ; for width=4
{
	int B1,B2,B3,B4; 
	B1 = convert(bid,0); 
	B2 = convert(bid,1); 
	B3 = convert(bid,2); 
	B4 = convert(bid,3); 

	int *L1,*L2,*L3,*L4; 
	L1 =table[0]; 
	L2 =table[1]; 
	L3 =table[2]; 
	L4 =table[3]; 
	int *Frq = this->freq; 
		
	int *hamdist=HamDist; 

	int i=nsingles; 
	Frq+=i; 
	while(i!=0)
	{
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)+HAMDIST((*L3)^B3)+HAMDIST((*L4)^B4)]++; 
		L1++; L2++; L3++; L4++; 
		i--; 
	}

	i=nmultiples; 
	while(i!=0)
	{
		cnt[HAMDIST((*L1)^B1)+HAMDIST((*L2)^B2)+HAMDIST((*L3)^B3)+HAMDIST((*L4)^B4)]+=*Frq; 
		L1++; L2++; L3++; L4++; Frq++; 
		i--; 
	}
}

void CLList::mismatchCount(int *bid, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; 
{
	int i;
	for(i=0;i<=L;i++)
	{
		cnt[i]=0;
	}
	
	if (this->twidth==2) return mismatchCount2(bid, cnt);
	if (this->twidth==3) return mismatchCount3(bid, cnt);
	if (this->twidth==4) return mismatchCount4(bid, cnt);

	int *B = tmpinttw; 
	for(i=0;i<twidth;i++)
	{
		B[i] = convert(bid,i); 
	}
		
	int *hamdist=HamDist; 

	int sz = nsingles+nmultiples; 
	for(i=0;i<sz;i++)
	{
		int dist =0; 
		for(int j=0;j<twidth;j++)
		{
			dist += HAMDIST((table[j][i])^B[j]); 
		}
		cnt[dist]+=freq[i]; 
	}

}


void CLList::mismatchCount(int col, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; 
{
	int i;
	for(i=0;i<=L;i++)
	{
		cnt[i]=0;
	}
	
	int *B = tmpinttw; 
	for(i=0;i<twidth;i++)
	{
		B[i] = convert(col,i); 
	}
		
	int *hamdist=HamDist; 

	int sz = nsingles+nmultiples; 
	for(i=0;i<sz;i++)
	{
		int dist =0; 
		for(int j=0;j<twidth;j++)
		{
			dist += HAMDIST((table[j][i])^B[j]); 
		}
		cnt[dist]+=freq[i]; 
	}

}

double CLList::calcScore(int *bid, double *kernel, int *mmcnt) //calculates the score. tmpcnt is int[L+1] 
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


	int i;
	this->mismatchCount(bid, mmcnt); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*mmcnt[i];
	}
	
	
	if (UseLookupTable)
	{
		if (lookuptable.size()<LookupTableMaxSize)
		{
			lookuptable.insert(Mymap::value_type(key, res));
		}
	}
	return res; 
}

void CLList::calcScoreAllLmers(double *kernel, int *mmcnt, double n0, double *cntEst) // calcs freqs for all L-mers and fills cntEst[0..4^L-1]
{  /// check for compatibility with maxwidth !!! -- maybe problematic for L>6 eg
	int *cnt = mmcnt;
	int nLmers = 1<<(2*L);
	for (int B1=0;B1<nLmers;B1++)
	{
		int i; 
		for(i=0;i<=L;i++) cnt[i]=0;

		int *L1=table[0]; 
		int *Frq = this->freq; 
		int *hamdist=HamDist; 

		i=this->nsingles; 
		Frq+=i; 
		while(i!=0)
		{	
			cnt[HAMDIST((*L1)^B1)]++; 
			L1++; 
			i--; 
		}

		i=this->nmultiples; 
		while(i!=0)
		{
			cnt[HAMDIST((*L1)^B1)]+=*Frq; 
			L1++; Frq++; 
			i--; 
		}

		double res = n0; 
		for(i=0;i<=L;i++)
		{
			res+=kernel[i]*cnt[i];
		}

		cntEst[B1]=res;
	}
}


void CLList::calcScoreAllLmersAdd(double *kernel, int *mmcnt, double n0, double *cntEst, double weight) // calcs freqs for all L-mers and adds to cntEst[0..4^L-1]
{
	int *cnt = mmcnt;
	int nLmers = 1<<(2*L);
	for (int B1=0;B1<nLmers;B1++)
	{
		int i; 
		for(i=0;i<=L;i++) cnt[i]=0;
		mismatchCount(B1,cnt); 
		double res = n0; 
		for(i=0;i<=L;i++)
		{
			res+=kernel[i]*cnt[i];
		}
		cntEst[B1]+=res*weight;
	}
}


void CLList::addFromLTree(class CLTree *tree)
{
	int *tmpbid = new int[L]; 
	tree->addToList(this,0, L-1, 1, tmpbid);  // adds L-mers with cnt=1
	tree->addToList(this,0, L-1, 0, tmpbid); // adds L-mers with cnt>1
	delete []tmpbid; 
}

void CLList::addSeq(int *bid, int cnt) // adds one sequence : for efficiency, first add sequences with cnt=1, then add those with cnt>1
{
	int pos = nmultiples+nsingles;

	for(int i=0;i<twidth; i++)
	{
		int si = convert(bid, i);
		table[i][pos]= si; 
	}
	freq[pos]= cnt; 

	if ((nmultiples>0)||(cnt>1))
	{
		nmultiples++; 
	}
	else
	{
		nsingles++; 
	}	
}


int CLList::convert(int *bid, int i)
{
/*
    int si= i*L12THR; 
	int ei= si+L12THR; 
	if (ei>L)
	{
		ei = L; 
	}
	int res = 0; 
	while(si<ei)
	{
		res = (res<<2)+bid[si];
		si++; 
	}
	return res; 
 */
    
    int ei= L-i*L12THR; 
	int si= ei-L12THR; 
	if (si<0)
	{
		si = 0; 
	}
	int res = 0; 
	while(si<ei)
	{
		res = (res<<2)+bid[si];
		si++; 
	}
	return res; 
}

int CLList::convert(int col, int i)
{
	return (col>>(2*i*L12THR))&((1<<(2*L12THR))-1);
}

char *CLList::convertInt2Str(int col, char *str, int L) // returns L-mer for idx=col
{
    str[(L-1)-0] = ::globalConverter.icidx[col%4]; 
    col >>= 2; 
    
    for(int i=1;i<L;i++)
    {
        str[(L-1)-i] = ::globalConverter.icidx[col%4]; 
        col >>= 2; 
    }
    str[L]=0; 
    return str; 
}


void CLList::clear()
{
	nmultiples= nsingles= 0; 
}

double CLList::calcInnerProd(CLList *L2, double *c, int *mmcnt)
{
	if (L<=L12THR) return calcInnerProd1(L2, c, mmcnt); 
	if (this->twidth==2)  return calcInnerProd2(L2, c, mmcnt);

	int N1 = nmultiples+nsingles; 
	int N2 = L2->nmultiples+L2->nsingles; 

	int *freq2 = L2->freq; 
	int **table2 = L2->table;
	
	int *hamdist=HamDist; 


	int i,j; 
	
	for (i=0;i<=L;i++)
	{
		mmcnt[i]=0;
	}
	
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			int dist =0; 
			for(int k=0;k<twidth;k++)
			{
				dist += HAMDIST((table[k][i])^(table2[k][j])); 
			}
			mmcnt[dist]+=freq[i]*freq2[j]; 
		}
	}

	double s = 0; 
	for (i=0;i<=L;i++)
	{
		s+=mmcnt[i]*c[i];
	}

	return(s); 

}





double CLList::calcInnerProd1(CLList *L2, double *c, int *mmcnt)
{
	int N1 = nmultiples+nsingles; 
	int N2 = L2->nmultiples+L2->nsingles; 

	int N1s = nsingles; 
	int N2s = L2->nsingles; 

	int *freq2 = L2->freq; 
	
	int *hamdist=HamDist; 


	int i,j; 
	
	for (i=0;i<=L;i++)
	{
		mmcnt[i]=0;
	}
	
	int *table0 = table[0];
	for(i=0;i<N1s;i++)
	{
		int table0i =*table0; table0++; 
		int *table20 = L2->table[0];
		for(j=0;j<N2s;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))]++; 
			table20++; 
		}

		for(j=N2s;j<N2;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))]+=freq2[j]; 
			table20++; 
		}

	}
	for(i=N1s;i<N1;i++)
	{
		int table0i =*table0; table0++; 
		int *table20 = L2->table[0];
		int freqi = freq[i]; 
		for(j=0;j<N2s;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))]+=freqi; 
			table20++; 
		}

		for(j=N2s;j<N2;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))]+=freqi*freq2[j]; 
			table20++; 
		}

	}


	double s = 0; 
	for (i=0;i<=L;i++)
	{
		s+=mmcnt[i]*c[i];
	}

	return(s); 

}



double CLList::calcInnerProd2(CLList *L2, double *c, int *mmcnt)
{
	int N1 = nmultiples+nsingles; 
	int N2 = L2->nmultiples+L2->nsingles; 

	int N1s = nsingles; 
	int N2s = L2->nsingles; 

	int *freq2 = L2->freq; 
	
	int *hamdist=HamDist; 


	int i,j; 
	
	for (i=0;i<=L;i++)
	{
		mmcnt[i]=0;
	}
	
	int *table0 = table[0];
	int *table1 = table[1];
	for(i=0;i<N1s;i++)
	{
		int table0i =*table0; table0++; 
		int table1i =*table1; table1++; 
		int *table20 = L2->table[0];
		int *table21 = L2->table[1];
		for(j=0;j<N2s;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))+HAMDIST(table1i^(*table21))]++; 
			table20++; 
			table21++; 
		}

		for(j=N2s;j<N2;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))+HAMDIST(table1i^(*table21))]+=freq2[j]; 
			table20++; 
			table21++; 
		}

	}
	for(i=N1s;i<N1;i++)
	{
		int table0i =*table0; table0++; 
		int table1i =*table1; table1++; 
		int *table20 = L2->table[0];
		int *table21 = L2->table[1];
		int freqi = freq[i]; 
		for(j=0;j<N2s;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))+HAMDIST(table1i^(*table21))]+=freqi; 
			table20++; 
			table21++; 
		}

		for(j=N2s;j<N2;j++)
		{
			mmcnt[HAMDIST(table0i^(*table20))+HAMDIST(table1i^(*table21))]+=freqi*freq2[j]; 
			table20++; 			
			table21++; 
		}

	}


	double s = 0; 
	for (i=0;i<=L;i++)
	{
		s+=mmcnt[i]*c[i];
	}

	return(s); 

}

int cmmcnt(int v) // return number of mismatches between x and y. (v = x^y)
{
	v = (v | (v>>1))& 1431655765; //1431655765=(01010101010101010101010101010101) in base 2 // I used this to reduce it to set bit count problem
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	return( ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24); // count
}
