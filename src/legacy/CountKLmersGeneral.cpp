/* CountKLmersGeneral.cpp
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

#include "CountKLmersGeneral.h"
#include "global.h"
#include "KLmer.h"


CCountKLmersGeneral::CCountKLmersGeneral(int L, int K)
{

	this->K = K; 
	this->L = L; 

	ncol = 1<<(2*K); // nrow = C(L,K), ncol= 4^K
	nrow = Combinations(L,K); 

	wdata = new int[nrow*L]; 
	w = new int *[nrow]; //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	
	table = new int *[nrow]; 
	
	int i,j; 
	
	for(i=0;i<nrow;i++)
	{
		table[i] = new int[ncol]; 
	
		for(j=0;j<ncol;j++)
		{
			table[i][j]=0; 
		}

		w[i]=wdata+i*L; 
	
		for(j=0;j<L;j++)
		{
			w[i][j]=0; 
		}
	}


	int *partialw = new int[L]; 
	
	i = fillwij(0,0,0, partialw);
	
	delete []partialw; 


	
}

void CCountKLmersGeneral::addSequence(int *seqBID, int size)
{

	int *si = seqBID; 
	for(int i=0;i<=(size-L); i++)
	{
		for(int irow =0; irow<nrow; irow++)
		{
			int jcol = 0;
			int *wi = w[irow]; 
			for(int j=0;j<L;j++)
			{
				jcol += wi[j]*si[j]; 
			}
			table[irow][jcol]++; 
		}
		si++;  
	}
}

char *CCountKLmersGeneral::convertCol2KmerString(int col, char *sKmer) // returns K-mer for idx=col
{		
	for(int i=0;i<K;i++)
	{
		sKmer[i] = ::globalConverter.icidx[col%4]; 
		col >>= 2; 
	}
	sKmer[K]=0; 
	return sKmer; 
}

int *CCountKLmersGeneral::convertCol2bid(int col, int *bid) // returns K-mer for idx=col
{		
	for(int i=0;i<K;i++)
	{
		bid[i] = col%4; 
		col >>= 2; 
	}
	return bid; 
}

char *CCountKLmersGeneral::convertRow2KLmerString(int row, char *sKmer, char *sKLmer) // maps Kmer to KLmer for idx=row
{
	int k=0; 
	for(int j=0;j<L;j++)
	{
		if (w[row][j]==0)
		{
			sKLmer[j]='.'; 
		}
		else
		{
			sKLmer[j]=sKmer[k]; 
			k++; 
		}
	}
	sKLmer[L]=0; 
	return sKLmer; 
}


int CCountKLmersGeneral::fillwij(int pos, int nfilled, int idx, int *partial)
{
	if (pos==L)
	{
		for(int j=0;j<L;j++)
		{
			w[idx][j]=partial[j]; 
		}
		return idx+1; 
	}

	if ((L-pos)>(K-nfilled)) // gap
	{
		partial[pos] = 0;
		idx = fillwij(pos+1, nfilled, idx, partial); 
	}
	if (nfilled<K) // known bit
	{
		partial[pos] = 1<<(2*nfilled); 
		idx = fillwij(pos+1, nfilled+1, idx, partial); 
	}
	return idx; 
}


CCountKLmersGeneral::~CCountKLmersGeneral(void)
{
	int i; 
	delete []wdata; 
	delete []w; 
	for(i=0;i<nrow;i++)
	{
		delete []table[i];
	}
	delete []table;

}


CLKTree *CCountKLmersGeneral::generateFreqTree() // generates a tree structure for gapped-kmers frequencies
{
	CLKTree *tree = new CLKTree(); 
	int *bid = new int[L]; 
	int *bidj = new int[K]; 
	int i,j; 
	for(j=0;j<ncol;j++)
	{
		this->convertCol2bid(j, bidj); 
		for(i=0;i<nrow;i++)
		{
			int k=0; 
			for(int t=0;t<L;t++)
			{
				if (w[i][t]==0)
				{
					bid[t]=4; //gap
				}
				else
				{
					bid[t]=bidj[k]; 
					k++; 
				}
			}
			tree->addSeq(bid, L, table[i][j]); 
		}
	}
	delete []bid;
	delete []bidj;

	return tree; 
}



void CCountKLmersGeneral::calcPosNegRatio(int **negTable) // replaces the counts by floor(1000000*(Np-Nn)/(Np+Nn)). neg table should be of exactly same size
{
	int i,j; 
	for(i=0;i<nrow;i++)
	{
		for(j=0;j<ncol;j++)
		{
			int Np = table[i][j];
			int Nn = negTable[i][j];
			if ((Np+Nn)==0)
			{
				table[i][j]=0;
			}
			else
			{
				table[i][j]=(int)(1000000*(Np-Nn)/(1.0*(Np+Nn))); 
			}
		}
	}
}
