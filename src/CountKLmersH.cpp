/* CountKLmersH.cpp
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

#include "CountKLmersH.h"


int mypower(int m, int n)
{
	int r = 1; 
	for(int i=0;i<n; i++)
	{
		r*=m; 
	}
	return r; 
}

CCountKLmersH::CCountKLmersH(int L, int K)
{

	this->K = K; 
	this->L = L; 

	// initialize table:  table [nh=L-K,..,L][posidx=0,..,C(L,nh)-1][validx=0,..,(b-1)^(L-nh)-1]
	table = new int **[L+1];
	ntable = new double **[L+1];
	for(int nh=L-K; nh<=L;nh++)
	{
		int mxPosIdxi = Combinations(L, nh); 
		table[nh] = new int *[mxPosIdxi]; 
		ntable[nh] = new double *[mxPosIdxi]; 

		int mxValIdxi = mypower(3, L-nh); //3 = b-1
		for(int pidxi=0;pidxi<mxPosIdxi; pidxi++)
		{
			table[nh][pidxi] = new int [mxValIdxi];
			ntable[nh][pidxi] = new double [mxValIdxi];

			for(int vidxi=0;vidxi<mxValIdxi;vidxi++)
			{
				table[nh][pidxi][vidxi]=0; 
			}
		}
	}
	
	for(int n=0;n<=L;n++)
		for(int r=0;r<=L;r++)
			nCr[n][r]=Combinations(n,r); 
	
}

void CCountKLmersH::btadd(int k, int vpar, int nh, int posidx, int validx ,int *seq)
{
	if (k==this->L)
	{
		table[nh][posidx][validx]+=vpar; 
		return;
	}
	//v'[k]=h
	btadd(k+1,vpar,nh+1,posidx+tcombinations(k, nh+1),validx,seq); 


	if ((L-K)<=(L-k-1 +nh)) // if there is enough positions available for 'h's down the road 
	{
		validx*=3; //(b-1) 

		//posidx+=Combinations(k, nh); //tabulate this later for speed
		
		//v'[k]=j
		for(int j=1;j<=3;j++)
		{
			if (seq[k]<j) btadd(k+1,vpar,nh, posidx, validx, seq);
			if (seq[k]==j) btadd(k+1,-seq[k]*vpar,nh,posidx, validx, seq); 
			validx++;
		}
	}
}

double CCountKLmersH::btest(int k, int vpar, int nh, int posidx, int validx ,int *seq)
{
	double res =0; 
	if (k==this->L)
	{
		return ntable[nh][posidx][validx]*vpar; 
	}
	//v'[k]=h
	res+=btest(k+1,vpar,nh+1,posidx+tcombinations(k, nh+1),validx,seq); 

	if (k<(K +nh)) //((L-K)<=(L-k-1 +nh)) // if there is enough positions available for 'h's down the road 
	{
		validx*=3; //(b-1) 

		//posidx+=Combinations(k, nh); //tabulate this later for speed
		//v'[k]=j
		for(int j=1;j<=3;j++)
		{
			if (seq[k]<j) res+=btest(k+1,vpar,nh, posidx, validx,seq);
			if (seq[k]==j) res+=btest(k+1,-seq[k]*vpar,nh,posidx, validx, seq); 
			validx++;
		}

	}
	return res; 
}


void CCountKLmersH::addSequence(int *seqBID, int size)
{
	int *si = seqBID; 
	for(int i=0;i<=(size-L); i++)
	{
		btadd(0,1,0,0,0,si);
		si++;  
	}
}

void CCountKLmersH::normalize()
{
	for (int n=0;n<=K; n++)
	{
		btnorm(0,0,1,n); 
	}
}


void CCountKLmersH::print(FILE *f)
{
	for (int n=0;n<=K; n++)
	{
		//btprint(0,0,n,f); 
		btprint(0,0,1,n,f); // normalized
		
	}
}

void CCountKLmersH::btnorm(int i, int valpar, int norm2par, int n) // calc normalized table
{
	if (i==n) 
	{
		int nh = L-n; 
		int mxPosIdxi = Combinations(L, nh);  
		double bln = 1.0*(1<<(2*(K-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		//double bln = 1.0*(1<<(2*(L-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		for(int j=0;j<mxPosIdxi;j++)
		{
			ntable[nh][j][valpar] = table[nh][j][valpar]/ (bln * norm2par); 
		}
		return;
	}

	valpar = valpar*3; //(b-1)
	for(int vj=1; vj<=3; vj++) //(b-1)
	{
		btnorm(i+1, valpar+vj-1, norm2par*vj*(vj+1), n); 
	}
	
}

void CCountKLmersH::btprint(int i, int valpar, int norm2par, int n, FILE *f)
{
	if (i==n) 
	{
		int nh = L-n; 
		int mxPosIdxi = Combinations(L, nh);  
		double bln = 1.0*(1<<(2*(K-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		//double bln = 1.0*(1<<(2*(L-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		for(int j=0;j<mxPosIdxi;j++)
		{
			//ntable[nh][j][valpar] = table[nh][j][valpar]/ (bln * norm2par); 
			//fprintf(f, "\t%d", table[nh][j][valpar]);
			fprintf(f, "\t%e", table[nh][j][valpar]/sqrt(bln * norm2par));
		}
		return;
	}

	valpar = valpar*3; //(b-1)
	for(int vj=1; vj<=3; vj++) //(b-1)
	{
		btprint(i+1, valpar+vj-1, norm2par*vj*(vj+1),n,  f); 
	}
}



void CCountKLmersH::btprint(int i, int valpar,  int n, FILE *f)
{
	if (i==n) 
	{
		int nh = L-n; 
		int mxPosIdxi = Combinations(L, nh);  
		//double bln = 1.0*(1<<(2*(K-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		//double bln = 1.0*(1<<(2*(L-n))); // b^(L-n)  or alternatively use b^(K-n) to scale everything by b^(L-K)
		for(int j=0;j<mxPosIdxi;j++)
		{
			//ntable[nh][j][valpar] = table[nh][j][valpar]/ (bln * norm2par); 
			fprintf(f, "\t%d", table[nh][j][valpar]);
		}
		return;
	}

	valpar = valpar*3; //(b-1)
	for(int vj=1; vj<=3; vj++) //(b-1)
	{
		btprint(i+1, valpar+vj-1, n, f); 
	}
}



double CCountKLmersH::estimate(int *bid) // given the sequence of length L, what is the ML estimate ? 
{
	return btest(0,1,0,0,0,bid);
}

CCountKLmersH::~CCountKLmersH(void)
{

	// delete table:  table [nh=L-K,..,L][posidx=0,..,C(L,nh)-1][validx=0,..,(b-1)^(L-nh)-1]
	for(int nh=L-K; nh<=L;nh++)
	{
		int mxPosIdxi = Combinations(L, nh); 

		for(int pidxi=0;pidxi<mxPosIdxi; pidxi++)
		{
			delete []table[nh][pidxi]; 
			delete []ntable[nh][pidxi]; 
		}
		delete []table[nh]; 
		delete []ntable[nh]; 
	}
	delete []table; 
	delete []ntable; 
}


int CCountKLmersH::tcombinations(int n, int r)
{
	if (n<r) return 0; 
	return nCr[n][r];
}
