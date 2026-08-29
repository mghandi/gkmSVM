/* CountKLmersH.h
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
class CCountKLmersH
{
public:
	CCountKLmersH(int L, int K);
	~CCountKLmersH(void);
	
	int ***table; // keeps the counts: table [nh=L-K,..,L][posidx=0,..,C(L,nh)-1][validx=0,..,(b-1)^(L-nh)-1]
	double ***ntable; // normalized table 

	int K,L; 
	
	void addSequence(int *seqBID, int size); 
	void normalize(); // calc normalized table // call this after training and before estimate 
	
	void print(FILE *f); // prints counts to file 

	double estimate(int *bid); // given the sequence of length L, what is the ML estimate ? 
	
	

/*	int nrow,ncol; // nrow = C(L,K), ncol= 4^K
	int **w; //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	


	
	
	char *convertCol2KmerString(int col, char *sKmer); // returns K-mer for idx=col
	char *convertRow2KLmerString(int row, char *sKmer, char *sKLmer); // maps Kmer to KLmer for idx=row
	int *convertCol2bid(int col, int *bid); // returns K-mer for idx=col
 
	void calcPosNegRatio(int **negTable); // replaces the counts by 1000000*(Np-Nn)/(Np+Nn). neg table should be of exactly same size
*/
private: 
//	int *wdata; //allocated memory for w
//	int fillwij(int pos, int nfilled, int idx, int *partial);

	void btadd(int k, int vpar, int nh, int posidx, int validx ,int *seq); // calcs z = P^{\top}x
	double btest(int k, int vpar, int nh, int posidx, int validx ,int *seq); // calcs x^{\hat} = Pz
	void btnorm(int i, int valpar, int norm2par, int n); // calc normalized table
	void btprint(int i, int valpar,  int n, FILE *f); // prints table to file 
	void btprint(int i, int valpar, int norm2par, int n, FILE *f); // prints normalized table to file 


	int nCr[50][50];
	int tcombinations(int n, int r); 
};



