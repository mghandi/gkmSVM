/* MLEstimKLmers.h
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
#include "LKTree.h"

//typedef std::tr1::unordered_map<int, double> Mymap;
class CMLEstimKLmers
{
public:

	int K,L; 
	int nrow,ncol; // nrow = C(L,K), ncol= 4^K

	double mu_x, mu_y, s2y,s2x, N0;  // mean and variance of KLmers (muy, s2y) and Kmers (mux,s2x), and normalization factor (N0)

	double *wm;  // w[0]...w[k] : weight for mismatches using the ML alg. 

	int **table; //KLmers frequencies table  
	
	int **wt;  //matrix for weights sum(wij*sj, 0<=j<L) gives the col index

	int *Nmm; ///Nmm[0..K] number of KLmers with distance m of any L-mer
	int *hamdist; //hamming distance for dist(i,j) = hamdit[i^j]
	
	int UseLookupTable; 
	Mymap lookuptable; 

	int UseTree; 
	CLKTree *tree; 

	double estimate(int *bid, double *ocntm=NULL); 
	int convert2int(int *bid, int L);

	CMLEstimKLmers(int L, int K, double *wm, int **KLmersfreq,int **wKLmersTable,  CLKTree *tree=NULL);
	~CMLEstimKLmers(void);

	void calcMean(); 


};

