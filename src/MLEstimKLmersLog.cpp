/* MLEstimKLmersLog.cpp
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

#include "MLEstimKLmersLog.h"
#include "global.h"



CMLEstimKLmersLog::CMLEstimKLmersLog(int L, int K, double *wm, int **KLmersfreq, int **wKLmersTable)
{
	this->L = L; 
	this->K = K; 
	this->wt = wKLmersTable;  //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	ncol = 1<<(2*K); // nrow = C(L,K), ncol= 4^K
	nrow = Combinations(L,K); 
	
	this->wm = new double [K+1]; 

	int i,j; 
	for(i=0;i<=K;i++)
	{
		this->wm[i] = wm[i]; 
	}

	//this->table = KLmersfreq; 
	table = new double *[nrow]; 
	
	for(i=0;i<nrow;i++)
	{
		table[i] = new double[ncol]; 
	
		for(j=0;j<ncol;j++)
		{
			table[i][j]=log(1.0+KLmersfreq[i][j]); 
		}
	}


	N0 = 0; 
	Nmm = new int[K+1]; 
	int p3 = 1; // 3^m
	for(int m = 0; m<=K; m++)
	{

		Nmm[m] = Combinations(L, K-m)*Combinations(L-K+m, m)*p3;
		p3 = 3*p3; 

		N0+=wm[m]*Nmm[m]*pow(4.0, L-K); 
	}

	/// hamdist 
	hamdist = new int[2*ncol]; 
	for(i=0;i<2*ncol;i++) 
	{
		hamdist[i] = 0; 
		int mask=3; 
		for(int k=0;k<L;k++) 
		{
			if ((mask&i)!=0) hamdist[i]++; 				
			mask <<=2; 
		}
	}

	///

	this->calcMean(); 
}

CMLEstimKLmersLog::~CMLEstimKLmersLog(void)
{
	delete []wm; 
	delete []Nmm; 
	delete []hamdist; 

	for(int i=0;i<nrow;i++)
	{
		delete []table[i];
	}
	delete []table;
}


void CMLEstimKLmersLog::calcMean()
{
	double m = 0; 
	double **table = this->table;
	int i,j;
	for( i=0;i<nrow;i++)
	{
		for( j=0;j<ncol;j++)
		{
			m +=table[i][j];
		}
	}

	m = m/(nrow*ncol); 
	this->mu_y = m; 
	this->mu_x = m/(1<<(2*(L-K))); 

	for( i=0;i<nrow;i++)
	{
		for( j=0;j<ncol;j++)
		{
			table[i][j]-=m;
		}
	}

}

double CMLEstimKLmersLog::estimate(int *bid) // given the sequence of length L, what is the ML estimate ? 
{


	int *hamdist = this->hamdist; 
	int i,j,m;
	int dif = 0; 

	//int *cntm = new int[K+1]; // sum yi for hdist=m
	double *dcntm = new double [K+1]; // sum yi for hdist=m
	for(i=0;i<=K;i++) 
	{
		//cntm[i] = 0; 
		dcntm[i] = 0; 
	}

	for(i=0; i<nrow; i++) 
	{

		int ci = 0; 
		for(j=0;j<L;j++) 
		{
			ci+=wt[i][j]*bid[j]; // binary representation of the selected K bits 
		}
		double *tablei = table[i]; 
		for(j=0; j<ncol; j++) 
		{
			dif = hamdist[j^ci]; 
			dcntm[dif]+=tablei[j]; 
		}

	}

	double wsum = 0; 
	for(m=0;m<=K;m++) 
	{
		//dcntm[m] +=	cntm[m]; 
		//wsum+=wm[m] * (dcntm[m] -this->mu_y * Nmm[m]); 
		wsum+=wm[m] * dcntm[m]; 
	}

	//double mle = mu_x+(wsum/N0); 
	double mle = wsum/N0; 

	
	//delete []cntm; 
	delete []dcntm; 
	

	return mle; 
}
