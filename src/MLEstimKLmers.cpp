/* MLEstimKLmers.cpp
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

#include "MLEstimKLmers.h"
#include "global.h"


CMLEstimKLmers::CMLEstimKLmers(int L, int K, double *wm, int **KLmersfreq, int **wKLmersTable, CLKTree *tree)
{
	this->L = L; 
	this->K = K; 
	this->wt = wKLmersTable;  //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	ncol = 1<<(2*K); // nrow = C(L,K), ncol= 4^K
	nrow = Combinations(L,K); 
	
	UseLookupTable= 0;
	if (L<=15) 
	{
		UseLookupTable = 1; 
	}

	this->tree = tree;
	this->UseTree = (tree!=NULL); 

	this->wm = new double [K+1]; 

	int i; 
	for(i=0;i<=K;i++)
	{
		this->wm[i] = wm[i]; 
	}

	this->table = KLmersfreq; 

	N0 = 0; 
	Nmm = new int[K+1]; 
	int p3 = 1; // 3^m
	for(int m = 0; m<=K; m++)
	{

		//Nmm[m] = Combinations(L, K-m)*Combinations(L-K+m, m)*p3;
		Nmm[m] = Combinations(L, K)*Combinations(K, m)*p3;
		p3 = 3*p3; 

		N0+=wm[m]*Nmm[m]; 
	}
	N0 *=pow(4.0, L-K);  // if comment this line, it will give xj, assuming yi 's are averages, if uncomment this line, it will give xj's assuming yi 's are sum of aij*xj. 


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

CMLEstimKLmers::~CMLEstimKLmers(void)
{
	delete []wm; 
	delete []Nmm; 
	delete []hamdist; 
}


void CMLEstimKLmers::calcMean()
{
	int s = 0; 
	double m = 0; 
	double mm=0;
	int **table = this->table;

	for(int i=0;i<nrow;i++)
	{
		for(int j=0;j<ncol;j++)
		{
			s +=table[i][j];
			mm +=table[i][j]*table[i][j];
		}
		if (s>200000000)
		{
			s -= 200000000; 
			m += 200000000.0; 
		}
	}

	m = (m+s)/(nrow*ncol); 
	this->mu_y = m; 
	this->s2y = mm/(nrow*ncol) -m*m;
	this->mu_x = m/(1<<(2*(L-K))); 
	this->s2x =0.5*s2y/(1<<(2*(L-K)));  // just a rough estimate 
	
    sprintf(globtmpstr,"\n mu_y= %lf\n s2y= %lf\nmu_x= %lf\n", mu_y,s2y,mu_x);Printf(globtmpstr);
}

int CMLEstimKLmers::convert2int(int *bid, int L)
{
	int s = 0; 
	for(int i=0;i<L;i++)
	{
		s<<=2; 
		s+=*bid; 
		bid++;
	}
	return s; 
}


double CMLEstimKLmers::estimate(int *bid, double *ocntm) // given the sequence of length L, what is the ML estimate ? 
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
/*
	 convert
	 c1.insert(Mymap::value_type(76251,3.16));
	 c1.insert(Mymap::value_type(76211,2.6));
	 c1.insert(Mymap::value_type(751, 23.16));
	 
	 printf("%d",c1.find(76211)->second);
	 */

	int *hamdist = this->hamdist; 
	int i,j,m;
	int dif = 0; 

	double *dcntm = ocntm; // sum yi for hdist=m
	if (ocntm==NULL) dcntm=new double[K+1]; // sum yi for hdist=m
	intptr_t *cntm=new intptr_t[K+1]; // sum yi for hdist=m
	
	for(i=0;i<=K;i++) 
	{
		cntm[i] = 0; 
		dcntm[i] = 0; 
	}


	if (this->UseTree)
	{
		tree->mismatchCount(bid,L,cntm, dcntm);
	}
	else
	{

		for(i=0; i<nrow; i++) 
		{

			int ci = 0; 
			for(j=0;j<L;j++) 
			{
				ci+=wt[i][j]*bid[j]; // binary representation of the selected K bits 
			}
			int *tablei = table[i]; 
			for(j=0; j<ncol; j++) 
			{
				dif = hamdist[j^ci]; 
				cntm[dif]+=tablei[j]; 
			}

			for(m=0;m<=K;m++) 
			{
				if (cntm[m]>10000000)
				{
					dcntm[m] += 10000000.0; 
					 cntm[m] -= 10000000; 
					 // you may ant to put a counter here to verify Nmm

					// printf("cnt>100000000:");
					// for(int hh=0;hh<L;hh++)printf("%d", bid[hh]);
				}
			}
		}
	}
	double wsum = 0; 
	for(m=0;m<=K;m++) 
	{
		dcntm[m] +=	cntm[m]; 
//		wsum+=wm[m] * (dcntm[m] -this->mu_y * Nmm[m]); 
		wsum+=wm[m] * dcntm[m] ; 
	}

	//double mle = mu_x+(wsum/N0); 
	double mle = wsum/N0; 

	
	if (ocntm==NULL) delete []dcntm; 
	delete []cntm;
	//ncol = 1<<(2*K); // nrow = C(L,K), ncol= 4^K
	//nrow = Combinations(L,K); 
	//w = new int *[nrow]; //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	
	if (UseLookupTable)
	{
		if (lookuptable.size()<262144)
		{
			lookuptable.insert(Mymap::value_type(key, mle));
		}
	}
	return mle; 
}
