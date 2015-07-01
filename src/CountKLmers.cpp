/* CountKLmers.cpp
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

#include "CountKLmers.h"
#include "global.h"
#include "KLmer.h"


CCountKLmers::CCountKLmers(int L, int K, int halfMem) // set RCSymmetric=1 to limit the first base of the motif to A,C
{

	this->K = K; 
	this->L = L; 
	this->halfMem = halfMem; 

	ncol = halfMem?(1<<(2*K-1)):(1<<(2*K)); // nrow = C(L-1,K-1), ncol=RCSym? (4^K)/2 : 4^K
	nrow = Combinations(L-1,K-1); 

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

void CCountKLmers::addSequence(int *seqBID, int size)
{
	int halfMem = this->halfMem;

	int *si = seqBID; 
	for(int i=0;i<=(size-L); i++)
	{
		if (!halfMem || (*si<2)) // first base can only be A or C. 
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
		}
		si++;  
	}
}

char *CCountKLmers::convertCol2KmerString(int col, char *sKmer) // returns K-mer for idx=col
{
	if (this->halfMem)
	{
		sKmer[0] = ::globalConverter.icidx[col%2]; 
		col >>= 1; 
	}
	else
	{
		sKmer[0] = ::globalConverter.icidx[col%4]; 
		col >>= 2; 
	}
	for(int i=1;i<K;i++)
	{
		sKmer[i] = ::globalConverter.icidx[col%4]; 
		col >>= 2; 
	}
	sKmer[K]=0; 
	return sKmer; 
}

char *CCountKLmers::convertRow2KLmerString(int row, char *sKmer, char *sKLmer) // maps Kmer to KLmer for idx=row
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


int CCountKLmers::fillwij(int pos, int nfilled, int idx, int *partial)
{
	if (pos==L)
	{
		for(int j=0;j<L;j++)
		{
			w[idx][j]=partial[j]; 
		}
		return idx+1; 
	}

	if ((pos>0)&&((L-pos)>(K-nfilled))) // first position always selected
	{
		partial[pos] = 0;
		idx = fillwij(pos+1, nfilled, idx, partial); 
	}
	if (nfilled<K)
	{
		if (nfilled==0)
		{
			partial[pos] = 1; 
		}
		else
		{
			partial[pos] = halfMem?(1<<(2*nfilled-1)):(1<<(2*nfilled)); // if RCSym, the first base can  be A or C (0 or 1)
		}
		idx = fillwij(pos+1, nfilled+1, idx, partial); 
	}
	return idx; 
}


CCountKLmers::~CCountKLmers(void)
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

double *matrixMultiply(int **rij, int **nij, int *deg, int n, double *x, double *y, double scale)
{

	for(int i=0; i<n; i++) 
	{
			double mfxi=0; 
			int *crij = rij[i]; 
			int *cnij = nij[i]; 
			for(int j=0; j<deg[i]; j++)
			{
				int jklmer = cnij[j]; 
				mfxi+=crij[jklmer]*x[jklmer]; 
			}
			y[i] = mfxi*scale; 
	}

	return y; 
}

//
//out.15.7.17555.10k
// /Win Debug: c:\data\test.klm 0 2 1 200 100 0 0.2 c:\data\outfn
// /Win Debug: c:\data\out.13.7.17555.10k 1 13 7 11000 100 0 0.2 c:\data\outfn

//out.15.7.17555.t1.2k
/*
int mainKLmerLinApprox(int argc, char * argv[]) // mainKLmerLinApprox : given a list of KL-mers and the corresponding value, it calculates wi's that best approximates the score. 
{
	int OUTPUTA = 1; 
	if (argc!=10)
	{
		if (OUTPUTA)
		{
			printf ("\n **this version, also outputs matrix A\n "); 
		}
		printf(" given a list of top KL-mers and the corresponding values, it calculates wi's inorder to approximate the score function.\n");

		printf("usage: ./KLmerLinApprox <KLmerFN> <takeLog> <L> <K> <maxnKLmers> <niter> <rseed> <learningRate> <outputfn>\n");
		printf("for example: ./KLmerLinApprox out.15.7.17555.10k 0 15 7 10000 100 0 0.2 outfn\n");
		return 0; 
	} 

	char *infn = argv[1]; 
	int takeLog = atoi(argv[2]); 
	int L = atoi(argv[3]); 
	int K = atoi(argv[4]); 
	int nmax = atoi(argv[5]); 
	int niter = atoi(argv[6]); 
	int rseed = atoi(argv[7]); 	srand(rseed);
	double learnRate = atof(argv[8]);
	char *outfn = argv[9]; 

	int i,j; 
	CKLmer **klmers = new CKLmer*[nmax]; // 
	for(i=0;i<nmax;i++)
	{
		klmers[i] = new CKLmer(L,K); 
	}

	int n=0; 
	char s[1000]; char tmps[1000]; 
	double escore; 
	int tmpi; 
	FILE *fi = fopen(infn, "r"); 
	while (!feof(fi))
	{
		fgets(s, 1000, fi); 
		if (length(s)>L)
		{
			klmers[n]->readKLmer(s); 
			sscanf(s, "%s%d%d%lf", tmps, &tmpi, &tmpi, &escore);
			klmers[n]->score = (takeLog?log(escore):escore); 
			klmers[n]->w = 0; 
			n++; 
			if (n>=nmax) 
			{
				printf("\n reached maximum n %d limit. increase nmax to avoid data truncation.\n ",nmax); 
				break; 
			}
		}
	}
	fclose(fi); 

	int **rij = new int*[n]; // number of common Kmers
	int *deg = new int[n];  // number of of indexes with rij>0
	int **nij = new int*[n];  // list of indexes with rij>0
	for(i=0;i<n;i++)
	{
		rij[i] = new int[n]; 
		nij[i] = new int[n]; 
		deg[i] = 0; 
		for(j=0;j<n;j++)
		{
			rij[i][j] = klmers[i]->commonKMerCnt(klmers[j]); 
			if (rij[i][j]>0) 
			{
				nij[i][deg[i]]=j; 
				deg[i]++; 
			}
		}
	}


	// output matrix A
	if (OUTPUTA)
	{
		FILE *foa = fopen("Amatrix.txt", "w"); 
		for(i=0;i<n;i++)
		{
			fprintf(foa, "%d",rij[i][0]); 
			for(j=1;j<n;j++)
			{
				fprintf(foa, "\t%d",rij[i][j]); 
			}
			fprintf(foa, "\n"); 
		}
		fclose(foa); 
	}
	// 
	double *fbg = new double[n]; 
	double *w = new double[n]; 

	for(i=0;i<n;i++)
	{
		fbg[i] = 0;//klmers[i]->calcfbg(); 
	}

	for(i=0;i<n;i++)
	{
		klmers[i]->w=0; 
		w[i] = 0; 
	}

	int M=1<<(2*(L-K)); // number of Kmers in a KLmer

	int *crij; 
	int *cnij;
	int jklmer;

	
	// iterate
	int iter; 
	double sserr; 
	
	int IterativeMethod=0; 
	if (!IterativeMethod)
	{
	
	for (iter=0; iter<niter; iter++)
	{
		double alpha = learnRate *(1.0/(iter+1.0)); 
		alpha = learnRate *pow(0.95,iter)/n; 
		
		//learnRate
		sserr = 0; 
		for(int jter=0;jter<n;jter++)
		{
			i = myrandom(n); 
			double mfxi=0; 
			crij = rij[i]; 
			cnij = nij[i]; 
			for(j=0; j<deg[i]; j++)
			{
				jklmer = cnij[j]; 
				mfxi+=crij[jklmer]*w[jklmer]; 
			}
			mfxi = mfxi/M + fbg[i];
			double erri=klmers[i]->score - mfxi; 
			w[i] = w[i]+alpha*(erri); 
			sserr+=erri*erri;
		}
		printf("%d\t%lf\n",iter, sqrt(sserr /n)); 
	}

	}

	double *y = new double[n];
	double *Pw = new double[n];
	double *P2w = new double[n];

	if(IterativeMethod)
	{
 


	// using the iterative method: 
	int N=Combinations(L, K)*1<<(2*K); // number of KLmers covering any kmer


	double MN= M; MN*=N; 

	
	
	double *scores = new double[n]; 
	for(i=0;i<n;i++)
	{
		scores[i] = klmers[i]->score; 
	}

	// y = P*s
	MN = (M/1.0)*n; 
	matrixMultiply(rij, nij, deg, n, scores, y, 1/MN); 

	// init w; 
	for(i=0;i<n;i++) 
	{
		w[i] = 0; 
	}

	double lambda = 10.0; 
	lambda = learnRate; 
	for (iter=0;iter<niter; iter++)
	{
		sserr=0; 
		matrixMultiply(rij, nij, deg, n, w, Pw, 1/MN); 
		matrixMultiply(rij, nij, deg, n, Pw, P2w, 1/MN); 
		for(i=0;i<n;i++)
		{
			w[i] = w[i] + lambda*(y[i]-P2w[i]); 
			sserr+=(y[i]-P2w[i])*(y[i]-P2w[i]); 
		}
		printf("itm- %d\t%lf\n",iter, sqrt(sserr /n)); 
	}
	// end iter 


	}

	// output
	FILE *fo = fopen(outfn, "w"); 
	fprintf(fo, "KLmer\tscore\twLinApprx\n"); 
	for(i=0;i<n;i++) 
	{
		fprintf(fo, "%s\t%lf\t%lf\n",klmers[i]->seq,  klmers[i]->score, klmers[i]->w=w[i]); 
	}
	fclose(fo);
	

	// 
	delete Pw; 
	delete P2w; 

	delete []y; 
	delete []fbg; 

	for(i=0;i<n;i++)
	{
		delete []rij[i]; 
		delete []nij[i]; 
	}
	delete []rij; 
	delete []nij; 
	delete []deg; 
	delete []w; 

	for(i=0;i<nmax;i++)
	{
		delete klmers[i]; 
	}
	delete []klmers;
}



*/

	

