/* CalcWmML.cpp
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

#include "CalcWmML.h"

#include "global.h"


CCalcWmML::CCalcWmML(int L, int K, int b)
{
	this->L = L; 
	this->K = K; 
	this->b = b; 

	wm = new double [K+1]; 
	kernel = new double [L+1]; 
	kernelTruncated = new double [L+1]; 


	c = new double [L+1]; //c_m : the inner product of all L-merEsts for two L-mer with m mismatches 
	cTr = new double [L+1]; // c using truncated g. 
	h = new double [L+1]; // h_m : for inner product of all gapped-kmer counts . 

	calcwm(); 
	calcKernel(); 
	calcc(); 

}

double *CCalcWmML::calcwm()
{
	double **wL = new double*[K+1]; 
	double **wLp = new double*[K+1]; 
	double **hv; 

	int i,j; 
	for(i=0;i<=K;i++)
	{
		wL[i]= new double[K+1]; 
		wLp[i]= new double[K+1]; 

		for(j=0;j<=K;j++)
		{
			wL[i][j] = wLp[i][j] = 1.0; 
		}
	}


	for (int iL=1; iL<=L; iL++)
	{
		for(int iK=1; iK<=K; iK++)
		{
			wL[iK][0] = wLp[iK][0] + (b-1)* wLp[iK-1][0]; 
		
			for (int jM=1; jM<=iK; jM++)
			{
				wL[iK][jM] = (wL[iK-1][jM-1] * (iK-iL))/iK;  
			}
		}

		hv = wLp; wLp=wL; wL=hv; 
	}


	double nnorm = dCombinations(L,K)*pow(b,1.0*L); 

	for (i=0;i<=K;i++)
	{
		wm[i] = wLp[K][i]/nnorm; 
		//	wm[i] = wLp[K][i]; 
	}


	for(i=0;i<=K;i++)
	{
		delete []wL[i];
		delete []wLp[i];
	}
	delete []wL;
	delete []wLp;

	// h[m] 
	for (i=0; i<=L; i++)
	  {
	    if ((L-i) >= K)
	      {
		h[i] = dCombinations(L-i,K); 
	      }
	  }

	/*

	if ((L==15)&&(K==7)&&(b==4))
	{
		wm[0] =  18575560.00; 
		wm[1] =  -5145388.57; 
		wm[2] =   1463262.86; 
		wm[3] =   -424597.71;
		wm[4] =    125173.71;
		wm[5] =    -37374.86;
		wm[6] =     11276.57;
		wm[7] =     -3432.00;
	}

	if ((L==13)&&(K==7)&&(b==4))
	{ 
		wm[0] =   5382976.0;
		wm[1] =  -1397214.9;
		wm[2] =    379120.0;
		wm[3] =   -106206.4;
		wm[4] =     30470.4;
		wm[5] =     -8904.0;
		wm[6] =      2640.0;
		wm[7] =      -792.0;
	}

	if ((L==2)&&(K==1)&&(b==4))
	{
		wm[0] =  7; 
		wm[1] =  -1;
	}
	*/

	return wm; 
}


double *CCalcWmML::calcKernel(void)
{
	for(int m=0;m<=L;m++)
	{
		kernel[m]=0;
		for(int i=0;i<=m;i++)
		{
			kernel[m]+=wm[i]*dCombinations(L-m,K-i)*dCombinations(m,i);
		}
	}

	kernelTruncatedLength = 0; 
	int hn=1; 
	int i; 
	for(i=0;i<=L;i++)
	{
		if (kernel[i]<1e-50) hn=0; 
		if (hn)
		{
			kernelTruncated[i]=kernel[i]; 
			n0=kernel[i]/2;
			kernelTruncatedLength=i+1; 
		}
		else
		{
			kernelTruncated[i]=0.0;
		}		
	}


	return kernel;
}

CCalcWmML::~CCalcWmML(void)
{
	delete []wm; 
	delete []kernel; 
	delete []kernelTruncated;
	delete []h; 
	delete []c; 
	delete []cTr; 
}



double *CCalcWmML::calcc()
{
	for (int m=0;m<=L;m++)
	{
		c[L-m] = 0; 
		cTr[L-m] = 0; 
		for(int m1 = 0;m1<=L; m1++)
			for(int m2 = 0;m2<=L; m2++)
				for(int t = 0;t<=L; t++)
				{
					int r= m1+m2-2*t-L+m; 
					if ((t<=m)&&((m1-t)<=(L-m))&&(r<=(m1-t))&&(r>=0))
					{
						double cc = dCombinations(m,t)*dCombinations(L-m,m1-t)*dCombinations(m1-t,r)*pow(b-1, 1.0*t)*pow(b-2, 1.0*r); 
						c[L-m]+= cc*kernel[m1]*kernel[m2]; 
						cTr[L-m]+= cc*kernelTruncated[m1]*kernelTruncated[m2]; 
					}
				}

	}
	return c; 
}



double CCalcWmML::calcWildcardKernelWeightsm(int L, int M, int b, double lambda, int m) //weights corresponding to wildcard kernel of LK2004 
{
	double w = 0; 
	for(int k=(L-M); k<=L; k++){
		if (L-m >= k){
			w = w+ pow(lambda,1.0*(L-k)) * dCombinations(L-m, k);
		}
	}
	return w;
}

double *CCalcWmML::calcWildcardKernelWeights(int L, int M, int b, double lambda, double* res) //weights corresponding to wildcard kernel of LK2004 
{
    for(int m=0; m<=L; m++){
        res[m]=calcWildcardKernelWeightsm(L,M,b,lambda, m);
    }
    return res;
}

double CCalcWmML::calcMismatchKernelWeightsm(int L, int M, int b, int m)
{
    double w = 0; 
	for(int m1=0; m1<=M; m1++){
        for(int m2=0; m2<=M; m2++){
			for(int t=0; t<=M; t++){
				int r = m2+m1-m-2*t; 
				w = w+ dCombinations(L-m,t)*pow(b-1, 1.0*t)*dCombinations(m,r)*pow(b-2,1.0*r)* dCombinations(m-r, m1-t-r);
			}
		}
	}
	return w;
}

double *CCalcWmML::calcMismatchKernelWeights(int L, int M, int b, double* res) //weights corresponding to mismatch kernel of LK2004 
{
    for(int m=0; m<=L; m++){
        res[m]=calcMismatchKernelWeightsm(L,M,b, m);
    }
    return res;
}

