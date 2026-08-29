/* EstimLogRatio.cpp
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

#include "EstimLogRatio.h"
#include "global.h"

double Sqr(double x) {return x*x;}



double CEstimLogRatio::estimateLogRatio2(double u, double v, double *kernel, int L)  
{
    double eps = 1E-90; 
    double n0=1; 
    for (int i=0;i<=L;i++)
    {
        if (fabs(kernel[i])<n0)
        {
            if (fabs(kernel[i])>eps)
            {
                n0 = fabs(kernel[i]);
            }
        }
    }
    n0 = n0/2; 
    if (u<0) u = 0; 
    if (v<0) v = 0; 
    return log((u+n0)/(v+n0))/log(10.0);
}


double CEstimLogRatio::estimateLogRatio(double u, double v, double s2n, double s2lr, double mur)
{
	int bi=0; 
	for(int i=0;i<Ndatapoints;i++)
	{
		//double z2 = z[i]*z[i];
		q[i] = Sqr(r[i]-mur)/(2*s2lr);
		q[i]+= Sqr(u-v*z[i])/(2*s2n*z2p1[i]);
		q[i]+=lz2p1[i];
		if (q[i]<q[bi]) 
		{
			bi = i; 
		}
	}

	//printf("\n u=%lf\tv=%lf\ts2n=%lf\ts2lr=%lf\tmur=%lf\trest=%lf\n",u,v,s2n,s2lr,mur, r[bi]);

	return r[bi];
}

CEstimLogRatio::CEstimLogRatio(void)
{
	double c= 4.0 / (0.5*Ndatapoints);  // to cover range 4>log[z]>-4
	for(int i=0;i<Ndatapoints;i++)
	{
		this->r[i] = c*(i-Ndatapoints/2);
		this->z[i] = exp(r[i]);
		z2p1[i]=1+z[i]*z[i];
		lz2p1[i]=0.5*log(z2p1[i]);
	
	
	}
}


CEstimLogRatio::~CEstimLogRatio(void)
{
}
