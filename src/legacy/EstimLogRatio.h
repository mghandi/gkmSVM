/* EstimLogRatio.h
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

#define Ndatapoints 5000
class CEstimLogRatio
{
public:
	CEstimLogRatio(void);
	~CEstimLogRatio(void);
	double estimateLogRatio(double u, double v, double s2n, double s2lr,double mur=0);  //estimates log(x/y), where u=x+n1, v=y+n2 are given. s2n is variance of noise and s2lr is variance of log(x/y) and mur is the mean of log(x/y)
    
    double estimateLogRatio2(double u, double v, double *kernel, int L) ; 


private:
	double r[Ndatapoints];  //log ratio 
	double z[Ndatapoints];  //ratio = exp(r) 

	double z2p1[Ndatapoints];//=1+z[i]*z[i];
	double lz2p1[Ndatapoints]; //=0.5*log(z2p1[i]);
	
	double q[Ndatapoints]; //f(r|u,v) 

};

