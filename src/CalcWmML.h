/* CalcWmML.h
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
class CCalcWmML
{
public:

	int K,L,b; 

	double *wm; 
	double *kernel; 
	double *kernelTruncated; //only first positive coeffs 

	double *c; //c_m : the inner product of all L-merEsts for two L-mer with m mismatches 
	double *cTr; // c using truncated g. 
	double *h; //h_m : the elements of A^{\top}A -- used for gkm count vector inner product 
	
	double n0; //smallest coeff divided by 2 
	int kernelTruncatedLength; // number of first positive coeffs
	CCalcWmML(int L, int K, int b);
	~CCalcWmML(void);

    
     static double *calcWildcardKernelWeights(int L, int M, int b, double lambda, double* res); //weights corresponding to wildcard kernel of LK2004 
     static double *calcMismatchKernelWeights(int L, int M, int b, double* res); //weights corresponding to mismatch kernel of LK2004 

    
private:
	double *calcwm(); 
	double *calcKernel(); // calculates the impulse response. 
	double *calcc(); 


     static double calcWildcardKernelWeightsm(int L, int M, int b, double lambda, int m);
     static double calcMismatchKernelWeightsm(int L, int M, int b, int m); 

};
