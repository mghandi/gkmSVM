/* LList.h
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
#include "LTree.h"

#define L12THRM 7  /* number of bases per integer, it should be <=15, but smaller value may result in faster performance because of the cpu chache limit.*/

#define USEMMLOOKUPTABLE 1  /* defines whether use a lookup table or bitwise operations to count number of mismatches */
#define HAMDIST(x) hamdist[x]

//#define USEMMLOOKUPTABLE 0  /* defines whether use a lookup table or bitwise operations to count number of mismatches */
//#define HAMDIST(x) cmmcnt(x)

class CLList
{
public:
	CLList(int L, int maxnLmers, int *extHdist=0);
	~CLList(void);
	
	int **table; // bp 1-L12THR are in table[1], bp[L12THR+1-2*L12THR] in table[2] , etc 
	int *freq; //count of Lmeri 

	int L; 
	int twidth; // (L+L12THR-1)/L12THR

	int nmultiples,nsingles; //number of L-mers with multiple frequencies and single frequencies. 

	int *HamDist; 
	int *extHamdist; // if extHamdist is available it does not allocate its own Hamdist
	
	void mismatchCount(int *bid, int *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; 
	void mismatchCount(int col, int *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ; 
	double calcScore(int *bid, double *kernel, int *mmcnt); //calculates the score. tmpcnt is int[L+1] 


	void addFromLTree(class CLTree *tree); 

	void addSeq(int *bid, int cnt); // adds one sequence : for efficiency, first add sequences with cnt=1, then add those with cnt>1

	void clear(); //clears the list 
	
	void calcScoreAllLmers(double *kernel, int *mmcnt, double n0, double *cntEst); // calcs freqs for all L-mers and fills cntEst[0..4^L-1]
	void calcScoreAllLmersAdd(double *kernel, int *mmcnt, double n0, double *cntEst, double weight); // calcs freqs for all L-mers and adds fills cntEst[0..4^L-1] weighted by weight

	int *tmpinttw; 

	int UseLookupTable; 
	int LookupTableMaxSize; 
	Mymap lookuptable; 

	
	double calcInnerProd(CLList *L2, double *c, int *mmcnt);
	double calcInnerProd1(CLList *L2, double *c, int *mmcnt);
	double calcInnerProd2(CLList *L2, double *c, int *mmcnt);

    
    char *convertInt2Str(int col, char *str, int L); // returns L-mer for idx=col

	int L12THR;

private:
	int convert(int *bid, int i);
	int convert(int col, int i);

	void mismatchCount2(int *bid, int *cnt); // optimized for width=2
	void mismatchCount3(int *bid, int *cnt); // optimized for width=3
	void mismatchCount4(int *bid, int *cnt); // optimized for width=4

};

