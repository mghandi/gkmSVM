/* LTree.h
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
#include "LList.h"
class CLTree
{
public:
	CLTree *daughter[MAX_ALPHABET_SIZE];
	void deleteTree(int n); //call with n=L from outside 
	void mismatchCount(int *bid, int n, int *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void mismatchCount(int *bid, int n, int *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	int addSequence(int *bid, int n, int L);  //adds all the L-subseqs 
	void addSequences(char *FsaFileName, int L, int maxSequenceLength,int addrevcompl=0, int numberOfCVPartitions=0, int selectPartitionNumber=0);
	
	double calcScore(int *bid, int L, double *kernel, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
	double calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
	double calcScore(int *bid,int *bidrc, int L, int slen, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 

	int count(int *bid, int n); //returns the number of times the sequence is found in the tree. n is length of the sequence //call with n=L from outside
	int leavesCount(int withMultiplicity, int n);  //returns the number of sequences in the tree.  //call with n=L from outside

	void addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid); // used by LList. adds sequences to lis from tree 

	CLTree(void);
	~CLTree(void);

private:
	void addSeq(int *bid, int n, int cnt);  //call with n=L from outside

};

class CLTreeMemorize
{	
public: 
	CLTreeMemorize(int UseLookupTable, int LookupTableMaxSize, class CLTree *tree);
	~CLTreeMemorize(void);

	class CLTree *tree; 
	double calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 

	int UseLookupTable; 
	int LookupTableMaxSize; 
	Mymap lookuptable; 
};
