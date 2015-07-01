/* LTreef.h
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
#include "LTreeS.h"



class CLTreef
{
public:
	fintptr_t daughter[MAX_ALPHABET_SIZE];
	int nonEmptyDaughterIdxs[MAX_ALPHABET_SIZE]; int nonEmptyDaughterCnt; //keeps the list of non empty daughters. this is good for sparser trees
	void deleteTree(int n, int alphabetSize); //call with n=L from outside
	void imismatchCount(int *bid, int n, int *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void fmismatchCount(int *bid, int n, myFlt *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void imismatchCount(int *bid, int n, int *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void fmismatchCount(int *bid, int n, myFlt *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void fimismatchCount(class CLTreef *iTree, int n, myFlt *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void iimismatchCount(class CLTreef *iTree, int n, int *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	void iimismatchCountGeneral(class CLTreef *iTree, int n, int *cnt, int maxmm, int alphabetSize); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
	int addSequence(int *bid, int n, int L);  //adds all the L-subseqs 
	int addSequence(int *bid, int n, int L, myFlt w);  //adds all the L-subseqs w times 
	void addSequences(char *FsaFileName, int L, int maxSequenceLength,int addrevcompl=0, int numberOfCVPartitions=0, int selectPartitionNumber=0);

	double calciScore(int *bid, int L, double *kernel, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
	double calciScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
	double calciScore(class CLTreef *iTree, int L, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1]. iTree is the suffix tree containing all the sequences. sequence count are assumed to be integer
	double calcfScore(int *bid, int L, double *kernel, myFlt *tmpcnt); //calculates the score. tmpcnt is myFlt[L+1] 
	double calcfScore(int *bid, int L, double *kernel, int maxmm,  myFlt *tmpcnt); //calculates the score. tmpcnt is myFlt[L+1] 
	double calcfScore(class CLTreef *iTree, int L, double *kernel, int maxmm, myFlt *tmpcnt); //calculates the score. tmpcnt is myFlt[L+1]. iTree is the suffix tree containing all the sequences. sequence count are assumed to be integer


	void addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid, int alphabetSize); // used by LList. adds sequences to lis from tree

	int icount(int *bid, int n); //returns the number of times the sequence is in found in the tree. n is length of the sequence //call with n=L from outside
	myFlt fcount(int *bid, int n); //returns the number of times the sequence is in found in the tree. n is length of the sequence //call with n=L from outside
	int leavesCount(int withMultiplicity, int n, int alphabetSize);  //returns the number of sequences in the tree.  //call with n=L from outside


//	void DFST( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize); //calculate svmScoreunorm for all sequences in given tree // without nonEmptyDaughterCnt
//	void DFSTn(CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize);// without nonEmptyDaughterCnt
	void DFST( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize); //calculate svmScoreunorm for all sequences in given tree // with nonEmptyDaughterCnt
	void DFSTn(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize);// with nonEmptyDaughterCnt

	CLTreef(void);
	~CLTreef(void);

private:
	void addSeq(int *bid, int n, int cnt);  //call with n=L from outside
	void addSeq(int *bid, int n, myFlt cnt);  //call with n=L from outside

};


