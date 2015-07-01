/* LTreeS.h
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
#include "GTree2.h"
#include "CbinMMtree.h"

union CLTreeSptr {
    class CLTreeS *t;
	LTreeSnodeData *node;
};

extern	int gLM1; //L-1
extern	int gMAXMM; //MaxMismatch
extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]
extern	myFlt **gMMProfile0; //mismatchprofile[seqidi][mm][seqidj]
extern	LTreeSnodeData ** gDFSlist[1000]; 
//extern CLTreeSptr **gDFSlistT[1000]; // without nonEmptyDaughterCnt
extern CLTreeS **gDFSlistT[1000]; // with nonEmptyDaughterCnt
extern	int *gDFSMMlist[1000]; 
extern CbinMMtree **gDFSMMtree[1000]; // for the iDL bound



class CLTreeS
{
public:
	CLTreeSptr daughter[MAX_ALPHABET_SIZE]; int maxSeqID;  int minSeqID;
	int nonEmptyDaughterIdxs[MAX_ALPHABET_SIZE]; int nonEmptyDaughterCnt; //keeps the list of non empty daughters. this is good for sparser trees
	#ifdef FAST_TRACK
		int *FT_seq; int FT_cnt;  int FT_seqID; // these are for fast track. FT_seq is a link to the last sequence added to this node. FT_cnt is the number of downstream nodes.
	#endif
	int addSequence(int *bid, int n, int L, int seqID);  //adds all the L-subseqs 
//	int addSequence(int *bid, int n, int L);  //adds all the L-subseqs 
    void addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int mnSeqID, int mxSeqID); // similar to addseq (but adds multiple seqs at once)
    
	void deleteTree(int n, int alphabetSize,int DontDeleteNodeData); //call with n=L from outside
	void initTree(); //initialize the tree
	
    int addToList(LTreeSnodeData **list, int n, int single, int listlen, int alphabetSize);
    void addToGTree(GTree2 *gtree, int n,int *tmpArray,int alphabetSize,int L);

    
	//void DFS( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt, int pos);
//	int DFSn0(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, LTreeSnodeData *nodei, int bid);
//	void DFSn1(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, LTreeSnodeData *nodei, int bid, int nmulti);

//	void DFSsingle( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt, int pos); //single matching L-mers
//	void DFSmulti( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt, int pos);//multi matching L-mers



	//	void DFS( /*LTreeSnodeData **matchingLmers,*/ int listlen,  /*int *curMismatchCnt, */int pos);
//	void DFSn( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt);
//	void DFSnSingle( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt);
//	void DFSnMulti( LTreeSnodeData **matchingLmers, int listlen,  int *curMismatchCnt);
	
	//void DFSnf(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt);

//	void DFST( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize);// without nonEmptyDaughterCnt
//	void DFSTn(CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize);// without nonEmptyDaughterCnt
    void DFST( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize);

    void DFSTiDL( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt,CbinMMtree **curMMtree, int pos, int alphabetSize); // this version has iDL (or more generally MMtree bound)

    
    
	void DFSTn(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize);
    void DFSTnIDL(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, CbinMMtree **curMMtree, int alphabetSize); // this version has iDL (or more generally MMtree bound)

	//void DFSTf( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos);  // this is for fast version -- yet to be implemented
	//void DFSTnf(CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt);

	int leavesCount(int withMultiplicity, int n, int alphabetSize, int *nodesAtDepth);  //returns the number of sequences in the tree.  //call with n=L from outside. // it also counts the number of nodes at each depth

    void cloneReorder(CLTreeS *newTree, int *order, int n, int L, int alphabetSize, int *tmpArray,int *tmpArray2); // reorders and clones to the new tree (without replicating the data nodes) // used by iDL bound

    int *reorder(int *lmer, int *order, int L, int *output); // reorders the L-mer

//	void mismatchCount(int *bid, int n, int *cnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
//	void mismatchCount(int *bid, int n, int *cnt, int maxmm); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) ;  //call with n=L from outside
//	void addSequences(char *FsaFileName, int L, int maxSequenceLength,int addrevcompl=0, int numberOfCVPartitions=0, int selectPartitionNumber=0);
	
//	double calcScore(int *bid, int L, double *kernel, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
//	double calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 
//	double calcScore(int *bid,int *bidrc, int L, int slen, double *kernel, int maxmm, int *tmpcnt); //calculates the score. tmpcnt is int[L+1] 

//	int count(int *bid, int n); //returns the number of times the sequence is found in the tree. n is length of the sequence //call with n=L from outside
//	int leavesCount(int withMultiplicity, int n);  //returns the number of sequences in the tree.  //call with n=L from outside

//	void addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid); // used by LList. adds sequences to lis from tree 

	//void DFST( class CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos);


	CLTreeS(void);
	~CLTreeS(void);

private:
	void addSeq(int *bid, int n, int *lmerbid, int seqID);  //call with n=L from outside
//	void addSeq(int *bid, int n, int cnt);  //call with n=L from outside

};
