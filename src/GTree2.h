//
//  GTree2.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 6/10/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

// gapped k-mer tree

#ifndef __gkmsvmXC__GTree2__
#define __gkmsvmXC__GTree2__

#include <stdio.h>

#include "global.h"
#include "GTreeLeafData.h"
#include "GTreeLeafData2.h"

/*union GTreePtr {
    class GTree *t;
    class GTreeLeafData *node;
};*/

union GTreePtr2 {
    class GTree2 *t;
    class GTreeLeafData2 *node;
};

extern	int gLM1; //L-1
extern	int gMAXMM; //MaxMismatch
extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]
extern	myFlt **gMMProfile0; //mismatchprofile[seqidi][mm][seqidj]
extern	LTreeSnodeData ** gDFSlist[1000];
//extern  GTreeLeafData *gGTreeLeaves; // list of all the leaf nodes
extern  GTreeLeafData2 *gGTreeLeaves2; // list of all the leaf nodes // new format replacing gGTree2Leaves
extern  int gGTreeLeavesCnt; // number of all leaf nodes


///extern CLTreeSptr **gDFSlistT[1000]; // without nonEmptyDaughterCnt
//extern CLTreeS **gDFSlistT[1000]; // with nonEmptyDaughterCnt
//extern	int *gDFSMMlist[1000];
//extern CbinMMtree **gDFSMMtree[1000]; // for the iDL bound



class GTree2
{
public:
    GTreePtr2 daughter[MAX_ALPHABET_SIZE+1];
    //int maxSeqID;  int minSeqID;
    //int nonEmptyDaughterIdxs[MAX_ALPHABET_SIZE]; int nonEmptyDaughterCnt; //keeps the list of non empty daughters. this is good for sparser trees
    
  //  int addSequence(int *bid, int n, int L, int seqID);  //adds all the L-subseqs
    //	int addSequence(int *bid, int n, int L);  //adds all the L-subseqs
    //void addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int mnSeqID, int mxSeqID); // similar to addseq (but adds multiple seqs at once)
    
    void deleteTree(int n, int alphabetSize); //call with n=L from outside
    void initTree(); //initialize the tree
    
    //int addToList(LTreeSnodeData **list, int n, int single, int listlen, int alphabetSize);
    // void DFST( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize);
    
    //void DFSTiDL( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt,CbinMMtree **curMMtree, int pos, int alphabetSize); // this version has iDL (or more generally MMtree bound)
    
    
    
    //void DFSTn(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize);
    //void DFSTnIDL(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, CbinMMtree **curMMtree, int alphabetSize); // this version has iDL (or more generally MMtree bound)
    
    //int leavesCount(int withMultiplicity, int n, int alphabetSize, int *nodesAtDepth);  //returns the number of sequences in the tree.  //call with n=L from outside. // it also counts the number of nodes at each depth
    
    //void cloneReorder(CLTreeS *newTree, int *order, int n, int L, int alphabetSize, int *tmpArray,int *tmpArray2); // reorders and clones to the new tree (without replicating the data nodes) // used by iDL bound
    
    //int *reorder(int *lmer, int *order, int L, int *output); // reorders the L-mer
    
    
    GTree2(void);
    ~GTree2(void);
    
    //private:
   // void addSeq(int *bid, int n, int *lmerbid, int seqID, int nGapsRemained, int curGapPosSeq);  //call with n=L from outside
     void addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int nGapsRemained, int curGapPosSeq); // similar to addseq (but adds multiple seqs at once)
    
};

#endif /* defined(__gkmsvmXC__GTree2__) */
