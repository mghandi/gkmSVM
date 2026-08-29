//
//  GTree.cpp
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 6/10/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

// gapped k-mer tree
#include "GTree.h"




void GTree::deleteTree(int n, int alphabetSize)
{
    if (n>1)
    {
        for (int i=0;i<=alphabetSize;i++)
        {
            if (daughter[i].t!=NULL)
            {
                daughter[i].t->deleteTree(n-1,alphabetSize);
                delete daughter[i].t;
                daughter[i].t=NULL;
            }
        }
    }
    
//    for(int i=0;i<=MAX_ALPHABET_SIZE;i++){
//        daughter[i].t=NULL;
//    }
//    maxSeqID=0;
//    minSeqID=0;
}


int GTree::addSequence(int *bid, int n, int L, int seqID)  //adds all the L-subseqs
{
    n = n-L+1;
    if (n<0) n=0;
    for(int i=0;i<n;i++)
    {
        addSeq(bid,L,bid, seqID, gMAXMM, 0);
        bid++;
    }
    return n; 
}


GTree::GTree(void)
{
    //daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
    for(int i=0;i<=MAX_ALPHABET_SIZE;i++){
        daughter[i].t=NULL;
    }
//  maxSeqID=0;
//  minSeqID=0;
    
//    nonEmptyDaughterCnt=0;
}


GTree::~GTree(void)
{
}

void GTree::initTree()
{
    //daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
    for(int i=0;i<=MAX_ALPHABET_SIZE;i++){
        daughter[i].t=NULL;
    }
//    maxSeqID=0;
//    minSeqID=0;
    
//    nonEmptyDaughterCnt=0;
}
void GTree::addSeq(int *bid, int n, int *lmerbid, int seqID, int nGapsRemained, int curGapPosSeq)  //call with n=L from outside
//nGapsRemained: number of remaining gaps. (initially call with Dmax)
//curGapPosSeq: contains the bases in the gapped positions up to this depth. these are coded in 1 integer. Assumes that: Dmax * ceiling(log2(MAX_ALPHABET_SIZE)) < 32
{
//    if (seqID>maxSeqID) maxSeqID=seqID;
//    if (seqID<minSeqID) minSeqID=seqID;
    if (n==1)
    {
        int tbid = *bid; // target base id
        if(nGapsRemained==1){
            // this last pos should be gap
            curGapPosSeq=(*bid)+(curGapPosSeq<<NBITS);
            tbid =MAX_ALPHABET_SIZE; // Gap
        }
        
        if (this->daughter[tbid].node==NULL)
        {
            //gGTreeLeaves[gGTreeLeavesCnt] = new GTreeLeafData();
            daughter[tbid].node = &gGTreeLeaves[gGTreeLeavesCnt];
            gGTreeLeavesCnt++;
        }
        daughter[tbid].node->add(seqID, curGapPosSeq);
    }
    else
    {
        if(n>nGapsRemained){
            if (this->daughter[*bid].t == NULL)
            {
                this->daughter[*bid].t = new GTree();
    //            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
            }
            daughter[*bid].t->addSeq(bid+1, n-1, lmerbid, seqID, nGapsRemained, curGapPosSeq);
        }
        if (nGapsRemained>0){
            if (this->daughter[MAX_ALPHABET_SIZE].t == NULL)
            {
                this->daughter[MAX_ALPHABET_SIZE].t = new GTree();
                //            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=MAX_ALPHABET_SIZE;
            }
            daughter[MAX_ALPHABET_SIZE].t->addSeq(bid+1, n-1, lmerbid, seqID, nGapsRemained-1, (*bid)+(curGapPosSeq<<NBITS));
            
        }
    }
}

/*
void GTree::addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int nGapsRemained, int curGapPosSeq){ // similar to addseq (but adds multiple seqs at once)
    
    if (n==1)
    {
        int tbid = *bid; // target base id
        if(nGapsRemained==1){
            // this last pos should be gap
            curGapPosSeq=(*bid)+(curGapPosSeq<<NBITS);
            tbid =MAX_ALPHABET_SIZE; // Gap
        }
        
        if (this->daughter[tbid].node==NULL)
        {
            //gGTreeLeaves[gGTreeLeavesCnt] = new GTreeLeafData();
            daughter[tbid].node = &gGTreeLeaves2[gGTreeLeavesCnt];
            gGTreeLeavesCnt++;
        }
        daughter[tbid].node->addLTreeSnodeData(nodeData, curGapPosSeq);
    }
    else
    {
        if(n>nGapsRemained){
            if (this->daughter[*bid].t == NULL)
            {
                this->daughter[*bid].t = new GTree();
                //            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
            }
            daughter[*bid].t->addLTreeSnodeData(bid+1, n-1, nodeData, nGapsRemained, curGapPosSeq);

        }
        if (nGapsRemained>0){
            if (this->daughter[MAX_ALPHABET_SIZE].t == NULL)
            {
                this->daughter[MAX_ALPHABET_SIZE].t = new GTree();
                //            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=MAX_ALPHABET_SIZE;
            }
            daughter[MAX_ALPHABET_SIZE].t->addLTreeSnodeData(bid+1, n-1, nodeData, nGapsRemained-1, (*bid)+(curGapPosSeq<<NBITS));

        }
    }
}
*/