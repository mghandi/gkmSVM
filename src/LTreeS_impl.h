/* LTreeS.cpp
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

#include "LTreeS.h"
#include <vector>
#include "Sequence.h"
#include "globalvar.h"

namespace GKM_NS {

CLTreeS::CLTreeS(void)
{
  //daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
  for(int i=0;i<MAX_ALPHABET_SIZE;i++){
    daughter[i].t=NULL;
  }
  maxSeqID=0;
  minSeqID=0;
  
  nonEmptyDaughterCnt=0;
}

CLTreeS::~CLTreeS(void)
{
}

void CLTreeS::initTree()
{
  //daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
  for(int i=0;i<MAX_ALPHABET_SIZE;i++){
    daughter[i].t=NULL;
  }
  maxSeqID=0;
  minSeqID=0;
  
  nonEmptyDaughterCnt=0;
}

void CLTreeS::addSeq(int *bid, int n, int *lmerbid, int seqID)  //call with n=L from outside
{
  if (seqID>maxSeqID) maxSeqID=seqID; 
  if (seqID<minSeqID) minSeqID=seqID; 
  
  
  if (n==1)
  {
    if (this->daughter[*bid].t==NULL)
    {
      LTreeSnodeData* nodeData = new LTreeSnodeData; 
      nodeData->n = 1;
      nodeData->seqIDs.i= seqID;
      //			nodeData->baseID = lmerbid;
      this->daughter[*bid].node=nodeData;
      this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
    }
    else
    {
      LTreeSnodeData* curnodeData =this->daughter[*bid].node;
      if (curnodeData->n==1)
      {
        intintptr newseqids; 
        newseqids.p= new int[2]; 
        newseqids.p[0]= curnodeData->seqIDs.i; 
        newseqids.p[1]= seqID;
        curnodeData->seqIDs.p = newseqids.p;
        curnodeData->n=2;
      }
      else
      {
        if ((curnodeData->n & (curnodeData->n-1))==0) // i.e. n is power of 2
        {
          // expand memory 
          intintptr newseqids; 
          newseqids.p= new int[(curnodeData->n)<<1]; 
          for(int j=0;j< curnodeData->n; j++)
          {
            newseqids.p[j]= curnodeData->seqIDs.p[j]; 
          }		
          delete []curnodeData->seqIDs.p; 
          curnodeData->seqIDs.p = newseqids.p; 
        }
        curnodeData->seqIDs.p[curnodeData->n] = seqID;
        curnodeData->n += 1;
      }
    }
  }
  else
  {
    if (this->daughter[*bid].t == NULL)
    {
      this->daughter[*bid].t = new CLTreeS(); 
      this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
    }
    daughter[*bid].t->addSeq(bid+1, n-1, lmerbid, seqID); 
  }
}

void CLTreeS::addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int mnSeqID, int mxSeqID)
  // similar to add seq, but adds multiple seqs together, used by iDLbound only to generate reordered trees
  //call with n=L from outside
{
  if (mxSeqID>maxSeqID) maxSeqID=mxSeqID;
  if (mnSeqID<minSeqID) minSeqID=mnSeqID;
  
  
  if (n==1)
  {
    if (this->daughter[*bid].t==NULL)
    {
      //            LTreeSnodeData* nodeData = new LTreeSnodeData;
      //            nodeData->n = 1;
      //            nodeData->seqIDs.i= seqID;
      //            nodeData->baseID = lmerbid;
      this->daughter[*bid].node=nodeData;
      this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
    }
    else
    {
      Printf(" nonempty node not expected Error !\n");
      return; //exit(1);
      /*
      LTreeSnodeData* curnodeData =this->daughter[*bid].node;
      if (curnodeData->n==1)
      {
      intintptr newseqids;
      newseqids.p= new int[2];
      newseqids.p[0]= curnodeData->seqIDs.i;
      newseqids.p[1]= seqID;
      curnodeData->seqIDs.p = newseqids.p;
      curnodeData->n=2;
      }
      else
      {
      if ((curnodeData->n & (curnodeData->n-1))==0) // i.e. n is power of 2
      {
      // expand memory
      intintptr newseqids;
      newseqids.p= new int[(curnodeData->n)<<1];
      for(int j=0;j< curnodeData->n; j++)
      {
      newseqids.p[j]= curnodeData->seqIDs.p[j];
      }
      delete []curnodeData->seqIDs.p;
      curnodeData->seqIDs.p = newseqids.p; 
      }
      curnodeData->seqIDs.p[curnodeData->n] = seqID;
      curnodeData->n += 1;
      }
      */
    }
  }
  else
  {
    if (this->daughter[*bid].t == NULL)
    {
      this->daughter[*bid].t = new CLTreeS(); 
      this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
    }
    daughter[*bid].t->addLTreeSnodeData(bid+1, n-1, nodeData,  mnSeqID,  mxSeqID );
  }
}

int CLTreeS::addToList(LTreeSnodeData **list, int n, int single, int listlen, int alphabetSize) // // adds all the Lmers (that are (not) present in single sequence) to a list
{
  if (n==1)
  {
    for(int bid=0;bid<alphabetSize;bid++)
    {
      
      if (this->daughter[bid].t==NULL) continue; 
      
      if ((this->daughter[bid].node->n==1)==single)
      {
        list[listlen]= this->daughter[bid].node; 
        listlen++;
      }
    }
  }
  else
  {
    for(int bid=0;bid<alphabetSize;bid++)
    {
      if (this->daughter[bid].t != NULL)
      {
        listlen=daughter[bid].t->addToList(list, n-1, single, listlen, alphabetSize);
      }
    }
  }
  
  return listlen; 
}

int *CLTreeS::reorder(int *lmer, int *order, int L, int *output){
  for(int i=0;i<L;i++){
    output[i]=lmer[order[i]];
  }
  return(output);
}

void CLTreeS::cloneReorder(CLTreeS *newTree, int *order, int n, int L,int alphabetSize, int *tmpArray, int *tmpArray2){ // reorders and clones
  
  
  if (n==1)
  {
    for(int bid=0;bid<alphabetSize;bid++)
    {
      
      if (this->daughter[bid].t==NULL) continue;
      tmpArray[L-n]=bid;
      tmpArray2=reorder(tmpArray, order, L, tmpArray2);
      newTree->addLTreeSnodeData(tmpArray2, L, this->daughter[bid].node, this->minSeqID, this->maxSeqID);
      
      //for(int k=0;k<L;k++){printf("%c",::globalConverter.icidx[tmpArray[k]]);}
      //printf("->");
      //for(int k=0;k<L;k++){printf("%c",::globalConverter.icidx[tmpArray2[k]]);}
      //printf("\n");
      
      
    }
  }
  else
  {
    for(int bid=0;bid<alphabetSize;bid++)
    {
      if (this->daughter[bid].t != NULL)
      {
        tmpArray[L-n]=bid;
        daughter[bid].t->cloneReorder(newTree, order, n-1, L,alphabetSize, tmpArray, tmpArray2);
      }
    }
  }
  
}

void addmmprof(aint *mmprofile_i,int *nodej_seqIDs_p,int nn, int curnodeid)
{
  for(int j=0;j<nn;j++)
  {
    if (*nodej_seqIDs_p>curnodeid) return;
    ++mmprofile_i[*nodej_seqIDs_p++];
  }
}

void CLTreeS::DFSTnIDL(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, CbinMMtree **curMMtree, int alphabetSize, const KernelContext *ctx)
{
  
  // note: similar to DFST, we should first check if mismatch not allowed, don't iterate over mismatch places and go directly to match
  
  
  int bid;
  for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
  {
    bid = this->nonEmptyDaughterIdxs[ibid];
    //		if(daughter[bid].node==NULL) continue;
    LTreeSnodeData *nodei=daughter[bid].node;
    int nodei_n=nodei->n;
    if (nodei_n==1)
    {
      int curnodeid = nodei->seqIDs.i;
      aint **mmprofile=ctx->mmProfile[curnodeid];
      
      for(int i=0;i<listlen;i++)
      {
        CLTreeS *imatchingLmer =matchingLmers[i];
        CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
        CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;
        
        int fbid;
        for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
        {
          fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
          
          if (bid==fbid)
          {
            if(icurMMtreeMatch!=NULL){
              
              LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
              if (nodej->n==1)
              {
                { int jj = nodej->seqIDs.i; int in = (jj <= curnodeid); mmprofile[curMismatchCnt[i]][in ? jj : curnodeid] += in; } // Phase 6: rows hold j <= i only; branchless (a data-dependent branch here cost 17%)
                
              }
              else
              {
                
                /* for(int j=0;j<nodej->n;j++)
                {
                if (nodej->seqIDs.p[j]>curnodeid) break;
                mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                }*/
                addmmprof(mmprofile[curMismatchCnt[i]],nodej->seqIDs.p,nodej->n, curnodeid);
                
              }
          }
        }
          else
          {
            //                        if (curMismatchCnt[i]<gMAXMM){
            if(icurMMtreeMisMatch!=NULL){
              //it tolerates one more mismatch
              
              LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
              if (nodej->n==1)
              {
                { int jj = nodej->seqIDs.i; int in = (jj <= curnodeid); mmprofile[1+curMismatchCnt[i]][in ? jj : curnodeid] += in; }
                
              }
              else
              {
                /*
                for(int j=0;j<nodej->n;j++)
                {
                if (nodej->seqIDs.p[j]>curnodeid) break;
                mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                }*/
                addmmprof(mmprofile[1+curMismatchCnt[i]],nodej->seqIDs.p,nodej->n, curnodeid);
                
              }
              
            }
            
      }
      }
    }
    }
    
    else
    {
      
      for(int i=0;i<listlen;i++)
      {
        CLTreeS *imatchingLmer =matchingLmers[i];
        CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
        CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;
        
        
        int fbid;
        //      if(icurMMtreeMisMatch!=NULL){
        //it tolerates one more mismatch
        for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
        {
          fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
          
          if (bid==fbid)
          {
            if(icurMMtreeMatch!=NULL){
              
              
              LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
              if (nodej->n==1)
              {
                
                
                for(int k=0;k<nodei_n;k++)
                {
                  int curnodeid = nodei->seqIDs.p[k];
                  aint **mmprofile=ctx->mmProfile[curnodeid];
                  
                  //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                  aint *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                  
                  { int jj = nodej->seqIDs.i; int in = (jj <= curnodeid); mmprofile_curMismatchCnt_i[in ? jj : curnodeid] += in; }
                  
                }
                
                
              }
              else
              {
                /*
                for(int j=0;j<nodej->n;j++)
                {
                if (nodej->seqIDs.p[j]>curnodeid) break;
                mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                }*/
                
                
                for(int k=0;k<nodei_n;k++)
                {
                  int curnodeid = nodei->seqIDs.p[k];
                  aint **mmprofile=ctx->mmProfile[curnodeid];
                  
                  //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                  aint *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                  
                  addmmprof(mmprofile_curMismatchCnt_i,nodej->seqIDs.p,nodej->n, curnodeid);
                }
                
              }
            }
        }
          else
          {
            //                            if (curMismatchCnt[i]<gMAXMM){
            if(icurMMtreeMisMatch!=NULL){
              LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
              int nodej_n=nodej->n;
              if (nodej_n==1)
              {
                
                for(int k=0;k<nodei->n;k++)
                {
                  int curnodeid = nodei->seqIDs.p[k];
                  aint **mmprofile=ctx->mmProfile[curnodeid];
                  
                  aint *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                  
                  { int jj = nodej->seqIDs.i; int in = (jj <= curnodeid); mmprofile_curMismatchCntP1_i[in ? jj : curnodeid] += in; }
                  
                }
              }
              else
              {
                for(int k=0;k<nodei->n;k++)
                {
                  int curnodeid = nodei->seqIDs.p[k];
                  //int **mmprofile=ctx->mmProfile[curnodeid];
                  
                  //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                  //int *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                  addmmprof(ctx->mmProfile[curnodeid][1+curMismatchCnt[i]],nodej->seqIDs.p,nodej_n, curnodeid);
                }
              }
              
            }
            
          }
      }
        /*        }else{
        // no more mismatches. just use the match
        //for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
        LTreeSnodeData *nodej=imatchingLmer->daughter[bid].node;
        
        if (nodej!=NULL){  // if the neighbor has daughter with matching base
        if (nodej->n==1)
        {
        mmprofile_curMismatchCnt_i[nodej->seqIDs.i]++;
        }
        else
        {
        
        
        addmmprof(mmprofile_curMismatchCnt_i,nodej->seqIDs.p,nodej->n, curnodeid);
        
        }
        }
        
        
        
    }*/
        //              }
        
    }
  }
}
  }

void CLTreeS::DFSTiDL( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt,CbinMMtree **curMMtree, int pos, int alphabetSize, const KernelContext *ctx)
{
  if(pos==ctx->LM1) //LM1 is L-1
  {
    DFSTnIDL(matchingLmers, listlen, curMismatchCnt, curMMtree, alphabetSize, ctx); // process the node.
  }
  else
  {
    
    
    int n = listlen * alphabetSize;
    // Phase 6: the three scratch arrays of this level come from a per-thread, per-depth buffer that
    // is reused across visits (the recursion is a DFS: one live frame per depth), instead of three
    // heap allocations per internal node visit.
    static thread_local std::vector<std::vector<char> > dfsBufs;
    if ((int)dfsBufs.size() <= pos) dfsBufs.resize(pos + 1);
    std::vector<char> &buf = dfsBufs[pos];
    size_t need = (size_t)n * (sizeof(CLTreeS*) + sizeof(CbinMMtree*) + sizeof(int)) + 16;
    if (buf.size() < need) buf.resize(need);
    CLTreeS **newlist = (CLTreeS **)buf.data();
    CbinMMtree **newMMtree = (CbinMMtree **)(newlist + n);
    int *newMismatchCnt = (int *)(newMMtree + n);
    
    int newlistlen = 0;
    CLTreeS **newlistnewlistlen = newlist; //&newlist[newlistlen]
    int *newMismatchCntnewlistlen = newMismatchCnt;
    CbinMMtree **newMMtreenewlistlen =newMMtree;
    //int alphabetSize = ::globalConverter.b; // for DNA it is 4
    //		for(int bid=0;bid<alphabetSize;bid++)
    //		{
    int bid;
    for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
    {
      bid = this->nonEmptyDaughterIdxs[ibid];
      //if(daughter[bid].t==NULL) continue;
      newlistlen = 0;
      newlistnewlistlen = newlist;
      newMismatchCntnewlistlen = newMismatchCnt;
      newMMtreenewlistlen =newMMtree;
      int daughter_maxSeqID = daughter[bid].t->maxSeqID;
      for(int i=0;i<listlen;i++)
      {
        //CLTreeSptr *imatchingLmer =matchingLmers[i];
        CLTreeS *imatchingLmer =matchingLmers[i];
        
        CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
        CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;
        
        int fbid;
        
        if(icurMMtreeMisMatch!=NULL){
          //it tolerates one more mismatch
          
          for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
          {
            fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
            
            if (bid==fbid){  // Match
              if(icurMMtreeMatch!=NULL){ // it allows match
                CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
                if (newnode->minSeqID >daughter_maxSeqID) continue;
                
                *newlistnewlistlen=newnode;
                newlistnewlistlen++;
                *newMismatchCntnewlistlen=curMismatchCnt[i];
                newMismatchCntnewlistlen++;
                *newMMtreenewlistlen=icurMMtreeMatch;
                newMMtreenewlistlen++;
                newlistlen++;
              }
            } else {
              CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
              if (newnode->minSeqID >daughter_maxSeqID) continue;
              
              *newlistnewlistlen=newnode;
              newlistnewlistlen++;
              *newMismatchCntnewlistlen=curMismatchCnt[i]+1;
              newMismatchCntnewlistlen++;
              *newMMtreenewlistlen=icurMMtreeMisMatch;
              newMMtreenewlistlen++;
              newlistlen++;
              // note : infact mismatchCnt is redundant, we can use the tree for final mismatch cnt. and we don't need to add new var, use intptr for a child
            }
            
          }
          
        }else{
          //no more mismatches tolerated at this depth
          //if(icurMMtreeMatch!=NULL){  not needed because we assume it is not possible both child0 and 1 be null at internal nodes (depth<L)
          CLTreeS *newnode =imatchingLmer->daughter[bid].t;
          
          if (newnode!=NULL){  // if the neighbor has daughter with matching base
            if (newnode->minSeqID >daughter_maxSeqID) continue;
            
            *newlistnewlistlen=newnode;
            newlistnewlistlen++;
            *newMismatchCntnewlistlen=curMismatchCnt[i];
            newMismatchCntnewlistlen++;
            newlistlen++;
            *newMMtreenewlistlen=icurMMtreeMatch;
            newMMtreenewlistlen++;
          }
        }
      }
      
      if (newlistlen!=0)
      {
        daughter[bid].t->DFSTiDL(newlist,newlistlen,newMismatchCnt,newMMtree, pos+1,alphabetSize, ctx);
      }
    }
    //delete []newlist;
    //delete []newMismatchCnt;
    
    // (scratch arrays are per-thread, per-depth buffers now: nothing to free)


    
  }
}

void CLTreeS::deleteTree(int n, int alphabetSize, int DontDeleteNodeData)
{
  if (n>1)
  {
    for (int i=0;i<alphabetSize;i++)
    {
      if (daughter[i].t!=NULL)
      {
        daughter[i].t->deleteTree(n-1,alphabetSize, DontDeleteNodeData);
        delete daughter[i].t;
      }
    }
  }
  if (n==1)
  {
    if (!DontDeleteNodeData){
      for (int i=0;i<alphabetSize;i++)
      {
        if (daughter[i].node!=NULL)
        {
          if (daughter[i].node->n>1)
          {
            delete []daughter[i].node->seqIDs.p;
          }
          delete daughter[i].node;
        }
      }
    }
  }
  
  for(int i=0;i<MAX_ALPHABET_SIZE;i++){
    daughter[i].t=NULL;
  }
  maxSeqID=0;
  minSeqID=0;
}

int CLTreeS::addSequence(int *bid, int n, int L, int seqID)  //adds all the L-subseqs 
{
  n = n-L+1;
  if (n<0) n=0;
  for(int i=0;i<n;i++)
  {
    addSeq(bid,L,bid, seqID);
    bid++;
  }
  return n; 
}

int CLTreeS::leavesCount(int withMultiplicity, int n, int alphabetSize, int *nodesAtDepth)  //returns the number of sequences in the tree.  //call with n=L from outside
{
  if(nodesAtDepth!=NULL){
    (*nodesAtDepth)++;
    nodesAtDepth++;
  }
  int nleaves = 0; 
  for (int i=0;i<alphabetSize;i++)
  {
    if (daughter[i].t!=NULL)
    {
      if (n==1)
      {
        if (withMultiplicity)
        {
          nleaves +=daughter[i].node->n; 
        }
        else 
        {
          nleaves++; 					
        }
      }
      else
      {
        nleaves+=daughter[i].t->leavesCount(withMultiplicity, n-1,alphabetSize, nodesAtDepth/*+1*/);
      }
    }
  }
  return nleaves; 
}

} // namespace GKM_NS
