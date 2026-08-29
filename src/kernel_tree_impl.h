/* kernel_tree_impl.h : the tree (suffix-trie) kernel driver, compiled once per alphabet size
 * (see trie_b4.cpp / trie_b32.cpp and the note in global.h). Was the second half of mainGkmKernel.cpp.
 */
#include <iostream>
#include <thread>
#include <atomic>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "globalvar.h"
#include "gkmOptions.h"
#include "gkmMainHelpers.h"
#include "SequenceSet.h"
#include "KernelFile.h"
#include "Sequence.h"
#include "CalcWmML.h"
#include "LTreeS.h"
#include "CiDLPasses.h"

namespace GKM_NS {

double calcinnerprod(int i, int j, double *c, const KernelContext &kc)
{
	double res = 0; 
	for(int m=0;m<=kc.maxmm;m++)
	{
		res+=kc.mmProfile[i][m][j]*c[m]; 
	}
	return(res); 
}


double calcinnerprod(int i, int j, double *c, double n0, double C, int nA, int nB, double btL, const KernelContext &kc) // gives inner prodict of the pseudo-counts . nA is the number of L-mers in A and is equal to length(A)-L+1, btL is b^L
{
	double res = 0; 
	for(int m=0;m<=kc.maxmm;m++)
	{
		res+=kc.mmProfile[i][m][j]*c[m]; 
	}

	res = res+(nA+nB)*n0*C+btL*n0*n0; 
	return(res); 
}

// One worker: takes passes from the shared counter until none are left (Phase 6: passes have
// unequal cost, so a static j % nThreads assignment left threads idle).
void task1(int L, CiDLPasses *iDL, CLTreeS *seqsTS, int M, std::atomic<int> *nextPass, const KernelContext *kc){
  for(;;){
    int j = nextPass->fetch_add(1);
    if (j >= M) break;
    
    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    //for(int j=0;j<iDL->M;j++){
    CLTreeS *seqsTSj= new CLTreeS();
    seqsTS->cloneReorder(seqsTSj, iDL->passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    //gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    //gDFSMMlist[0][0] = 0;
    //if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
    //    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL->passTrees+j, 0, globalConverter.b);
    int zero=0;
    if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
      seqsTSj->DFSTiDL(&seqsTSj,1, &zero, iDL->passTrees+j, 0, globalConverter.b, kc);
    seqsTSj->deleteTree(L, globalConverter.b, 1);
    delete seqsTSj;
    
    // print mismatch profile:
    /* for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=kc.maxmm; dd++){
    printf("%d ",kc.mmProfile[si][dd][sj]);
    }
    }
  }
    */
    
    
    //}
    delete []tmpArray1;
    delete []tmpArray2;
    //snprintf(globtmpstr, GKM_TMPSTR_LEN,"ended pass %d out of %d.\n",j+1,iDL->M);Printf(globtmpstr);
    }
  
  }

int gkmKernelSuffixTree(OptsGkmKernel &opt)  //maingKernel
{
  
  int L = opt.L;
  int K = opt.K;
  int maxseqlen =	opt.maxseqlen;
  int useTgkm = opt.useTgkm;
  int maxnmm = opt.maxnmm; //auto 
  int nMAXSEQUENCES = opt.maxnumseq; 
  bool addRC = opt.addRC;
  bool usePseudocnt= opt.usePseudocnt; 
  
  const char *posSeqsFN = opt.posfile.c_str();
  const char *negSeqsFN = opt.negfile.c_str();
  const char *outFN = opt.outfile.c_str();
  
  int i = 1; 
  
  //char tmps[1000];
  
  char **seqname = new char *[nMAXSEQUENCES];
  
  CCalcWmML wmc(L, K, globalConverter.b);
  //double *kernel = wmc.kernelTruncated;
  if (maxnmm==-1)
  { 
    maxnmm=L;
    if (useTgkm==1)
    {
      maxnmm = 2*(wmc.kernelTruncatedLength-1);
      if (maxnmm>L) 
      {
        maxnmm=L;
      }
    }
    if (useTgkm==2)
    {
      maxnmm = L-K; 
    }		
    if (useTgkm==3)  //wildcard kernel 
    {
      maxnmm = opt.wildcardMismatchM;
    }
    if (useTgkm==4)  //mismatch kernel 
    {
      maxnmm = 2*opt.wildcardMismatchM;
    }
    if (maxnmm>L)
    {
      maxnmm=L;
    }
  }
  double n0 = wmc.n0; 
  (void)n0; // computed for symmetry with the other paths; not used here
  double *c = wmc.cTr; 
  
  n0 = c[maxnmm]/2; 
  
  if (!useTgkm)
  {
    n0 = 0; 
    //kernel = wmc.kernel; 
    c = wmc.c; // same as kernel
  }
  if (useTgkm==2)
  {
    //	n0 = 0; 
    //	kernel = wmc.kernel; 
    c = wmc.h; 
    n0 = c[maxnmm]/2;
    
  }
  if (useTgkm==3)  //wildcard kernel 
  {
    c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, opt.wildcardLambda, c);
    n0 = c[maxnmm]/2;
    
  }
  if (useTgkm==4)  //mismatch kernel 
  {
    c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, c);
    n0 = c[maxnmm]/2;
    
  }
  
  snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);
  for(int ii=0;ii<=maxnmm;ii++) {
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n c[%d] = %e",ii,c[ii] ); 	Printf(globtmpstr);
  }
  Printf("\n");
  
  int npos=0; 
  int nneg=0;
  
  KernelContext kc; // per-call state (was the kc.mmProfile/kc.maxmm/kc.LM1 globals)
  kc.maxmm=maxnmm; //MaxMismatch
  
  CLTreeS *seqsTS= new CLTreeS();
  
  
  
  int **seqsB = new int *[nMAXSEQUENCES]; 
  int **seqsBrc  = new int *[nMAXSEQUENCES]; 
  
  int *LmersCnt = new int [nMAXSEQUENCES]; 
  
  CSequence sgii(maxseqlen+3);
  CSequence *sgi = &sgii;
  
  int ntotal = 0; //number of lmers
  int nseqs=0;
  
  char**seqname2= NULL; 
  
  seqname2 = new char *[nMAXSEQUENCES];
  std::vector<SeqRecord> records; // row identity, written to <outfile>.index
  
  //read positive sequence file
  FILE *sfi = fopen(posSeqsFN, "r"); 
  if (sfi == NULL) return gkmCannotOpen(posSeqsFN);
  
  char *tmpSeq=new char[maxseqlen+3];
  int  *tmpSeqB=new int[maxseqlen+3];
  
  int multiseq_allowed = opt.mergeByName ? 1 : 0; // Phase 3: default one record = one row
  
  
  while (!feof(sfi))
  {
    
    if (sgi->readFsa(sfi) < 0) return 1; 
    
    if(sgi->getLength()>0)
    {
      int dupl_seq_idx= -1;
      if(multiseq_allowed){
        dupl_seq_idx = find_str(seqname,nseqs, sgi->getName());
      }
      if (dupl_seq_idx>=0){
        records[dupl_seq_idx].length += sgi->getLength();
        sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
        ntotal-=LmersCnt[dupl_seq_idx];
        LmersCnt[dupl_seq_idx]+= seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        ntotal+=LmersCnt[dupl_seq_idx];
        if(addRC)
        {
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
          LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        }
      }else{
        if (nseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
        seqname2[nseqs] = new char[strlen(sgi->getName())+1]; // was a fixed 100 bytes: names >= 100 chars overflowed it
        strcpy(seqname2[nseqs], sgi->getName()); // buffer sized strlen+1 above
        seqname[nseqs]=seqname2[nseqs]; 
        
        seqsB[nseqs] = new int [sgi->getLength()]; 
        sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]); 
        LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
        
        if(addRC)
        {
          seqsBrc[nseqs] = new int [sgi->getLength()];
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]); 
          LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
        }
        else
        {
          seqsBrc[nseqs]=NULL; 
        }
        
        ntotal = ntotal + LmersCnt[nseqs]; 
        records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), 0});
        nseqs++;
      }
    }
  }
  fclose(sfi);
  npos = nseqs;
  
  //read negative sequence file
  sfi = fopen(negSeqsFN, "r"); 
  if (sfi == NULL) return gkmCannotOpen(negSeqsFN);
  while (!feof(sfi))
  {
    if (sgi->readFsa(sfi) < 0) return 1; 
    
    if(sgi->getLength()>0)
    {
      // check if segName already exists
      int dupl_seq_idx= -1;
      if(multiseq_allowed){
        dupl_seq_idx = find_str(seqname+npos,nseqs-npos, sgi->getName());
      }
      if (dupl_seq_idx>=0){
        dupl_seq_idx+=npos;
        records[dupl_seq_idx].length += sgi->getLength();
        sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
        ntotal-=LmersCnt[dupl_seq_idx];
        LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        ntotal+=LmersCnt[dupl_seq_idx];
        if(addRC)
        {
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
          LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        }
      }else{
        if (nseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
        seqname2[nseqs] = new char[strlen(sgi->getName())+1];
        strcpy(seqname2[nseqs], sgi->getName()); // buffer sized strlen+1 above
        seqname[nseqs]=seqname2[nseqs];
        
        seqsB[nseqs] = new int [sgi->getLength()];
        sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]);
        LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
        if(addRC)
        {
          seqsBrc[nseqs] = new int [sgi->getLength()];
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]);
          LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
          
        }
        else
        {
          seqsBrc[nseqs]=NULL; 
        }
        
        ntotal = ntotal + LmersCnt[nseqs]; 
        records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), 0});
        nseqs++; 
        
      }
    }
  }
  fclose(sfi);
  
  delete []tmpSeq;
  delete []tmpSeqB;
  
  nneg = nseqs - npos;
  for(int i=0;i<nseqs;i++) { records[i].id = seqRecordId(i, npos); records[i].label = (i<npos) ? 1 : -1; records[i].nlmers = LmersCnt[i]; }
  
  // global vars init: 
  kc.LM1=L-1;
  kc.maxmm=maxnmm; //MaxMismatch
  // Phase 6: the mismatch profile is only ever read for j <= i (calcinnerprod is called for i >= j and
  // the leaf code only records pairs with j <= i), so row i holds i+1 entries: half the memory.
  // Rows are allocated per tile below (Phase 6: peak memory bounded by tileRows / tileMemoryMB).
  kc.mmProfile=new aint **[nseqs];
  for(int i=0;i<nseqs;i++) kc.mmProfile[i] = NULL;
  
  int *nodesAtDepthCnt = new int[L];
  for(int i=0;i<L; i++){
    nodesAtDepthCnt[i]=0;
  }
  
  int uniqueLmerCnt = seqsTS->leavesCount(0,L, globalConverter.b, nodesAtDepthCnt);
  snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n npos %d \n nneg %d \n  ntotal %d \n nunique %d\n",npos,nneg,ntotal,uniqueLmerCnt);Printf(globtmpstr);
    // if no IDL bound
    /*
    seqsTS->DFST(gDFSlistT[0],1, gDFSMMlist[0], 0, globalConverter.b);
    
    for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=kc.maxmm; dd++){
    printf("%d ",kc.mmProfile[si][dd][sj]);
    }
    }
    }
    */
    // else if IDL bound then
    for(int i=0;i<L; i++){
      //     snprintf(globtmpstr, GKM_TMPSTR_LEN,"d%d , %d\n", i, nodesAtDepthCnt[i]);Printf(globtmpstr);
    }
    
    CiDLPasses iDL;
    //iDL.newIDLPasses(L, kc.maxmm);
    double p=1.0/::globalConverter.b;
    
    //iDL.initPassOrderAll(L, kc.maxmm);
    //iDL.newGreedyIDLPasses(L,iDL.M,  kc.maxmm, nodesAtDepthCnt, p);
    
    iDL.newGreedyIDLPasses(L,2*L,  kc.maxmm, nodesAtDepthCnt, p);
    
    //iDL.newGreedy2IDLPasses(L,2*L,  kc.maxmm, nodesAtDepthCnt, p);
    //iDL.newGreedy2IDLPasses(L,L,  kc.maxmm, nodesAtDepthCnt, p);
    
    //iDL.newPassOrderDesignCover( L, kc.maxmm, 3);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
    //iDL.newGreedyIDLPasses(L,iDL.M,  kc.maxmm, nodesAtDepthCnt, p);
    
    
    /*    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    for(int j=0;j<iDL.M;j++){
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"pass %d out of %d.\n",j+1,iDL.M);Printf(globtmpstr);
    CLTreeS *seqsTSj= new CLTreeS();
    seqsTS->cloneReorder(seqsTSj, iDL.passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    gDFSMMlist[0][0] = 0;
    if(!((iDL.passTrees[j]->child0==NULL)&&(iDL.passTrees[j]->child1==NULL))) // i.e. if not empty tree
    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    seqsTSj->deleteTree(L, globalConverter.b, 1);
    delete seqsTSj;
    
    // print mismatch profile:
    / * for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=kc.maxmm; dd++){
    printf("%d ",kc.mmProfile[si][dd][sj]);
    }
    }
    }
    * /
    }
    delete []tmpArray1;
    delete []tmpArray2;
    */
    
    //myThreads[0] = std::thread(task1, L, 0, &iDL, seqsTS);
    //myThreads[1] = std::thread(task1, L, 1, &iDL, seqsTS);
    //myThreads[0].join();
    //myThreads[1].join();
    
  delete []nodesAtDepthCnt;
  
  /*
  for(int i=0;i<nseqs;i++)
  {
  for(int k=0;k<nseqs;k++)
  {
  for (int j=0;j<=kc.maxmm;j++)
  {
  printf("(%d,%d)[%d] = %d\n",i,k,j, kc.mmProfile[i][j][k]);
  }
  }
  }
  */	
  
  /// calc C 
  double C =0; 
  for(int m=0;m<=L;m++)
  {
    C+=dCombinations(L,m)*pow(1.0*globalConverter.b-1,m)*wmc.kernelTruncated[m];
    //		C+=dCombinations(L,m)*pow(3.0,m)*wmc.kernelTruncated[m];
  }
  
  //	double btL=pow(4.0,L);
  double btL=pow(1.0*globalConverter.b,L);
  double *norm = new double [nseqs]; 
  
  GkmkWriter bin;
  FILE *fo = NULL;
  if (opt.OutputBinary) { gkmFillGkmkHeader(bin, opt, maxnmm, nseqs, npos); bin.values.reserve((size_t)nseqs*(nseqs+1)/2); }
  else { fo = fopen(outFN, "w"); if (fo == NULL) return gkmCannotOpen(outFN); }
  /*
  if (OutputMismatchProfileOnly)
  {
  fprintf(fo, "%d\tL (length)\n", L); 
  fprintf(fo, "%d\td (maximum number of mismatches)\n", kc.maxmm); 
  fprintf(fo, "%d\tNp (number of sequences in positive class)\n", npos);
  fprintf(fo, "%d\tNn (number of sequences in negative class)\n", nneg); 
  
  for (int nmm=0;nmm<=kc.maxmm;nmm++)
  {
  fprintf(fo, "d=%d\n",nmm);
  
  
  for(i=0;i<nseqs;i++)
  {
  
  if (outputClassLabel)
  {
  fprintf(fo, "%d\t", (i<npos)?1:-1);
  }
  if (OutputSeqNames)
  {
  fprintf(fo, "%s\t", seqname[i]);
  }
  
  for(int j=0;j<nseqs;j++)
  {
  if(i>=j)
  {
  fprintf(fo, "%d\t",kc.mmProfile[i][nmm][j]);
  }
  }
  fprintf(fo, "\n"); 
  }
  }
  }
  else 
  */

  // ---- Phase 6: tiles of rows. The mismatch profile of rows [lo, hi] is built by a full set of
  // passes with the DFS pruned to that band (node id ranges), the rows are written, the tile is freed.
  // Counts are integers, so the result is identical to the untiled computation.
  int tileRows = opt.tileRows;
  if (tileRows <= 0) {
    double bytesPerRow = (double)(kc.maxmm + 1) * sizeof(aint) * nseqs; // upper bound (the last row is the longest)
    double budget = (double)opt.tileMemoryMB * 1048576.0;
    tileRows = (int)(budget / (bytesPerRow > 0 ? bytesPerRow : 1));
    if (tileRows < 1) tileRows = 1;
  }
  if (tileRows > nseqs) tileRows = nseqs;
  int ntiles = (nseqs + tileRows - 1) / tileRows;
  if (ntiles > 1) { snprintf(globtmpstr, GKM_TMPSTR_LEN, "Computing the kernel in %d tiles of %d rows.\n", ntiles, tileRows); Printf(globtmpstr); }
  // One contiguous block, sized for the largest (last) tile and reused by every tile: freeing and
  // reallocating a block per tile does not lower the resident set (the allocator keeps freed large
  // blocks, and each tile's block is bigger than the previous one), a single reused block does.
  size_t maxTileCounters = 0;
  for (int tile = 0; tile < ntiles; tile++) {
    int lo = tile * tileRows, hi = lo + tileRows - 1; if (hi > nseqs - 1) hi = nseqs - 1;
    size_t c = 0; for(int i=lo;i<=hi;i++) c += (size_t)(kc.maxmm+1) * (i+1);
    if (c > maxTileCounters) maxTileCounters = c;
  }
  aint *tileBlock = new aint[maxTileCounters > 0 ? maxTileCounters : 1];
  for (int tile = 0; tile < ntiles; tile++)
  {
    int lo = tile * tileRows, hi = lo + tileRows - 1; if (hi > nseqs - 1) hi = nseqs - 1;
    kc.rowLo = lo; kc.rowHi = hi;
    size_t tileCounters = 0;
    for(int i=lo;i<=hi;i++) tileCounters += (size_t)(kc.maxmm+1) * (i+1);
    for(size_t k=0;k<tileCounters;k++) tileBlock[k]=0;
    {
      size_t off = 0;
      for(int i=lo;i<=hi;i++)
      {
        kc.mmProfile[i] = new aint*[kc.maxmm+1];
        for (int j=0;j<=kc.maxmm;j++) { kc.mmProfile[i][j] = tileBlock + off; off += (size_t)(i+1); }
      }
    }
    int nThreads=std::thread::hardware_concurrency();
    if(nThreads==0){nThreads = iDL.M;}
    if(nThreads>iDL.M){nThreads =iDL.M;}
    if(opt.maxnThread<nThreads){nThreads=opt.maxnThread;}
    if(nThreads<1){nThreads=1;}
    std::atomic<int> nextPass(0);
    
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"Running %d passes on %d thread%s.\n", iDL.M, nThreads, (nThreads==1)?"":"s"); Printf(globtmpstr);
    if (nThreads<=1){
      task1( L, &iDL, seqsTS, iDL.M, &nextPass, &kc);
    }else{
      
#ifndef MULTI_THREAD_SAFE
      Printf("Warning -- MULTI_THREAD_SAFE is not enabled (see src/global.h). Some values may be approximated.\n");
#endif
      
      
      std::thread *myThreads = new std::thread[nThreads];
      int j;
      
      for(j=0;j<nThreads;j++){
        myThreads[j] = std::thread(task1, L, &iDL, seqsTS, iDL.M, &nextPass, &kc);
        // myThreads[j].join();
      }
      for(j=0;j<nThreads;j++){
        myThreads[j].join();
      }
      delete []myThreads;
    }
    
    for(i=lo;i<=hi;i++)
    {
      if (usePseudocnt)
      {
        norm[i] = sqrt(calcinnerprod(i,i,c,n0,C,LmersCnt[i], LmersCnt[i], btL, kc));
      }
      else
      {
        norm[i] = sqrt(calcinnerprod(i,i,c,kc));
      }
    }
    
    for(i=lo;i<=hi;i++)
    {
      
      //if (outputClassLabel)
      //{
      //	fprintf(fo, "%d\t", (i<npos)?1:-1);
      //}
      //if (OutputSeqNames)
      //{
      //	fprintf(fo, "%s\t", seqname[i]);
      //}
      
      for(int j=0;j<=i;j++)
      {
        double v;
        if (i==j) v = 1.0;
        else if (usePseudocnt) v = (norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c, n0,C,LmersCnt[i], LmersCnt[j], btL, kc)/(norm[i]*norm[j]);
        else v = (norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c,kc)/(norm[i]*norm[j]);
        if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
        else bin.add(v);
      }
      if (fo) fprintf(fo, "\n"); 
    }
    for(int i=lo;i<=hi;i++)
    {
      delete []kc.mmProfile[i];
      kc.mmProfile[i] = NULL;
    }
  }
  delete []tileBlock;
  if (fo) fclose(fo); 
  else if (bin.write(opt.outfile, records) != 0) return gkmCannotOpen(outFN);
  if (writeIndexSidecar(opt.outfile, records) != 0) { snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n WARNING: could not write %s.index\n", outFN); Printf(globtmpstr); }
  
  delete []norm;
  delete []LmersCnt;
  seqsTS->deleteTree(L, globalConverter.b, 0);
  delete seqsTS;
  //delete []curmmcnt; 
  
  for(int i=0;i<nseqs;i++)
  {//printf("\n4 %d\n",i);
    delete []seqsB[i]; 
    if (seqsBrc[i]!=NULL) delete []seqsBrc[i]; 
  }
  delete []kc.mmProfile;
  
  delete []seqname; 
  if (seqname2!=NULL)
  {
    for(i=0;i<nseqs;i++)
    {
      delete []seqname2[i]; 
    }
    delete []seqname2; 
  }
  
  return 0; 
}

} // namespace GKM_NS
