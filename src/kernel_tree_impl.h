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
#include "GeneralB.h"
#include "MultiTrack.h"
#include <unordered_map>

namespace GKM_NS {

double calcinnerprod(int i, int j, const double *c, const KernelContext &kc)
{
	double res = 0; 
	for(int m=0;m<kc.nclasses;m++) // Phase 7: rows are mismatch classes; unreachable ones are NULL (for one alphabet: m = 0..maxmm, all present)
	{
		const aint *row = kc.mmProfile[i][m];
		if (row) res+=row[j]*c[m]; 
	}
	return(res); 
}


double calcinnerprod(int i, int j, const double *c, double n0, double C, int nA, int nB, double btL, const KernelContext &kc) // gives inner prodict of the pseudo-counts . nA is the number of L-mers in A and is equal to length(A)-L+1, btL is b^L
{
	double res = 0; 
	for(int m=0;m<kc.nclasses;m++) // Phase 7: rows are mismatch classes; unreachable ones are NULL (for one alphabet: m = 0..maxmm, all present)
	{
		const aint *row = kc.mmProfile[i][m];
		if (row) res+=row[j]*c[m]; 
	}

	res = res+(nA+nB)*n0*C+btL*n0*n0; 
	return(res); 
}

// One worker: takes passes from the shared counter until none are left (Phase 6: passes have
// unequal cost, so a static j % nThreads assignment left threads idle).
void task1(int L, CiDLPasses *iDL, CLTreeS *seqsTS, int M, std::atomic<int> *nextPass, const KernelContext *kc, int b){
  for(;;){
    int j = nextPass->fetch_add(1);
    if (j >= M) break;
    
    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    //for(int j=0;j<iDL->M;j++){
    // Phase 7: prefix-split passes use the identity order -> traverse the shared trie itself (read-only)
    CLTreeS *seqsTSj = seqsTS;
    if (!iDL->identityOrder) { seqsTSj = new CLTreeS(); seqsTS->cloneReorder(seqsTSj, iDL->passOrder[j], L,L,b, tmpArray1, tmpArray2); }
    // Phase 7: class-index increment per depth of the reordered trie = that of the original position placed there
    std::vector<int> step(L);
    for(int d=0; d<L; d++) step[d] = kc->stepPos[iDL->passOrder[j][d]];
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, b);
    //gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    //gDFSMMlist[0][0] = 0;
    //if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
    //    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL->passTrees+j, 0, b);
    int zero=0;
    if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
      seqsTSj->DFSTiDL(&seqsTSj,1, &zero, iDL->passTrees+j, 0, b, kc, step.data());
    if (!iDL->identityOrder) { seqsTSj->deleteTree(L, b, 1); delete seqsTSj; }
    
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
    //gkmMsg("ended pass %d out of %d.\n",j+1,iDL->M);
    }
  
  }

// Everything after the trie is built (shared by the single-alphabet and the multi-track drivers):
// pass design, the tiled mismatch-class profile, norms, the kernel matrix (text or .gkmk) and the
// .index sidecar. kc must have LM1, maxmm, nclasses, stepPos set and mmProfile allocated (rows NULL);
// reach = the class indices whose rows are recorded; c = the coefficient per class index.
static int gkmKernelPassesAndOutput(OptsGkmKernel &opt, KernelContext &kc, CLTreeS *seqsTS, int L, int b, double p,
                                    int nseqs, int npos, int nneg, int ntotal, const int *LmersCnt,
                                    const std::vector<SeqRecord> &records, const double *c, const std::vector<int> &reach,
                                    bool usePseudocnt, double n0, double C, double btL, int maxnmm,
                                    int hdrB, const std::string &hdrAlphabet, double coefDen)
{
  int i;
  const size_t nrowsUsed = reach.size();
  int *nodesAtDepthCnt = new int[L];
  for(int i=0;i<L; i++){
    nodesAtDepthCnt[i]=0;
  }
  
  int uniqueLmerCnt = seqsTS->leavesCount(0,L, b, nodesAtDepthCnt);
  gkmMsg("\n npos %d \n nneg %d \n  ntotal %d \n nunique %d\n",npos,nneg,ntotal,uniqueLmerCnt);
    // if no IDL bound
    /*
    seqsTS->DFST(gDFSlistT[0],1, gDFSMMlist[0], 0, b);
    
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
      //     gkmMsg("d%d , %d\n", i, nodesAtDepthCnt[i]);
    }
    
    CiDLPasses iDL;
    //iDL.newIDLPasses(L, kc.maxmm);
    // p (probability of a match) is a parameter: 1/b for one alphabet
    
    //iDL.initPassOrderAll(L, kc.maxmm);
    //iDL.newGreedyIDLPasses(L,iDL.M,  kc.maxmm, nodesAtDepthCnt, p);
    
    // Phase 7: the greedy iDL design materialises every mismatch pattern with <= maxmm ones; for long
    // words (multi-track, large -d) the prefix-split design enumerates them implicitly instead
    // maxmm == 0 (exact word matching, e.g. K > L with the truncated filter): the greedy design divides by
    // Dmax in its pass orders (UB; SIGFPE on Linux), the prefix-split DAG handles it
    // The thread budget (computed once; also decides the pass design below).
    int nThreadsBudget = (int)std::thread::hardware_concurrency();
    if (nThreadsBudget <= 0) nThreadsBudget = 1;
    if (opt.maxnThread < nThreadsBudget) nThreadsBudget = opt.maxnThread;
    if (nThreadsBudget < 1) nThreadsBudget = 1;
    // -P 0 (automatic): prefix-split for long words (pattern table too large), for maxmm == 0, and -- measured
    // 2026-08-30/31, dev/PERFORMANCE_PLAN.md section 5 -- for words of >= 16 positions whenever the mismatch count
    // is capped (-d >= 0) and at least 4 threads are available: it traverses the shared trie (no per-pass clone)
    // and its 64 passes balance far better (2.8-3.9x faster wall on 28 threads at l = 20-24, d = 4-6, identical
    // kernel). Greedy stays the choice for short words (l = 10 DNA, d = 3: prefix-split was 28 % slower at 4 000
    // sequences), for 1-3 threads (20-35 % less total CPU) and for -d -1.
    bool prefixSplit = (opt.passDesign == 2) || (opt.passDesign == 0 && (kc.maxmm == 0 || CiDLPasses::patternCount(L, kc.maxmm) > GKM_MAX_PATTERN_TABLE
                                                                         || (opt.maxnmm >= 0 && L >= 16 && nThreadsBudget >= 4)));
    if (prefixSplit) {
      int q = 6; if (q > L) q = L;
      iDL.newPrefixSplitPasses(L, kc.maxmm, q);
      gkmMsg("Long words: %d prefix-split passes (%.3g mismatch patterns).\n", iDL.M, CiDLPasses::patternCount(L, kc.maxmm));
    } else if (opt.passDesign == 4) {
      iDL.newGreedyIDLPasses(L,2*L,  kc.maxmm, nodesAtDepthCnt, p);   // legacy cost model: one mean match probability, identity-order trie widths
    } else {
      // A2b: cost model from the actual l-mers (per-position match probabilities, per-order trie widths);
      // a stride sample of at most ~2M distinct l-mers keeps this to a few seconds
      std::vector<int> lmers; std::vector<int> prefix(L);
      long stride = (uniqueLmerCnt > 2000000) ? (uniqueLmerCnt + 1999999) / 2000000 : 1, counter = 0;
      lmers.reserve((size_t)(uniqueLmerCnt / stride + 1) * L);
      seqsTS->collectLmers(L, L, b, prefix.data(), lmers, stride, counter);
      iDL.newGreedyIDLPassesB(L, 2*L, kc.maxmm, lmers, b);
    }
    
    //iDL.newGreedy2IDLPasses(L,2*L,  kc.maxmm, nodesAtDepthCnt, p);
    //iDL.newGreedy2IDLPasses(L,L,  kc.maxmm, nodesAtDepthCnt, p);
    
    //iDL.newPassOrderDesignCover( L, kc.maxmm, 3);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
    //iDL.newGreedyIDLPasses(L,iDL.M,  kc.maxmm, nodesAtDepthCnt, p);
    
    
    /*    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    for(int j=0;j<iDL.M;j++){
    gkmMsg("pass %d out of %d.\n",j+1,iDL.M);
    CLTreeS *seqsTSj= new CLTreeS();
    seqsTS->cloneReorder(seqsTSj, iDL.passOrder[j], L,L,b, tmpArray1, tmpArray2);
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, b);
    gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    gDFSMMlist[0][0] = 0;
    if(!((iDL.passTrees[j]->child0==NULL)&&(iDL.passTrees[j]->child1==NULL))) // i.e. if not empty tree
    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, b);
    seqsTSj->deleteTree(L, b, 1);
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
  
  double *norm = new double [nseqs]; 
  
  GkmkWriter bin;
  FILE *fo = NULL;
  if (opt.OutputBinary) { gkmFillGkmkHeader(bin, opt, maxnmm, nseqs, npos, hdrB, hdrAlphabet); bin.values.reserve((size_t)nseqs*(nseqs+1)/2); }
  else { fo = fopen(opt.outfile.c_str(), "w"); if (fo == NULL) return gkmCannotOpen(opt.outfile.c_str()); }
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

  // ---- Accumulation mode (Plan A1, dev/PERFORMANCE_PLAN.md section 2). The kernel is a linear function of the
  // mismatch-class counts, K[i][j] = sum_class c[class] N[i][class][j], and for the full filter and the gkm counts
  // the coefficients scaled by coefDen (prod b_i) are integers: the leaf can add c~[class] straight into one int64
  // lower triangle per thread -- no per-class counters, hence no tiles and no repeated traversal -- and the result
  // is exact (integer) and independent of the thread count. Not applicable with pseudocounts, or when the scaled
  // coefficients are not integers (truncated filter) or could overflow int64.
  bool kernelMode = false;
  std::vector<long long> kcoef;
  {
    const char *why = NULL;
    if (opt.accMode == 1) why = "-K 1";
    else if (usePseudocnt) why = "pseudocounts (-p)";
    else if (!(coefDen > 0)) why = "no coefficient scale";
    else {
      kcoef.assign(kc.nclasses, 0);
      double maxAbs = 0;
      for (size_t q = 0; q < nrowsUsed; q++) {
        double v = c[reach[q]] * coefDen, r = (double)std::llround(v);
        // exact integers only (relative 1e-12, i.e. double rounding): the truncated filter's least-squares coefficients
        // are near-integers at ~1e-8 and would change the output in the 7th digit if rounded
        if (fabs(v - r) > 1e-12 * (fabs(r) > 1.0 ? fabs(r) : 1.0) || fabs(r) >= 9007199254740992.0) { why = "the scaled coefficients are not exact integers (truncated filter or wildcard/mismatch weights)"; break; }
        kcoef[reach[q]] = (long long)r;
        if (fabs(r) > maxAbs) maxAbs = fabs(r);
      }
      double maxw = 0; for (int i = 0; i < nseqs; i++) if (LmersCnt[i] > maxw) maxw = LmersCnt[i];
      if (!why && maxAbs * maxw * maxw > 4.0e18) why = "the accumulators could overflow int64 (too many windows per record for these coefficients)";
    }
    kernelMode = (why == NULL);
    if (!kernelMode && opt.accMode == 2) gkmMsg("Kernel-mode accumulation (-K 2) is not applicable here: %s; using the mismatch-class profile.\n", why);
  }

  if (kernelMode)
  {
    const size_t tri = (size_t)nseqs * (nseqs + 1) / 2;
    int nThreads = nThreadsBudget;
    if (nThreads > iDL.M) nThreads = iDL.M;
    if (nThreads < 1) nThreads = 1;
    double triBytes = (double)tri * sizeof(std::atomic<long long>);
    int copies = (int)((double)opt.accMemoryMB * 1048576.0 / (triBytes > 0 ? triBytes : 1));
    if (copies < 1) copies = 1;
    if (copies > nThreads) copies = nThreads;
    std::atomic<long long> *acc = new std::atomic<long long>[tri * copies];
    for (size_t x = 0; x < tri * copies; x++) acc[x].store(0, std::memory_order_relaxed);
    std::vector<KernelContext> kcs(nThreads, kc);
    for (int t = 0; t < nThreads; t++) {
      kcs[t].rowLo = 0; kcs[t].rowHi = nseqs - 1;
      kcs[t].kacc = acc + (size_t)(t % copies) * tri;
      kcs[t].kcoef = kcoef.data();
      kcs[t].kaccAtomic = (copies < nThreads);
    }
    gkmMsg("Kernel-mode accumulation: %d int64 triangle%s of %.2f GB (%s adds; coefficients x %g).\n", copies, (copies == 1) ? "" : "s",
           triBytes * copies / 1073741824.0, (copies < nThreads) ? "shared, atomic" : "private", coefDen);
    std::atomic<int> nextPass(0);
    gkmMsg("Running %d passes on %d thread%s.\n", iDL.M, nThreads, (nThreads==1)?"":"s");
    if (nThreads <= 1) {
      task1(L, &iDL, seqsTS, iDL.M, &nextPass, &kcs[0], b);
    } else {
      std::thread *myThreads = new std::thread[nThreads];
      for (int j = 0; j < nThreads; j++) myThreads[j] = std::thread(task1, L, &iDL, seqsTS, iDL.M, &nextPass, &kcs[j], b);
      for (int j = 0; j < nThreads; j++) myThreads[j].join();
      delete []myThreads;
    }
    for (int cp = 1; cp < copies; cp++) {
      const std::atomic<long long> *src = acc + (size_t)cp * tri;
      for (size_t x = 0; x < tri; x++) acc[x].store(acc[x].load(std::memory_order_relaxed) + src[x].load(std::memory_order_relaxed), std::memory_order_relaxed);
    }
    for (i = 0; i < nseqs; i++) norm[i] = sqrt((double)acc[(size_t)i * (i + 1) / 2 + i].load(std::memory_order_relaxed));
    for (i = 0; i < nseqs; i++)
    {
      const std::atomic<long long> *row = acc + (size_t)i * (i + 1) / 2;
      for (int j = 0; j <= i; j++)
      {
        double v;
        if (i == j) v = 1.0;
        else v = (norm[i]*norm[j]<1E-50) ? 0.0 : (double)row[j].load(std::memory_order_relaxed) / (norm[i]*norm[j]);
        if (v == 0.0) v = 0.0;
        if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
        else bin.add(gkmCanon(v));
      }
      if (fo) fprintf(fo, "\n");
    }
    delete []acc;
  }
  else
  {
    // ---- Phase 6: tiles of rows. The mismatch profile of rows [lo, hi] is built by a full set of
    // passes with the DFS pruned to that band (node id ranges), the rows are written, the tile is freed.
    // Counts are integers, so the result is identical to the untiled computation.
    int tileRows = opt.tileRows;
    if (tileRows <= 0) {
      double bytesPerRow = (double)nrowsUsed * sizeof(aint) * nseqs; // upper bound (the last row is the longest)
      double budget = (double)opt.tileMemoryMB * 1048576.0;
      tileRows = (int)(budget / (bytesPerRow > 0 ? bytesPerRow : 1));
      if (tileRows < 1) tileRows = 1;
    }
    if (tileRows > nseqs) tileRows = nseqs;
    int ntiles = (nseqs + tileRows - 1) / tileRows;
    if (ntiles > 1) { gkmMsg("Computing the kernel in %d tiles of %d rows.\n", ntiles, tileRows); }
    // One contiguous block, sized for the largest (last) tile and reused by every tile: freeing and
    // reallocating a block per tile does not lower the resident set (the allocator keeps freed large
    // blocks, and each tile's block is bigger than the previous one), a single reused block does.
    size_t maxTileCounters = 0;
    for (int tile = 0; tile < ntiles; tile++) {
      int lo = tile * tileRows, hi = lo + tileRows - 1; if (hi > nseqs - 1) hi = nseqs - 1;
      size_t c = 0; for(int i=lo;i<=hi;i++) c += nrowsUsed * (i+1);
      if (c > maxTileCounters) maxTileCounters = c;
    }
    aint *tileBlock = new aint[maxTileCounters > 0 ? maxTileCounters : 1];
    for (int tile = 0; tile < ntiles; tile++)
    {
      int lo = tile * tileRows, hi = lo + tileRows - 1; if (hi > nseqs - 1) hi = nseqs - 1;
      kc.rowLo = lo; kc.rowHi = hi;
      size_t tileCounters = 0;
      for(int i=lo;i<=hi;i++) tileCounters += nrowsUsed * (i+1);
      for(size_t k=0;k<tileCounters;k++) tileBlock[k]=0;
      {
        size_t off = 0;
        for(int i=lo;i<=hi;i++)
        {
          kc.mmProfile[i] = new aint*[kc.nclasses];
          for (int j=0;j<kc.nclasses;j++) kc.mmProfile[i][j] = NULL;
          for (size_t q=0;q<nrowsUsed;q++) { kc.mmProfile[i][reach[q]] = tileBlock + off; off += (size_t)(i+1); }
        }
      }
      int nThreads = nThreadsBudget;
      if(nThreads>iDL.M){nThreads =iDL.M;}
      if(nThreads<1){nThreads=1;}
      std::atomic<int> nextPass(0);
    
      gkmMsg("Running %d passes on %d thread%s.\n", iDL.M, nThreads, (nThreads==1)?"":"s");
      if (nThreads<=1){
        task1( L, &iDL, seqsTS, iDL.M, &nextPass, &kc, b);
      }else{
      
  #ifndef MULTI_THREAD_SAFE
        Printf("Warning -- MULTI_THREAD_SAFE is not enabled (see src/global.h). Some values may be approximated.\n");
  #endif
      
      
        std::thread *myThreads = new std::thread[nThreads];
        int j;
      
        for(j=0;j<nThreads;j++){
          myThreads[j] = std::thread(task1, L, &iDL, seqsTS, iDL.M, &nextPass, &kc, b);
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
          if (v == 0.0) v = 0.0; // canonical +0 (an exact zero can be -0.0 on one compiler and +0.0 on another)
          if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
          else bin.add(gkmCanon(v));
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
  }
  if (fo) fclose(fo); 
  else if (bin.write(opt.outfile, records) != 0) return gkmCannotOpen(opt.outfile.c_str());
  if (writeIndexSidecar(opt.outfile, records) != 0) { gkmMsg("\n WARNING: could not write %s.index\n", opt.outfile.c_str()); }
  
  delete []norm;
  return 0;
}

int gkmKernelSuffixTree(OptsGkmKernel &opt, const CConverter &conv)  //maingKernel
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
  
  CCalcWmML wmc(L, K, conv.b);
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
    c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, conv.b, opt.wildcardLambda, c);
    n0 = c[maxnmm]/2;
    
  }
  if (useTgkm==4)  //mismatch kernel 
  {
    c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, conv.b, c);
    n0 = c[maxnmm]/2;
    
  }
  
  gkmMsg("\n maximumMismatch = %d\n", maxnmm);
  for(int ii=0;ii<=maxnmm;ii++) {
    gkmMsg("\n c[%d] = %e",ii,c[ii] );
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
  
  CSequence sgii(maxseqlen+3, &conv);
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
  // Phase 7: one alphabet -> the class index is the mismatch count: rows 0..maxmm, step 1 everywhere
  std::vector<int> stepPos(L, 1);
  kc.stepPos = stepPos.data();
  kc.nclasses = maxnmm+1;
  std::vector<int> reach(kc.nclasses);
  for(int m=0;m<kc.nclasses;m++) reach[m]=m;
  // Phase 6: the mismatch profile is only ever read for j <= i (calcinnerprod is called for i >= j and
  // the leaf code only records pairs with j <= i), so row i holds i+1 entries: half the memory.
  // Rows are allocated per tile below (Phase 6: peak memory bounded by tileRows / tileMemoryMB).
  kc.mmProfile=new aint **[nseqs];
  for(int i=0;i<nseqs;i++) kc.mmProfile[i] = NULL;
  
  /// calc C 
  double C =0; 
  for(int m=0;m<=L;m++)
  {
    C+=dCombinations(L,m)*pow(1.0*conv.b-1,m)*wmc.kernelTruncated[m];
    //		C+=dCombinations(L,m)*pow(3.0,m)*wmc.kernelTruncated[m];
  }
  
  //	double btL=pow(4.0,L);
  double btL=pow(1.0*conv.b,L);
  {
    int rc = gkmKernelPassesAndOutput(opt, kc, seqsTS, L, conv.b, 1.0/conv.b, nseqs, npos, nneg, ntotal, LmersCnt, records, c, reach,
                                      usePseudocnt, n0, C, btL, maxnmm, conv.b, std::string(conv.alphabet, conv.alphabet + conv.b), btL);
    if (rc != 0) return rc;
  }
  delete []LmersCnt;
  seqsTS->deleteTree(L, conv.b, 0);
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


// Phase 7: kernel for multi-track records (T >= 2 aligned tracks, each with its own alphabet). The
// l-mers are the track-major windows (l = T*L), the mismatch profile is indexed by the block class
// (dev/PHASE7_PLAN.md D1-D7) and the coefficients come from GeneralB.
int gkmKernelMultiTrack(OptsGkmKernel &opt, const TrackAlphabets &ta)
{
  const int T = ta.T();
  const int L = opt.L, K = opt.K, ell = T * L, b = ta.bmax();
  const int nMAXSEQUENCES = opt.maxnumseq;
  const bool addRC = opt.addRC;
  AlphabetVector av(ta.alphabetVector(L));
  GeneralBTables tab(av, K);
  int maxnmm = opt.maxnmm;
  if (maxnmm == -1) maxnmm = tab.autoMaxmm(opt.useTgkm);
  if (maxnmm > ell) maxnmm = ell;
  const double *c = tab.table(opt.useTgkm);
  if (c == NULL) { Printf("\n ERROR: filter types 3 and 4 (wildcard, mismatch kernels) are only available for a single alphabet.\n"); return 1; }
  std::vector<int> reach = av.reachable(maxnmm);
  gkmMsg("\n %d tracks, window %d -> l-mers of length %d; %s\n", T, L, ell, av.describe().c_str());
  gkmMsg("\n maximumMismatch = %d (%d of %d classes reachable)\n", maxnmm, (int)reach.size(), av.nclasses);
  for (size_t q = 0; q < reach.size(); q++) gkmMsg("\n c[%s] = %e", av.classLabel(reach[q]).c_str(), c[reach[q]]);
  Printf("\n");

  CLTreeS *seqsTS = new CLTreeS();
  std::vector<int> LmersCnt;
  std::vector<SeqRecord> records;
  std::vector<int> win;
  int ntotal = 0, nseqs = 0, npos = 0;
  const char *fns[2] = { opt.posfile.c_str(), opt.negfile.c_str() };
  for (int which = 0; which < 2; which++) {
    FILE *sfi = fopen(fns[which], "r");
    if (sfi == NULL) { delete seqsTS; return gkmCannotOpen(fns[which]); }
    std::unordered_map<std::string, int> byName; // -N: same name within a file -> one row
    std::string pending;
    MultiTrackRecord rec;
    int r;
    while ((r = readMfaRecord(sfi, T, rec, pending, opt.maxseqlen)) == 1) {
      int idx = -1;
      if (opt.mergeByName) { std::unordered_map<std::string, int>::iterator it = byName.find(rec.name); if (it != byName.end()) idx = it->second; }
      if (idx < 0) {
        if (nseqs >= nMAXSEQUENCES) { fclose(sfi); delete seqsTS; return gkmTooManySequences(nMAXSEQUENCES); }
        idx = nseqs++;
        LmersCnt.push_back(0);
        records.push_back(SeqRecord{idx, "", rec.name, 0, 0, 0});
        if (opt.mergeByName) byName[rec.name] = idx;
      }
      records[idx].length += rec.length();
      win.clear();
      int nw = encodeWindows(rec, ta, L, addRC, win);
      for (int w = 0; w < nw; w++) seqsTS->addSeq(&win[(size_t)w * ell], ell, &win[(size_t)w * ell], idx);
      LmersCnt[idx] += nw; ntotal += nw;
    }
    fclose(sfi);
    if (r < 0) { delete seqsTS; return 1; }
    if (which == 0) npos = nseqs;
  }
  int nneg = nseqs - npos;
  for (int i = 0; i < nseqs; i++) { records[i].id = seqRecordId(i, npos); records[i].label = (i < npos) ? 1 : -1; records[i].nlmers = LmersCnt[i]; }
  {
    long maxw = 0; for (int i = 0; i < nseqs; i++) if (LmersCnt[i] > maxw) maxw = LmersCnt[i];
    if ((double)maxw * (double)maxw >= 2147483647.0) { gkmMsg("\n ERROR: a record has %ld windows; the profile counters hold at most 2^31 window pairs per sequence pair.\n", maxw); delete seqsTS; return 1; }
  }

  KernelContext kc;
  kc.LM1 = ell - 1;
  kc.maxmm = maxnmm;
  kc.stepPos = av.step.data();
  kc.nclasses = av.nclasses;
  kc.mmProfile = new aint **[nseqs > 0 ? nseqs : 1];
  for (int i = 0; i < nseqs; i++) kc.mmProfile[i] = NULL;
  double p = 0; for (int i = 0; i < ell; i++) p += 1.0 / av.B[i]; p /= ell; // mean match probability (pass-design heuristic only)
  double coefDen = 1.0; for (int i = 0; i < ell; i++) coefDen *= av.B[i]; // prod b_i: scales H (and the gkm counts) to integers
  int rc = gkmKernelPassesAndOutput(opt, kc, seqsTS, ell, b, p, nseqs, npos, nneg, ntotal, LmersCnt.data(), records, c, reach,
                                    false, 0.0, 0.0, 0.0, maxnmm, b, ta.canonical(), coefDen);
  seqsTS->deleteTree(ell, b, 0);
  delete seqsTS;
  delete []kc.mmProfile;
  return rc;
}

} // namespace GKM_NS
