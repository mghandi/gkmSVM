//
//  CiDLPasses.cpp
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#include "CiDLPasses.h"
#include <algorithm>
#include "CbinMMtable.h"
#include "global.h"

CiDLPasses::CiDLPasses(void){
  L=0; // Length
  M=0; // number of passes
  identityOrder=0;
  //Dmax=0; // maximum number of mismatches
  passOrder=NULL; // order of the bases to be traversed. each of passOrder[0..M-1] is a permutation of 0..L-1
  passTrees=NULL; // for each of the M passes, there is a tree that holds all the mismatch (binary) L-mers that
}

void CiDLPasses::initPassOrderIDL(int L){
  deletePassOrder();
  this->L = L;
  M = L;
  this->passOrder = new int *[M];
  for(int i=0;i<M;i++){
    passOrder[i]= new int[L];
    for(int j=0;j<L;j++){
      passOrder[i][j]=(i+j)%L;
    }
  }
}

void CiDLPasses::initPassOrderAll(int L, int Dmax){ // generates choose(L, Dmax) permutations
  deletePassOrder();
  this->L = L;
  this->M=Combinations(L, Dmax);
  
  class CbinMMtable tt;
  tt.createTable(L, Dmax);
  
  this->passOrder = new int *[M];
  
  for(int i=0;i<M;i++){
    passOrder[i]= new int[L];
    for(int j=0;j<L;j++){
      passOrder[i][j]=(i+j)%L;
    }
  }
  
  int n=0;
  for(int i=0;i<tt.nrow;i++){
    Printf("\n");
    int n1 =0;
    for(int j=0;j<L;j++){
      n1+=tt.table[i][j];
    }
    if(n1==Dmax){
      int n0i=0;
      int n1i=0;
      for(int j=0;j<L;j++){
        int ii;
        if(tt.table[i][j]==1){
          n1i++;
          ii = L-n1i;
        }else{
          ii = n0i;
          n0i++;
        }
        passOrder[n][ii]=j;
      }
      n++;
    }
  }
  if(n!=M){
    Printf("\nError number mismatch!\n");
    return;//exit(10);
  }
  
  for(int i=0;i< this->M;i++){
    Printf("\n");
    for(int j=0;j<this->L;j++){
      //            printf("%d ", this->passOrder[i][j]);
      gkmMsg("%d ", this->passOrder[i][j]);
      
    }
  }
}
void CiDLPasses::newPassOrderDesignCover(int L, int Dmax, int k){// generates M passes, that gaurantee that the first k places are matches (in all the trees)
  
  // source: http://www.ccrwest.org/cover.html
  // A (v,k,t)-covering design is a collection of k-element subsets, called blocks, of {1,2,...,v}, such that any t-element subset is contained in at least one block.  This site contains a collection of good (v,k,t)-coverings. Each of these coverings gives an upper bound for the corresponding C(v,k,t), the smallest possible number of blocks in such a covering design.
  
  // we need (v=L, t=Dm, k=L-k)-covering design
  // eg. L=10, Dm=3, L-k=10-3=7
  // https://www.ccrwest.org/cover/t_pages/t3/k7/C_10_7_3.html
  /*
  1  2  3  4  5  6  7
  1  2  3  4  8  9 10
  1  5  6  7  8  9 10
  2  3  4  5  6  8  9
  2  3  4  5  6  7 10
  1  2  3  4  7  8  9
  
  8 9 10   1  2  3  4  5  6  7
  5 6 7    1  2  3  4  8  9 10
  2 3 4    1  5  6  7  8  9 10
  1 7 10   2  3  4  5  6  8  9
  1 8 9    2  3  4  5  6  7 10
  5 6 10   1  2  3  4  7  8  9
  
  */
  
  this->L = L;
  
  if((L==10)&&(Dmax==3)&&(k==3)){
    this->M = 6;
    this->passOrder = new int *[M];
    for(int i=0;i<M;i++){
      passOrder[i]= new int[L];
      for(int j=0;j<L;j++){
        passOrder[i][j]=(i+j)%L;
      }
    }
    int DC[6][10]=  {{8,9,10,1,2,3,4,5,6,7},
    {5,6,7,1,2,3,4,8,9,10},
    {2,3,4,1,5,6,7,8,9,10},
    {1,7,10,2,3,4,5,6,8,9},
    {1,8,9,2,3,4,5,6,7,10},
    {5,6,10,1,2,3,4,7,8,9}};
    
    
    for(int i=0;i<M;i++){
      for(int j=0;j<L;j++){
        passOrder[i][j]=DC[i][j]-1;
      }
    }
  }
  
  
  if((L==10)&&(Dmax==3)&&(k==4)){
    this->M = 10;
    this->passOrder = new int *[M];
    for(int i=0;i<M;i++){
      passOrder[i]= new int[L];
      for(int j=0;j<L;j++){
        passOrder[i][j]=(i+j)%L;
      }
    }
    int DC[10][10]={
      {6, 8 , 9 , 10 , 1,2,3,4,5,7},
      {1, 7, 9, 10,    2,3,4,5,6,8},
      {10, 1, 2, 8,    3,4,5,6,7,9},
      {3, 9, 1, 2,     4,5,6,7,8,10},
      {3, 4, 10,2,     9, 1,5,6,7,8},
      {1, 3, 4, 5,     9,10,2,6,7,8},
      {4, 5, 6, 2,    8,9,10,1,3,7},
      {6, 7,3,5,     4,8,9,10,1,2},
      {8,4, 6, 7,     2,3,5,9,10,1},
      {5, 7, 8,9,      1,2,3,4,6,10}};
    
    
    for(int i=0;i<M;i++){
      for(int j=0;j<L;j++){
        passOrder[i][j]=DC[i][j]-1;
      }
    }
  }
  
  
}
int CiDLPasses::isCoprime(int a, int b){
  int r=1;
  for(int i=2;i<a;i++){
    if((a%i==0)&&((i+b)%i==0)){
      r = 0;
    }
  }
  return(r);
}

void CiDLPasses::initPassOrderIDL(int L, int M, int Dmax){
  deletePassOrder();
  this->L = L;
  //M = L;
  this->M = M;
  
  if(L==M){
    this->passOrder = new int *[M];
    for(int i=0;i<M;i++){
      passOrder[i]= new int[L];
      for(int j=0;j<L;j++){
        passOrder[i][j]=(i+j)%L;
      }
    }
  }
  if(2*L==M){
    this->passOrder = new int *[M];
    for(int i=0;i<L;i++){
      passOrder[i]= new int[L];
      for(int j=0;j<L;j++){
        passOrder[i][j]=(i+j)%L;
      }
    }
    /* reverse orders
    for(int i=0;i<L;i++){
    passOrder[i+L]= new int[L];
    for(int j=0;j<L;j++){
    passOrder[i+L][j]=(L+i-j)%L;
    }
    }
    */
    
    /* noori method: 
    Let r=integer part of $\ell / (\ell-m)$;
    If \ell and r are not coprime choose a minumum number \ell'>\ell such that \ell' and r are coprime
    Find a permutation pi on $\{1,\cdots,\ell'\}$ with  pi(a)-pi(b)=r(a-b) (mod \ell')
    Modify pi to a permutation pi' over \{1,\cdots,\ell\} (This should be explained more precisely)
    */
    //int L2 = L;
    //int r = 2;
    //::gMAXMM
    //while(!this->isCoprime(L2,r)){L2++;}
    
    //1,5,8,2,6,9,3,7,10,4
    passOrder[L]= new int[L];
    /*   passOrder[L][0]=1-1;
    passOrder[L][1]=5-1;
    passOrder[L][2]=8-1;
    passOrder[L][3]=2-1;
    passOrder[L][4]=6-1;
    passOrder[L][5]=9-1;
    passOrder[L][6]=3-1;
    passOrder[L][7]=7-1;
    passOrder[L][8]=10-1;
    passOrder[L][9]=4-1;
    
    //1,4,7,9,2,5,8,10,3,6
    
    passOrder[L][0]=1-1;
    passOrder[L][1]=4-1;
    passOrder[L][2]=7-1;
    passOrder[L][3]=9-1;
    passOrder[L][4]=2-1;
    passOrder[L][5]=5-1;
    passOrder[L][6]=8-1;
    passOrder[L][7]=10-1;
    passOrder[L][8]=3-1;
    passOrder[L][9]=6-1;
    */
    //for(int j=0;j<L;j++){
    //    passOrder[L][j]=((j*L)/3)%L;
    //}
    //        int q=3; int r=1;
    //        int q=2; int r=2;
    int q = (L)/Dmax;
    int r = L-q*Dmax;
    //        q=4;r=2;
    //        passOrder[L][0]=0;
    passOrder[L][L-1]=0;
    int pii=0;
    for(int j=1;j<L;j++){
      if((pii<(r*q+r)%L)||(pii>=(L-q))){
        pii+=q+1;
      }else{
        pii+=q;
      }
      pii%=L;
      //          passOrder[L][j]=pii;
      passOrder[L][L-j-1]=pii;
    }
    
    
    for(int i=1;i<L;i++){
      passOrder[i+L]= new int[L];
      for(int j=0;j<L;j++){
        passOrder[i+L][j]=passOrder[L][(j+i)%L];
      }
    }
    /*  print passes:
    Printf("\n");
    for(int i=0;i< this->M;i++){
    for(int j=0;j<this->L;j++){
    gkmMsg("%d ", this->passOrder[i][j]);
    }
    Printf("\n");
    
    }
    */
    }
  
  // notye : it does not work for other Ms . use random pass for that
  
  }

double CiDLPasses::calcSlope(int *lmer, int *order, int L){
  // calculates slope (max #1s / total)
  double res = 0;
  double n1=0;
  for(int i=0;i<L;i++){
    n1 = n1 + lmer[order[i]];
    if(n1/(i+1)>res){
      res = n1/(i+1);
    }
  }
  return(res);
}


int *CiDLPasses::reorder(int *lmer, int *order, int L, int *output){
  for(int i=0;i<L;i++){
    output[i]=lmer[order[i]];
  }
  return(output);
}

void CiDLPasses::newIDLPasses(int L, int Dmax){
  this->L = L;
  
  if(this->passOrder==NULL){
    initPassOrderIDL(L);
  }
  
  deletePassTrees();
  passTrees = new CbinMMtree *[M];
  for(int i=0;i<M;i++){
    passTrees[i] = new CbinMMtree();
  }
  
  CbinMMtable mmtable;
  mmtable.createTable(L, Dmax);
  
  int *jlmer = new int[L];
  for(int i=0;i<mmtable.nrow; i++){
    int *lmer = mmtable.table[i];
    double sMin = 1.5;
    int jMin = 0;
    for(int j=0;j<M; j++){
      double sj = calcSlope(lmer, passOrder[j], L);
      if(sj<sMin){
        sMin = sj;
        jMin = j;
      }
    }
    jlmer = reorder(lmer, passOrder[jMin], L, jlmer);
    passTrees[jMin]->addSeq(jlmer, L);
    /*
    printf(" T%d, ", jMin); 
    for(int kk=0;kk<L;kk++){
    printf("%d", lmer[kk]);
    }
    printf(" --> ");
    for(int kk=0;kk<L;kk++){
    printf("%d", jlmer[kk]);
    }
    printf("\n");
    */
  }
  
  delete []jlmer;
  mmtable.deleteTable();
}

void CiDLPasses::newGreedyIDLPasses(int L, int M,  int Dmax, int *nodesAtDepthCnt, double p){
  
  this->L = L;
  
  if(this->passOrder==NULL){
    initPassOrderIDL(L, M, Dmax);
  }
  
  double *w = new double[L];
  for(int i=0;i<L;i++){
    w[i]=(1.0*nodesAtDepthCnt[i])*nodesAtDepthCnt[i];
    // w[i]=1.0;
  }
  
  deletePassTrees();
  passTrees = new CbinMMtree *[M];
  for(int i=0;i<M;i++){
    passTrees[i] = new CbinMMtree();
  }
  
  CbinMMtable mmtable;
  mmtable.createTable(L, Dmax);
  
  // A2a (dev/PERFORMANCE_PLAN.md section 5, 2026-08-31): the original rule put every pattern into the pass where
  // its own estimated DFS cost is minimal -- the objective was total work, with no term for the heaviest pass,
  // which is what the threads wait for. Now: patterns in decreasing order of their minimum cost (LPT), each
  // assigned to the least-loaded pass among those whose cost is within (1 + GKM_GREEDY_EPS) of its minimum.
  // Any assignment gives the identical kernel; only the balance of the passes changes. Measured 2026-08-31 on
  // 28 threads: -10 % wall for l = 20 at 13 398 sequences (d = 4), -4 % for l = 24 (d = 6), but +10 % for DNA
  // l = 10 (d = 3) -- so the load-aware rule is used for words of >= 16 positions only (eps = 0 reproduces the
  // min-cost rule up to ties). The cost model itself (A2b) is the larger lever.
  const double GKM_GREEDY_EPS = (L >= 16) ? 0.2 : 0.0;
  int n = mmtable.nrow;
  double *minCost = new double[n];
  int *order = new int[n];
  for(int i=0;i<n;i++){
    double sMin = 1.0E300;
    for(int j=0;j<M; j++){ double sj = calcCost(mmtable.table[i], passOrder[j], w, p, L); if(sj<sMin) sMin = sj; }
    minCost[i] = sMin; order[i] = i;
  }
  std::sort(order, order+n, [&](int a, int b){ return minCost[a] > minCost[b] || (minCost[a] == minCost[b] && a < b); });
  double *load = new double[M];
  for(int j=0;j<M;j++) load[j]=0;
  int *jlmer = new int[L];
  for(int ri=0;ri<n; ri++){
    int i = order[ri];
    int *lmer = mmtable.table[i];
    double bound = minCost[i]*(1.0+GKM_GREEDY_EPS);
    int jMin = -1; double bestCost = 0;
    for(int j=0;j<M; j++){
      double sj = calcCost(lmer, passOrder[j], w, p, L);
      if(sj<=bound && (jMin<0 || load[j]<load[jMin])){ jMin = j; bestCost = sj; }
    }
    load[jMin] += bestCost;
    jlmer = reorder(lmer, passOrder[jMin], L, jlmer);
    passTrees[jMin]->addSeq(jlmer, L);
  }
  delete []minCost; delete []order; delete []load;
  
  delete []w;
  delete []jlmer;
  mmtable.deleteTable();
}


double CiDLPasses::calcCost(int *lmer, int *order, double *w, double p, int L){
  double res = 0;
  double pj=1.0;
  for(int i=0;i<L;i++){
    if(lmer[order[i]]==0){
      pj *= p;
    }else{
      pj *= (1-p);
    }
    res += w[i]*pj;
  }
  return res;
}

void CiDLPasses::newGreedy2IDLPasses(int L, int M,  int Dmax, int *nodesAtDepthCnt, double p){
  
  deletePassOrder();
  this->L = L;
  this->M = M;
  
  this->passOrder = new int *[M];
  for(int i=0;i<M;i++){
    passOrder[i]= new int[L];
    for(int j=0;j<L;j++){
      passOrder[i][j]=j;
    }
  }
  
  double *w = new double[L];
  for(int i=0;i<L;i++){
    w[i]=(1.0*nodesAtDepthCnt[i])*nodesAtDepthCnt[i];
  }
  
  deletePassTrees();
  passTrees = new CbinMMtree *[M];
  for(int i=0;i<M;i++){
    passTrees[i] = new CbinMMtree();
  }
  
  CbinMMtable mmtable;
  mmtable.createTable(L, Dmax);
  
  
  //int *rndi = new int[mmtable.nrow];
  //for(int i=0;i<mmtable.nrow;i++){
  //    rndi[i]=i;
  //}
  //randomPermute(rndi, mmtable.nrow);
  int n = mmtable.nrow;
  int ** table = mmtable.table;
  double totalCost = 0;
  double *minCost = new double[n];
  int *minTree = new int[n];
  for(int i=0;i<n;i++){
    minTree[i]=0;
    minCost[i]=calcCost(table[i], passOrder[0], w, p, L);
    totalCost +=minCost[i];
  }
  
  double *sumcost=new double[L];
  
  for(int m=1;m<M; m++){
    gkmMsg("  %d total cost = %lf\n", m, totalCost);
    // calc avgCost
    for(int j=0;j<L;j++){sumcost[j]=0;}
    for(int i=0;i<n;i++){
      for(int j=0;j<L;j++){
        sumcost[j]+=table[i][j]*minCost[i];
      }
    }
    // sort the sumcost to get the new pass
    int *passi = passOrder[m];
    for(int i=0;i<L;i++){
      for(int j=0;j<i;j++){
        if(sumcost[passi[i]]<sumcost[passi[j]]){
          int h= passi[i];
          passi[i]=passi[j];
          passi[j]=h;
        }
      }
    }
    
    for(int i=0;i<L;i++){
      gkmMsg(" %d ", passi[i]);
    }
    Printf("\n");
    //for(int i=0;i<L;i++){
    //    printf(" %lf ", sumcost[passi[i]]);
    //}
    //printf("\n");
    
    // update minCosts
    totalCost=0;
    for(int i=0;i<n;i++){
      double costm=calcCost(table[i], passOrder[m], w, p, L);
      if(costm<minCost[i]){
        minTree[i]=m;
        minCost[i]=costm;
      }
      totalCost +=minCost[i];
    }
  }
  
  // add sequences to the trees
  int *jlmer = new int[L];
  
  for(int i=0;i<n;i++){
    jlmer = reorder(table[i], passOrder[minTree[i]], L, jlmer);
    
    passTrees[minTree[i]]->addSeq(jlmer, L);
    
  }
  
  delete []w;
  delete []jlmer;
  delete []sumcost;
  delete []minCost;
  delete []minTree;
  mmtable.deleteTable();
}

void CiDLPasses::deletePassOrder(){
  if(passOrder!=NULL){
    
    for(int i=0;i<M;i++){
      delete []passOrder[i];
    }
    delete[]passOrder;
    passOrder=NULL;
  }
}


void CiDLPasses::deletePassTrees(){
  if(passTrees!=NULL){
    if(arena.empty()){
      for(int i=0;i<M;i++){
        if(passTrees[i]!=NULL){ passTrees[i]->deleteTree(); delete passTrees[i]; }
      }
    }
    delete []passTrees;
    passTrees=NULL;
  }
  arena.clear();
  identityOrder=0;
}

double CiDLPasses::patternCount(int L, int Dmax){
  double n = 0;
  for(int i=0;i<=Dmax && i<=L;i++) n += dCombinations(L, i);
  return n;
}

void CiDLPasses::newPrefixSplitPasses(int L, int Dmax, int q){
  if(q>L) q=L;
  if(q<0) q=0;
  deletePassTrees();
  deletePassOrder();
  this->L = L;
  // the passes: every q-bit prefix with at most Dmax ones
  std::vector<int> prefixes;
  for(int pat=0; pat<(1<<q); pat++){
    int ones=0; for(int i=0;i<q;i++) ones += (pat>>i)&1;
    if(ones<=Dmax) prefixes.push_back(pat);
  }
  M = (int)prefixes.size();
  passOrder = new int *[M];
  for(int j=0;j<M;j++){ passOrder[j] = new int[L]; for(int d=0;d<L;d++) passOrder[j][d]=d; }
  identityOrder = 1;
  // arena: chain nodes (M*q) + DAG nodes for depths q..L and 0..Dmax mismatches used
  int nd = L-q+1;
  arena.assign((size_t)M*q + (size_t)nd*(Dmax+1), CbinMMtree());
  CbinMMtree *dag = arena.data() + (size_t)M*q; // dag[(d-q)*(Dmax+1)+u]
  for(int d=q; d<=L; d++){
    for(int u=0;u<=Dmax;u++){
      CbinMMtree &node = dag[(size_t)(d-q)*(Dmax+1)+u];
      if(d==L){ node.child0=node.child1=NULL; continue; } // leaf level: never dereferenced by the DFS
      node.child0 = &dag[(size_t)(d+1-q)*(Dmax+1)+u];
      node.child1 = (u<Dmax) ? &dag[(size_t)(d+1-q)*(Dmax+1)+u+1] : NULL;
    }
  }
  passTrees = new CbinMMtree *[M];
  for(int j=0;j<M;j++){
    int pat = prefixes[j], ones=0;
    CbinMMtree *chain = arena.data() + (size_t)j*q;
    for(int i=0;i<q;i++){
      int bit = (pat>>i)&1; ones += bit;
      CbinMMtree *next = (i+1<q) ? &chain[i+1] : &dag[(size_t)0*(Dmax+1)+ones];
      chain[i].child0 = bit ? NULL : next;
      chain[i].child1 = bit ? next : NULL;
    }
    passTrees[j] = (q>0) ? &chain[0] : &dag[0];
  }
}

CiDLPasses::~CiDLPasses(void){
  deletePassTrees();
  deletePassOrder();
}
