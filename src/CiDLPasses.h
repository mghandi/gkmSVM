//
//  CiDLPasses.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#ifndef __gkmsvmXC__CiDLPasses__
#define __gkmsvmXC__CiDLPasses__

#include <stdio.h>
#include <vector>
#include "CbinMMtree.h"

class CiDLPasses
{
public:
    CiDLPasses();
    int L; // Length
    int M; // number of passes
    //int Dmax; // maximum number of mismatches
    int **passOrder; // order of the bases to be traversed. each of passOrder[0..M-1] is a permutation of 0..L-1
    CbinMMtree **passTrees; // for each of the M passes, there is a tree that holds all the mismatch (binary) L-mers that are checked in that pass. each L-mer should belong to one pass (tree) . the ordering is based on pass[i][] order
    
    void initPassOrderIDL(int L); // generates M=L circular pass orders
    void initPassOrderIDL(int L, int M, int Dmax); // generates M=L or 2L circular pass orders
    void initPassOrderAll(int L, int Dmax); // generates choose(L, Dmax) permutations
    void newPassOrderDesignCover(int L, int Dmax, int k);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
        
    void newIDLPasses(int L, int Dmax); // generates M circular passes and assigns L-mers to the trees based on min iDL rule.
    void newGreedyIDLPasses(int L, int M, int Dmax, int *nodesAtDepthCnt, double p); // uses the greedty alg to assign to trees. p =1/b prob(match)
    // A2b (2026-08-31): the same greedy assignment with a cost model measured on the actual l-mers: per-position
    // match probabilities p_j = sum_s f_j(s)^2 and, for every pass order, the number of distinct prefixes at each
    // depth of the trie reordered by that pass (exact, by sorting packed codes; l * bits <= 64, lexicographic
    // otherwise). lmers: U distinct l-mers of L symbols each (a stride sample of them for very large U).
    void newGreedyIDLPassesB(int L, int M, int Dmax, const std::vector<int> &lmers, int b);
    double calcCostB(int *lmer, int *order, const double *w, const double *pPos, int L);
    void newGreedy2IDLPasses(int L, int M, int Dmax, int *nodesAtDepthCnt, double p); // uses the greedty 2alg to assign to trees. p =1/b prob(match)
    // Phase 7 (long words): one pass per match/mismatch pattern of the first q positions (with <= Dmax
    // mismatches), the remaining positions bounded by an implicit (depth, mismatches used) DAG that plugs
    // into the DFS through the same child0/child1 pointers. Every pattern with <= Dmax mismatches is
    // enumerated exactly once, no pattern table is materialised, and all passes use the identity order
    // (identityOrder = 1: the driver traverses the trie itself instead of a reordered clone).
    void newPrefixSplitPasses(int L, int Dmax, int q);
    int identityOrder;      // 1 when every passOrder is the identity (prefix-split passes)
    static double patternCount(int L, int Dmax); // number of binary L-mers with <= Dmax ones
private:
    std::vector<CbinMMtree> arena; // nodes of the prefix-split passes (shared DAG; not deleted through deleteTree)
public:

    void deletePassOrder();
    double calcSlope(int *lmer, int *order, int L); // calculates slope (max #1s / total)
    double calcCost(int *lmer, int *order, double *w, double p, int L); // calcs a cost for traversing a lmer in specific order. for Greedy2alg
    
    int *reorder(int *lmer, int *order, int L, int *output); // reorders the L-mer
    
    int isCoprime(int a, int b); // checks if a and b are co-prime
    
    ~CiDLPasses(void);
    void deletePassTrees();
};


#endif /* defined(__gkmsvmXC__CiDLPasses__) */
