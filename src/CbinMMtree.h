//
//  CbinMMtree.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#ifndef __gkmsvmXC__CbinMMtree__
#define __gkmsvmXC__CbinMMtree__

#include <stdio.h>
#include "global.h"

class CbinMMtree
{
public:
    CbinMMtree();
    CbinMMtree *child0;     //match
    CbinMMtree *child1;     //mismatch
    int addSeq(int *seq, int n); //adds a binary sequence of length n to the tree
    int deleteTree(); // deletes the tree

    int addtree(int n0, int n1); //makes a tree consisting all the possible combinations of n0 zeros and n1 ones. returns number of leaves
    int addLDtree(int L, int Dmax); // makes a tree consisting of all L-mers with at most Dmax 1s. returns number of leaves
    int addTreeToTable(int ** table, int frompos, int n, int *tmpArray); // copies all the n-mers to a table. frompos should be 0 if called from outside 
    
    /*
    int makeLDtable(int **table, int L, int Dmax); //fills in the elements in a table all the possible combinations of n0 zeros and n1 ones. returns number of rows added
    int makeTable(int **table, int n0, int n1);    //fills in the elements in a table all the possible L-mers with at most Dmax 1s. returns number of rows added
     */
    double calcAddCost(int *lmer, double *w,  int L, double p); // calculates the cost of the additional edges // p = 1/b (prob of match)

    ~CbinMMtree(void);
};




#endif /* defined(__gkmsvmXC__CbinMMtree__) */
