//
//  GTreeLeafData.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 6/10/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#ifndef __gkmsvmXC__GTreeLeafData__
#define __gkmsvmXC__GTreeLeafData__

#include <stdio.h>
#include "global.h"

class GTreeLeafData
{
public:

    int n; // number of L-mers in the list
    intintptr seqIDs_gbits; //if n==1, it is int and contains the ID, otherwise it is int* and is the array of IDs;
    //  LPTr Lmers; //pointer to the starts of the sequences
    int first_gbits; // gbits for the case n==1, otherwise, seqIDs and gapped_bits are both written in seqIDs_gbits (2 numbers for each L-mer)

    GTreeLeafData(void);
    ~GTreeLeafData(void);
    
    void add(int seqID, int gbits);
    void addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq);
    
    
    void process();// calculates the mismatch profiles
    int calcdist(int difx); // calculates number of mismatches
private:
    
};

extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]


#endif /* defined(__gkmsvmXC__GTreeLeafData__) */
