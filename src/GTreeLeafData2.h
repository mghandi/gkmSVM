//
//  GTreeLeafData2.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 6/10/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#ifndef __gkmsvmXC__GTreeLeafData2__
#define __gkmsvmXC__GTreeLeafData2__

#include <stdio.h>
#include "global.h"

class GTreeLeafData2  // handles addLTreeSnodeData instead of individual seqIDs , hence hopefully faster and less memory req.
{
public:
    
    int n; // number of L-mers in the list
    LTreeSnodeDataptr seqIDsets; //if n==1, it is LTreeSnodeData* and contains one seqIDset, otherwise it is LTreeSnodeData** and is the array of seqIDsets;  // each seqIDset is the list of seqIDs containing one l-mer
    intintptr gbits; // gbits for the case n==1, or array of gbits for n>1
    
    GTreeLeafData2(void);
    ~GTreeLeafData2(void);
    
    //void add(int seqID, int gbits);
    void addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq);
    
    
    void process();// calculates the mismatch profiles
    int calcdist(int difx); // calculates number of mismatches
private:
    
};

extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]


#endif /* defined(__gkmsvmXC__GTreeLeafData2__) */
