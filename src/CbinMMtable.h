//
//  CbinMMtable.h
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#ifndef __gkmsvmXC__CbinMMtable__
#define __gkmsvmXC__CbinMMtable__

#include <stdio.h>
#include "global.h"

class CbinMMtable
{
public:
    CbinMMtable();

    int **table;
    int *dat;
    int L, Dmax;
    int nrow;
    
    void deleteTable(); // deletes the table
    int createTable(int L, int Dmax); // creates a table with all l-mer with max Dmax 1s.
    
    ~CbinMMtable(void);
};

#endif /* defined(__gkmsvmXC__CbinMMtable__) */
