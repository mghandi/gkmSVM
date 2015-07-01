//
//  CbinMMtable.cpp
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#include "CbinMMtable.h"
#include "CbinMMtree.h"



CbinMMtable::CbinMMtable(void){
    table=nullptr;
    dat=nullptr;
    L=Dmax=nrow=0;
}

int  CbinMMtable::createTable(int L, int Dmax){
    CbinMMtree *mmtree = new CbinMMtree();
    nrow = mmtree->addLDtree(L,Dmax);
    this->L = L;
    this->Dmax = Dmax;
    dat = new int[nrow*L];
    table = new int*[nrow];
    for(int i=0;i<nrow;i++){
        table[i]= dat+i*L;
    }
    int *tmpArray = new int[L];
    mmtree->addTreeToTable(table, 0, L, tmpArray);
    delete []tmpArray;
    mmtree->deleteTree();
    delete mmtree;
    return nrow;
}
    
void CbinMMtable::deleteTable(){
    if(dat!=nullptr){
        delete []dat;
        delete []table;
        table=nullptr;
        dat=nullptr;
        L=Dmax=nrow=0;
    }
}

CbinMMtable::~CbinMMtable(void){
    this->deleteTable(); 
}

/*
 class CbinMMtable tt;
 tt.createTable(4, 2);
 for(int i=0;i<tt.nrow;i++){
   printf("\n");
   for(int j=0;j<4;j++){
     printf("%d", tt.table[i][j]);
   }
 } 
 */
