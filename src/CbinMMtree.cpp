//
//  CbinMMtree.cpp
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 5/21/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

// binary mismatch tree
#include "CbinMMtree.h"



CbinMMtree::CbinMMtree(void){
    child0 = child1 = nullptr;
}

CbinMMtree::~CbinMMtree(void){
}

int CbinMMtree::addSeq(int *seq, int n){
    if(n<=0){
        return(0);
    }
    if(*seq){
        if(child1==nullptr){
            child1 = new CbinMMtree();
        }
        child1->addSeq(seq+1,n-1);
    }else{
        if(child0==nullptr){
            child0= new CbinMMtree();
        }
        child0->addSeq(seq+1,n-1);
    }
    return(0);
}

int CbinMMtree::addLDtree(int L, int Dmax){
    int res = 0;
    for(int i=0; i<=Dmax; i++){
        res+=this->addtree(L-i, i);
    }
    return(res);
}

int CbinMMtree::addtree(int n0, int n1){
    if((n0==0)&&(n1==0)){
        return(1);
    }
    int res = 0;
    if(n0>0){
        if(child0==nullptr){
            child0 = new CbinMMtree();
        }
        res+=child0->addtree(n0-1, n1);
    }
    if(n1>0){
        if(child1==nullptr){
            child1 = new CbinMMtree();
        }
        res+=child1->addtree(n0, n1-1);
    }

    return(res);
}


int CbinMMtree::addTreeToTable(int ** table, int frompos, int n, int *tmpArray){
 
    if(frompos==n){
        for(int i=0;i<n;i++){
            table[0][i]=tmpArray[i];
          //  printf("%d",tmpArray[i]);
        }
        //printf("\n");
        return(1);
    }
    int nadded = 0;
    int res = 0;
    if(child0!=nullptr){
        tmpArray[frompos]=0;
        nadded = child0->addTreeToTable(table, frompos+1, n, tmpArray);
        res+=nadded;
    }
    if(child1!=nullptr){
        tmpArray[frompos]=1;
        nadded = child1->addTreeToTable(table+res, frompos+1, n, tmpArray);
        res+=nadded;
    }
    return(res);
}
/*

int CbinMMtree::makeLDtable(int **table, int L, int Dmax){
    int res = 0;
    for(int i=0; i<=Dmax; i++){
     res+=this->makeTable(table+res, L-i, i);
    }
    return(res);
 }
 
int CbinMMtree::makeTable(int **table, int n0, int n1){
    if((n0==0)&&(n1==0)){
        return(1);
    }
    int nadded = 0;
    int res = 0;
    if(n0>0){
        nadded=makeTable(table+res, n0-1, n1);
        for(int i=0;i<nadded;i++){
            table[i][n0+n1-1]=0;
        }
        res+=nadded;
    }
    if(n1>0){
        nadded=makeTable(table+res, n0, n1-1);
        for(int i=0;i<nadded;i++){
            table[i][n0+n1-1]=1;
        }
        res+=nadded;
    }
    return(res);
 }

 */


double CbinMMtree::calcAddCost(int *lmer, double *w, int L, double p){ // calculates the cost of the additional edges
    // w is the weights . nodesatDepth^2 for example
    CbinMMtree *node = this;
    CbinMMtree *son = nullptr;
    double pj = 1;
    
    for(int i=0;i<L;i++){
        if (lmer[i]==0){
            son = node->child0;
            pj = pj *p;
        }else{
            son = node->child1;
            pj = pj *(1-p);
        }
        if (son==nullptr){
            double res = w[i]*pj;
            for(int j=i+1; j<L;j++){
                if (lmer[j]==0){
                    pj = pj *p;
                }else{
                    pj = pj *(1-p);
                }
                res += w[j]*pj;
            }
            return res;
        }
    }
    return 0;
}

int CbinMMtree::deleteTree(){
    if(child0!=nullptr){
        child0->deleteTree();
        delete child0;
    }
    if(child1!=nullptr){
        child1->deleteTree();
        delete child1;
    }
    return(0);
}

