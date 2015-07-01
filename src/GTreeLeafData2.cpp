//
//  GTreeLeafData2.cpp
//  gkmsvmXC
//
//  Created by Mahmoud Ghandi on 6/10/15.
//  Copyright (c) 2015 Broad. All rights reserved.
//

#include "GTreeLeafData2.h"

GTreeLeafData2::GTreeLeafData2(void){
    this->n=0;
}

GTreeLeafData2::~GTreeLeafData2(void){
}

/*void GTreeLeafData2::add(int seqID, int gbits){
    
    if(n==0){
        n = 1;
        this->seqIDs_gbits.i = seqID;
        this->first_gbits = gbits;
    }else{
        
        if (n==1){
            intintptr newseqids;
            newseqids.p= new int[2*2];
            newseqids.p[0]= seqIDs_gbits.i;
            newseqids.p[1]= first_gbits;
            
            newseqids.p[2]= seqID;
            newseqids.p[3]= gbits;
            
            this->seqIDs_gbits = newseqids;
            n = 2;
        }else{
            // n>1
            
            
            if ((n & (n-1))==0) // i.e. n is power of 2
            {
                // expand memory
                intintptr newseqids;
                newseqids.p= new int[n<<2];
                for(int j=0;j< (2*n); j++)
                {
                    newseqids.p[j]= seqIDs_gbits.p[j];
                }
                delete []seqIDs_gbits.p;
                seqIDs_gbits.p = newseqids.p;
            }
            seqIDs_gbits.p[2*n] = seqID;
            seqIDs_gbits.p[2*n+1] = gbits;
            n++;
        }
    }
    
}
*/

void GTreeLeafData2::addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq){
    // adds many sequences all at once
     if(n==0){
         n = 1;
         this->seqIDsets.p= nodeData;
         this->gbits.i=curGapPosSeq;
     }else{
         if (n==1){
             
             LTreeSnodeDataptr newseqIDsets; //if n==1, it is LTreeSnodeData* and contains one seqIDset, otherwise it is LTreeSnodeData** and is the array of seqIDsets;  // each seqIDset is the list of seqIDs containing one l-mer
             intintptr newgbits; // gbits for the case n==1, or array of gbits for n>1
             
             newgbits.p=new int[2];
             newgbits.p[0]=gbits.i;
             newgbits.p[1]=curGapPosSeq;
             gbits.p=newgbits.p;
             
             newseqIDsets.pp=new LTreeSnodeData*[2];
             newseqIDsets.pp[0]=seqIDsets.p;
             newseqIDsets.pp[1]=nodeData;
             seqIDsets.pp=newseqIDsets.pp;
             
             n = 2;
         }else{
         // n>1
         
             
             if ((n & (n-1))==0) // i.e. n is power of 2
             {
                 // expand memory
                 
                 LTreeSnodeDataptr newseqIDsets; //if n==1, it is LTreeSnodeData* and contains one seqIDset, otherwise it is LTreeSnodeData** and is the array of seqIDsets;  // each seqIDset is the list of seqIDs containing one l-mer
                 intintptr newgbits; // gbits for the case n==1, or array of gbits for n>1
                 
                 newgbits.p=new int[2*n];
                 newseqIDsets.pp=new LTreeSnodeData*[2*n];
                 for(int j=0;j< n; j++)
                 {
                     newgbits.p[j]=gbits.p[j];
                     newseqIDsets.pp[j]=seqIDsets.pp[j];
                 }
                 delete []gbits.p;
                 delete []seqIDsets.pp;
                 gbits.p=newgbits.p;
                 seqIDsets.pp=newseqIDsets.pp;
                 
             }
             gbits.p[n]=curGapPosSeq;
             seqIDsets.pp[n]=nodeData;

             n++;
        }
    }
    
}

void GTreeLeafData2::process(){// calculates the mismatch profiles
    if(n==1){
        LTreeSnodeData *nodeData= this->seqIDsets.p;
        if(nodeData->n==1){
          //  gMMProfile[nodeData->seqIDs.i][0][nodeData->seqIDs.i]++; // it is not needed as it can be computed based on the sequence length !


        }else{
            for(int i=1;i<nodeData->n;i++){
                for(int j=0;j<i;j++){
                    gMMProfile[nodeData->seqIDs.p[i]][0][nodeData->seqIDs.p[j]]++;
                }
            }
       }
    }else{

        for(int ni=0;ni<n;ni++){

            LTreeSnodeData *nodeDatai= this->seqIDsets.pp[ni];
            // first within ni
            if(nodeDatai->n==1){
                //  gMMProfile[nodeData->seqIDs.i][0][nodeData->seqIDs.i]++; // it is not needed as it can be computed based on the sequence length !
            }else{
                for(int i=1;i<nodeDatai->n;i++){
                    for(int j=0;j<i;j++){
                        gMMProfile[nodeDatai->seqIDs.p[i]][0][nodeDatai->seqIDs.p[j]]++;
                    }
                }
            }
            
            
            //now ni x nj
            int gbitsni = this->gbits.p[ni];
            
            if(nodeDatai->n==1){
                
                for(int nj=0;nj<ni;nj++){
                    LTreeSnodeData *nodeDataj= this->seqIDsets.pp[nj];
                    int gbitsnj = this->gbits.p[nj];
                    int difij =calcdist(gbitsni^gbitsnj);
                    
                    if(nodeDataj->n==1){
                        gMMProfile[nodeDatai->seqIDs.i][difij][nodeDataj->seqIDs.i]++;
                        
                    }else{
                        for(int j=0;j<nodeDataj->n;j++){
                            gMMProfile[nodeDatai->seqIDs.i][difij][nodeDataj->seqIDs.p[j]]++;
                        }
                    }
                }
                
            }else{
                
                
                for(int nj=0;nj<ni;nj++){
                    LTreeSnodeData *nodeDataj= this->seqIDsets.pp[nj];
                    int gbitsnj = this->gbits.p[nj];
                    int difij =calcdist(gbitsni^gbitsnj);
                    
                    if(nodeDataj->n==1){
                        
                        for(int i=0;i<nodeDatai->n;i++){
                            gMMProfile[nodeDatai->seqIDs.p[i]][difij][nodeDataj->seqIDs.i]++;
                        }
                        
                    }else{
                        for(int i=0;i<nodeDatai->n;i++){

                            for(int j=0;j<nodeDataj->n;j++){
                                gMMProfile[nodeDatai->seqIDs.p[i]][difij][nodeDataj->seqIDs.p[j]]++;
                            }
                        }
                    }
                }
                
           }
        }
    }
}

int GTreeLeafData2::calcdist(int difx){
    int res = 0;
    while(difx>0){
        res+=((difx&NBITSONES)!=0);
        difx>>=NBITS;
    }
    return(res);
}