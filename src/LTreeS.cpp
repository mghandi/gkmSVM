/* LTreeS.cpp
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "LTreeS.h"
#include "Sequence.h"
#include "globalvar.h"

CLTreeS::CLTreeS(void)
{
	//daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
	for(int i=0;i<MAX_ALPHABET_SIZE;i++){
		daughter[i].t=NULL;
	}
	maxSeqID=0;
	minSeqID=0;

	nonEmptyDaughterCnt=0;
	#ifdef FAST_TRACK
	FT_seq=NULL; FT_cnt=0;
	#endif
}

CLTreeS::~CLTreeS(void)
{
}

void CLTreeS::initTree()
{
	//daughter[0].t=daughter[1].t=daughter[2].t=daughter[3].t=NULL;
	for(int i=0;i<MAX_ALPHABET_SIZE;i++){
		daughter[i].t=NULL;
	}
	maxSeqID=0;
	minSeqID=0;

	nonEmptyDaughterCnt=0;
	#ifdef FAST_TRACK
	FT_seq=NULL; FT_cnt=0;
	#endif
}

void CLTreeS::addSeq(int *bid, int n, int *lmerbid, int seqID)  //call with n=L from outside
{
	if (seqID>maxSeqID) maxSeqID=seqID; 
	if (seqID<minSeqID) minSeqID=seqID; 

	#ifdef FAST_TRACK
		FT_seq=bid; FT_cnt++; FT_seqID = seqID;
	#endif

	if (n==1)
	{
		if (this->daughter[*bid].t==NULL)
		{
			LTreeSnodeData* nodeData = new LTreeSnodeData; 
			nodeData->n = 1;
			nodeData->seqIDs.i= seqID;
//			nodeData->baseID = lmerbid;
			this->daughter[*bid].node=nodeData;
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		else
		{
			LTreeSnodeData* curnodeData =this->daughter[*bid].node;
			if (curnodeData->n==1)
			{
				intintptr newseqids; 
				newseqids.p= new int[2]; 
				newseqids.p[0]= curnodeData->seqIDs.i; 
				newseqids.p[1]= seqID;
				curnodeData->seqIDs.p = newseqids.p;
				curnodeData->n=2;
			}
			else
			{
				if ((curnodeData->n & (curnodeData->n-1))==0) // i.e. n is power of 2
				{
					// expand memory 
					intintptr newseqids; 
					newseqids.p= new int[(curnodeData->n)<<1]; 
					for(int j=0;j< curnodeData->n; j++)
					{
						newseqids.p[j]= curnodeData->seqIDs.p[j]; 
					}		
					delete []curnodeData->seqIDs.p; 
					curnodeData->seqIDs.p = newseqids.p; 
				}
				curnodeData->seqIDs.p[curnodeData->n] = seqID;
				curnodeData->n += 1;
			}
		}
	}
	else
	{
		if (this->daughter[*bid].t == NULL)
		{
			this->daughter[*bid].t = new CLTreeS(); 
			this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
		}
		daughter[*bid].t->addSeq(bid+1, n-1, lmerbid, seqID); 
	}
}

void CLTreeS::addLTreeSnodeData(int *bid, int n, LTreeSnodeData* nodeData, int mnSeqID, int mxSeqID)
// similar to add seq, but adds multiple seqs together, used by iDLbound only to generate reordered trees
//call with n=L from outside
{
    if (mxSeqID>maxSeqID) maxSeqID=mxSeqID;
    if (mnSeqID<minSeqID) minSeqID=mnSeqID;
    
#ifdef FAST_TRACK
//    FT_seq=bid; FT_cnt++; FT_seqID = seqID;
    Printf(" fast track not supported with iDL bound! ERROR \n");
    return;//exit(1);
#endif
    
    if (n==1)
    {
        if (this->daughter[*bid].t==NULL)
        {
//            LTreeSnodeData* nodeData = new LTreeSnodeData;
//            nodeData->n = 1;
//            nodeData->seqIDs.i= seqID;
//            nodeData->baseID = lmerbid;
            this->daughter[*bid].node=nodeData;
            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
        }
        else
        {
            Printf(" nonempty node not expected Error !\n");
            return; //exit(1);
            /*
            LTreeSnodeData* curnodeData =this->daughter[*bid].node;
            if (curnodeData->n==1)
            {
                intintptr newseqids;
                newseqids.p= new int[2];
                newseqids.p[0]= curnodeData->seqIDs.i;
                newseqids.p[1]= seqID;
                curnodeData->seqIDs.p = newseqids.p;
                curnodeData->n=2;
            }
            else
            {
                if ((curnodeData->n & (curnodeData->n-1))==0) // i.e. n is power of 2
                {
                    // expand memory
                    intintptr newseqids;
                    newseqids.p= new int[(curnodeData->n)<<1];
                    for(int j=0;j< curnodeData->n; j++)
                    {
                        newseqids.p[j]= curnodeData->seqIDs.p[j];
                    }
                    delete []curnodeData->seqIDs.p;
                    curnodeData->seqIDs.p = newseqids.p; 
                }
                curnodeData->seqIDs.p[curnodeData->n] = seqID;
                curnodeData->n += 1;
            }
            */
        }
    }
    else
    {
        if (this->daughter[*bid].t == NULL)
        {
            this->daughter[*bid].t = new CLTreeS(); 
            this->nonEmptyDaughterIdxs[this->nonEmptyDaughterCnt++]=*bid;
        }
        daughter[*bid].t->addLTreeSnodeData(bid+1, n-1, nodeData,  mnSeqID,  mxSeqID );
    }
}

void CLTreeS::addToGTree(GTree2 *gtree, int n,int *tmpArray,int alphabetSize, int L){
    if (n==1)
    {
        for(int bid=0;bid<alphabetSize;bid++)
        {
            
            if (this->daughter[bid].t==NULL) continue;
            tmpArray[L-n]=bid;
            gtree->addLTreeSnodeData(tmpArray, L, daughter[bid].node, gMAXMM, 0);
        }
    }
    else
    {
        for(int bid=0;bid<alphabetSize;bid++)
        {
            if (this->daughter[bid].t != NULL)
            {
                tmpArray[L-n]=bid;
                daughter[bid].t->addToGTree(gtree, n-1,tmpArray, alphabetSize,L);
            }
        }
    }
    
}


int CLTreeS::addToList(LTreeSnodeData **list, int n, int single, int listlen, int alphabetSize) // // adds all the Lmers (that are (not) present in single sequence) to a list
{
	if (n==1)
	{
		for(int bid=0;bid<alphabetSize;bid++)
		{

			if (this->daughter[bid].t==NULL) continue; 

			if ((this->daughter[bid].node->n==1)==single)
			{
				list[listlen]= this->daughter[bid].node; 
				listlen++;
			}
		}
	}
	else
	{
		for(int bid=0;bid<alphabetSize;bid++)
		{
			if (this->daughter[bid].t != NULL)
			{
				listlen=daughter[bid].t->addToList(list, n-1, single, listlen, alphabetSize);
			}
		}
	}

	return listlen; 
}

int *CLTreeS::reorder(int *lmer, int *order, int L, int *output){
    for(int i=0;i<L;i++){
        output[i]=lmer[order[i]];
    }
    return(output);
}
void CLTreeS::cloneReorder(CLTreeS *newTree, int *order, int n, int L,int alphabetSize, int *tmpArray, int *tmpArray2){ // reorders and clones
    

    if (n==1)
    {
        for(int bid=0;bid<alphabetSize;bid++)
        {
            
            if (this->daughter[bid].t==NULL) continue;
            tmpArray[L-n]=bid;
            tmpArray2=reorder(tmpArray, order, L, tmpArray2);
            newTree->addLTreeSnodeData(tmpArray2, L, this->daughter[bid].node, this->minSeqID, this->maxSeqID);
            
            //for(int k=0;k<L;k++){printf("%c",::globalConverter.icidx[tmpArray[k]]);}
            //printf("->");
            //for(int k=0;k<L;k++){printf("%c",::globalConverter.icidx[tmpArray2[k]]);}
            //printf("\n");
            
            
        }
    }
    else
    {
        for(int bid=0;bid<alphabetSize;bid++)
        {
            if (this->daughter[bid].t != NULL)
            {
                tmpArray[L-n]=bid;
                daughter[bid].t->cloneReorder(newTree, order, n-1, L,alphabetSize, tmpArray, tmpArray2);
            }
        }
    }
    
}
/*
int CLTreeS::DFSn0(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, LTreeSnodeData *nodei, int bid)
{
	int nmulti =0; 
//	for(nmulti = 0; nmulti<listlen;nmulti++)
//	{
//		if (matchingLmers[nmulti]->n==1) break; 
//	}

	int i,k; 

//	for(int bid=0;bid<4;bid++)
//	{
//		if(daughter[bid].node==NULL) continue; 
//		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int nodei_seqIDs_i=nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[nodei_seqIDs_i];

			int *curMismatchCnt_j = curMismatchCnt; 
			LTreeSnodeData **matchingLmers_j = matchingLmers;

			for(nmulti=0;nmulti<listlen;nmulti++)
			{
				if (matchingLmers[nmulti]->n==1) break; 
				if ((*matchingLmers_j)->baseID[gLM1]==bid)
				{
					int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
					int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
					int matchingLmersj_n=(*matchingLmers_j)->n; 
					for(k=0;k<matchingLmersj_n;k++)
					{
						if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
						mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
						matchingLmers_j_seqIDs_p++;
					}
				}
				else
				{
					if ((*curMismatchCnt_j)<gMAXMM)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
				}
				curMismatchCnt_j++;
				matchingLmers_j++;

			}

			int j; 
			for(j=nmulti;j<listlen;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					mmprofile[*curMismatchCnt_j][(*matchingLmers_j)->seqIDs.i]++;
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
				}
				curMismatchCnt_j++;
				matchingLmers_j++;
			}
		}
		else
		{

			for(nmulti = 0; nmulti<listlen;nmulti++)
			{
				if (matchingLmers[nmulti]->n==1) break; 
			}
			int j; 

			for (int i=0;i<nodei->n;i++)
			{
				int nodei_seqIDs_pi=nodei->seqIDs.p[i]; 
				int **mmprofile=gMMProfile[nodei_seqIDs_pi];
				//int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];

				int *curMismatchCnt_j = curMismatchCnt; 
				LTreeSnodeData **matchingLmers_j = matchingLmers;

				for(j=0;j<nmulti;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						{
							int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
							int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
							int matchingLmersj_n=(*matchingLmers_j)->n; 
							for(k=0;k<(*matchingLmers_j)->n;k++)
							{
//								mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.p[k]]++;
								if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
								mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
								matchingLmers_j_seqIDs_p++;

							}
						}
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}

				for(j=nmulti;j<listlen;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.i]++;
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}
			}
		}
		return (nmulti); 
	}
	
//}


void CLTreeS::DFSn1(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, LTreeSnodeData *nodei, int bid, int nmulti)
{

	int i,j,k; 

//	for(int bid=0;bid<4;bid++)
//	{
//		if(daughter[bid].node==NULL) continue; 
//		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int nodei_seqIDs_i=nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[nodei_seqIDs_i];

			int *curMismatchCnt_j = curMismatchCnt; 
			LTreeSnodeData **matchingLmers_j = matchingLmers;

			for(j=0;j<nmulti;j++)
			{
				if ((*matchingLmers_j)->baseID[gLM1]==bid)
				{
					int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
					int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
					int matchingLmersj_n=(*matchingLmers_j)->n; 
					for(k=0;k<matchingLmersj_n;k++)
					{
						if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
						mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
						matchingLmers_j_seqIDs_p++;
					}
				}
				else
				{
					if ((*curMismatchCnt_j)<gMAXMM)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
				}
				curMismatchCnt_j++;
				matchingLmers_j++;

			}

			for(j=nmulti;j<listlen;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					mmprofile[*curMismatchCnt_j][(*matchingLmers_j)->seqIDs.i]++;
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
				}
				curMismatchCnt_j++;
				matchingLmers_j++;
			}
		}
		else
		{
			for (int i=0;i<nodei->n;i++)
			{
				int nodei_seqIDs_pi=nodei->seqIDs.p[i]; 
				int **mmprofile=gMMProfile[nodei_seqIDs_pi];
				//int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];

				int *curMismatchCnt_j = curMismatchCnt; 
				LTreeSnodeData **matchingLmers_j = matchingLmers;

				for(j=0;j<nmulti;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						{
							int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
							int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
							int matchingLmersj_n=(*matchingLmers_j)->n; 
							for(k=0;k<(*matchingLmers_j)->n;k++)
							{
//								mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.p[k]]++;
								if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
								mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
								matchingLmers_j_seqIDs_p++;

							}
						}
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}

				for(j=nmulti;j<listlen;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.i]++;
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}
			}
		}
	}
//}

void CLTreeS::DFSn(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt)
{
	int nmulti =0; 
	for(nmulti = 0; nmulti<listlen;nmulti++)
	{
		if (matchingLmers[nmulti]->n==1) break; 
	}

	int i,j,k; 

	for(int bid=0;bid<4;bid++)
	{
		if(daughter[bid].node==NULL) continue; 
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int nodei_seqIDs_i=nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[nodei_seqIDs_i];

			int *curMismatchCnt_j = curMismatchCnt; 
			LTreeSnodeData **matchingLmers_j = matchingLmers;

			for(j=0;j<nmulti;j++)
			{
				if ((*matchingLmers_j)->baseID[gLM1]==bid)
				{
					int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
					int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
					int matchingLmersj_n=(*matchingLmers_j)->n; 
					for(k=0;k<matchingLmersj_n;k++)
					{
						if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
						mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
						matchingLmers_j_seqIDs_p++;
					}
				}
				else
				{
					if ((*curMismatchCnt_j)<gMAXMM)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
				}
				curMismatchCnt_j++;
				matchingLmers_j++;

			}

			for(j=nmulti;j<listlen;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					mmprofile[*curMismatchCnt_j][(*matchingLmers_j)->seqIDs.i]++;
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
				}
				curMismatchCnt_j++;
				matchingLmers_j++;
			}
		}
		else
		{
			for (int i=0;i<nodei->n;i++)
			{
				int nodei_seqIDs_pi=nodei->seqIDs.p[i]; 
				int **mmprofile=gMMProfile[nodei_seqIDs_pi];
				//int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];

				int *curMismatchCnt_j = curMismatchCnt; 
				LTreeSnodeData **matchingLmers_j = matchingLmers;

				for(j=0;j<nmulti;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						{
							int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
							int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
							int matchingLmersj_n=(*matchingLmers_j)->n; 
							for(k=0;k<(*matchingLmers_j)->n;k++)
							{
//								mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.p[k]]++;
								if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
								mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
								matchingLmers_j_seqIDs_p++;

							}
						}
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}

				for(j=nmulti;j<listlen;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.i]++;
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}
			}
		}
	}
}

void CLTreeS::DFSnMulti(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt)
{
	int i,j,k; 

	for(int bid=0;bid<4;bid++)
	{
		if(daughter[bid].node==NULL) continue; 
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int nodei_seqIDs_i=nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[nodei_seqIDs_i];

			int *curMismatchCnt_j = curMismatchCnt; 
			LTreeSnodeData **matchingLmers_j = matchingLmers;

			for(j=0;j<listlen;j++)
			{
				if ((*matchingLmers_j)->baseID[gLM1]==bid)
				{
					int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
					int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
					int matchingLmersj_n=(*matchingLmers_j)->n; 
					for(k=0;k<matchingLmersj_n;k++)
					{
						if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
						mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
						matchingLmers_j_seqIDs_p++;
					}
				}
				else
				{
					if ((*curMismatchCnt_j)<gMAXMM)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_i) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
				}
				curMismatchCnt_j++;
				matchingLmers_j++;

			}

		}
		else
		{
			for (int i=0;i<nodei->n;i++)
			{
				int nodei_seqIDs_pi=nodei->seqIDs.p[i]; 
				int **mmprofile=gMMProfile[nodei_seqIDs_pi];
				//int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];

				int *curMismatchCnt_j = curMismatchCnt; 
				LTreeSnodeData **matchingLmers_j = matchingLmers;

				for(j=0;j<listlen;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)]; 
						int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
						int matchingLmersj_n=(*matchingLmers_j)->n; 
						for(k=0;k<matchingLmersj_n;k++)
						{
//							mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.p[k]]++;
							if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
							mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
							matchingLmers_j_seqIDs_p++;
						}
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						{
							int *mmprofile_curMismatchCnt_j = mmprofile[(*curMismatchCnt_j)+1]; 
							int *matchingLmers_j_seqIDs_p = (*matchingLmers_j)->seqIDs.p; 
							int matchingLmersj_n=(*matchingLmers_j)->n; 
							for(k=0;k<(*matchingLmers_j)->n;k++)
							{
//								mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.p[k]]++;
								if (*matchingLmers_j_seqIDs_p>nodei_seqIDs_pi) break; 
								mmprofile_curMismatchCnt_j[*matchingLmers_j_seqIDs_p]++;
								matchingLmers_j_seqIDs_p++;

							}
						}
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}

			}
		}
	}
}

void CLTreeS::DFSnSingle(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt)
{
	int i,j,k; 

	for(int bid=0;bid<4;bid++)
	{
		if(daughter[bid].node==NULL) continue; 
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int nodei_seqIDs_i=nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[nodei_seqIDs_i];

			int *curMismatchCnt_j = curMismatchCnt; 
			LTreeSnodeData **matchingLmers_j = matchingLmers;

			for(j=0;j<listlen;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					mmprofile[*curMismatchCnt_j][(*matchingLmers_j)->seqIDs.i]++;
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
				}
				curMismatchCnt_j++;
				matchingLmers_j++;
			}
		}
		else
		{
			for (int i=0;i<nodei->n;i++)
			{
				int nodei_seqIDs_pi=nodei->seqIDs.p[i]; 
				int **mmprofile=gMMProfile[nodei_seqIDs_pi];
				//int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];

				int *curMismatchCnt_j = curMismatchCnt; 
				LTreeSnodeData **matchingLmers_j = matchingLmers;

				for(j=0;j<listlen;j++)
				{
					if ((*matchingLmers_j)->baseID[gLM1]==bid)
					{
						mmprofile[(*curMismatchCnt_j)][(*matchingLmers_j)->seqIDs.i]++;
					}
					else
					{
						if ((*curMismatchCnt_j)<gMAXMM)
						mmprofile[(*curMismatchCnt_j)+1][(*matchingLmers_j)->seqIDs.i]++;
					}
					curMismatchCnt_j++;
					matchingLmers_j++;
				}
			}
		}
	}
}
*/
/*
void CLTreeS::DFSnf(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt)
{
	int nmulti =0; 

	if(daughter[0].node==NULL) 
	{
		if(daughter[1].node==NULL) 
		{
			if(daughter[2].node==NULL) 
			{
				DFSn0(matchingLmers, listlen, curMismatchCnt, daughter[3].node, 3); 
			}
			else
			{
				nmulti = DFSn0(matchingLmers, listlen, curMismatchCnt, daughter[2].node, 2); 
				if(daughter[3].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[3].node, 3, nmulti); 
			}
		}
		else
		{
			nmulti = DFSn0(matchingLmers, listlen, curMismatchCnt, daughter[1].node, 1); 
			if(daughter[2].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[2].node, 2, nmulti); 
			if(daughter[3].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[3].node, 3, nmulti); 
		}
	}
	else
	{
		nmulti = DFSn0(matchingLmers, listlen, curMismatchCnt, daughter[0].node, 0); 
		if(daughter[1].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[1].node, 1, nmulti); 
		if(daughter[2].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[2].node, 2, nmulti); 
		if(daughter[3].node!=NULL) DFSn1(matchingLmers, listlen, curMismatchCnt, daughter[3].node, 3, nmulti); 

	}

}

void CLTreeS::DFS( LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, int pos)
{
	//LTreeSnodeData **matchingLmers = gDFSlist[pos]; 
	//int *curMismatchCnt= gDFSMMlist[pos];
		
	if(pos==gLM1) //LM1 is L-1
	{
		DFSnf(matchingLmers, listlen, curMismatchCnt); // process the node. 
	}
	else
	{
		LTreeSnodeData **newlist = gDFSlist[pos+1]; 
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		LTreeSnodeData **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		for(int bid=0;bid<4;bid++)
		{
			if(daughter[bid].t==NULL) continue; 
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;

			for(int i=0;i<listlen;i++)
			{
				if(matchingLmers[i]->baseID[pos]==bid)
				{
					*newlistnewlistlen=matchingLmers[i];
					newlistnewlistlen++;
					*newMismatchCntnewlistlen=curMismatchCnt[i];
					newMismatchCntnewlistlen++; 
					newlistlen++;
				}
				else
				{
					if (curMismatchCnt[i]<gMAXMM)
					{
						*newlistnewlistlen=matchingLmers[i];
						newlistnewlistlen++;
						*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
						newMismatchCntnewlistlen++; 
						newlistlen++;
					}
				}
			}
			if (newlistlen!=0)
			{
				daughter[bid].t->DFS(newlist,newlistlen,newMismatchCnt,pos+1); 
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
*/
/*
// without nonEmptyDaughterCnt
void CLTreeS::DFSTn(CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize)
{
	int i,j,k; 

	for(int bid=0;bid<alphabetSize;bid++)
	{
		bid = this->nonEmptyDaughterIdxs[ibid];
		if(daughter[bid].node==NULL) continue; 
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int curnodeid = nodei->seqIDs.i; 
			int **mmprofile=gMMProfile[curnodeid];

			for(int fbid=0;fbid<alphabetSize;fbid++)
			{

				if (bid==fbid)
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].node!=NULL)
						{
							if (matchingLmers[i][fbid].node->n==1)
							{
								mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.i]++;
							}
							else
							{
								for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
								{
									if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break; 
									mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.p[j]]++;
								}
							}
						}
					}
				}
				else
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].node!=NULL)
						{
							if (curMismatchCnt[i]<gMAXMM)
							{
								if (matchingLmers[i][fbid].node->n==1)
								{
									mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.i]++;
								}
								else
								{
									for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
									{
										if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break; 
										mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.p[j]]++;
									}
								}
							}
						}
					}

				}
			}
		}
		else
		{
			for(int k=0;k<nodei->n;k++)
			{
				int curnodeid = nodei->seqIDs.p[k]; 
				int **mmprofile=gMMProfile[curnodeid];

				for(int fbid=0;fbid<alphabetSize;fbid++)
				{

					if (bid==fbid)
					{
						for(int i=0;i<listlen;i++)
						{
							if(matchingLmers[i][fbid].node!=NULL)
							{
								if (matchingLmers[i][fbid].node->n==1)
								{
									mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.i]++;
								}
								else
								{
									for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
									{
										if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break; 
										mmprofile[curMismatchCnt[i]][matchingLmers[i][fbid].node->seqIDs.p[j]]++;
									}
								}
							}
						}
					}
					else
					{
						for(int i=0;i<listlen;i++)
						{
							if(matchingLmers[i][fbid].node!=NULL)
							{
								if (curMismatchCnt[i]<gMAXMM)
								{
									if (matchingLmers[i][fbid].node->n==1)
									{
										mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.i]++;
									}
									else
									{
										for(int j=0;j<matchingLmers[i][fbid].node->n;j++)
										{
											if (matchingLmers[i][fbid].node->seqIDs.p[j]>curnodeid) break; 
											mmprofile[curMismatchCnt[i]+1][matchingLmers[i][fbid].node->seqIDs.p[j]]++;
										}
									}
								}
							}
						}

					}
				}
			}
		}
	}

}
// without nonEmptyDaughterCnt

void CLTreeS::DFST( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize)
{
	//LTreeSnodeData **matchingLmers = gDFSlist[pos]; 
	//int *curMismatchCnt= gDFSMMlist[pos];
		
	if(pos==gLM1) //LM1 is L-1
	{
		DFSTn(matchingLmers, listlen, curMismatchCnt, alphabetSize); // process the node.
	}
	else
	{
		CLTreeSptr **newlist = gDFSlistT[pos+1]; 
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		CLTreeSptr **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		//int alphabetSize = ::globalConverter.b; // for DNA it is 4
		for(int bid=0;bid<alphabetSize;bid++)
		{
			if(daughter[bid].t==NULL) continue; 
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;
			int daughter_maxSeqID = daughter[bid].t->maxSeqID; 
			for(int fbid=0;fbid<alphabetSize;fbid++) //  foreign bid
			{
				if (bid==fbid)
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].t!=NULL)
						{
							if (matchingLmers[i][fbid].t->minSeqID >daughter_maxSeqID) continue; 
							*newlistnewlistlen=matchingLmers[i][fbid].t->daughter;
							newlistnewlistlen++;
							*newMismatchCntnewlistlen=curMismatchCnt[i];
							newMismatchCntnewlistlen++; 
							newlistlen++;
						}
					}
				}
				else
				{
					for(int i=0;i<listlen;i++)
					{
						if(matchingLmers[i][fbid].t!=NULL)
						{
							if (curMismatchCnt[i]<gMAXMM)
							{
								if (matchingLmers[i][fbid].t->minSeqID >daughter_maxSeqID) continue; 

								*newlistnewlistlen=matchingLmers[i][fbid].t->daughter;
								newlistnewlistlen++;
								*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
								newMismatchCntnewlistlen++; 
								newlistlen++;
							}
						}
					}
				}
			}

			if (newlistlen!=0)
			{
				daughter[bid].t->DFST(newlist,newlistlen,newMismatchCnt,pos+1,alphabetSize);
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
*/

// with nonEmptyDaughterCnt
void CLTreeS::DFSTn(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int alphabetSize)
{
//	for(int bid=0;bid<alphabetSize;bid++)
//	{
	int bid;
	for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
	{
		bid = this->nonEmptyDaughterIdxs[ibid];
//		if(daughter[bid].node==NULL) continue;
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int curnodeid = nodei->seqIDs.i;
			int **mmprofile=gMMProfile[curnodeid];

			for(int i=0;i<listlen;i++)
			{
				CLTreeS *imatchingLmer =matchingLmers[i];
				int fbid;
				for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
				{
					fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];

					if (bid==fbid)
					{

						LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
						if (nodej->n==1)
						{
							mmprofile[curMismatchCnt[i]][nodej->seqIDs.i]++;
						}
						else
						{
							for(int j=0;j<nodej->n;j++)
							{
								if (nodej->seqIDs.p[j]>curnodeid) break;
								mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
							}
						}
					}
					else
					{
						if (curMismatchCnt[i]<gMAXMM){
							LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
							if (nodej->n==1)
							{
								mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.i]++;
							}
							else
							{
								for(int j=0;j<nodej->n;j++)
								{
									if (nodej->seqIDs.p[j]>curnodeid) break;
									mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
								}
							}

						}

					}
				}
			}
		}

		else
		{
			for(int k=0;k<nodei->n;k++)
			{
				int curnodeid = nodei->seqIDs.p[k];
				int **mmprofile=gMMProfile[curnodeid];

				for(int i=0;i<listlen;i++)
				{
					CLTreeS *imatchingLmer =matchingLmers[i];
					int fbid;
					for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
					{
						fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];

						if (bid==fbid)
						{

							LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
							if (nodej->n==1)
							{
								mmprofile[curMismatchCnt[i]][nodej->seqIDs.i]++;
							}
							else
							{
								for(int j=0;j<nodej->n;j++)
								{
									if (nodej->seqIDs.p[j]>curnodeid) break;
									mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
								}
							}
						}
						else
						{
							if (curMismatchCnt[i]<gMAXMM){
								LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
								if (nodej->n==1)
								{
									mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.i]++;
								}
								else
								{
									for(int j=0;j<nodej->n;j++)
									{
										if (nodej->seqIDs.p[j]>curnodeid) break;
										mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
									}
								}

							}

						}
					}
				}

			}
		}
	}
}

void addmmprof(int *mmprofile_i,int *nodej_seqIDs_p,int nn, int curnodeid)
{
    for(int j=0;j<nn;j++)
    //for(int j=0, length=nn; j<length; j++)

   // while(nn--)
    {
        if (*nodej_seqIDs_p>curnodeid) return;
        mmprofile_i[*nodej_seqIDs_p++]++;
  }
}


// with nonEmptyDaughterCnt
// with iDL bound
void CLTreeS::DFSTnIDL(CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, CbinMMtree **curMMtree, int alphabetSize)
{
    
    // note: similar to DFST, we should first check if mismatch not allowed, don't iterate over mismatch places and go directly to match


    int bid;
    for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
    {
        bid = this->nonEmptyDaughterIdxs[ibid];
        //		if(daughter[bid].node==NULL) continue;
        LTreeSnodeData *nodei=daughter[bid].node;
        int nodei_n=nodei->n;
        if (nodei_n==1)
        {
            int curnodeid = nodei->seqIDs.i;
            int **mmprofile=gMMProfile[curnodeid];
            
            for(int i=0;i<listlen;i++)
            {
                CLTreeS *imatchingLmer =matchingLmers[i];
                CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
                CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;

                int fbid;
                for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
                {
                    fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
                    
                    if (bid==fbid)
                    {
                        if(icurMMtreeMatch!=nullptr){
                            
                            LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
                            if (nodej->n==1)
                            {
                                mmprofile[curMismatchCnt[i]][nodej->seqIDs.i]++;
                            }
                            else
                            {

                               /* for(int j=0;j<nodej->n;j++)
                                {
                                    if (nodej->seqIDs.p[j]>curnodeid) break;
                                    mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                                }*/
                                addmmprof(mmprofile[curMismatchCnt[i]],nodej->seqIDs.p,nodej->n, curnodeid);

                            }
                        }
                    }
                    else
                    {
//                        if (curMismatchCnt[i]<gMAXMM){
                        if(icurMMtreeMisMatch!=nullptr){
                            //it tolerates one more mismatch

                            LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
                            if (nodej->n==1)
                            {
                                mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.i]++;
                            }
                            else
                            {
                                /*
                                for(int j=0;j<nodej->n;j++)
                                {
                                    if (nodej->seqIDs.p[j]>curnodeid) break;
                                    mmprofile[1+curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                                }*/
                                addmmprof(mmprofile[1+curMismatchCnt[i]],nodej->seqIDs.p,nodej->n, curnodeid);

                            }
                            
                        }
                        
                    }
                }
            }
        }
        
        else
        {
            
                for(int i=0;i<listlen;i++)
                {
                    CLTreeS *imatchingLmer =matchingLmers[i];
                    CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
                    CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;
               
                    
                    int fbid;
              //      if(icurMMtreeMisMatch!=nullptr){
                            //it tolerates one more mismatch
                        for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
                        {
                            fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
                            
                            if (bid==fbid)
                            {
                                if(icurMMtreeMatch!=nullptr){

                                
                                    LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
                                    if (nodej->n==1)
                                    {
                                        
                                        
                                        for(int k=0;k<nodei_n;k++)
                                        {
                                            int curnodeid = nodei->seqIDs.p[k];
                                            int **mmprofile=gMMProfile[curnodeid];
                                            
                                            //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                                            int *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                                            
                                            mmprofile_curMismatchCnt_i[nodej->seqIDs.i]++;

                                        }
                                        
                                        
                                    }
                                    else
                                    {
                                        /*
                                        for(int j=0;j<nodej->n;j++)
                                        {
                                            if (nodej->seqIDs.p[j]>curnodeid) break;
                                            mmprofile[curMismatchCnt[i]][nodej->seqIDs.p[j]]++;
                                        }*/

                                        
                                        for(int k=0;k<nodei_n;k++)
                                        {
                                            int curnodeid = nodei->seqIDs.p[k];
                                            int **mmprofile=gMMProfile[curnodeid];
                                            
                                            //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                                            int *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                                            
                                            addmmprof(mmprofile_curMismatchCnt_i,nodej->seqIDs.p,nodej->n, curnodeid);
                                        }

                                    }
                                }
                            }
                            else
                            {
    //                            if (curMismatchCnt[i]<gMAXMM){
                                if(icurMMtreeMisMatch!=nullptr){
                                    LTreeSnodeData *nodej=imatchingLmer->daughter[fbid].node;
                                    int nodej_n=nodej->n;
                                    if (nodej_n==1)
                                    {
                                        
                                        for(int k=0;k<nodei->n;k++)
                                        {
                                            int curnodeid = nodei->seqIDs.p[k];
                                            int **mmprofile=gMMProfile[curnodeid];
                                            
                                            int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                                            
                                            mmprofile_curMismatchCntP1_i[nodej->seqIDs.i]++;

                                        }
                                    }
                                    else
                                    {
                                        for(int k=0;k<nodei->n;k++)
                                        {
                                            int curnodeid = nodei->seqIDs.p[k];
                                            //int **mmprofile=gMMProfile[curnodeid];
                                            
                                            //int *mmprofile_curMismatchCntP1_i=mmprofile[1+curMismatchCnt[i]];
                                            //int *mmprofile_curMismatchCnt_i=mmprofile[curMismatchCnt[i]];
                                            addmmprof(gMMProfile[curnodeid][1+curMismatchCnt[i]],nodej->seqIDs.p,nodej_n, curnodeid);
                                        }
                                    }
                                    
                                }
                                
                            }
                        }
            /*        }else{
                            // no more mismatches. just use the match
                        //for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
                        LTreeSnodeData *nodej=imatchingLmer->daughter[bid].node;
                        
                        if (nodej!=NULL){  // if the neighbor has daughter with matching base
                            if (nodej->n==1)
                            {
                                mmprofile_curMismatchCnt_i[nodej->seqIDs.i]++;
                            }
                            else
                            {
             
                                
                                addmmprof(mmprofile_curMismatchCnt_i,nodej->seqIDs.p,nodej->n, curnodeid);
                                
                            }
                        }
                        
                        
                        
                    }*/
  //              }
                
            }
        }
    }
}

// with nonEmptyDaughterCnt
void CLTreeS::DFST( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt, int pos, int alphabetSize)
{
    //static int cntt=0;
    //cntt++; if(cntt%1000==0){printf(" cnt %d \n",cntt);}
	//LTreeSnodeData **matchingLmers = gDFSlist[pos];
	//int *curMismatchCnt= gDFSMMlist[pos];
	#ifdef FAST_TRACK
		if((this->FT_cnt==1)&&(listlen==1)&&((*matchingLmers)->FT_cnt==1)){
			// fast track: this is for speed up only.
			int *bid1 = this->FT_seq;
			int *bid2 = (*matchingLmers)->FT_seq;

			for(int k=0;k<pos;k++){
				if (bid1[k]!=bid2[k]){
					(*curMismatchCnt)++;
					if ((*curMismatchCnt) > gMAXMM){
						return;
					}
				}
			}
			gMMProfile[this->FT_seqID][*curMismatchCnt][(*matchingLmers)->FT_seqID]++;
			gMMProfile[(*matchingLmers)->FT_seqID][*curMismatchCnt][this->FT_seqID]++;
			return;
		}
	#endif
	if(pos==gLM1) //LM1 is L-1
	{
		DFSTn(matchingLmers, listlen, curMismatchCnt, alphabetSize); // process the node.
	}
	else
	{
//		CLTreeS **newlist = new CLTreeS *[alphabetSize*listlen];
//		int *newMismatchCnt= new int[alphabetSize*listlen];
		CLTreeS **newlist = gDFSlistT[pos+1];
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0;
		CLTreeS **newlistnewlistlen = newlist; //&newlist[newlistlen]
		int *newMismatchCntnewlistlen = newMismatchCnt;
		//int alphabetSize = ::globalConverter.b; // for DNA it is 4
//		for(int bid=0;bid<alphabetSize;bid++)
//		{
		int bid;
		for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
		{
			bid = this->nonEmptyDaughterIdxs[ibid];
			//if(daughter[bid].t==NULL) continue;
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;
			int daughter_maxSeqID = daughter[bid].t->maxSeqID;
			for(int i=0;i<listlen;i++)
			{
				//CLTreeSptr *imatchingLmer =matchingLmers[i];
				CLTreeS *imatchingLmer =matchingLmers[i];
				int fbid;
				
                
                if ((curMismatchCnt[i]<gMAXMM)
                         &&((pos>=5)||(curMismatchCnt[i]<(gMAXMM*(1+pos))/(gLM1+1)))
                        //                            &&(curMismatchCnt[i]<((gMAXMM)*(1+pos))/(gLM1+1))
                        //                            &&((curMismatchCnt[i]+1)<=((1+gMAXMM)*(1+pos))/(gLM1+1))
                        //                            &&((curMismatchCnt[i]+1)<=(1.2*gMAXMM*(1+pos))/(gLM1+1))
                        //                            &&((pos>gLM1/3)||(curMismatchCnt[i]<gMAXMM/3))
                        //                            &&((pos>gLM1/2)||(curMismatchCnt[i]<gMAXMM/2))
                        )
                {  //it tolerates one more mismatch

                    for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
                    {
                        fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
                        
                        if (bid==fbid){
                            CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
                            if (newnode->minSeqID >daughter_maxSeqID) continue;
                            
                            *newlistnewlistlen=newnode;
                            newlistnewlistlen++;
                            *newMismatchCntnewlistlen=curMismatchCnt[i];
                            newMismatchCntnewlistlen++;
                            newlistlen++;
                        } else {
                            CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
                            if (newnode->minSeqID >daughter_maxSeqID) continue;
                            
                            *newlistnewlistlen=newnode;
                            newlistnewlistlen++;
                            *newMismatchCntnewlistlen=curMismatchCnt[i]+1;
                            newMismatchCntnewlistlen++;
                            newlistlen++;
                        }
                        
                    }
                    
                }else{
                   //no more mismatches tolerated at this depth
                    
                    CLTreeS *newnode =imatchingLmer->daughter[bid].t;
                    
                    if (newnode!=NULL){  // if the neighbor has daughter with matching base
                        if (newnode->minSeqID >daughter_maxSeqID) continue;
                        
                        *newlistnewlistlen=newnode;
                        newlistnewlistlen++;
                        *newMismatchCntnewlistlen=curMismatchCnt[i];
                        newMismatchCntnewlistlen++;
                        newlistlen++;
                    }
                }
			}

			if (newlistlen!=0)
			{
				daughter[bid].t->DFST(newlist,newlistlen,newMismatchCnt,pos+1,alphabetSize);
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}




// with nonEmptyDaughterCnt
void CLTreeS::DFSTiDL( CLTreeS **matchingLmers, int listlen, int *curMismatchCnt,CbinMMtree **curMMtree, int pos, int alphabetSize)
{
    if(pos==gLM1) //LM1 is L-1
    {
        DFSTnIDL(matchingLmers, listlen, curMismatchCnt, curMMtree, alphabetSize); // process the node.
    }
    else
    {
        CLTreeS **newlist = gDFSlistT[pos+1];
        int *newMismatchCnt= gDFSMMlist[pos+1];
        CbinMMtree **newMMtree= gDFSMMtree[pos+1];
        
        int newlistlen = 0;
        CLTreeS **newlistnewlistlen = newlist; //&newlist[newlistlen]
        int *newMismatchCntnewlistlen = newMismatchCnt;
        CbinMMtree **newMMtreenewlistlen =newMMtree;
        //int alphabetSize = ::globalConverter.b; // for DNA it is 4
        //		for(int bid=0;bid<alphabetSize;bid++)
        //		{
        int bid;
        for(int ibid=0;ibid<this->nonEmptyDaughterCnt;ibid++)
        {
            bid = this->nonEmptyDaughterIdxs[ibid];
            //if(daughter[bid].t==NULL) continue;
            newlistlen = 0;
            newlistnewlistlen = newlist;
            newMismatchCntnewlistlen = newMismatchCnt;
            newMMtreenewlistlen =newMMtree;
            int daughter_maxSeqID = daughter[bid].t->maxSeqID;
            for(int i=0;i<listlen;i++)
            {
                //CLTreeSptr *imatchingLmer =matchingLmers[i];
                CLTreeS *imatchingLmer =matchingLmers[i];
                
                CbinMMtree *icurMMtreeMatch = curMMtree[i]->child0;
                CbinMMtree *icurMMtreeMisMatch = curMMtree[i]->child1;
                
                int fbid;
                
                if(icurMMtreeMisMatch!=nullptr){
                    //it tolerates one more mismatch
                    
                    for(int jbid=0;jbid<imatchingLmer->nonEmptyDaughterCnt;jbid++)
                    {
                        fbid = imatchingLmer->nonEmptyDaughterIdxs[jbid];
                        
                        if (bid==fbid){  // Match
                            if(icurMMtreeMatch!=nullptr){ // it allows match
                                CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
                                if (newnode->minSeqID >daughter_maxSeqID) continue;
                            
                                *newlistnewlistlen=newnode;
                                newlistnewlistlen++;
                                *newMismatchCntnewlistlen=curMismatchCnt[i];
                                newMismatchCntnewlistlen++;
                                *newMMtreenewlistlen=icurMMtreeMatch;
                                newMMtreenewlistlen++;
                                newlistlen++;
                            }
                        } else {
                            CLTreeS *newnode =imatchingLmer->daughter[fbid].t;
                            if (newnode->minSeqID >daughter_maxSeqID) continue;
                            
                            *newlistnewlistlen=newnode;
                            newlistnewlistlen++;
                            *newMismatchCntnewlistlen=curMismatchCnt[i]+1;
                            newMismatchCntnewlistlen++;
                            *newMMtreenewlistlen=icurMMtreeMisMatch;
                            newMMtreenewlistlen++;
                            newlistlen++;
                            // note : infact mismatchCnt is redundant, we can use the tree for final mismatch cnt. and we don't need to add new var, use intptr for a child
                        }
                        
                    }
                    
                }else{
                    //no more mismatches tolerated at this depth
                //if(icurMMtreeMatch!=nullptr){  not needed because we assume it is not possible both child0 and 1 be null at internal nodes (depth<L)
                    CLTreeS *newnode =imatchingLmer->daughter[bid].t;
                    
                    if (newnode!=NULL){  // if the neighbor has daughter with matching base
                        if (newnode->minSeqID >daughter_maxSeqID) continue;
                        
                        *newlistnewlistlen=newnode;
                        newlistnewlistlen++;
                        *newMismatchCntnewlistlen=curMismatchCnt[i];
                        newMismatchCntnewlistlen++;
                        newlistlen++;
                        *newMMtreenewlistlen=icurMMtreeMatch;
                        newMMtreenewlistlen++;
                    }
                }
            }
            
            if (newlistlen!=0)
            {
                daughter[bid].t->DFSTiDL(newlist,newlistlen,newMismatchCnt,newMMtree, pos+1,alphabetSize);
            }
        }
        //delete []newlist;
        //delete []newMismatchCnt;
    }
}
/*
void CLTreeS::DFSTf( CLTreeSptr **matchingLmers, int listlen, int *curMismatchCnt, int pos)
{
	//LTreeSnodeData **matchingLmers = gDFSlist[pos]; 
	//int *curMismatchCnt= gDFSMMlist[pos];
		
	if(pos==gLM1) //LM1 is L-1
	{
		DFSTn(matchingLmers, listlen, curMismatchCnt); // process the node. 
	}
	else
	{
		CLTreeSptr **newlist = gDFSlistT[pos+1]; 
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		CLTreeSptr **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		for(int bid=0;bid<4;bid++)
		{
			if(daughter[bid].t==NULL) continue; 
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;
			int daughter_maxSeqID = daughter[bid].t->maxSeqID; 
			for(int fbid=0;fbid<4;fbid++) //  foreign bid
			{
				for(int i=0;i<listlen;i++)
				{
					class CLTreeS *matchingLmers_i_fbid_t = matchingLmers[i][fbid].t; 
					if(matchingLmers_i_fbid_t!=NULL)
					{
						if (bid==fbid)
						{
							if (matchingLmers_i_fbid_t->minSeqID >daughter_maxSeqID) continue; 
							*newlistnewlistlen=matchingLmers[i][fbid].t->daughter;
							newlistnewlistlen++;
							*newMismatchCntnewlistlen=curMismatchCnt[i];
							newMismatchCntnewlistlen++; 
							newlistlen++;
						}
						else
						{
							if (curMismatchCnt[i]<gMAXMM)
							{
								if (matchingLmers_i_fbid_t->minSeqID >daughter_maxSeqID) continue; 

								*newlistnewlistlen=matchingLmers_i_fbid_t->daughter;
								newlistnewlistlen++;
								*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
								newMismatchCntnewlistlen++; 
								newlistlen++;
							}
						}
					}
				}
			}

			if (newlistlen!=0)
			{
				daughter[bid].t->DFSTf(newlist,newlistlen,newMismatchCnt,pos+1); 
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}

*/
/*
void CLTreeS::DFSsingle( LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, int pos)
{
	if(pos==gLM1) //LM1 is L-1
	{
		DFSnSingle(matchingLmers, listlen, curMismatchCnt); // process the node. 
	}
	else
	{
		LTreeSnodeData **newlist = gDFSlist[pos+1]; 
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		LTreeSnodeData **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		for(int bid=0;bid<4;bid++)
		{
			class CLTreeS *daughter_bid = daughter[bid].t; 
			if(daughter_bid==NULL) continue; 
			int daughter_maxSeqID = daughter_bid->maxSeqID; 
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;

			for(int i=0;i<listlen;i++)
			{
				if (matchingLmers[i]->seqIDs.i > daughter_maxSeqID) break; 
				if(matchingLmers[i]->baseID[pos]==bid)
				{
					*newlistnewlistlen=matchingLmers[i];
					newlistnewlistlen++;
					*newMismatchCntnewlistlen=curMismatchCnt[i];
					newMismatchCntnewlistlen++; 
					newlistlen++;
				}
				else
				{
					if (curMismatchCnt[i]<gMAXMM)
					{
						*newlistnewlistlen=matchingLmers[i];
						newlistnewlistlen++;
						*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
						newMismatchCntnewlistlen++; 
						newlistlen++;
					}
				}
			}
			if (newlistlen!=0)
			{
				daughter[bid].t->DFSsingle(newlist,newlistlen,newMismatchCnt,pos+1); 
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}

void CLTreeS::DFSmulti( LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, int pos)
{
	if(pos==gLM1) //LM1 is L-1
	{
//		DFSnf(matchingLmers, listlen, curMismatchCnt); // process the node. 
		DFSnMulti(matchingLmers, listlen, curMismatchCnt); // process the node. 
	}
	else
	{
		LTreeSnodeData **newlist = gDFSlist[pos+1]; 
		int *newMismatchCnt= gDFSMMlist[pos+1];

		int newlistlen = 0; 
		LTreeSnodeData **newlistnewlistlen = newlist;
		int *newMismatchCntnewlistlen = newMismatchCnt;
		for(int bid=0;bid<4;bid++)
		{
			class CLTreeS *daughter_bid = daughter[bid].t; 
			if(daughter_bid==NULL) continue; 
			int daughter_maxSeqID = daughter_bid->maxSeqID; 
			newlistlen = 0;
			newlistnewlistlen = newlist;
			newMismatchCntnewlistlen = newMismatchCnt;

			for(int i=0;i<listlen;i++)
			{
				if (*(matchingLmers[i]->seqIDs.p) > daughter_maxSeqID) break; 
				if(matchingLmers[i]->baseID[pos]==bid)
				{
					*newlistnewlistlen=matchingLmers[i];
					newlistnewlistlen++;
					*newMismatchCntnewlistlen=curMismatchCnt[i];
					newMismatchCntnewlistlen++; 
					newlistlen++;
				}
				else
				{
					if (curMismatchCnt[i]<gMAXMM)
					{
						*newlistnewlistlen=matchingLmers[i];
						newlistnewlistlen++;
						*newMismatchCntnewlistlen=curMismatchCnt[i]+1;
						newMismatchCntnewlistlen++; 
						newlistlen++;
					}
				}
			}
			if (newlistlen!=0)
			{
				daughter[bid].t->DFSmulti(newlist,newlistlen,newMismatchCnt,pos+1); 
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
*/
/*

void CLTreeS::DFSn(LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt)
{
	int nmulti =0; 
	for(nmulti = 0; nmulti<listlen;nmulti++)
	{
		if (matchingLmers[nmulti]->n==1) break; 
	}

	int i,j,k; 

	for(int bid=0;bid<4;bid++)
	{
		if(daughter[bid].node==NULL) continue; 
		LTreeSnodeData *nodei=daughter[bid].node;
		if (nodei->n==1)
		{
			int **mmprofile=gMMProfile[nodei->seqIDs.i];

			for(j=0;j<nmulti;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					for(k=0;k<matchingLmers[j]->n;k++)
					{
						mmprofile[curMismatchCnt[j]][matchingLmers[j]->seqIDs.p[k]]++;
					}
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					for(k=0;k<matchingLmers[j]->n;k++)
					{
						mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
					}
				}
			}

			for(j=nmulti;j<listlen;j++)
			{
				if (matchingLmers[j]->baseID[gLM1]==bid)
				{
					mmprofile[curMismatchCnt[j]][matchingLmers[j]->seqIDs.i]++;
				}
				else
				{
					if (curMismatchCnt[j]<gMAXMM)
					mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.i]++;
				}
			}
		}
		else
		{
			for (int i=0;i<nodei->n;i++)
			{
				int **mmprofile=gMMProfile[nodei->seqIDs.p[i]];
				for(j=0;j<nmulti;j++)
				{
					if (matchingLmers[j]->baseID[gLM1]==bid)
					{
						for(k=0;k<matchingLmers[j]->n;k++)
						{
							mmprofile[curMismatchCnt[j]][matchingLmers[j]->seqIDs.p[k]]++;
						}
					}
					else
					{
						if (curMismatchCnt[j]<gMAXMM)
						for(k=0;k<matchingLmers[j]->n;k++)
						{
							mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.p[k]]++;
						}
					}
				}

				for(j=nmulti;j<listlen;j++)
				{
					if (matchingLmers[j]->baseID[gLM1]==bid)
					{
						mmprofile[curMismatchCnt[j]][matchingLmers[j]->seqIDs.i]++;
					}
					else
					{
						if (curMismatchCnt[j]<gMAXMM)
						mmprofile[curMismatchCnt[j]+1][matchingLmers[j]->seqIDs.i]++;
					}
				}
			}
		}
	}
}
*/
/*
void CLTreeS::DFS( LTreeSnodeData **matchingLmers, int listlen, int *curMismatchCnt, int pos)
{
	if(pos==gLM1) //LM1 is L-1
	{
		DFSn(matchingLmers, listlen, curMismatchCnt); // process the node. 
	}
	else
	{
		//LTreeSnodeData **newlist = new LTreeSnodeData*[listlen]; 
		LTreeSnodeData **newlist = gDFSlist[pos]; 
		//int *newMismatchCnt= new int[listlen];
		int *newMismatchCnt= gDFSMMlist[pos];
		int newlistlen = 0; 
		for(int bid=0;bid<4;bid++)
		{
			if(daughter[bid].t==NULL) continue; 
			newlistlen = 0;
			for(int i=0;i<listlen;i++)
			{
				if(matchingLmers[i]->baseID[pos]==bid)
				{
					newlist[newlistlen]=matchingLmers[i];
					newMismatchCnt[newlistlen]=curMismatchCnt[i];
					newlistlen++;
				}
				else
				{
					if (curMismatchCnt[i]<gMAXMM)
					{
						newlist[newlistlen]=matchingLmers[i];
						newMismatchCnt[newlistlen]=curMismatchCnt[i]+1;
						newlistlen++;
					}
				}
			}
			if (newlistlen!=0)
			{
				daughter[bid].t->DFS(newlist,newlistlen,newMismatchCnt,pos+1); 
			}
		}
		//delete []newlist;
		//delete []newMismatchCnt;
	}
}
*/
/*
void CLTreeS::addSeq(int *bid, int n, int cnt)
{
	if (n==1)
	{
		this->daughter[*bid] = (CLTreeS *)((intptr_t)(daughter[*bid])+cnt); 
	}
	else
	{
		if (this->daughter[*bid] == NULL)
		{
			this->daughter[*bid] = new CLTreeS(); 
		}
		daughter[*bid]->addSeq(bid+1, n-1, cnt); 
	}
}
*/


void CLTreeS::deleteTree(int n, int alphabetSize, int DontDeleteNodeData)
{
	if (n>1)
	{
		for (int i=0;i<alphabetSize;i++)
		{
			if (daughter[i].t!=NULL)
			{
				daughter[i].t->deleteTree(n-1,alphabetSize, DontDeleteNodeData);
				delete daughter[i].t;
			}
		}
	}
	if (n==1)
	{
        if (!DontDeleteNodeData){
            for (int i=0;i<alphabetSize;i++)
            {
                if (daughter[i].node!=NULL)
                {
                    if (daughter[i].node->n>1)
                    {
                        delete daughter[i].node->seqIDs.p;
                    }
                    delete daughter[i].node;
                }
            }
        }
	}

    for(int i=0;i<MAX_ALPHABET_SIZE;i++){
        daughter[i].t=NULL;
    }
    maxSeqID=0;
    minSeqID=0;
}


int CLTreeS::addSequence(int *bid, int n, int L, int seqID)  //adds all the L-subseqs 
{
	n = n-L+1;
	if (n<0) n=0;
	for(int i=0;i<n;i++)
	{
		addSeq(bid,L,bid, seqID);
		bid++;
	}
	return n; 
}

int CLTreeS::leavesCount(int withMultiplicity, int n, int alphabetSize, int *nodesAtDepth)  //returns the number of sequences in the tree.  //call with n=L from outside
{
    if(nodesAtDepth!=nullptr){
        (*nodesAtDepth)++;
        nodesAtDepth++;
    }
	int nleaves = 0; 
	for (int i=0;i<alphabetSize;i++)
	{
		if (daughter[i].t!=NULL)
		{
			if (n==1)
			{
				if (withMultiplicity)
				{
					nleaves +=daughter[i].node->n; 
				}
				else 
				{
					nleaves++; 					
				}
			}
			else
			{
				nleaves+=daughter[i].t->leavesCount(withMultiplicity, n-1,alphabetSize, nodesAtDepth/*+1*/);
			}
		}
	}
	return nleaves; 
}



/*
int CLTreeS::count(int *bid, int n) //returns the number of times the sequence is in found in the tree
{
	if (n==1)
	{
		return (intptr_t)(daughter[*bid]); 
	}
	else
	{
		if (this->daughter[*bid] == NULL)
		{
			return 0; 
		}
		return daughter[*bid]->count(bid+1, n-1); 
	}
}


void  CLTreeS::addSequences(char *FsaFileName, int L, int maxSequenceLength, int addrevcompl, int numberOfCVPartitions, int selectPartitionNumber) // adds all the sequence from a FSA file 
{
	if (numberOfCVPartitions==0)
	{
		numberOfCVPartitions = 1; 
	}
	selectPartitionNumber = selectPartitionNumber % numberOfCVPartitions; 

	FILE *fi;

	fi = fopen(FsaFileName,"r"); 
		
	CSequence *s = new CSequence(maxSequenceLength+3);

	int counter = 0; 

	while (!feof(fi))
	{
		s->readFsa(fi); 
		if(s->getLength()>0)
		{
			counter++; 
			if (counter % numberOfCVPartitions==selectPartitionNumber)
			{ 
				this->addSequence(s->getSeqBaseId(), s->getLength(), L); 

				if (addrevcompl)
				{
					this->addSequence(s->getReverseComplement()->getSeqBaseId(), s->getLength(), L); 
				}
			}
		}
	}

	fclose(fi); 

	delete s; 
}

void CLTreeS::mismatchCount(int *bid, int n, int *cnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=(intptr_t)(daughter[i]);
		cnt++; 
		i=(i+1)&3;  //(i+1)%4
		*cnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*cnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*cnt+=(intptr_t)(daughter[i]);
	}
	else
	{
		n--; 
		bid++; 
		
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);

		cnt++; 
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt);

	}
}

void CLTreeS::mismatchCount(int *bid, int n, int *cnt, int maxmm) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*cnt+=(intptr_t)(daughter[i]);
		if (maxmm!=0)
		{
			cnt++; 
			i=(i+1)&3;  //(i+1)%4
			*cnt+=(intptr_t)(daughter[i]);
			i=(i+1)&3; 
			*cnt+=(intptr_t)(daughter[i]);
			i=(i+1)&3; 
			*cnt+=(intptr_t)(daughter[i]);
		}
	}
	else
	{
		if (maxmm!=0)
		{
			n--; 
			bid++; 
		
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);

			maxmm--;
			cnt++; 
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
			i=(i+1)&3; 
			if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,cnt,maxmm);
		}
		else
		{   // fast track for maxmm==0
			CLTreeS *cur = this;  
			while (--n)
			{
				if ((cur = cur->daughter[*bid++])==NULL) return;  
			}
			*cnt+=(intptr_t)(cur->daughter[*bid]); 
		}
	}
}

double CLTreeS::calcScore(int *bid, int L, double *kernel, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->mismatchCount(bid, L, tmpcnt); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}

double CLTreeS::calcScore(int *bid, int L, double *kernel, int maxmm, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	this->mismatchCount(bid, L, tmpcnt,maxmm); 
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}

double CLTreeS::calcScore(int *bid,int *bidrc, int L, int slen, double *kernel, int maxmm, int *tmpcnt)//calculates the score. tmpcnt is int[L+1]
{
	int i;
	for(i=0;i<=L;i++)
	{
		tmpcnt[i]=0;
	}
	slen = slen-L+1; 
	for(i=0;i<slen;i++)
	{
		this->mismatchCount(bid, L, tmpcnt,maxmm); 
		bid++; 
	}
	if (bidrc!=NULL)
	{
		for(i=0;i<slen;i++)
		{
			this->mismatchCount(bidrc, L, tmpcnt,maxmm); 
			bidrc++; 
		}
	}
	double res = 0; 
	for(i=0;i<=L;i++)
	{
		res+=kernel[i]*tmpcnt[i];
	}
	return res; 
}


void CLTreeS::addToList(class CLList *list, int n, int Lm1, int single, int *tmpbid) // used by LList. adds sequences to lis from tree 
{
	for (int i=0;i<4;i++)
	{
		if (daughter[i]!=NULL)
		{
			tmpbid[n] =i; 
			if (n==Lm1)
			{
				int frq= (intptr_t)(daughter[i]); 
				if ((frq==1)==single)
				{
					list->addSeq(tmpbid, frq); 
				}
			}
			else
			{
				daughter[i]->addToList(list, n+1, Lm1, single, tmpbid); 
			}
		}
	}
}
*/
