/* LKTree.cpp
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

#include "LKTree.h"


CLKTree::CLKTree(void)
{
	GAPdaughter = daughter[0]=daughter[1]=daughter[2]=daughter[3]=NULL;
}

CLKTree::~CLKTree(void)
{
}

void CLKTree::addSeq(int *bid, int n, int cnt)
{//printf("\n %d %d %d",*bid, n,cnt);
	if (n==1)
	{
		this->daughter[*bid] = (CLKTree *)((intptr_t)(daughter[*bid])+cnt); 
	}
	else
	{
		if (this->daughter[*bid] == NULL)
		{
			this->daughter[*bid] = new CLKTree(); 
		}
		daughter[*bid]->addSeq(bid+1, n-1, cnt); 
	}
}

void CLKTree::deleteTree(int n)
{
	if (n>1)
	{
		for (int i=0;i<5;i++)
		{
			if (daughter[i]!=NULL)
			{
				daughter[i]->deleteTree(n-1);
				delete daughter[i];
			}
		}
	}
}

void CLKTree::mismatchCount(int *bid, int n, intptr_t *mcnt, double *dcnt) // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 
{
	int i = *bid; 
	if (n==1)
	{
		*mcnt+=(intptr_t)(GAPdaughter);
		*mcnt+=(intptr_t)(daughter[i]);
		if (*mcnt>10000000) {*mcnt-=10000000; *dcnt+=10000000.0; }
		mcnt++; dcnt++;
		i=(i+1)&3;  //(i+1)%4
		*mcnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*mcnt+=(intptr_t)(daughter[i]);
		i=(i+1)&3; 
		*mcnt+=(intptr_t)(daughter[i]);
		if (*mcnt>10000000) {*mcnt-=10000000; *dcnt+=10000000.0; }
	}
	else
	{
		n--; 
		bid++; 
		// assuming the tree is full (so I don't check whether each daughter is null or not) // this assumption would speed up a litle bit, however, I don't make that to make it feasible for partial KLmers like G3 !
		if (GAPdaughter!=NULL) GAPdaughter->mismatchCount(bid,n,mcnt,dcnt);
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,mcnt,dcnt);

		mcnt++; dcnt++;
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,mcnt, dcnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,mcnt, dcnt);
		i=(i+1)&3; 
		if (daughter[i]!=NULL) daughter[i]->mismatchCount(bid,n,mcnt, dcnt);

	}
}

