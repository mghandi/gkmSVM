/* gkmCommonLib.cpp
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<iostream>

#include "global.h"

#include "LTreeS.h"
#include "GTree.h"
//#include "GTreeLeafData.h"

int gLM1; //L-1
int gMAXMM; //MaxMismatch
int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]
myFlt **gMMProfile0; //gMMProfile[0]

GTreeLeafData *gGTreeLeaves; // list of all the leaf nodes
GTreeLeafData2 *gGTreeLeaves2; // list of all the leaf nodes
int gGTreeLeavesCnt; // number of all leaf nodes

LTreeSnodeData ** gDFSlist[1000];
//CLTreeSptr **gDFSlistT[1000]; // without nonEmptyDaughterCnt
CLTreeS **gDFSlistT[1000]; // with nonEmptyDaughterCnt
int *gDFSMMlist[1000];
CbinMMtree **gDFSMMtree[1000]; // for the iDL bound 


int strlength(char *s)
{
	int i =0;
	while ((*s!=0)&&(*s!=10)&&(*s!=13)&&(*s!=EOF)) {i++; s++;}
	return i;
}




void randomPermute(double *x, int N)
{
	int i,j; 
	double h; 
	for(i=1; i<N; i++)
	{
		j = myrandom(i+1); 
		h = x[i]; 
		x[i]=x[j]; 
		x[j]=h; 
	}
}

void randomPermute(int *x, int N)
{
	int i,j; 
	int h; 
	for(i=1; i<N; i++)
	{
		j = myrandom(i+1); 
		h = x[i]; 
		x[i]=x[j]; 
		x[j]=h; 
	}
}


int Combinations(int n, int r)
{
	//if ((log(1.0*n)*r/log(2.0) > 64)&&(log(1.0*n)*(n-r)/log(2.0) > 64)) {printf("solution might be too big. Use dCombinations instead!");}

	if (r<0) return 0; 
	if (n<0) return Combinations(r-n-1, r)*((r%2==0)?1:-1); 
	if (n<r) return 0; 
    if ((n==0)&&(r==0)) return 1.0; 

	int i,j; 

	int *nn,*no, *h; 
	nn = new int[r+1]; 
	no = new int[r+1]; 

	for(i=0;i<=r;i++)
	{
		nn[i]=no[i]=0; 
	}
	nn[0]=no[0]=1;
	
	for(i=1;i<=n;i++)
	{
		h = no; no = nn; nn=h; 
		for(j=1;j<=r;j++)
		{
			nn[j] = no[j]+no[j-1]; 
		}
	}
	int res = nn[r]; 
	delete []nn; 
	delete []no; 

	return res; 
}


double dCombinations(int n, int r)
{
	if (r<0) return 0; 
	if (n<0) return dCombinations(r-n-1, r)*((r%2==0)?1:-1); 
	if (n<r) return 0; 
    if ((n==0)&&(r==0)) return 1.0; 

	int i,j; 

	double *nn,*no, *h; 
	nn = new double[r+1]; 
	no = new double[r+1]; 

	for(i=0;i<=r;i++)
	{
		nn[i]=no[i]=0; 
	}
	nn[0]=no[0]=1;
	
	for(i=1;i<=n;i++)
	{
		h = no; no = nn; nn=h; 
		for(j=1;j<=r;j++)
		{
			nn[j] = no[j]+no[j-1]; 
		}
	}
	double res = nn[r]; 
	delete []nn; 
	delete []no; 

	return res; 
}


int stringcompare(char *s1, char*s2, int maxlength)
{
	for (int i=0; i<maxlength; i++) 
	{
		if (s2[i]!=s1[i]) return 0; 
		if (((s1[i]==0)||(s1[i]==13))&&((s2[i]==0)||(s2[i]==13))) return 1; 
	}
	return 1; 
}


int convert2int(int *bid, int L)
{
	int s = 0; 
	for(int i=0;i<L;i++)
	{
		s<<=2; 
		s+=*bid; 
		bid++;
	}
	return s; 
}


int convertint2intRC(int x, int L)
{
	int s =0; 
	for(int i=0;i<L;i++)
	{
		s<<=2;
		s+=3-(x%4);
		x>>=2;
	}
	return s; 
}

char globtmpstr[10000]; // global temp string;
void Printf(const char *str){ // this to replace printf
    sprintf(globtmpstr, "%s", str);
    Printf(globtmpstr);
}

#define RPACKAGE
#ifdef RPACKAGE
  #include "R_ext/Print.h"
  void Printf(char *str){ // this to replace printf
      Rprintf("%s", str);
  }
  unsigned int locrseed = 123456789;
  int myrandom(int M) {
      locrseed = (1103515245 * locrseed + 12345) % (1<<31);
      return locrseed%M;
  }
#else
  void Printf(char *str){ // this to replace printf
      printf("%s", str);
  }
  int myrandom(int M) // generates uniform random integer between 0..M-1  (M<<2^31)
  {
      int b30r;
      b30r = rand()%1024;
      b30r = (b30r<<10)+(rand()%1024);
      b30r = (b30r<<11)+(rand()%2048);
      return b30r%M;
  }
#endif

