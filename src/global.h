/* global.h
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

#pragma once

#ifndef __gkmsvmXC__global_h__
#define __gkmsvmXC__global_h__

//default values
#define MULTI_THREAD_SAFE  /*if remove (or comment) this line, it will become (about 2 times) faster, but may give approx results when using multithread */
#define MAX_ALPHABET_SIZE 4 /*for DNA set this to 4 (for faster speen and less memory usage); */
//#define MAX_ALPHABET_SIZE 24 /*for protein set this to 24 (or the maximu alphabet size); */
#define NBITS 2 /*ceiling log2 MAX_ALPHABET_SIZE */
#define DEF_L 10
#define DEF_K 6
#define DEF_D 3
#define DEF_MAXSEQLEN 10000
#define DEF_MAXNUMSEQ 1000000
#define DEF_TGKM 1
#define DEF_BATCHSIZE 100000
#ifndef MULTI_THREAD_SAFE
//    #define USE_GLOBAL /* may improve speed in some cases. but may make is slower in other cases.  */
#endif

#define TALK 0
//#define DEBUG 0 /*0 1*/
//#define TALK DEBUG /*1*/
//#define FAST_TRACK  // if enabled, it may make the algorithm faster for sparse trees.
#define PI 3.141593
#include <stdio.h>
#include <stdlib.h>

#include "Converter.h"
#include <iostream>
#include <math.h>

//	#include <tr1/unordered_map>  // this line for gcc
//  typedef std::tr1::unordered_map<int, double> Mymap;

#include <unordered_map>
typedef std::unordered_map<int, double> Mymap;

#ifdef MULTI_THREAD_SAFE
#include <atomic>
typedef std::atomic<int> aint;
#else 
typedef int aint;
#endif

int stringcompare(char *s1, char*s2, int maxlength) ; 
int strlength(char *s);
#define MYABS(x) (((x)<0)?-(x):x)

int Combinations(int n, int r);//
double dCombinations(int n, int r);//
int convert2int(int *bid, int L);
int convertint2intRC(int x, int L);

char *convertInt2Str(int col, char *str, int L); // returns L-mer for idx=col

extern char globtmpstr[]; // global temp string;
void Printf(char *str); // this to replace printf
void Printf(const char *str); // this to replace printf


int myrandom(int M);
void randomPermute(double *x, int N); 
void randomPermute(int *x, int N);
int find_str(char **strs, int N, char *str);

#define YSTMAXCHRPOS 1600000

const int YSTCHRSIZE[]={230208,  
                        813178,
                        316617,
                        1531919,
                        576869,	
                        270148,
                        1090947,
                        562643,	
                        439885,
                        745741,
                        666454,	
                        1078175,
                        924429,	
                        784334,
                        1091289,
                        948062}; 

//static CConverter globalConverter;
extern CConverter globalConverter;

#define freeMem(x) if(x!=NULL) delete []x

#define pi 3.14159265
#define sqr(x) ((x)*(x))
#define Epsilon 0.0000000000001
#define MAX_LINE_WIDTH 10000	/* maximum line width */

#define min(x,y) ((x<y)?x:y)
#define max(x,y) ((x>y)?x:y)
#define lcase(c) ((c>='a')?c:c-'A'+'a')
#define ucase(c) ((c>='a')?c-'a'+'A':c)



/*
struct Lmer{
  int seqID; 
int *baseID; 
};

union LPTr {
Lmer **all;
Lmer *one;
};
*/



union intintptr {
  int i;
  int *p;
};

struct LTreeSnodeData {
  int n;
  intintptr  seqIDs; //if n==1, it is int and contains the ID, otherwise it is int* and is the array of IDs;
  //  LPTr Lmers; //pointer to the starts of the sequences
  // int *baseID;
};


union LTreeSnodeDataptr {
  LTreeSnodeData *p;
  LTreeSnodeData **pp;
};


/*
 struct GTreeLeafData {
 int n;
 intintptr seqIDs_gbits; //if n==1, it is int and contains the ID, otherwise it is int* and is the array of IDs;
 //  LPTr Lmers; //pointer to the starts of the sequences
 int first_gbits; // gbits for the case n==1, otherwise, seqIDs and gapped_bits are both written in seqIDs_gbits (2 numbers for each L-mer)
 };
 */


#define myFlt double
union fintptr_t {
  myFlt f;
  class CLTreef *p;
  unsigned int i; 
};

#define NBITSONES ((1<<NBITS)-1)

#endif
