/* Converter.cpp: implementation of the CConverter class.
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

//#include "stdafx.h"
#include "Converter.h"
#include "global.h"
//#include <stdio.h>
//#include <ctype.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CConverter::CConverter()
{
	b=4;
	alphabet[0]='A';
	alphabet[1]='C';
	alphabet[2]='G';
	alphabet[3]='T';

	init();

}

void CConverter::init(){

	//char ci;
	int ici;
	//for (ci='A';ci<'z';ci++) cidx[ci]=0;
	for (ici=0;ici<256;ici++) cidx[ici]=0;

	//	cidx['a'] = 0; 		cidx['c'] = 1; 		cidx['g'] = 2; 		cidx['t'] = 3;
	//	cidx['A'] = 0; 		cidx['C'] = 1; 		cidx['G'] = 2; 		cidx['T'] = 3;
	for(ici=0;ici<b;ici++){
		cidx[toupper(alphabet[ici])]=ici;
		cidx[tolower(alphabet[ici])]=ici;
	}

	icidx = new char[b];
	icidxL = new char[b];

	//	icidx[0]  = 'A';		icidx[1]  = 'C';		icidx[2]  = 'G';		icidx[3]  = 'T';
	//	icidxL[0] = 'a';		icidxL[1] = 'c';		icidxL[2] = 'g';		icidxL[3] = 't';
	for(ici=0;ici<b;ici++){
		icidx[ici]=toupper(alphabet[ici]);
		icidxL[ici]=tolower(alphabet[ici]);

		bidcompl[ici]=b-ici-1; //only good for DNA
		if (b==16){
			// special case for dinucleotides
			bidcompl[ici]= ((3-(ici&3))<<2) + (3-((ici&12)>>2)) ; //only good for DNA
		}

	}

	//	bidcompl[0] = 3;
	//	bidcompl[1] = 2;
	//	bidcompl[2] = 1;
	//	bidcompl[3] = 0;

	for (ici=0;ici<256;ici++) //(ci='A';ci<'z';ci++)
	{
		//bcompl[ci] = icidx[bidcompl[cidx[ci]]];
		bcompl[ici] = icidx[bidcompl[cidx[ici]]];
	}

	for (ici=0;ici<256;ici++)
	{
		isACGT[ici] = 0;
		isInAlphabet[ici] = 0;
	}
	isACGT['a'] = 1; 		isACGT['c'] = 1; 		isACGT['g'] = 1; 		isACGT['t'] = 1;
	isACGT['A'] = 1; 		isACGT['C'] = 1; 		isACGT['G'] = 1; 		isACGT['T'] = 1;

	for(ici=0;ici<b;ici++){
		isInAlphabet[toupper(alphabet[ici])]=1;
		isInAlphabet[tolower(alphabet[ici])]=1;
	}
	

}
CConverter::~CConverter()
{
	delete []icidx; 
	delete []icidxL; 

}


/*

aa	0
ca	1
ga	2
ta	3
ac	4
cc	5
gc	6
tc	7
ag	8
cg	9
gg	10
tg	11
at	12
ct	13
gt	14
tt	15
*/


dinuclId CConverter::dnidx(char *dn)
{
	return cidx[dn[0]]+b*cidx[dn[1]];
}




void CConverter::convertBasetoDinucl(baseId x[], dinuclId y[], int N) // x[0..N], y[0..N-1]
{
	for(int i=0;i<N;i++)
		y[i] = x[i]+b*x[i+1];
}

void CConverter::convertBasetoDinucl(char x[], dinuclId y[], int N) // x[0..N], y[0..N-1]
{
	for(int i=0;i<N;i++)
		y[i] = cidx[x[i]]+b*cidx[x[i+1]];
}


void CConverter::readAlphabetFile(char *FN, int MAX_ALPHABET_SIZE_copy){
	FILE *f= fopen(FN,"r");
	static char sline[1000+3];

	b=0;
	fgets(sline, 1000, f);
	while(!feof(f)){
		alphabet[b++]=sline[0];
		fgets(sline, 1000, f);
	}
    sprintf(globtmpstr,"Alphabet Size = %d\n",b);Printf(globtmpstr);
	if(b>MAX_ALPHABET_SIZE_copy){
		Printf("ERROR: alphabet size greater than #MAX_ALPHABET_SIZE. Redefine #MAX_ALPHABET_SIZE in global.h\n \n");
		return;//exit(-1);
	}

	delete []icidx;
	delete []icidxL;
	init();
}

