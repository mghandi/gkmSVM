/* Converter.h: interface for the CConverter class.
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

//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONVERTER_H__109528B8_5F1B_4654_82F8_B8A43B5E9210__INCLUDED_)
#define AFX_CONVERTER_H__109528B8_5F1B_4654_82F8_B8A43B5E9210__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

//#include "global.h"

//typedef int contextId; 
//typedef int baseId; // {0=a,1=c,2=g,3=t} 
//typedef int dinuclId; // {0=aa,1=ca,..,15=tt} 

#include "globalTypeDefs.h"

class CConverter  
{
public:

	CConverter();
//	CConverter(char *alphabet=NULL,int alphabetSize=0);
	void init();

	virtual ~CConverter();

	int cidx[256];
	char *icidx; // 'ACGT'= 0123
	char *icidxL;  // 'acgt'= 0123
	char bcompl[256]; 
	baseId bidcompl[256];

	
	dinuclId dnidx(char *dn); 

	void convertBasetoDinucl(baseId x[], dinuclId y[], int N);// x[0..N], y[0..N-1]
	void convertBasetoDinucl(char x[], dinuclId y[], int N); // x[0..N], y[0..N-1]

	int isInAlphabet[256];
	int isACGT[256];

	char alphabet[256];
	int b; //alphabetSize;
	void readAlphabetFile(char *FN, int MAX_ALPHABET_SIZE_copy);

};

#endif // !defined(AFX_CONVERTER_H__109528B8_5F1B_4654_82F8_B8A43B5E9210__INCLUDED_)
