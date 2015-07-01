/* Sequence.h: interface for the CSequence class. 
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

#if !defined(AFX_SEQUENCE_H__49B2B093_380F_4384_AA74_60169F38FFD7__INCLUDED_)
#define AFX_SEQUENCE_H__49B2B093_380F_4384_AA74_60169F38FFD7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "global.h"
//#include "BendingData.h"
//#include "OccupancyData.h"

class CSequence  
{

private: 

	char *seq, *subseq; 
	char *seqName; 
	char *seqLabel; 
	int length;
	int maxLength; 
	double weight; 
	char *NameLink;
	dinuclId *dinucl; 
	baseId *seqBaseId; // 0123 instead of acgt
	CSequence *reverseComplement; 
//	char *seqName; //sequence name or ID (optional) // upto 99 chars

public:
	void mutateOneBase(int pos, baseId nwbid);
	char *getSubseq();
	CSequence *getReverseComplement(); 

	char *getSubseq(int p1,int p2);
	int *getSubseqBaseId(int p1,int p2, int *obid);

//	class CBendingData *allocBendingData();
//	class COccupancyData *allocOccupancyData();


	char *getSeq(); 
	int getLength(); 
	dinuclId *getDinucl(); 
	baseId *getSeqBaseId(); 
	void setName(char *newname); 
	char *getName(); 
	char *getLabel(); 
	double getWeight(); 
	void setWeight(double w); 
	char *getNameLink(); 
	void setNameLink(char *sp);  



//	class CBendingData *bendingData; 
//	class COccupancyData *occupancyData; 


	int readFsa(FILE *f,int SkipAlphabetCheck=false);  // reads one sequence from already opened file f;
	void writeFsa(FILE *f); 	 
	int readBasic(FILE *f);  // reads one sequence from already opened file f; 
	void writeBasic(FILE *f); 	 
	int readString(char *s);  // reads sequence from a string (coverts chat * to sequence)

	CSequence(int maxLength, CSequence *sCopyFrom = NULL );  // makes a replicate of s;
	
	virtual ~CSequence();

};

#endif // !defined(AFX_SEQUENCE_H__49B2B093_380F_4384_AA74_60169F38FFD7__INCLUDED_)
