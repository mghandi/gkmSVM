/* KLmer.h
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
class CKLmer
{
public:
	CKLmer(int L, int K);
	int K,L; 
	char *seq;
	int *iseq; //A:1 C:2 G:4 T:8 .:15
	double score, w; // optional 

	void readKLmer(char *s);
	int countHits(char *s, int size);

	int commonKMerCnt(CKLmer *klmerj); // returns number of distinct KMers that both this KLmer and KLmerj cover. 
	double calcfbg(); 
	
	~CKLmer(void);
};

