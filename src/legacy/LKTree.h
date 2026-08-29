/* LKTree.h
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

#define GAPdaughter daughter[4]
#include "global.h"
class CLKTree
{
public:
	CLKTree * daughter[5]; 
//	CLKTree * gdaughter; //gap daughter  
	
	void deleteTree(int n); 
	void addSeq(int *bid, int n, int cnt); 
	void mismatchCount(int *bid, int n, intptr_t *mcnt, double *dcnt); // fills the mismatch count array (howmany sequences with m mismatches to a given sequence exist) 

	CLKTree(void);
	~CLKTree(void);
};

