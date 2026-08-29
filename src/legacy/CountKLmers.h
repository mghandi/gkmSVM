/* CountKLmers.h
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
class CCountKLmers
{
public:
	CCountKLmers(int L, int K, int halfMem);
	~CCountKLmers(void);

	int **table; // keeps the count
	int K,L; 
	int halfMem; // set halfMem=1 to limit the first base of the motif to A,C
	int nrow,ncol; // nrow = C(L-1,K-1), ncol=RCSym? (4^K)/2 : 4^K
	int **w; //matrix for weights sum(wij*sj, 0<=j<L) gives the col index
	 

	void addSequence(int *seqBID, int size); 
	
	char *convertCol2KmerString(int col, char *sKmer); // returns K-mer for idx=col
	char *convertRow2KLmerString(int row, char *sKmer, char *sKLmer); // maps Kmer to KLmer for idx=row

private: 
	int *wdata; //allocated memory for w
	int fillwij(int pos, int nfilled, int idx, int *partial);
};

