/* mainGkmKernel.cpp : gkmKernel program
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

#include <iostream>       // std::cout
#include <thread>         // std::thread

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include<iostream>

#include "global.h"
#include "globalvar.h"
#include "gkmOptions.h"
#include "MultiTrack.h"
#include "SequenceSet.h"
#include "KernelFile.h"
#include "gkmMainHelpers.h"
#include <vector>

#include "Sequence.h"
#include "CalcWmML.h"
#include "SequenceNames.h"
#include "LTree.h"
#include "LTreef.h"
#include "LTreeS.h"
#include "LList.h"
#include "CiDLPasses.h"

using namespace std;

//#include "stdafx.h"

int gkmKernelSimple(OptsGkmKernel &opt, const CConverter &conv);
// the tree algorithm is instantiated per alphabet size (trie_b4.cpp / trie_b32.cpp, see global.h)
namespace gkm_b4  { int gkmKernelSuffixTree(OptsGkmKernel &opt, const CConverter &conv); int gkmKernelMultiTrack(OptsGkmKernel &opt, const TrackAlphabets &ta); }
namespace gkm_b32 { int gkmKernelSuffixTree(OptsGkmKernel &opt, const CConverter &conv); int gkmKernelMultiTrack(OptsGkmKernel &opt, const TrackAlphabets &ta); }
static int gkmKernelSuffixTree(OptsGkmKernel &opt, const CConverter &conv)
{
	if (conv.b <= 4) return gkm_b4::gkmKernelSuffixTree(opt, conv);
	return gkm_b32::gkmKernelSuffixTree(opt, conv);
}

void print_usage_and_exit_gkmKernel(const char *prog)
{
    Printf("\n");
    gkmMsg(" Usage: %s [options] <pos_seqfile> <neg_seqfile> <outfile>\n",prog );
    Printf("\n");
	Printf("  generates a lower triangle of kernel matrix (i.e. pairwise similarities)\n");
	Printf("  between the sequences.\n");
    Printf("\n");
    Printf(" Arguments:\n");
	Printf("  pos_seqfile: positive sequence file name (fasta format)\n");
	Printf("  neg_seqfile: negative sequence file name (fasta format)\n");
	Printf("  outfile: output file name\n");
    Printf(" \n");
    Printf(" Options:\n");
    gkmMsg("  -l L           set word length, default= %d\n", DEF_L);
	gkmMsg("  -k K           set number of informative columns, default= %d \n",DEF_K);
	gkmMsg("  -d maxMismatch set maximum number of mismatches to consider, default= %d \n", DEF_D);
    gkmMsg("%s","  -m maxSeqLen   set maximum sequence length in the sequence files,\n");
	gkmMsg("                 default= %d \n", DEF_MAXSEQLEN);
	gkmMsg("%s", "  -n maxNumSeq   set maximum number of sequences in the sequence files,\n");
	gkmMsg("                 default= %d\n", DEF_MAXNUMSEQ);
	gkmMsg("%s", "  -t filterType  set filter type: 0(use full filter), 1(use truncated filter:\n" );
	gkmMsg("%s", "                 this gaurantees non-negative counts for all L-mers), 2(use h[m],\n");
	gkmMsg("                 gkm count vector), 3(wildcard), 4(mismatch), default=%d\n",DEF_TGKM );
	gkmMsg("%s", "  -a algorithm   set algorithm type: 0(auto), 1(XOR Hashtable), 2(tree),\n");
	gkmMsg("%s", "                 default=0\n");
	gkmMsg("%s", "  -R             if set, reverse complement sequences will NOT be considered\n");
	gkmMsg("%s", "  -p             if set, a constant to will be added to the count estimates\n");
	gkmMsg("%s", "  -M             max mismatch for Mismatch kernel or wildcard kernel, default=2\n");
	gkmMsg("%s", "  -L             lambda for wildcard kernel, defaul=1.0\n");
    gkmMsg("%s", "  -A             alphabet: dna (default), rna, protein, an alphabet file (one symbol per line,\n");
    gkmMsg("%s", "                 optional complement), or =SYMBOLS (literal, e.g. =01). A comma-separated list\n");
    gkmMsg("%s", "                 gives one alphabet per track of a multi-track FASTA input (header line\n");
    gkmMsg("%s", "                 followed by one line per track), e.g. -A dna,=01 for DNA + methylation flags\n");
    gkmMsg("%s", "  -T             maximum number of threads, defaul=2*l\n");
    gkmMsg("%s", "  -N             merge records that share a name within a file into one row (default: one record = one row)\n");
    gkmMsg("%s", "  -b             write the kernel in the binary .gkmk format (float32) instead of text\n");
    gkmMsg("%s", "  -r rows        compute the kernel in tiles of this many rows (bounded memory; 0 = automatic, 1 GB budget)\n");
    gkmMsg("%s", "  -P design      pass design: 0 automatic (default), 1 greedy iDL table, 2 prefix-split (long words; same result)\n");
    
    Printf(" \n");
}

// Parse argv into opt. Returns 0 on success, 1 if usage was printed (nothing to run).
static int gkmKernelParseArgs(int argc, char** argv, OptsGkmKernel &opt)
{
    ::optind=1; // reset getopt()
    int c;
    if (argc == 1) { print_usage_and_exit_gkmKernel(argv[0]); return 1;}

	while ((c = getopt (argc, argv, "l:k:d:m:n:t:a:L:M:A:T:r:P:RpbN")) != -1)
	{
		switch (c) 
		{
			case 'l': opt.L = atoi(optarg); break;
			case 'k': opt.K = atoi(optarg); break;
			case 'd': opt.maxnmm = atoi(optarg); break;
			case 'm': opt.maxseqlen = atoi(optarg); break;
			case 'n': opt.maxnumseq = atoi(optarg); break;
			case 't': opt.useTgkm = atoi(optarg); break;
			case 'a': opt.alg = atoi(optarg); break;
			case 'R': opt.addRC = false; break;
			case 'p': opt.usePseudocnt = true; break;
			case 'b': opt.OutputBinary = true; break;
			case 'M': opt.wildcardMismatchM = atoi(optarg); break;
			case 'L': opt.wildcardLambda = atof(optarg); break;
			case 'A': opt.alphabetFN = optarg; break;
			case 'T': opt.maxnThread = atoi(optarg); break;
			case 'N': opt.mergeByName = true; break;
			case 'r': opt.tileRows = atoi(optarg); break;
			case 'P': opt.passDesign = atoi(optarg); break;
			default: print_usage_and_exit_gkmKernel(argv[0]); return 1;
		}
	}

	if (argc-optind != 3) { print_usage_and_exit_gkmKernel(argv[0]); return 1;}

	int index = optind;
	opt.posfile = argv[index++];
	opt.negfile = argv[index++];
	opt.outfile = argv[index++];
	return 0;
}

int mainGkmKernel(int argc, char** argv) //mainGkmKernel
{
	OptsGkmKernel opt;
	if (gkmKernelParseArgs(argc, argv, opt) != 0) return 1;
	return gkmKernelRun(opt);
}

// Validate, load the alphabet, dispatch. Shared by the CLI and the R binding.
int gkmKernelRun(OptsGkmKernel &opt)
{
	//check parameters
	if (opt.L < 1 || opt.K < 1)
	{
		Printf("\n ERROR: L and K must be positive.\n"); return 1;
	}

	if (opt.useTgkm < 0 || opt.useTgkm > 4)
	{
		Printf("\n ERROR: filter type (-t) must be between 0 and 4.\n"); return 1;
	}
	if (opt.alg < 0 || opt.alg > 2)
	{
		Printf("\n ERROR: algorithm (-a) must be 0, 1 or 2.\n"); return 1;
	}
	if (opt.maxseqlen < opt.L || opt.maxnumseq < 1)
	{
		Printf("\n ERROR: maxSeqLen must be at least L and maxNumSeq at least 1.\n"); return 1;
	}

	TrackAlphabets ta; // the alphabet(s) of this call (DNA unless -A); Phase 7: one per track
	if (ta.parse(opt.alphabetFN, GKM_MAX_ALPHABET) != 0) return 1;
	// K counts informative positions of the whole word: L for one track, T*L for multi-track input (Phase 7)
	if ((opt.K > ta.T() * opt.L) && (opt.useTgkm < 3))
	{
		if (ta.T() == 1) Printf("\n ERROR: K must be less than or equal to L!\n");
		else Printf("\n ERROR: K must be less than or equal to the number of tracks times L!\n");
		return 1;
	}
	if ((opt.maxnmm > 0) && (ta.T() * opt.L < opt.maxnmm))
	{
		if (ta.T() == 1) Printf("\n ERROR: maxMismatch must be less than or equal to L!\n");
		else Printf("\n ERROR: maxMismatch must be less than or equal to the number of tracks times L!\n");
		return 1;
	}
	if (ta.T() > 1) {
		if (opt.useTgkm > 2) { Printf("\n ERROR: filter types 3 and 4 (wildcard, mismatch kernels) are only available for a single alphabet.\n"); return 1; }
		if (opt.usePseudocnt) { Printf("\n ERROR: pseudocounts (-p) are only available for a single alphabet.\n"); return 1; }
		if (opt.alg == 1) Printf("\nAlgorithm is set to 2 (Tree) for multi-track input.\n");
		if (opt.addRC && !ta.conv[0]->hasComplement) { opt.addRC = false; Printf("\nAdd Reverse Complement option is turned off (the alphabet of track 1 declares no complement pairs).\n"); }
		if (ta.bmax() <= 4) return gkm_b4::gkmKernelMultiTrack(opt, ta);
		return gkm_b32::gkmKernelMultiTrack(opt, ta);
	}
	const CConverter &conv = *ta.conv[0];
	int alg = opt.alg;
	if (!opt.alphabetFN.empty()){
		if (opt.addRC&&!conv.hasComplement){
			opt.addRC=false;
			Printf("\nAdd Reverse Complement option is turned off (the alphabet declares no complement pairs).\n");
		}
		if ((alg!=2)&&(conv.b!=4)){
			alg=2;
			Printf("\nAlgorithm is set to 2 (Tree) to support alphabet size other than 4.\n");
		}
	}

	switch (alg) 
	{
		case 0:
			if ((opt.L-opt.K <= 4) || (opt.maxnmm >= 0 && opt.maxnmm <= 4))
			{
				return gkmKernelSuffixTree(opt, conv);
			}
			else
			{
				return gkmKernelSimple(opt, conv);
			}
		case 1:
			return gkmKernelSimple(opt, conv);
		default:
			return gkmKernelSuffixTree(opt, conv);
	}
}

int gkmKernelSimple(OptsGkmKernel &opt, const CConverter &conv)  //Use XOR precomputed hash table
{

	int L = opt.L;
	int K = opt.K;
	int maxseqlen =	opt.maxseqlen;
	int useTgkm = opt.useTgkm;
	int maxnmm = opt.maxnmm; //auto 
	int nMAXSEQUENCES = opt.maxnumseq;
	bool addRC = opt.addRC;
	//bool usePseudocnt= opt.usePseudocnt; 

	const char *posSeqsFN = opt.posfile.c_str();
	const char *negSeqsFN = opt.negfile.c_str();
	const char *outFN = opt.outfile.c_str();

	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];
	std::vector<SeqRecord> records; // row identity, written to <outfile>.index

	int i; 
	CSequence *sgi= new CSequence(maxseqlen+3, &conv);

	CCalcWmML wmc(L, K, conv.b);

	//for curiosity
	/*
	int jj;
	cout << "wmc.n0: " << wmc.n0 << endl;
	cout << "wmc.wm:" << endl;
	for (jj=0;jj<=K;jj++) {
		cout << jj << ": " << wmc.wm[jj] << endl;
	}
	cout << "jj\tkernel\tTruncated\tc\tcTr\th" << endl;
	for (jj=0;jj<=L;jj++) {
		cout << jj << "\t" << wmc.kernel[jj];
		cout << "\t" << wmc.kernelTruncated[jj];
		cout << "\t" << wmc.c[jj];
		cout << "\t" << wmc.cTr[jj];
		cout << "\t" << wmc.h[jj] << endl;
	}
	*/

	//double *kernel = wmc.kernelTruncated;
	if (maxnmm==-1)
	{ 
		maxnmm=L;
		if (useTgkm==1)
		{
			maxnmm = 2*(wmc.kernelTruncatedLength-1);
		}
		if (useTgkm==2)  //gapped kmer kernel 
		{
			maxnmm = L-K; 
		}
		if (useTgkm==3)  //wildcard kernel 
        {
            maxnmm = opt.wildcardMismatchM;
        }
		if (useTgkm==4)  //mismatch kernel 
        {
            maxnmm = 2*opt.wildcardMismatchM;
        }
        if (maxnmm>L)
        {
            maxnmm=L;
        }
	}
	double n0 = wmc.n0; 
	double *c = wmc.cTr; 
	
	n0 = c[maxnmm]/2; 

	if (useTgkm==0)
	{
		n0 = 0; 
		//kernel = wmc.kernel; 
		c = wmc.c; // same as kernel
	}
	if (useTgkm==2)
	{
	  //	n0 = 0; 
	  //	kernel = wmc.kernel; 
		c = wmc.h; 
		n0 = c[maxnmm]/2;

	}
    if (useTgkm==3)  //wildcard kernel 
    {
        c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, conv.b, opt.wildcardLambda, c);
    	n0 = c[maxnmm]/2;

    }
    if (useTgkm==4)  //mismatch kernel 
    {
        c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, conv.b, c);
    	n0 = c[maxnmm]/2;

    }

	gkmMsg("\n maximumMismatch = %d\n", maxnmm);
	for(int ii=0;ii<=maxnmm;ii++) {
		gkmMsg("\n c[%d] = %e",ii,c[ii] );
	}
	Printf("\n");

	int npos=0; 
	int nneg=0; 

	//char *tmps = new char[maxseqlen+3]; 
	int *mmcnt = new int[L+1];  //mismatch count

	CLList psetL(L,2*maxseqlen+5/* psetT->leavesCount(0,L)+1*/); // keeps all the sequences of length L
	psetL.UseLookupTable =0;
	int *hdist = psetL.HamDist; 
	int nseqs = 0; 

	//read positive sequence file
	FILE *sfi = fopen(posSeqsFN, "r"); 
	if (sfi == NULL) return gkmCannotOpen(posSeqsFN);
	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) return 1; 
		if(sgi->getLength()>0)
		{
			if (nseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
			seqsL[nseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			CLTree *psetT = new CLTree();// keeps all the sequences of length L
			psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L); 
			if(addRC)
			{
				psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L); 
			}			
			seqsL[nseqs]->addFromLTree(psetT); 
			psetT->deleteTree(L); 
			delete psetT; 
			records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), (long)(sgi->getLength()-L+1 > 0 ? sgi->getLength()-L+1 : 0) * (addRC ? 2 : 1)});
			nseqs++; 
		}
	}
	fclose(sfi); 

	npos = nseqs;

	//read negative sequence file
	sfi = fopen(negSeqsFN, "r"); 
	if (sfi == NULL) return gkmCannotOpen(negSeqsFN);
	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) return 1; 
		if(sgi->getLength()>0)
		{
			if (nseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
			seqsL[nseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			CLTree *psetT = new CLTree();// keeps all the sequences of length L
			psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L); 
			if(addRC)
			{
				psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L); 
			}			
			seqsL[nseqs]->addFromLTree(psetT); 
			psetT->deleteTree(L); 
			delete psetT; 
			records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), (long)(sgi->getLength()-L+1 > 0 ? sgi->getLength()-L+1 : 0) * (addRC ? 2 : 1)});
			nseqs++; 
		}
	}
	fclose(sfi); 

	nneg = nseqs - npos;
	for(i=0;i<nseqs;i++) { records[i].id = seqRecordId(i, npos); records[i].label = (i<npos) ? 1 : -1; }
	(void)nneg; (void)n0; // kept for symmetry with the tree path (pseudo-counts are not applied by the XOR algorithm)

	for(i=0;i<nseqs;i++)
	{
		norm[i] = sqrt(seqsL[i]->calcInnerProd(seqsL[i],c,mmcnt));
	}

	GkmkWriter bin;
	FILE *fo = NULL;
	if (opt.OutputBinary) { gkmFillGkmkHeader(bin, opt, maxnmm, nseqs, npos, conv); bin.values.reserve((size_t)nseqs*(nseqs+1)/2); }
	else { fo = fopen(outFN, "w"); if (fo == NULL) return gkmCannotOpen(outFN); }

	for(i=0;i<nseqs;i++)
	{
		for(int j=0;j<=i;j++)
		{
			double v = (i==j) ? 1.0 : ((norm[i]*norm[j]<1E-50)?0.0:seqsL[i]->calcInnerProd(seqsL[j],c,mmcnt)/(norm[i]*norm[j]));
			if (v == 0.0) v = 0.0; // canonical +0
			if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
			else bin.add(gkmCanon(v));
		}
		if (fo) fprintf(fo, "\n"); 
	}

	if (fo) fclose(fo); 
	else if (bin.write(opt.outfile, records) != 0) return gkmCannotOpen(outFN);
	if (writeIndexSidecar(opt.outfile, records) != 0) { gkmMsg("\n WARNING: could not write %s.index\n", outFN); }
	//delete []tmps; 
	delete []mmcnt;
	for(i=0;i<nseqs;i++)
	{
		delete seqsL[i];
	}
	delete []seqsL; 
	delete []norm; 
	delete sgi; 

	return 0; 
}


