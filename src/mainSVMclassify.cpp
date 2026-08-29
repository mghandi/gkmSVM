/* mainSVMclassify.cpp
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
#include<unistd.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<iostream>

#include "global.h"
#include "globalvar.h"
#include "gkmOptions.h"
#include "gkmMainHelpers.h"

#include "Sequence.h"
#include "CalcWmML.h"
#include "SequenceNames.h"
#include "LTree.h"
#include "LTreef.h"
#include "LTreeS.h"
#include "LList.h"

using namespace std;

int svmClassifySimple(OptsSVMClassify &opt);
// the tree algorithm is instantiated per alphabet size (trie_b4.cpp / trie_b32.cpp, see global.h)
namespace gkm_b4  { int svmClassifySuffixTree(OptsSVMClassify &opt); }
namespace gkm_b32 { int svmClassifySuffixTree(OptsSVMClassify &opt); }
static int svmClassifySuffixTree(OptsSVMClassify &opt)
{
	if (::globalConverter.b <= 4) return gkm_b4::svmClassifySuffixTree(opt);
	return gkm_b32::svmClassifySuffixTree(opt);
}

void print_usage_and_exit(const char *prog)
{
    Printf("\n");
    snprintf(globtmpstr, GKM_TMPSTR_LEN, " Usage: %s [options] <test_seqfile> <sv_seqfile> <sv_alphafile> <outfile>\n",prog );Printf(globtmpstr);
    Printf("\n");
	Printf("  given support vectors SVs and corresponding coefficients alphas and a set of \n");
	Printf("  sequences, calculates the SVM scores for the sequences.\n");
    Printf("\n");
	Printf(" Arguments:\n");
	Printf("  test_seqfile: sequence file name to test/score (fasta format)\n");
	Printf("  sv_seqfile: sequence file name containing support vectors (fasta format)\n");
	Printf("  sv_alphafile: coefficient file name containing alphas for support vectors\n");
	Printf("  outfile: output file name\n");
    Printf("\n");
    Printf(" Options:\n");
	snprintf(globtmpstr, GKM_TMPSTR_LEN,"  -l L           set word length, default=%d\n",DEF_L); Printf(globtmpstr);
	snprintf(globtmpstr, GKM_TMPSTR_LEN,"  -k K           set number of informative columns, default=%d\n", DEF_K); Printf(globtmpstr);
    snprintf(globtmpstr, GKM_TMPSTR_LEN, "  -d maxMismatch set maximum number of mismatches to consider, default=%d\n",DEF_D); Printf(globtmpstr);
	Printf("  -m maxSeqLen   set maximum sequence length in the sequence files,\n");
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"                 default=%d\n",DEF_MAXSEQLEN); Printf(globtmpstr);
	Printf("  -n maxNumSeq   set maximum number of sequences in the sequence files,\n");
    snprintf(globtmpstr, GKM_TMPSTR_LEN, "                 default=%d\n",DEF_MAXNUMSEQ); Printf(globtmpstr);
	Printf("  -t filterType  set filter type: 0(use full filter), 1(use truncated filter:\n");
	Printf("                 this gaurantees non-negative counts for all L-mers), 2(use h[m],\n");
	snprintf(globtmpstr, GKM_TMPSTR_LEN, "                 gkm count vector), 3(wildcard), 4(mismatch), default=%d\n",DEF_TGKM); Printf(globtmpstr);
	Printf("  -a algorithm   set algorithm type: 0(auto), 1(XOR Hashtable), 2(tree),\n");
	Printf("                 default=0\n");
	Printf("  -b             set number of sequences to compute scores for in batch, \n");
	snprintf(globtmpstr, GKM_TMPSTR_LEN, "                 default=%d\n", DEF_BATCHSIZE); Printf(globtmpstr);
	Printf("  -R             if set, reverse complement sequences will NOT be considered\n");
	Printf("  -p             if set, a constant to count estimates will be added\n");
	Printf("  -M             max mismatch for Mismatch kernel or wildcard kernel, default=2\n");
	Printf("  -L             lambda for wildcard kernel, defaul=0.9\n");
	Printf("  -A             alphabets file name, if not specified, it is assumed the inputs are DNA sequences\n");
    Printf("\n");
}

static int svmClassifyParseArgs(int argc, char** argv, OptsSVMClassify &opt)
{
    ::optind=1; // reset getopt()
	int c;
    if (argc == 1) { print_usage_and_exit(argv[0]); return 1; }

	while ((c = getopt (argc, argv, "l:k:d:m:n:t:a:b:M:L:A:Rp")) != -1)
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
			case 'b': opt.batchSize= atoi(optarg); break;
			case 'R': opt.addRC = false; break;
			case 'p': opt.usePseudocnt = true; break;
			case 'M': opt.wildcardMismatchM = atoi(optarg); break;
			case 'L': opt.wildcardLambda = atof(optarg); break;
			case 'A': opt.alphabetFN = optarg; break;
			default: print_usage_and_exit(argv[0]); return 1;
		}
	}

	if (argc-optind != 4) { print_usage_and_exit(argv[0]); return 1; }

	int index = optind;
	opt.seqfile = argv[index++];
	opt.svseqfile = argv[index++];
	opt.alphafile = argv[index++];
	opt.outfile = argv[index++];
	return 0;
}

int mainSVMclassify(int argc, char** argv) // mainSVMclassify
{
	OptsSVMClassify opt;
	if (svmClassifyParseArgs(argc, argv, opt) != 0) return 1;
	return svmClassifyRun(opt);
}

int svmClassifyRun(OptsSVMClassify &opt)
{
	if (opt.L < 1 || opt.K < 1) { Printf("\n ERROR: L and K must be positive.\n"); return 1; }
	if ((opt.K > opt.L) && (opt.useTgkm < 3)) { Printf("\n ERROR: K must be less than or equal to L!\n"); return 1; }
	if ((opt.maxnmm > 0) && (opt.L < opt.maxnmm)) { Printf("\n ERROR: maxMismatch must be less than or equal to L!\n"); return 1; }
	if (opt.useTgkm < 0 || opt.useTgkm > 4) { Printf("\n ERROR: filter type (-t) must be between 0 and 4.\n"); return 1; }
	if (opt.alg < 0 || opt.alg > 2) { Printf("\n ERROR: algorithm (-a) must be 0, 1 or 2.\n"); return 1; }
	if (opt.batchSize < 1) { Printf("\n ERROR: batch size must be at least 1.\n"); return 1; }
	if (opt.maxseqlen < opt.L || opt.maxnumseq < 1) { Printf("\n ERROR: maxSeqLen must be at least L and maxNumSeq at least 1.\n"); return 1; }

	::globalConverter.resetToDNA(); // the alphabet is a global that outlives this call
	int alg = opt.alg;
	if (!opt.alphabetFN.empty()){
		if (::globalConverter.readAlphabetFile(opt.alphabetFN.c_str(), GKM_MAX_ALPHABET)!=0) return 1;
		if (opt.addRC&&!::globalConverter.hasComplement){
			opt.addRC=false;
			Printf("\nAdd Reverse Complement option is turned off (the alphabet declares no complement pairs).\n");
		}
		if ((alg!=2)&&(::globalConverter.b!=4)){
			alg=2;
			Printf("\nAlgorithm is set to 2 (Tree) to support alphabet size other than 4.\n");
		}
	}

	switch (alg) 
	{
		case 0:
			if (opt.L<=10)
			{
				return svmClassifySuffixTree(opt);
			}
			else
			{
				if ((opt.L-opt.K <= 2) || (opt.maxnmm >= 0 && opt.maxnmm <= 2))
				{
					return svmClassifySuffixTree(opt);
				}
				else
				{
					return svmClassifySimple(opt);
				}
			}
		case 1:
			return svmClassifySimple(opt);
		default:
			return svmClassifySuffixTree(opt);
	}
}

int svmClassifySimple(OptsSVMClassify &opt)
{
	int i;
	int L = opt.L; 
	int K = opt.K; 
	int maxnmm = opt.maxnmm; //auto 
	int useTgkm = opt.useTgkm; 
	int maxseqlen =	opt.maxseqlen; 
	int nMAXSEQUENCES = opt.maxnumseq; 
	bool addRC = opt.addRC; 
	//int batchSize = opt.batchSize; //batch size 

	const char *SeqsFN = opt.seqfile.c_str();
	const char *SVSeqsFN = opt.svseqfile.c_str();
	const char *SVSeqIDsFN = opt.alphafile.c_str();
	const char *outFN = opt.outfile.c_str();

	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];
	char **seqname = new char *[nMAXSEQUENCES];

	CSequence *sgi = NULL; // the SV reader owns its own CSequence; the test-sequence one is allocated below

	CCalcWmML wmc(L, K, ::globalConverter.b);
	if (maxnmm==-1)
	{ 
		maxnmm=L;
		if (useTgkm==1)
		{
			maxnmm = 2*(wmc.kernelTruncatedLength-1);
			if (maxnmm>L) 
			{
				maxnmm=L;
			}
		}
		if (useTgkm==2)
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
	}
	double n0 = wmc.n0; 
	(void)n0; // computed for symmetry with the other paths; not used here
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
		//n0 = 0; 
		//kernel = wmc.kernel; 
		c = wmc.h; 
		n0 = c[maxnmm]/2;

	}
    if (useTgkm==3)  //wildcard kernel
    {
        c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, ::globalConverter.b, opt.wildcardLambda, c);
    	n0 = c[maxnmm]/2;

    }
    if (useTgkm==4)  //mismatch kernel
    {
        c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, ::globalConverter.b, c);
    	n0 = c[maxnmm]/2;

    }

    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);

	for(int ii=0;ii<=maxnmm;ii++) {
		snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n c[%d] = %e",ii,c[ii] );Printf(globtmpstr);
	}
	Printf("\n");

	int *mmcnt = new int[L+1];  //mismatch count
	CLList psetL(L,2*maxseqlen+5); // keeps current sequence, used for computation of norm. 
	psetL.UseLookupTable =0;
	int *hdist = psetL.HamDist; 
	int nsvseqs = 0; 

	CSequenceNames *svsn= new CSequenceNames(); 
	svsn->readSeqNamesandWeights(SVSeqIDsFN); 
    snprintf(globtmpstr, GKM_TMPSTR_LEN,"\n  %d SV ids read. \n",svsn->Nseqs);Printf(globtmpstr);


	svsn->openSeqFile(SVSeqsFN, maxseqlen);

	double *alphaovernorm = norm; 

	if (svsn->error) return 1;
	while ((sgi = svsn->nextSeq()) != NULL)
	{
		if(sgi->getLength()>0)
		{
			if (nsvseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
			seqsL[nsvseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			CLTree *psetT = new CLTree();// keeps all the sequences of length L
			psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(), L); 
			if(addRC)
			{
				psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(), L); 
			}			
			seqsL[nsvseqs]->addFromLTree(psetT); 
			psetT->deleteTree(L); 
			delete psetT; 

			alphaovernorm[nsvseqs] = sgi->getWeight()/sqrt(seqsL[nsvseqs]->calcInnerProd(seqsL[nsvseqs], c, mmcnt));

			nsvseqs++; 
		}
	}

	if (svsn->error) return 1;
	double rho = svsn->rho;
	snprintf(globtmpstr, GKM_TMPSTR_LEN,"  %d SV seqs read \n",nsvseqs);Printf(globtmpstr);

	FILE *sfi = fopen(SeqsFN, "r"); 
	if (sfi == NULL) return gkmCannotOpen(SeqsFN);

	int nseqs = nsvseqs; //add test sequences at the end of the svseqs

	sgi = new CSequence(maxseqlen+3);

	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) return 1; 
		if(sgi->getLength()>0)
		{
			if (nseqs>=nMAXSEQUENCES) return gkmTooManySequences(nMAXSEQUENCES);
			seqsL[nseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			CLTree *psetT = new CLTree();// keeps all the sequences of length L
			psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(), L); 
			if(addRC)
			{
				psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(), L); 
			}			
			seqsL[nseqs]->addFromLTree(psetT); 
			psetT->deleteTree(L); 
			delete psetT; 

			seqname[nseqs] = new char [strlength(sgi->getName())+1]; //XXX: should be freed...
			strcpy(seqname[nseqs], sgi->getName()); // buffer sized strlength+1 above

			//NOTE: no alpha for test sequence. 
			alphaovernorm[nseqs] = (1.0)/sqrt(seqsL[nseqs]->calcInnerProd(seqsL[nseqs], c, mmcnt));

			nseqs++; 
		}
	}

	fclose(sfi);
	FILE *fo = fopen(outFN, "w"); 
	if (fo == NULL) return gkmCannotOpen(outFN);

	//test sequences
	for(i=nsvseqs;i<nseqs;i++)
	{
		double svmscore = 0.0;

		//sv sequences
		for(int j=0;j<nsvseqs;j++)
		{
			svmscore += (seqsL[i]->calcInnerProd(seqsL[j],c,mmcnt)*alphaovernorm[i]*alphaovernorm[j]);
		}

		fprintf(fo, "%s\t%f\n",seqname[i], gkmCanon(svmscore - rho));
	}
	fclose(fo); // was never closed: in a long-lived R session the scores stayed in the stdio buffer

	for(i=0;i<nseqs;i++) delete seqsL[i];
	for(i=nsvseqs;i<nseqs;i++) delete []seqname[i];
	delete []seqsL;
	delete []norm;
	delete []seqname;
	delete []mmcnt;
	delete svsn;
	delete sgi;
	return 0;
}

