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
#include "SequenceSet.h"
#include "KernelFile.h"
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

int gkmKernelSimple(OptsGkmKernel &opt);
int gkmKernelSuffixTree(OptsGkmKernel &opt);

void print_usage_and_exit_gkmKernel(const char *prog)
{
    Printf("\n");
    sprintf(globtmpstr, " Usage: %s [options] <pos_seqfile> <neg_seqfile> <outfile>\n",prog );Printf(globtmpstr);
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
    sprintf(globtmpstr,"  -l L           set word length, default= %d\n", DEF_L); Printf(globtmpstr);
	sprintf(globtmpstr,"  -k K           set number of informative columns, default= %d \n",DEF_K); Printf(globtmpstr);
	sprintf(globtmpstr,"  -d maxMismatch set maximum number of mismatches to consider, default= %d \n", DEF_D); Printf(globtmpstr);
    sprintf(globtmpstr,"%s","  -m maxSeqLen   set maximum sequence length in the sequence files,\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"                 default= %d \n", DEF_MAXSEQLEN); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -n maxNumSeq   set maximum number of sequences in the sequence files,\n"); Printf(globtmpstr);
	sprintf(globtmpstr, "                 default= %d\n", DEF_MAXNUMSEQ); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -t filterType  set filter type: 0(use full filter), 1(use truncated filter:\n" ); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "                 this gaurantees non-negative counts for all L-mers), 2(use h[m],\n"); Printf(globtmpstr);
	sprintf(globtmpstr, "                 gkm count vector), 3(wildcard), 4(mismatch), default=%d\n",DEF_TGKM ); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -a algorithm   set algorithm type: 0(auto), 1(XOR Hashtable), 2(tree),\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "                 default=0\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -R             if set, reverse complement sequences will NOT be considered\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -p             if set, a constant to will be added to the count estimates\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -M             max mismatch for Mismatch kernel or wildcard kernel, default=2\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -L             lambda for wildcard kernel, defaul=1.0\n"); Printf(globtmpstr);
    sprintf(globtmpstr,"%s", "  -A             alphabets file name, if not specified, it is assumed the inputs are DNA sequences \n"); Printf(globtmpstr);
    sprintf(globtmpstr,"%s", "  -T             maximum number of threads, defaul=2*l\n"); Printf(globtmpstr);
    sprintf(globtmpstr,"%s", "  -N             merge records that share a name within a file into one row (default: one record = one row)\n"); Printf(globtmpstr);
    sprintf(globtmpstr,"%s", "  -b             write the kernel in the binary .gkmk format (float32) instead of text\n"); Printf(globtmpstr);
    
    Printf(" \n");
}

// Parse argv into opt. Returns 0 on success, 1 if usage was printed (nothing to run).
static int gkmKernelParseArgs(int argc, char** argv, OptsGkmKernel &opt)
{
    ::optind=1; // reset getopt()
    int c;
    if (argc == 1) { print_usage_and_exit_gkmKernel(argv[0]); return 1;}

	while ((c = getopt (argc, argv, "l:k:d:m:n:t:a:L:M:A:T:RpbN")) != -1)
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
	if ((opt.K > opt.L) &&(opt.useTgkm<3))
	{
		Printf("\n ERROR: K must be less than or equal to L!\n"); return 1;
	}

	if ((opt.maxnmm > 0) && (opt.L < opt.maxnmm))
	{
		Printf("\n ERROR: maxMismatch must be less than or equal to L!\n"); return 1;
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

	// the alphabet lives in a global that outlives this call: start every call from DNA
	globalConverter.resetToDNA();
	int alg = opt.alg;
	if (!opt.alphabetFN.empty()){
		if (globalConverter.readAlphabetFile(opt.alphabetFN.c_str(), MAX_ALPHABET_SIZE)!=0) return 1;
		if (opt.addRC&&(globalConverter.b!=4)&&(::globalConverter.b!=16)){
			opt.addRC=false;
			Printf("\nAdd Reverse Complement option is turned off.\n");
		}
		if ((alg!=2)&&(globalConverter.b!=4)){
			alg=2;
			Printf("\nAlgorithm is set to 2 (Tree) to support alphabet size other than 4.\n");
		}
	}

	switch (alg) 
	{
		case 0:
			if ((opt.L-opt.K <= 4) || (opt.maxnmm >= 0 && opt.maxnmm <= 4))
			{
				return gkmKernelSuffixTree(opt);
			}
			else
			{
				return gkmKernelSimple(opt);
			}
		case 1:
			return gkmKernelSimple(opt);
		default:
			return gkmKernelSuffixTree(opt);
	}
}

static void fillGkmkHeader(GkmkWriter &w, const OptsGkmKernel &opt, int maxnmm, int n, int npos)
{
	w.hdr.n = n; w.hdr.npos = npos;
	w.hdr.L = opt.L; w.hdr.K = opt.K; w.hdr.maxnmm = maxnmm; w.hdr.useTgkm = opt.useTgkm;
	w.hdr.b = globalConverter.b; w.hdr.addRC = opt.addRC ? 1 : 0; w.hdr.usePseudocnt = opt.usePseudocnt ? 1 : 0;
	w.hdr.alphabet.assign(globalConverter.alphabet, globalConverter.alphabet + globalConverter.b);
}

static int cannotOpen(const char *fn)
{
	sprintf(globtmpstr,"\n ERROR: cannot open file %s\n", fn); Printf(globtmpstr);
	return 1;
}

static int tooManySequences(int limit)
{
	sprintf(globtmpstr,"\n ERROR: more than %d sequences (set -n / maxnumseq to at least the number of sequences).\n", limit); Printf(globtmpstr);
	return 1;
}


int gkmKernelSimple(OptsGkmKernel &opt)  //Use XOR precomputed hash table
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
	CSequence *sgi= new CSequence(maxseqlen+3);

	CCalcWmML wmc(L, K, globalConverter.b);

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
        c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, opt.wildcardLambda, c);
    	n0 = c[maxnmm]/2;

    }
    if (useTgkm==4)  //mismatch kernel 
    {
        c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, c);
    	n0 = c[maxnmm]/2;

    }

	sprintf(globtmpstr,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);
	for(int ii=0;ii<=maxnmm;ii++) {
		sprintf(globtmpstr,"\n c[%d] = %e",ii,c[ii] ); Printf(globtmpstr);
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
	if (sfi == NULL) return cannotOpen(posSeqsFN);
	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) return 1; 
		if(sgi->getLength()>0)
		{
			if (nseqs>=nMAXSEQUENCES) return tooManySequences(nMAXSEQUENCES);
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
	if (sfi == NULL) return cannotOpen(negSeqsFN);
	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) return 1; 
		if(sgi->getLength()>0)
		{
			if (nseqs>=nMAXSEQUENCES) return tooManySequences(nMAXSEQUENCES);
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

	for(i=0;i<nseqs;i++)
	{
		norm[i] = sqrt(seqsL[i]->calcInnerProd(seqsL[i],c,mmcnt));
	}

	GkmkWriter bin;
	FILE *fo = NULL;
	if (opt.OutputBinary) { fillGkmkHeader(bin, opt, maxnmm, nseqs, npos); bin.values.reserve((size_t)nseqs*(nseqs+1)/2); }
	else { fo = fopen(outFN, "w"); if (fo == NULL) return cannotOpen(outFN); }

	for(i=0;i<nseqs;i++)
	{
		for(int j=0;j<=i;j++)
		{
			double v = (i==j) ? 1.0 : ((norm[i]*norm[j]<1E-50)?0.0:seqsL[i]->calcInnerProd(seqsL[j],c,mmcnt)/(norm[i]*norm[j]));
			if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
			else bin.add(gkmCanon(v));
		}
		if (fo) fprintf(fo, "\n"); 
	}

	if (fo) fclose(fo); 
	else if (bin.write(opt.outfile, records) != 0) return cannotOpen(outFN);
	if (writeIndexSidecar(opt.outfile, records) != 0) { sprintf(globtmpstr,"\n WARNING: could not write %s.index\n", outFN); Printf(globtmpstr); }
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


double calcinnerprod(int i, int j, double *c)
{
	double res = 0; 
	for(int m=0;m<=::gMAXMM;m++)
	{
		res+=::gMMProfile[i][m][j]*c[m]; 
	}
	return(res); 
}


double calcinnerprod(int i, int j, double *c, double n0, double C, int nA, int nB, double btL) // gives inner prodict of the pseudo-counts . nA is the number of L-mers in A and is equal to length(A)-L+1, btL is b^L
{
	double res = 0; 
	for(int m=0;m<=::gMAXMM;m++)
	{
		res+=::gMMProfile[i][m][j]*c[m]; 
	}

	res = res+(nA+nB)*n0*C+btL*n0*n0; 
	return(res); 
}

void task1(int L, int j0, CiDLPasses *iDL, CLTreeS *seqsTS, int M, int nThreads){
  //    sprintf(globtmpstr,"started pass %d out of %d.\n",j+1,iDL->M);Printf(globtmpstr);
  for(int j=0;j<M;j++) if((j%nThreads)==j0){
    
    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    //for(int j=0;j<iDL->M;j++){
    CLTreeS *seqsTSj= new CLTreeS();
    seqsTS->cloneReorder(seqsTSj, iDL->passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    //gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    //gDFSMMlist[0][0] = 0;
    //if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
    //    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL->passTrees+j, 0, globalConverter.b);
    int zero=0;
    if(!((iDL->passTrees[j]->child0==NULL)&&(iDL->passTrees[j]->child1==NULL))) // i.e. if not empty tree
      seqsTSj->DFSTiDL(&seqsTSj,1, &zero, iDL->passTrees+j, 0, globalConverter.b);
    seqsTSj->deleteTree(L, globalConverter.b, 1);
    delete seqsTSj;
    
    // print mismatch profile:
    /* for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=gMAXMM; dd++){
    printf("%d ",gMMProfile[si][dd][sj]);
    }
    }
  }
    */
    
    
    //}
    delete []tmpArray1;
    delete []tmpArray2;
    //sprintf(globtmpstr,"ended pass %d out of %d.\n",j+1,iDL->M);Printf(globtmpstr);
    }
  
  }

int gkmKernelSuffixTree(OptsGkmKernel &opt)  //maingKernel
{
  
  int L = opt.L;
  int K = opt.K;
  int maxseqlen =	opt.maxseqlen;
  int useTgkm = opt.useTgkm;
  int maxnmm = opt.maxnmm; //auto 
  int nMAXSEQUENCES = opt.maxnumseq; 
  bool addRC = opt.addRC;
  bool usePseudocnt= opt.usePseudocnt; 
  
  const char *posSeqsFN = opt.posfile.c_str();
  const char *negSeqsFN = opt.negfile.c_str();
  const char *outFN = opt.outfile.c_str();
  
  int i = 1; 
  
  //char tmps[1000];
  
  char **seqname = new char *[nMAXSEQUENCES];
  
  CCalcWmML wmc(L, K, globalConverter.b);
  //double *kernel = wmc.kernelTruncated;
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
    if (maxnmm>L)
    {
      maxnmm=L;
    }
  }
  double n0 = wmc.n0; 
  double *c = wmc.cTr; 
  
  n0 = c[maxnmm]/2; 
  
  if (!useTgkm)
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
    c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, opt.wildcardLambda, c);
    n0 = c[maxnmm]/2;
    
  }
  if (useTgkm==4)  //mismatch kernel 
  {
    c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, globalConverter.b, c);
    n0 = c[maxnmm]/2;
    
  }
  
  sprintf(globtmpstr,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);
  for(int ii=0;ii<=maxnmm;ii++) {
    sprintf(globtmpstr,"\n c[%d] = %e",ii,c[ii] ); 	Printf(globtmpstr);
  }
  Printf("\n");
  
  int npos=0; 
  int nneg=0;
  
  gMAXMM=maxnmm; //MaxMismatch
  
  CLTreeS *seqsTS= new CLTreeS();
  
  
  
  int **seqsB = new int *[nMAXSEQUENCES]; 
  int **seqsBrc  = new int *[nMAXSEQUENCES]; 
  
  int *LmersCnt = new int [nMAXSEQUENCES]; 
  
  CSequence sgii(maxseqlen+3);
  CSequence *sgi = &sgii;
  
  int ntotal = 0; //number of lmers
  int nseqs=0;
  
  char**seqname2= NULL; 
  
  seqname2 = new char *[nMAXSEQUENCES];
  std::vector<SeqRecord> records; // row identity, written to <outfile>.index
  
  //read positive sequence file
  FILE *sfi = fopen(posSeqsFN, "r"); 
  if (sfi == NULL) return cannotOpen(posSeqsFN);
  
  char *tmpSeq=new char[maxseqlen+3];
  int  *tmpSeqB=new int[maxseqlen+3];
  
  int multiseq_allowed = opt.mergeByName ? 1 : 0; // Phase 3: default one record = one row
  
  
  while (!feof(sfi))
  {
    
    if (sgi->readFsa(sfi) < 0) return 1; 
    
    if(sgi->getLength()>0)
    {
      int dupl_seq_idx= -1;
      if(multiseq_allowed){
        dupl_seq_idx = find_str(seqname,nseqs, sgi->getName());
      }
      if (dupl_seq_idx>=0){
        records[dupl_seq_idx].length += sgi->getLength();
        sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
        ntotal-=LmersCnt[dupl_seq_idx];
        LmersCnt[dupl_seq_idx]+= seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        ntotal+=LmersCnt[dupl_seq_idx];
        if(addRC)
        {
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
          LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        }
      }else{
        if (nseqs>=nMAXSEQUENCES) return tooManySequences(nMAXSEQUENCES);
        seqname2[nseqs] = new char[strlen(sgi->getName())+1]; // was a fixed 100 bytes: names >= 100 chars overflowed it
        sprintf(seqname2[nseqs],"%s", sgi->getName());
        seqname[nseqs]=seqname2[nseqs]; 
        
        seqsB[nseqs] = new int [sgi->getLength()]; 
        sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]); 
        LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
        
        if(addRC)
        {
          seqsBrc[nseqs] = new int [sgi->getLength()];
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]); 
          LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
        }
        else
        {
          seqsBrc[nseqs]=NULL; 
        }
        
        ntotal = ntotal + LmersCnt[nseqs]; 
        records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), 0});
        nseqs++;
      }
    }
  }
  fclose(sfi);
  npos = nseqs;
  
  //read negative sequence file
  sfi = fopen(negSeqsFN, "r"); 
  if (sfi == NULL) return cannotOpen(negSeqsFN);
  while (!feof(sfi))
  {
    if (sgi->readFsa(sfi) < 0) return 1; 
    
    if(sgi->getLength()>0)
    {
      // check if segName already exists
      int dupl_seq_idx= -1;
      if(multiseq_allowed){
        dupl_seq_idx = find_str(seqname+npos,nseqs-npos, sgi->getName());
      }
      if (dupl_seq_idx>=0){
        dupl_seq_idx+=npos;
        records[dupl_seq_idx].length += sgi->getLength();
        sgi->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
        ntotal-=LmersCnt[dupl_seq_idx];
        LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        ntotal+=LmersCnt[dupl_seq_idx];
        if(addRC)
        {
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, tmpSeqB);
          LmersCnt[dupl_seq_idx] += seqsTS->addSequence(tmpSeqB, sgi->getLength(),L, dupl_seq_idx);
        }
      }else{
        if (nseqs>=nMAXSEQUENCES) return tooManySequences(nMAXSEQUENCES);
        seqname2[nseqs] = new char[strlen(sgi->getName())+1];
        sprintf(seqname2[nseqs],"%s", sgi->getName());
        seqname[nseqs]=seqname2[nseqs];
        
        seqsB[nseqs] = new int [sgi->getLength()];
        sgi->getSubseqBaseId(0, sgi->getLength()-1, seqsB[nseqs]);
        LmersCnt[nseqs] = seqsTS->addSequence(seqsB[nseqs], sgi->getLength(),L, nseqs);
        if(addRC)
        {
          seqsBrc[nseqs] = new int [sgi->getLength()];
          sgi->getReverseComplement()->getSubseqBaseId(0, sgi->getLength()-1, seqsBrc[nseqs]);
          LmersCnt[nseqs] = LmersCnt[nseqs] + seqsTS->addSequence(seqsBrc[nseqs], sgi->getLength(),L, nseqs);
          
        }
        else
        {
          seqsBrc[nseqs]=NULL; 
        }
        
        ntotal = ntotal + LmersCnt[nseqs]; 
        records.push_back(SeqRecord{nseqs, "", sgi->getName(), 0, sgi->getLength(), 0});
        nseqs++; 
        
      }
    }
  }
  fclose(sfi);
  
  delete []tmpSeq;
  delete []tmpSeqB;
  
  nneg = nseqs - npos;
  for(int i=0;i<nseqs;i++) { records[i].id = seqRecordId(i, npos); records[i].label = (i<npos) ? 1 : -1; records[i].nlmers = LmersCnt[i]; }
  
  // global vars init: 
  gLM1=L-1;
  gMAXMM=maxnmm; //MaxMismatch
  gMMProfile=new aint **[nseqs];
  for(int i=0;i<nseqs;i++)
  {
    gMMProfile[i] = new aint*[gMAXMM+1];
    for (int j=0;j<=gMAXMM;j++)
    {
      gMMProfile[i][j]=new aint[nseqs];
      for(int k=0;k<nseqs;k++)
      {
        gMMProfile[i][j][k]=0;
      }
    }
  }
  
  int *nodesAtDepthCnt = new int[L];
  for(int i=0;i<L; i++){
    nodesAtDepthCnt[i]=0;
  }
  
  int uniqueLmerCnt = seqsTS->leavesCount(0,L, globalConverter.b, nodesAtDepthCnt);
  sprintf(globtmpstr,"\n npos %d \n nneg %d \n  ntotal %d \n nunique %d\n",npos,nneg,ntotal,uniqueLmerCnt);Printf(globtmpstr);
  {
    // if no IDL bound
    /*
    seqsTS->DFST(gDFSlistT[0],1, gDFSMMlist[0], 0, globalConverter.b);
    
    for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=gMAXMM; dd++){
    printf("%d ",gMMProfile[si][dd][sj]);
    }
    }
    }
    */
    // else if IDL bound then
    for(int i=0;i<L; i++){
      //     sprintf(globtmpstr,"d%d , %d\n", i, nodesAtDepthCnt[i]);Printf(globtmpstr);
    }
    
    CiDLPasses iDL;
    //iDL.newIDLPasses(L, gMAXMM);
    double p=1.0/::globalConverter.b;
    
    //iDL.initPassOrderAll(L, gMAXMM);
    //iDL.newGreedyIDLPasses(L,iDL.M,  gMAXMM, nodesAtDepthCnt, p);
    
    iDL.newGreedyIDLPasses(L,2*L,  gMAXMM, nodesAtDepthCnt, p);
    
    //iDL.newGreedy2IDLPasses(L,2*L,  gMAXMM, nodesAtDepthCnt, p);
    //iDL.newGreedy2IDLPasses(L,L,  gMAXMM, nodesAtDepthCnt, p);
    
    //iDL.newPassOrderDesignCover( L, gMAXMM, 3);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
    //iDL.newGreedyIDLPasses(L,iDL.M,  gMAXMM, nodesAtDepthCnt, p);
    
    
    /*    int *tmpArray1 = new int[L];
    int *tmpArray2 = new int[L];
    for(int j=0;j<iDL.M;j++){
    sprintf(globtmpstr,"pass %d out of %d.\n",j+1,iDL.M);Printf(globtmpstr);
    CLTreeS *seqsTSj= new CLTreeS();
    seqsTS->cloneReorder(seqsTSj, iDL.passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
    //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
    gDFSMMlist[0][0] = 0;
    if(!((iDL.passTrees[j]->child0==NULL)&&(iDL.passTrees[j]->child1==NULL))) // i.e. if not empty tree
    seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
    seqsTSj->deleteTree(L, globalConverter.b, 1);
    delete seqsTSj;
    
    // print mismatch profile:
    / * for(int si=0;si<nseqs; si++){
    for(int sj=0;sj<nseqs;sj++){
    printf("\n (s%d, s%d) = ",si,sj);
    for(int dd = 0; dd<=gMAXMM; dd++){
    printf("%d ",gMMProfile[si][dd][sj]);
    }
    }
    }
    * /
    }
    delete []tmpArray1;
    delete []tmpArray2;
    */
    
    //myThreads[0] = std::thread(task1, L, 0, &iDL, seqsTS);
    //myThreads[1] = std::thread(task1, L, 1, &iDL, seqsTS);
    //myThreads[0].join();
    //myThreads[1].join();
    
    int nThreads=std::thread::hardware_concurrency();
    if(nThreads==0){nThreads = iDL.M;}
    if(nThreads>iDL.M){nThreads =iDL.M;}
    nThreads = iDL.M/(int)(iDL.M/nThreads);
    if(opt.maxnThread<nThreads){nThreads=opt.maxnThread;}
    
    sprintf(globtmpstr,"Running %d passes on %d thread%s.\n", iDL.M, nThreads, (nThreads==1)?"":"s"); Printf(globtmpstr);
    if (nThreads<=1){
      task1( L, 0, &iDL, seqsTS, iDL.M, 1);
    }else{
      
#ifndef MULTI_THREAD_SAFE
      Printf("Warning -- MULTI_THREAD_SAFE is not enabled (see src/global.h). Some values may be approximated.\n");
#endif
      
      
      std::thread *myThreads = new std::thread[nThreads];
      int j;
      
      for(j=0;j<nThreads;j++){
        myThreads[j] = std::thread(task1, L, j, &iDL, seqsTS, iDL.M, nThreads);
        // myThreads[j].join();
      }
      for(j=0;j<nThreads;j++){
        myThreads[j].join();
      }
      delete []myThreads;
    }
    
  }
  
  delete []nodesAtDepthCnt;
  
  /*
  for(int i=0;i<nseqs;i++)
  {
  for(int k=0;k<nseqs;k++)
  {
  for (int j=0;j<=gMAXMM;j++)
  {
  printf("(%d,%d)[%d] = %d\n",i,k,j, gMMProfile[i][j][k]);
  }
  }
  }
  */	
  
  /// calc C 
  double C =0; 
  for(int m=0;m<=L;m++)
  {
    C+=dCombinations(L,m)*pow(1.0*globalConverter.b-1,m)*wmc.kernelTruncated[m];
    //		C+=dCombinations(L,m)*pow(3.0,m)*wmc.kernelTruncated[m];
  }
  
  //	double btL=pow(4.0,L);
  double btL=pow(1.0*globalConverter.b,L);
  double *norm = new double [nseqs]; 
  
  GkmkWriter bin;
  FILE *fo = NULL;
  if (opt.OutputBinary) { fillGkmkHeader(bin, opt, maxnmm, nseqs, npos); bin.values.reserve((size_t)nseqs*(nseqs+1)/2); }
  else { fo = fopen(outFN, "w"); if (fo == NULL) return cannotOpen(outFN); }
  /*
  if (OutputMismatchProfileOnly)
  {
  fprintf(fo, "%d\tL (length)\n", L); 
  fprintf(fo, "%d\td (maximum number of mismatches)\n", gMAXMM); 
  fprintf(fo, "%d\tNp (number of sequences in positive class)\n", npos);
  fprintf(fo, "%d\tNn (number of sequences in negative class)\n", nneg); 
  
  for (int nmm=0;nmm<=gMAXMM;nmm++)
  {
  fprintf(fo, "d=%d\n",nmm);
  
  
  for(i=0;i<nseqs;i++)
  {
  
  if (outputClassLabel)
  {
  fprintf(fo, "%d\t", (i<npos)?1:-1);
  }
  if (OutputSeqNames)
  {
  fprintf(fo, "%s\t", seqname[i]);
  }
  
  for(int j=0;j<nseqs;j++)
  {
  if(i>=j)
  {
  fprintf(fo, "%d\t",gMMProfile[i][nmm][j]);
  }
  }
  fprintf(fo, "\n"); 
  }
  }
  }
  else 
  */
  {
    for(i=0;i<nseqs;i++)
    {
      if (usePseudocnt)
      {
        norm[i] = sqrt(calcinnerprod(i,i,c,n0,C,LmersCnt[i], LmersCnt[i], btL));
      }
      else
      {
        norm[i] = sqrt(calcinnerprod(i,i,c));
      }
    }
    
    for(i=0;i<nseqs;i++)
    {
      
      //if (outputClassLabel)
      //{
      //	fprintf(fo, "%d\t", (i<npos)?1:-1);
      //}
      //if (OutputSeqNames)
      //{
      //	fprintf(fo, "%s\t", seqname[i]);
      //}
      
      for(int j=0;j<=i;j++)
      {
        double v;
        if (i==j) v = 1.0;
        else if (usePseudocnt) v = (norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c, n0,C,LmersCnt[i], LmersCnt[j], btL)/(norm[i]*norm[j]);
        else v = (norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c)/(norm[i]*norm[j]);
        if (fo) { if (i==j) fprintf(fo, "1.0\t"); else fprintf(fo, "%e\t", v); }
        else bin.add(gkmCanon(v));
      }
      if (fo) fprintf(fo, "\n"); 
    }
  }
  
  if (fo) fclose(fo); 
  else if (bin.write(opt.outfile, records) != 0) return cannotOpen(outFN);
  if (writeIndexSidecar(opt.outfile, records) != 0) { sprintf(globtmpstr,"\n WARNING: could not write %s.index\n", outFN); Printf(globtmpstr); }
  
  delete []norm;
  delete []LmersCnt;
  seqsTS->deleteTree(L, globalConverter.b, 0);
  delete seqsTS;
  //delete []curmmcnt; 
  
  for(int i=0;i<nseqs;i++)
  {//printf("\n4 %d\n",i);
    delete []seqsB[i]; 
    if (seqsBrc[i]!=NULL) delete []seqsBrc[i]; 
    for (int j=0;j<=gMAXMM;j++)
    {
      delete []gMMProfile[i][j];
    }
    delete []gMMProfile[i];
  }
  delete []gMMProfile;
  
  delete []seqname; 
  if (seqname2!=NULL)
  {
    for(i=0;i<nseqs;i++)
    {
      delete []seqname2[i]; 
    }
    delete []seqname2; 
  }
  
  return 0; 
}
