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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include<iostream>

#include "global.h"
#include "globalvar.h"

#include "Sequence.h"
#include "CountKLmers.h"
#include "CountKLmersGeneral.h"
#include "CountKLmersH.h"
#include "CalcWmML.h"
#include "MLEstimKLmers.h"
#include "MLEstimKLmersLog.h"
#include "KLmer.h"
#include "SequenceNames.h"
#include "EstimLogRatio.h"
#include "LTree.h"
#include "LTreef.h"
#include "LTreeS.h"
#include "LList.h"
#include "CiDLPasses.h"
#include "GTree2.h"

using namespace std;

//#include "stdafx.h"

typedef struct {
	int L;
	int K;
	int maxnmm;
	int maxseqlen;
	int maxnumseq;
	int useTgkm;
	bool addRC;
	bool usePseudocnt;
	bool OutputBinary;
	char *posfile;
	char *negfile;
	char *outfile;
    double wildcardLambda;  //parameter lambda for (LK2004)
    int wildcardMismatchM;  //parameter M for wildcard or mismatch kernels (LK2004)
    char *alphabetFN; //alphabets file name
} OptsGkmKernel;

int gkmKernelSimple(OptsGkmKernel &opt);
int gkmKernelSuffixTree(OptsGkmKernel &opt);

void print_usage_and_exit_gkmKernel(char *prog)
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
	//cout << "  -b             if set, the output matrix will be stored in a binary format" << endl;
	sprintf(globtmpstr,"%s", "  -M             max mismatch for Mismatch kernel or wildcard kernel, default=2\n"); Printf(globtmpstr);
	sprintf(globtmpstr,"%s", "  -L             lambda for wildcard kernel, defaul=1.0\n"); Printf(globtmpstr);
    sprintf(globtmpstr,"%s", "  -A             alphabets file name, if not specified, it is assumed the inputs are DNA sequences \n"); Printf(globtmpstr);
    
    Printf(" \n");

	//exit(0);
}
int mainGkmKernel(int argc, char** argv) //mainGkmKernel
{
    /*printf("welcome!\n");
    printf("\n argc = %d \n", argc);
    for(int i=0;i<argc;i++){
        printf("%s\n",argv[i]);
    }*/
    

	OptsGkmKernel opt;
    
    ::optind=1; // reset getopt()

    int c;

	opt.L = DEF_L; 
	opt.K = DEF_K; 
	opt.maxnmm = DEF_D;
	opt.maxseqlen = DEF_MAXSEQLEN; 
	opt.maxnumseq = DEF_MAXNUMSEQ; 
	opt.useTgkm = DEF_TGKM; 
	opt.addRC = true; 
	opt.usePseudocnt=false; 
	opt.OutputBinary=false; 
    
    opt.wildcardMismatchM= 2;
    opt.wildcardLambda = 1.0;
    opt.alphabetFN = NULL;
    
	int alg = 0;
    
    if (argc == 1) { print_usage_and_exit_gkmKernel(argv[0]); return 0;}

	while ((c = getopt (argc, argv, "l:k:d:m:n:t:a:L:M:A:Rpb")) != -1)
	{
		switch (c) 
		{
			case 'l':
                //printf("\ngetopt  l = %s = %d\n",optarg,atoi(optarg)  );
				opt.L = atoi(optarg);
				break;
			case 'k':
				opt.K = atoi(optarg);
				break;
			case 'd':
				opt.maxnmm = atoi(optarg);
				break;
			case 'm':
				opt.maxseqlen = atoi(optarg);
				break;
			case 'n':
				opt.maxnumseq = atoi(optarg);
				break;
			case 't':
				opt.useTgkm = atoi(optarg);
				break;
			case 'a':
				alg = atoi(optarg);
				break;
			case 'R':
				opt.addRC = false;
				break;
			case 'p':
				opt.usePseudocnt = true;
				break;
			case 'b':
				opt.OutputBinary = true;
				break;
			case 'M':
				opt.wildcardMismatchM = atoi(optarg);
				break;
			case 'L':
				opt.wildcardLambda = atof(optarg);
				break;
			case 'A':
				opt.alphabetFN = optarg;
				break;
			default:
                print_usage_and_exit_gkmKernel(argv[0]);return 0;
		}
	}

	if (argc-optind != 3) { print_usage_and_exit_gkmKernel(argv[0]); return 0;}

	int index = optind;
	opt.posfile = argv[index++];
	opt.negfile = argv[index++];
	opt.outfile = argv[index++];

	//check parameters
	if ((opt.K > opt.L) &&(opt.useTgkm<3))
	{
		Printf("K must be less than or equal to L!\n\n");
        print_usage_and_exit_gkmKernel(argv[0]); return 0;
	}

	if ((opt.maxnmm > 0) && (opt.L < opt.maxnmm))
	{
		Printf("maxMismatch must be less than or equal to L!\n\n");
        print_usage_and_exit_gkmKernel(argv[0]); return 0;
	}

	if (opt.alphabetFN!=NULL){
//		printf("\nj%da%dc%dg%dt%d ",globalConverter.isInAlphabet['j'],globalConverter.isInAlphabet['a'],globalConverter.isInAlphabet['c'],globalConverter.isInAlphabet['g'],globalConverter.isInAlphabet['t']);
		globalConverter.readAlphabetFile(opt.alphabetFN, MAX_ALPHABET_SIZE);
//		printf("\nj%da%dc%dg%dt%d ",globalConverter.isInAlphabet['j'],globalConverter.isInAlphabet['a'],globalConverter.isInAlphabet['c'],globalConverter.isInAlphabet['g'],globalConverter.isInAlphabet['t']);
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
				gkmKernelSuffixTree(opt);
			}
			else
			{
				gkmKernelSimple(opt);
			}
			break;
		case 1:
			gkmKernelSimple(opt);
			break;
		case 2:
			gkmKernelSuffixTree(opt);
			break;
		default:
            print_usage_and_exit_gkmKernel(argv[0]); return 0;
	}

	return 0;
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

	char *posSeqsFN = opt.posfile;
	char *negSeqsFN = opt.negfile;
	char *outFN = opt.outfile;

	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];

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
	if (sfi == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}
	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 
		if(sgi->getLength()>0)
		{

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
			nseqs++; 
		}
	}
	fclose(sfi); 

	npos = nseqs;

	//read negative sequence file
	sfi = fopen(negSeqsFN, "r"); 
	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 
		if(sgi->getLength()>0)
		{

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
			nseqs++; 
		}
	}
	fclose(sfi); 

	nneg = nseqs - npos;

	for(i=0;i<nseqs;i++)
	{
		norm[i] = sqrt(seqsL[i]->calcInnerProd(seqsL[i],c,mmcnt));
	}

	FILE *fo = fopen(outFN, "w"); 

	for(i=0;i<nseqs;i++)
	{
		for(int j=0;j<nseqs;j++)
		{
			if(i>j)
			{
				fprintf(fo, "%e\t",(norm[i]*norm[j]<1E-50)?0.0:seqsL[i]->calcInnerProd(seqsL[j],c,mmcnt)/(norm[i]*norm[j]));
			}
			else if (i==j)
			{
				fprintf(fo, "1.0\t");
			}
		}
		fprintf(fo, "\n"); 
	}

	fclose(fo); 
	//delete []tmps; 
	delete []mmcnt;
	for(i=0;i<nseqs;i++)
	{
		delete seqsL[i];
	}
	delete []seqsL; 
	delete []norm; 

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

	char *posSeqsFN = opt.posfile;
	char *negSeqsFN = opt.negfile;
	char *outFN = opt.outfile;

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
    int UseGTree = 0; // GTree algorithm for speed
    //if(UseGTree){
        //int maxn=(1<<(2*L)) * ::Combinations(L, gMAXMM);
       // if(maxn>100000000){maxn=100000000;};
       // int maxn=100000000;
       // gGTreeLeaves=new GTreeLeafData[maxn]; // list of all the leaf nodes
       // gGTreeLeavesCnt=0; // number of all leaf nodes
        
    //}


	CLTreeS *seqsTS= new CLTreeS();
//    GTree *seqsGTree= new GTree();
 //   GTree2 *seqsGTree2= new GTree2();

    
    
	int **seqsB = new int *[nMAXSEQUENCES]; 
	int **seqsBrc  = new int *[nMAXSEQUENCES]; 

	int *LmersCnt = new int [nMAXSEQUENCES]; 

	CSequence sgii(maxseqlen+3);
	CSequence *sgi = &sgii;

	int ntotal = 0; //number of lmers
	int nseqs=0;

	char**seqname2= NULL; 

	seqname2 = new char *[nMAXSEQUENCES];

	//read positive sequence file
	FILE *sfi = fopen(posSeqsFN, "r"); 
	if (sfi == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}

	while (!feof(sfi))
	{

		sgi->readFsa(sfi); 

		if(sgi->getLength()>0)
		{
			seqname2[nseqs] = new char[100]; 
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
			nseqs++; 
		}
	}
	fclose(sfi);
	npos = nseqs;

	//read negative sequence file
	sfi = fopen(negSeqsFN, "r"); 
	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 

		if(sgi->getLength()>0)
		{
			seqname2[nseqs] = new char[100]; 
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
			nseqs++; 
		}
	}
	fclose(sfi);

	nneg = nseqs - npos;

	// global vars init: 
	gLM1=L-1;
	gMAXMM=maxnmm; //MaxMismatch
	gMMProfile=new int **[nseqs]; 
	for(int i=0;i<nseqs;i++)
	{
		gMMProfile[i] = new int*[gMAXMM+1];
		for (int j=0;j<=gMAXMM;j++)
		{
			gMMProfile[i][j]=new int[nseqs];
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
	int minL2 = L; if (minL2<2) minL2 = 2; 
	for(int i=0;i<=minL2;i++)
	{
		//gDFSlist[i] = new LTreeSnodeData *[uniqueLmerCnt];
	//	gDFSlistT[i] = new CLTreeSptr *[uniqueLmerCnt];  // without nonEmptyDaughterCnt
		gDFSlistT[i] = new CLTreeS *[uniqueLmerCnt];  // with nonEmptyDaughterCnt
		gDFSMMlist[i] = new int[uniqueLmerCnt]; 
		gDFSMMtree[i] = new CbinMMtree *[uniqueLmerCnt];
  	}
	//int *curmmcnt = gDFSMMlist[0];
		
	//gDFSlistT[0][0] = seqsTS->daughter; // without nonEmptyDaughterCnt
	gDFSlistT[0][0] = seqsTS; // with nonEmptyDaughterCnt
	gDFSMMlist[0][0] = 0; 
	
//   int UseGTree = 1;
    if (!UseGTree){
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

        
        int *tmpArray1 = new int[L];
        int *tmpArray2 = new int[L];
        for(int j=0;j<iDL.M;j++){
            sprintf(globtmpstr,"pass %d out of %d.\n",j+1,iDL.M);Printf(globtmpstr);
            CLTreeS *seqsTSj= new CLTreeS();
            seqsTS->cloneReorder(seqsTSj, iDL.passOrder[j], L,L,globalConverter.b, tmpArray1, tmpArray2);
            //seqsTS->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
            gDFSlistT[0][0] = seqsTSj; // with nonEmptyDaughterCnt
            gDFSMMlist[0][0] = 0;
            if(!((iDL.passTrees[j]->child0==nullptr)&&(iDL.passTrees[j]->child1==nullptr))) // i.e. if not empty tree
                seqsTSj->DFSTiDL(gDFSlistT[0],1, gDFSMMlist[0], iDL.passTrees+j, 0, globalConverter.b);
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
            
            
        }
        delete []tmpArray1;
        delete []tmpArray2;
        
    
    }else{
        gGTreeLeaves2=new GTreeLeafData2[uniqueLmerCnt * ::Combinations(L, gMAXMM)]; // list of all the leaf nodes
        gGTreeLeavesCnt=0; // number of all leaf nodes
        GTree2 *seqsGTree2= new GTree2();
        int *tmpArray1 = new int[L];

        seqsTS->addToGTree(seqsGTree2, L,tmpArray1, MAX_ALPHABET_SIZE, L);
        delete []tmpArray1;
        sprintf(globtmpstr," gGTreeLeavesCnt = %d \n",gGTreeLeavesCnt);Printf(globtmpstr);
        
        for(int i=0;i<gGTreeLeavesCnt;i++){
            gGTreeLeaves2[i].process();
        }
        // normalize
        
        for(int i=0;i<nseqs;i++)
        {
            
            //            gMMProfile[i][0][i]/=Combinations(L, gMAXMM);;
            
            for (int j=0;j<=gMAXMM;j++)
            {
                //              int s=::Combinations(L-j, L-gMAXMM);
                for(int k=0;k<=i;k++)
                {
                    int mmijk = (gMMProfile[i][j][k]+gMMProfile[k][j][i]);
                    //                mmijk=mmijk/s;
                    gMMProfile[i][j][k]=gMMProfile[k][j][i]=mmijk;
                }
            }
        }

        for(int i=0;i<nseqs;i++)
        {
            
            //            gMMProfile[i][0][i]/=Combinations(L, gMAXMM);;
            
            for (int j=0;j<=gMAXMM;j++)
            {
                int s=::Combinations(L-j, L-gMAXMM);
//                int s=::Combinations(gMAXMM, j);
                for(int k=0;k<nseqs;k++)
                {
//                    gMMProfile[i][j][k]*=2;
                    gMMProfile[i][j][k]/=s;
                }
            }
        }


        for(int i=0;i<nseqs;i++)
        {
            
            gMMProfile[i][0][i]+=LmersCnt[i];//self;
        }
    }

    delete []nodesAtDepthCnt;
   
    //
	for(int i=0;i<=minL2;i++)
	{
		delete []gDFSlistT[i];
		delete []gDFSMMlist[i];
        delete []gDFSMMtree[i];
	}
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

	FILE *fo = fopen(outFN, "w"); 
	if (fo == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}
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

			for(int j=0;j<nseqs;j++)
			{
				if(i>j)
				{
					if (usePseudocnt)
					{
						fprintf(fo, "%e\t",(norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c, n0,C,LmersCnt[i], LmersCnt[j], btL)/(norm[i]*norm[j]));
					}
					else
					{
						fprintf(fo, "%e\t",(norm[i]*norm[j]<1E-50)?0.0:calcinnerprod(i,j,c)/(norm[i]*norm[j]));
					}

				}
				else if (i==j) 
				{
					fprintf(fo, "1.0\t");
				}
			}
			fprintf(fo, "\n"); 
		}
    }

	fclose(fo); 

	delete []norm;
	delete []LmersCnt;
	seqsTS->deleteTree(L, globalConverter.b, 0);
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
