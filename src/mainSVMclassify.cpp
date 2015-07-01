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

#include "Sequence.h"
#include "CountKLmers.h"
#include "CountKLmersGeneral.h"
#include "CountKLmersH.h"
//#include "CountKLmers3G.h"
#include "CalcWmML.h"
#include "MLEstimKLmers.h"
//#include "MLEstimKLmers3G.h"
#include "MLEstimKLmersLog.h"
#include "KLmer.h"
#include "SequenceNames.h"
#include "EstimLogRatio.h"
#include "LTree.h"
#include "LTreef.h"
#include "LTreeS.h"
#include "LList.h"
#include "SVMtrain.h"

using namespace std;

typedef struct {
	int L;
	int K;
	int maxnmm;
	int maxseqlen;
	int maxnumseq;
	int useTgkm;
	int batchSize;
	bool addRC;
	bool usePseudocnt;
	char *seqfile;
	char *svseqfile;
	char *alphafile;
	char *outfile;
    double wildcardLambda;  //parameter lambda for (LK2004)
    int wildcardMismatchM;  //parameter M for wildcard or mismatch kernels (LK2004)
    char *alphabetFN; //alphabets file name

} OptsSVMClassify;

double calcnorm(CSequence *sgi, int addRC, CLList *tmplist, double *c, int *mmcnt, int L, int maxmm);
double svmScoreunorm(int i, double *c); 

int svmClassifySuffixTree(OptsSVMClassify &opt);
int svmClassifySimple(OptsSVMClassify &opt);

void print_usage_and_exit(char *prog)
{
    Printf("\n");
    sprintf(globtmpstr, " Usage: %s [options] <test_seqfile> <sv_seqfile> <sv_alphafile> <outfile>\n",prog );Printf(globtmpstr);
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
	sprintf(globtmpstr,"  -l L           set word length, default=%d\n",DEF_L);
	sprintf(globtmpstr,"  -k K           set number of informative columns, default=%d\n", DEF_K);
    sprintf(globtmpstr, "  -d maxMismatch set maximum number of mismatches to consider, default=%d\n",DEF_D);
	Printf("  -m maxSeqLen   set maximum sequence length in the sequence files,\n");
    sprintf(globtmpstr,"                 default=%d\n",DEF_MAXSEQLEN);
	Printf("  -n maxNumSeq   set maximum number of sequences in the sequence files,\n");
	sprintf(globtmpstr, "                 default=%d\n",DEF_MAXNUMSEQ);
	Printf("  -t filterType  set filter type: 0(use full filter), 1(use truncated filter:\n");
	Printf("                 this gaurantees non-negative counts for all L-mers), 2(use h[m],\n");
	sprintf(globtmpstr, "                 gkm count vector), 3(wildcard), 4(mismatch), default=%d\n",DEF_TGKM);
	Printf("  -a algorithm   set algorithm type: 0(auto), 1(XOR Hashtable), 2(tree),\n");
	Printf("                 default=0\n");
	Printf("  -b             set number of sequences to compute scores for in batch, \n");
	sprintf(globtmpstr, "                 default=%d\n", DEF_BATCHSIZE);
	Printf("  -R             if set, reverse complement sequences will NOT be considered\n");
	Printf("  -p             if set, a constant to count estimates will be added\n");
	Printf("  -M             max mismatch for Mismatch kernel or wildcard kernel, default=2\n");
	Printf("  -L             lambda for wildcard kernel, defaul=0.9\n");
	Printf("  -A             alphabets file name, if not specified, it is assumed the inputs are DNA sequences\n");


    Printf("\n");

	//exit(0);
}

int mainSVMclassify(int argc, char** argv) // mainSVMclassify
{
	OptsSVMClassify opt;
    ::optind=1; // reset getopt()
	int c;

	opt.L = DEF_L; 
	opt.K = DEF_K; 
	opt.maxnmm = DEF_D;
	opt.maxseqlen = DEF_MAXSEQLEN; 
	opt.maxnumseq = DEF_MAXNUMSEQ; 
	opt.useTgkm = DEF_TGKM; 
	opt.batchSize = DEF_BATCHSIZE;
	opt.addRC = true; 
	opt.usePseudocnt=false;

    opt.wildcardMismatchM= 2;
    opt.wildcardLambda = 1.0;
    opt.alphabetFN = NULL;

	int alg = 0;

    if (argc == 1) { print_usage_and_exit(argv[0]); return 1; }

	while ((c = getopt (argc, argv, "l:k:d:m:n:t:a:b:M:L:A:Rp")) != -1)
	{
		switch (c) 
		{
			case 'l':
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
			case 'b':
				opt.batchSize= atoi(optarg);
				break;
			case 'R':
				opt.addRC = false;
				break;
			case 'p':
				opt.usePseudocnt = true;
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
				print_usage_and_exit(argv[0]);return 1;
		}
	}

	if (argc-optind != 4) { print_usage_and_exit(argv[0]);return 1; }

	int index = optind;
	opt.seqfile = argv[index++];
	opt.svseqfile = argv[index++];
	opt.alphafile = argv[index++];
	opt.outfile = argv[index++];

	if (opt.alphabetFN!=NULL){
		::globalConverter.readAlphabetFile(opt.alphabetFN, MAX_ALPHABET_SIZE);
		if (opt.addRC&&(::globalConverter.b!=4)&&(::globalConverter.b!=16)){
			opt.addRC=false;
			Printf("\nAdd Reverse Complement option is turned off.\n");
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
				svmClassifySuffixTree(opt);
			}
			else
			{
				if ((opt.L-opt.K <= 2) || (opt.maxnmm >= 0 && opt.maxnmm <= 2))
				{
					svmClassifySuffixTree(opt);
				}
				else
				{
					svmClassifySimple(opt);
				}
			}
			break;
		case 1:
			svmClassifySimple(opt);
			break;
		case 2:
			svmClassifySuffixTree(opt);
			break;
		default:
			print_usage_and_exit(argv[0]);return 1;
	}

	return 0;
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

	char *SeqsFN = opt.seqfile;
	char *SVSeqsFN = opt.svseqfile;
	char *SVSeqIDsFN = opt.alphafile;
	char *outFN = opt.outfile;

	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];
	char **seqname = new char *[nMAXSEQUENCES];

	CSequence *sgi= new CSequence(maxseqlen+3);

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

    sprintf(globtmpstr,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);

	for(int ii=0;ii<=maxnmm;ii++) {
		sprintf(globtmpstr,"\n c[%d] = %e",ii,c[ii] );Printf(globtmpstr);
	}
	Printf("\n");

	int *mmcnt = new int[L+1];  //mismatch count
	CLList psetL(L,2*maxseqlen+5); // keeps current sequence, used for computation of norm. 
	psetL.UseLookupTable =0;
	int *hdist = psetL.HamDist; 
	int nsvseqs = 0; 

	CSequenceNames *svsn= new CSequenceNames(); 
	svsn->readSeqNamesandWeights(SVSeqIDsFN); 
    sprintf(globtmpstr,"\n  %d SV ids read. \n",svsn->Nseqs);Printf(globtmpstr);


	svsn->openSeqFile(SVSeqsFN, maxseqlen);

	double *alphaovernorm = norm; 

	for(i=0;i<svsn->Nseqs;i++)
	{
		sgi = svsn->nextSeq(); 
		if(sgi==NULL)
		{
            sprintf(globtmpstr,"\n the sequences for only %d out of %d sequence names in SVs file (%s) were found. \n", i,svsn->Nseqs, SVSeqIDsFN);Printf(globtmpstr);

			break; 
		}

		if(sgi->getLength()>0)
		{
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

	sprintf(globtmpstr,"  %d SV seqs read \n",nsvseqs);Printf(globtmpstr);

	FILE *sfi = fopen(SeqsFN, "r"); 
	if (sfi == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}

	int nseqs = nsvseqs; //add test sequences at the end of the svseqs

	sgi = new CSequence(maxseqlen+3);

	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 
		if(sgi->getLength()>0)
		{
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
			sprintf(seqname[nseqs],"%s", sgi->getName());

			//NOTE: no alpha for test sequence. 
			alphaovernorm[nseqs] = (1.0)/sqrt(seqsL[nseqs]->calcInnerProd(seqsL[nseqs], c, mmcnt));

			nseqs++; 
		}
	}

	FILE *fo = fopen(outFN, "w"); 
	if (fo == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}

	//test sequences
	for(i=nsvseqs;i<nseqs;i++)
	{
		double svmscore = 0.0;

		//sv sequences
		for(int j=0;j<nsvseqs;j++)
		{
			svmscore += (seqsL[i]->calcInnerProd(seqsL[j],c,mmcnt)*alphaovernorm[i]*alphaovernorm[j]);
		}

		fprintf(fo, "%s\t%f\n",seqname[i], svmscore);
	}

	return 0;
}

//given fasta file for SVs and the corresponding weights, outputs and another file for the test sequences, gives the SVM score
int svmClassifySuffixTree(OptsSVMClassify &opt)
{
	int i;
	int L = opt.L; 
	int K = opt.K; 
	int maxnmm = opt.maxnmm; //auto 
	int useTgkm = opt.useTgkm; 
	int maxseqlen =	opt.maxseqlen; 
	int nMAXSEQUENCES = opt.maxnumseq; 
	bool addRC = opt.addRC; 
	int batchSize = opt.batchSize; //batch size 

	char *SeqsFN = opt.seqfile;
	char *SVSeqsFN = opt.svseqfile;
	char *SVSeqIDsFN = opt.alphafile;
	char *outFN = opt.outfile;
    
//	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];
	char **seqname = new char *[nMAXSEQUENCES];

	CSequence *sgi; //= new CSequence(maxseqlen+3);

    char *tmps = new char[maxseqlen+L+2+1000]; 

	CCalcWmML wmc(L, K, ::globalConverter.b);
	double *kernel = wmc.kernelTruncated;
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
	double *c = wmc.cTr; 
	if (!useTgkm)
	{
		n0 = 0; 
		kernel = wmc.kernel; 
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
        c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, 4, opt.wildcardLambda, c);
    	n0 = c[maxnmm]/2;

    }
    if (useTgkm==4)  //mismatch kernel
    {
        c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, 4, c);
    	n0 = c[maxnmm]/2;

    }
    
	sprintf(globtmpstr,"\n maximumMismatch = %d\n", maxnmm);Printf(globtmpstr);
	for(int ii=0;ii<=maxnmm;ii++) {
		sprintf(globtmpstr,"\n c[%d] = %e",ii,c[ii] ); 	Printf(globtmpstr);
	}
	Printf("\n");
	
//	char *tmps = new char[maxseqlen+L+2]; 
	int *mmcnt = new int[L+1];  //mismatch count

	CLList psetL(L,2*maxseqlen+5); // keeps current sequence, used for computation of norm. 
	psetL.UseLookupTable =0;
//	int *hdist = psetL.HamDist; 
	int nsvseqs = 0; 

	CSequenceNames *svsn= new CSequenceNames(); 
	svsn->readSeqNamesandWeights(SVSeqIDsFN); 
	sprintf(globtmpstr,"\n  %d SV ids read. \n",svsn->Nseqs);Printf(globtmpstr);

	svsn->openSeqFile(SVSeqsFN, maxseqlen);

	CLTreef *tSVs= new CLTreef(); //keeps all the sequences of length L in support vectors
	
	//double *alphaovernorm = norm; 

	for(i=0;i<svsn->Nseqs;i++)
	{
		sgi = svsn->nextSeq(); 
		if(sgi==NULL)
		{
            sprintf(globtmpstr,"\n the sequences for only %d out of %d sequence names in SVs file (%s) were found. \n", i,svsn->Nseqs, SVSeqIDsFN);Printf(globtmpstr);

			break; 
		}

		if(sgi->getLength()>0)
		{
//			printf(" %d %s\n", i,sgi->getName()); 
			//seqsL[nsvseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			
			double alphaovernorm = sgi->getWeight()/calcnorm(sgi, addRC, &psetL, c, mmcnt,L,maxnmm);

			tSVs->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L, alphaovernorm); 
			if(addRC)
			{
				tSVs->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L, alphaovernorm); 
			}	
			nsvseqs++; 
		}
	}
	
	sprintf(globtmpstr,"  %d SV seqs read \n",nsvseqs); Printf(globtmpstr);

	delete svsn; 

	CLTreeS *seqsTS= new CLTreeS();
//	int **seqsB = new int *[nMAXSEQUENCES]; 
//	int **seqsBrc  = new int *[nMAXSEQUENCES]; 
//	int *LmersCnt = new int [nMAXSEQUENCES]; 
	int **seqsB = new int *[batchSize+2]; 
	int **seqsBrc  = new int *[batchSize+2]; 
	int *LmersCnt = new int [batchSize+2]; 

		
//	CSequenceNames *seqsn= new CSequenceNames(); 
//	seqsn->readSeqNamesandWeights(SeqIDsFN); 
//	seqsn->openSeqFile(SeqsFN, maxseqlen);
	FILE *sfi = fopen(SeqsFN, "r"); 
	if (sfi == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}

	FILE *fo = fopen(outFN, "w"); 
	if (fo == NULL)
	{
		perror ("error occurred while opening a file");
		return 0;
	}

	int ntotal = 0; //number of lmers
	int nseqs = 0; 

	sgi = new CSequence(maxseqlen+3);

	while (!feof(sfi))
	{
		sgi->readFsa(sfi); 

//	for(int ii=0;ii<seqsn->Nseqs;ii++)
//	{
//		sgi = seqsn->nextSeq(); 
//		if(sgi==NULL)
//		{
//			printf("\n the sequences for only %d out of %d sequence names in SVs file (%s) were found. \n", ii,seqsn->Nseqs, SeqIDsFN); 
//			break; 
//		}


		if(sgi->getLength()>0)
		{
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

			norm[nseqs]=calcnorm(sgi, addRC, &psetL, c, mmcnt,L, maxnmm);
			seqname[nseqs] = new char [strlength(sgi->getName())+1]; //XXX: should be freed...
			sprintf(seqname[nseqs],"%s", sgi->getName());
			//seqname[nseqs] = sgi->getNameLink(); //seqsn->seqNames[ii]; 			

			ntotal = ntotal + LmersCnt[nseqs]; 
			nseqs++; 
		}
	
		if(((nseqs%batchSize)==0)||(feof(sfi)))
		{

			// global vars init: 
			gLM1=L-1;
			gMAXMM=maxnmm; //MaxMismatch
	
			gMMProfile0 = new myFlt*[gMAXMM+1];

			for (int j=0;j<=gMAXMM;j++)
			{
				gMMProfile0[j]=new myFlt[nseqs];
				for(int k=0;k<nseqs;k++)
				{
					gMMProfile0[j][k]=0;
				}
			}

			int uniqueLmerCnt = seqsTS->leavesCount(0,L, ::globalConverter.b, NULL);

			int minL2 = L; if (minL2<2) minL2 = 2; 
			for(i=0;i<=minL2;i++)
			{
				//gDFSlist[i] = new LTreeSnodeData *[uniqueLmerCnt];
//				gDFSlistT[i] = new CLTreeSptr *[uniqueLmerCnt];  // without nonEmptyDaughterCnt
				gDFSlistT[i] = new CLTreeS *[uniqueLmerCnt]; // with nonEmptyDaughterCnt
				gDFSMMlist[i] = new int[uniqueLmerCnt]; 
		
			}
			//int *curmmcnt = gDFSMMlist[0];
		
//			gDFSlistT[0][0] = seqsTS->daughter; // without nonEmptyDaughterCnt
			gDFSlistT[0][0] = seqsTS;// with nonEmptyDaughterCnt
			gDFSMMlist[0][0] = 0; 
			tSVs->DFST(gDFSlistT[0],1, gDFSMMlist[0], 0, ::globalConverter.b);

			for(i=0;i<=minL2;i++)
			{
				delete []gDFSlistT[i];
				delete []gDFSMMlist[i]; 
			}

			for(i=0;i<nseqs;i++)
			{
				fprintf(fo, "%s\t%f\n",seqname[i], svmScoreunorm(i,c)/norm[i]);
			}

			seqsTS->deleteTree(L, ::globalConverter.b, 0);
			seqsTS->initTree();

			for(int i=0;i<nseqs;i++)
			{//printf("\n4 %d\n",i);
				delete []seqsB[i]; 
				if (seqsBrc[i]!=NULL) delete []seqsBrc[i]; 
			}

			for (int j=0;j<=gMAXMM;j++)
			{
				delete []gMMProfile0[j];
			}
			delete []gMMProfile0;

			nseqs = 0; 
		}
	
	}

	fclose(fo); 

	delete []norm;
	delete []LmersCnt;

	delete []seqname; 
	//delete seqsn; 
	delete []mmcnt; 
	delete []tmps;
//	delete sgi; 
	tSVs->deleteTree(L,::globalConverter.b);
	delete tSVs; 
	return 0; 
}

double svmScoreunorm(int i, double *c)
{
	double res = 0; 
	for(int m=0;m<=::gMAXMM;m++)
	{
		res+=::gMMProfile0[m][i]*c[m]; 
	}
	return(res); 
}

double calcnorm(CSequence *sgi, int addRC, CLList *tmplist, double *c, int *mmcnt, int L, int maxnmm)
{
		if (::globalConverter.b==4){
			//calc norm 
			CLTree *psetT = new CLTree();// keeps all the sequences of length L
			psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L); 
			if(addRC)
			{
				psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L); 
			}	
			tmplist->clear(); 
			tmplist->addFromLTree(psetT); 
			psetT->deleteTree(L); 
			delete psetT; 
			return sqrt(tmplist->calcInnerProd(tmplist,c, mmcnt)); 
		}
		// using Tree based method to calc mismatch profile and calc norm
		CLTreef *psetT = new CLTreef();// keeps all the sequences of length L
		psetT->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L);
		if(addRC)
		{
			psetT->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L);
		}
		//tmplist->clear();
		//tmplist->addFromLTree(psetT);
		for(int i=0;i<=maxnmm;i++){
			mmcnt[i]=0;
		}
		psetT->iimismatchCountGeneral(psetT,L,mmcnt,maxnmm,::globalConverter.b);

		double s=0;
		for(int i=0;i<=maxnmm;i++){
			s += mmcnt[i]*c[i];
		}
		psetT->deleteTree(L, ::globalConverter.b);
		delete psetT;
		return(sqrt(s));
}
