/* classify_tree_impl.h : the tree (suffix-trie) classify driver, compiled once per alphabet size
 * (see trie_b4.cpp / trie_b32.cpp and the note in global.h). Was the second half of mainSVMclassify.cpp.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
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

namespace GKM_NS {

double calcnorm(CSequence *sgi, int addRC, CLList *tmplist, double *c, int *mmcnt, int L, int maxmm, int b, bool legacyNorm);
double svmScoreunorm(int i, double *c, const ScoreContext &sc); 

//given fasta file for SVs and the corresponding weights, outputs and another file for the test sequences, gives the SVM score
int svmClassifySuffixTree(OptsSVMClassify &opt, const CConverter &conv)
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

	const char *SeqsFN = opt.seqfile.c_str();
	const char *SVSeqsFN = opt.svseqfile.c_str();
	const char *SVSeqIDsFN = opt.alphafile.c_str();
	const char *outFN = opt.outfile.c_str();
    
//	CLList **seqsL = new CLList *[nMAXSEQUENCES];
	double *norm = new double [nMAXSEQUENCES];
	char **seqname = new char *[nMAXSEQUENCES];

	CSequence *sgi; //= new CSequence(maxseqlen+3);

    char *tmps = new char[maxseqlen+L+2+1000]; 

	CCalcWmML wmc(L, K, conv.b);
	double *kernel = wmc.kernelTruncated;
	(void)kernel; // computed for symmetry with the other paths; not used here
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
        c = wmc.calcWildcardKernelWeights(L,  opt.wildcardMismatchM, conv.b, opt.wildcardLambda, c); // was a literal 4: wrong for b != 4
    	n0 = c[maxnmm]/2;

    }
    if (useTgkm==4)  //mismatch kernel
    {
        c = wmc.calcMismatchKernelWeights(L,  opt.wildcardMismatchM, conv.b, c); // was a literal 4: wrong for b != 4
    	n0 = c[maxnmm]/2;

    }
    
	gkmMsg("\n maximumMismatch = %d\n", maxnmm);
	for(int ii=0;ii<=maxnmm;ii++) {
		gkmMsg("\n c[%d] = %e",ii,c[ii] );
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
	gkmMsg("\n  %d SV ids read. \n",svsn->Nseqs);

	svsn->openSeqFile(SVSeqsFN, maxseqlen, &conv);

	CLTreef *tSVs= new CLTreef(); //keeps all the sequences of length L in support vectors
	
	//double *alphaovernorm = norm; 

	if (svsn->error) return 1;
	while ((sgi = svsn->nextSeq()) != NULL)
	{
		if(sgi->getLength()>0)
		{
			//seqsL[nsvseqs] = new CLList(L, 2*maxseqlen+5, hdist);  
			
			double alphaovernorm = sgi->getWeight()/calcnorm(sgi, addRC, &psetL, c, mmcnt,L,maxnmm, conv.b, opt.legacyNorm);

			tSVs->addSequence(sgi->getSeqBaseId(), sgi->getLength(),L, alphaovernorm); 
			if(addRC)
			{
				tSVs->addSequence(sgi->getReverseComplement()->getSeqBaseId(), sgi->getLength(),L, alphaovernorm); 
			}	
			nsvseqs++; 
		}
	}
	
	if (svsn->error) return 1;
	double rho = svsn->rho;
	gkmMsg("  %d SV seqs read \n",nsvseqs);

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
	if (sfi == NULL) return gkmCannotOpen(SeqsFN);

	FILE *fo = fopen(outFN, "w"); 
	if (fo == NULL) return gkmCannotOpen(outFN);

	int ntotal = 0; //number of lmers
	int nseqs = 0; 

	sgi = new CSequence(maxseqlen+3, &conv);

	while (!feof(sfi))
	{
		if (sgi->readFsa(sfi) < 0) { fclose(fo); remove(outFN); return 1; } // the output was already opened: do not leave a partial file

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
			if (nseqs>=nMAXSEQUENCES) { fclose(fo); remove(outFN); return gkmTooManySequences(nMAXSEQUENCES); }
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

			norm[nseqs]=calcnorm(sgi, addRC, &psetL, c, mmcnt,L, maxnmm, conv.b, opt.legacyNorm);
			seqname[nseqs] = new char [strlength(sgi->getName())+1]; //XXX: should be freed...
			strcpy(seqname[nseqs], sgi->getName()); // buffer sized strlength+1 above
			//seqname[nseqs] = sgi->getNameLink(); //seqsn->seqNames[ii]; 			

			ntotal = ntotal + LmersCnt[nseqs]; 
			nseqs++; 
		}
	
		if((nseqs>0) && (((nseqs%batchSize)==0)||(feof(sfi)))) // an empty batch (EOF right after a full one) used to write gDFSlistT[0][0] into a zero-length array
		{

			// global vars init: 
			ScoreContext sctx; // per-batch state (was the sctx.mmProfile0/sctx.maxmm/gLM1 globals)
			sctx.LM1=L-1;
			sctx.maxmm=maxnmm; //MaxMismatch
	
			sctx.mmProfile0 = new myFlt*[sctx.maxmm+1];

			for (int j=0;j<=sctx.maxmm;j++)
			{
				sctx.mmProfile0[j]=new myFlt[nseqs];
				for(int k=0;k<nseqs;k++)
				{
					sctx.mmProfile0[j][k]=0;
				}
			}

			int uniqueLmerCnt = seqsTS->leavesCount(0,L, conv.b, NULL);

			int minL2 = L; if (minL2<2) minL2 = 2; 
			ScoreScratch scratch; // per-batch DFS lists, one per trie level (was gDFSlistT/gDFSMMlist)
			scratch.listT.resize(minL2+1); scratch.mmlist.resize(minL2+1);
			for(i=0;i<=minL2;i++)
			{
				scratch.listT[i] = new CLTreeS *[uniqueLmerCnt]; // with nonEmptyDaughterCnt
				scratch.mmlist[i] = new int[uniqueLmerCnt]; 
			}
		
			scratch.listT[0][0] = seqsTS;// with nonEmptyDaughterCnt
			scratch.mmlist[0][0] = 0; 
			tSVs->DFST(scratch.listT[0],1, scratch.mmlist[0], 0, conv.b, &sctx, &scratch);

			for(i=0;i<=minL2;i++)
			{
				delete []scratch.listT[i];
				delete []scratch.mmlist[i]; 
			}

			for(i=0;i<nseqs;i++)
			{
				{ double sc = svmScoreunorm(i,c,sctx)/norm[i] - rho; if (sc == 0.0) sc = 0.0; fprintf(fo, "%s\t%f\n",seqname[i], sc); }
			}

			seqsTS->deleteTree(L, conv.b, 0);
			seqsTS->initTree();

			for(int i=0;i<nseqs;i++)
			{
				delete []seqsB[i]; 
				if (seqsBrc[i]!=NULL) delete []seqsBrc[i]; 
				delete []seqname[i]; 
			}

			for (int j=0;j<=sctx.maxmm;j++)
			{
				delete []sctx.mmProfile0[j];
			}
			delete []sctx.mmProfile0;

			nseqs = 0; 
		}
	
	}

	fclose(fo); 
	fclose(sfi); 

	delete []norm;
	delete []LmersCnt;
	delete []seqsB;
	delete []seqsBrc;

	delete []seqname; 
	delete []mmcnt; 
	delete []tmps;
	delete sgi; 
	seqsTS->deleteTree(L, conv.b, 0);
	delete seqsTS;
	tSVs->deleteTree(L,conv.b);
	delete tSVs; 
	return 0; 
}

double svmScoreunorm(int i, double *c, const ScoreContext &sc)
{
	double res = 0; 
	for(int m=0;m<=sc.maxmm;m++)
	{
		res+=sc.mmProfile0[m][i]*c[m]; 
	}
	return(res); 
}

double calcnorm(CSequence *sgi, int addRC, CLList *tmplist, double *c, int *mmcnt, int L, int maxnmm, int b, bool legacyNorm)
{
		if (b==4){
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
			// The score only sums mismatch levels m <= maxnmm (svmScoreunorm), and so does the b!=4
			// branch below; this branch used to sum all m <= L, i.e. norm and score were inconsistent
			// whenever -d was smaller than the support of c[]. Truncate c[] the same way.
			// legacyNorm (-y / legacyNorm=TRUE) reproduces the pre-0.90 scores exactly: norm over all m <= L.
			if (legacyNorm) return sqrt(tmplist->calcInnerProd(tmplist, c, mmcnt)); 
			double *cTrunc = new double[L+1];
			for(int i=0;i<=L;i++){ cTrunc[i] = (i<=maxnmm) ? c[i] : 0.0; }
			double s = tmplist->calcInnerProd(tmplist,cTrunc, mmcnt); 
			delete []cTrunc;
			return sqrt(s); 
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
		psetT->iimismatchCountGeneral(psetT,L,mmcnt,maxnmm,b);

		double s=0;
		for(int i=0;i<=maxnmm;i++){
			s += mmcnt[i]*c[i];
		}
		psetT->deleteTree(L, b);
		delete psetT;
		return(sqrt(s));
}

} // namespace GKM_NS
