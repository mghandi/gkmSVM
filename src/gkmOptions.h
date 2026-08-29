/* gkmOptions.h : parameter structs shared by the command-line front ends and the R bindings.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The R bindings used to serialise their arguments into a fake argv and re-parse it with getopt
 * (silently truncating paths at 5000 characters and reporting invalid parameters as success).
 * Both front ends now fill these structs directly; the CLI keeps getopt.
 */
#pragma once
#include <string>

struct OptsGkmKernel {
	int L = 10;                 // word length
	int K = 6;                  // number of informative columns
	int maxnmm = 3;             // maximum mismatches, -1 = automatic
	int maxseqlen = 10000;
	int maxnumseq = 1000000;
	int useTgkm = 1;            // filter type 0..4
	int alg = 0;                // 0 auto, 1 XOR hash table, 2 tree
	bool addRC = true;
	bool usePseudocnt = false;
	bool OutputBinary = false;  // parsed for compatibility; implemented in Phase 4
	bool mergeByName = false;   // -N: records with the same FASTA name inside one file are merged into one kernel row (pre-Phase-3 behaviour)
	std::string posfile, negfile, outfile;
	double wildcardLambda = 1.0;
	int wildcardMismatchM = 2;
	std::string alphabetFN;     // empty = DNA
	int maxnThread = 1000;
	int tileRows = 0;          // -r: rows of the mismatch profile held at once (0 = automatic from tileMemoryMB)
	int tileMemoryMB = 1024;   // memory budget for the profile when tileRows is automatic
	int passDesign = 0;        // -P: 0 automatic (greedy iDL table up to GKM_MAX_PATTERN_TABLE patterns, prefix-split beyond), 1 greedy iDL, 2 prefix-split (Phase 7)
};

struct OptsSVMClassify {
	int L = 10;
	int K = 6;
	int maxnmm = 3;
	int maxseqlen = 10000;
	int maxnumseq = 1000000;
	int useTgkm = 1;
	int alg = 0;
	int batchSize = 100000;
	bool addRC = true;
	bool usePseudocnt = false;
	std::string seqfile, svseqfile, alphafile, outfile;
	double wildcardLambda = 1.0;
	int wildcardMismatchM = 2;
	std::string alphabetFN;
	bool legacyNorm = false;    // -y: normalise DNA test/SV sequences over all mismatch levels as versions <= 0.80 did (score still truncated at maxnmm); deprecated
};

struct OptsGkmTrain {
	std::string kernelfile, posfile, negfile, outprefix;
	double C = 1.0;             // C-SVC regularisation
	double posWeight = 1.0;     // class weight for the positive class (negative class weight = 1)
	double eps = 0.001;         // stopping tolerance
	bool shrinking = true;
	int nfold = 0;              // > 1: report nfold cross-validation accuracy before training
	int cacheMB = 200;
	int maxseqlen = 1000000;
	bool legacyPair = true;     // write {prefix}_svalpha.out + {prefix}_svseq.fa
	bool modelFile = true;      // write {prefix}.gkmmodel
	bool quiet = false;
	std::string alphabetFN;     // Phase 7: per-track alphabet specs ("dna,=01") when the inputs are multi-track FASTA; empty = single-track (any alphabet: the sequences are copied as they are)
	int maxseqlenPad = 0;       // (unused)
};
int gkmTrainRun(OptsGkmTrain &opt);
int mainSVMtrain(int argc, char *argv[]);

// Validate the options, load the alphabet, run. Return 0 on success, non-zero on error (a message
// has been printed through Printf). These are the entry points used by the R package.
int gkmKernelRun(OptsGkmKernel &opt);
int svmClassifyRun(OptsSVMClassify &opt);

// getopt front ends (the standalone binaries): parse argv into the struct and run.
int mainGkmKernel(int argc, char **argv);
int mainSVMclassify(int argc, char **argv);
