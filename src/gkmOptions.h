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
	std::string posfile, negfile, outfile;
	double wildcardLambda = 1.0;
	int wildcardMismatchM = 2;
	std::string alphabetFN;     // empty = DNA
	int maxnThread = 1000;
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
};

// Validate the options, load the alphabet, run. Return 0 on success, non-zero on error (a message
// has been printed through Printf). These are the entry points used by the R package.
int gkmKernelRun(OptsGkmKernel &opt);
int svmClassifyRun(OptsSVMClassify &opt);

// getopt front ends (the standalone binaries): parse argv into the struct and run.
int mainGkmKernel(int argc, char **argv);
int mainSVMclassify(int argc, char **argv);
