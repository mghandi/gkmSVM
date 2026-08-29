/* mainSVMtrain.cpp : gkmsvm_train - C-SVC on a precomputed gkm kernel matrix (LIBSVM backend).
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
 *
 * Phase 4b of dev/REFACTORING_PLAN.md: the previous solver (CSVMtrain, a fixed-iteration
 * Jaakkola-Diekhans-Haussler heuristic without a C parameter or a stopping criterion) is replaced by
 * LIBSVM (src/libsvm, BSD-3) with kernel_type = PRECOMPUTED. The output formats are unchanged:
 * {prefix}_svalpha.out (name TAB alpha, negative-class alphas negated) and {prefix}_svseq.fa, plus the
 * single-file {prefix}.gkmmodel of Phase 4, which also carries the bias rho.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <set>

#include "global.h"
#include "globalvar.h"
#include "gkmOptions.h"
#include "Sequence.h"
#include "KernelFile.h"
#include "libsvm/svm.h"

using namespace std;

static void print_usage_train()
{
	Printf("\n");
	Printf(" Usage: gkmsvm_train [options] <kernel_file> <pos_seqfile> <neg_seqfile> <out_prefix>\n");
	Printf("\n");
	Printf("  trains a C-SVC (LIBSVM) on a precomputed gkm kernel matrix.\n");
	Printf("\n");
	Printf(" Arguments:\n");
	Printf("  kernel_file: kernel matrix written by gkmsvm_kernel (text or binary .gkmk)\n");
	Printf("  pos_seqfile: positive sequence file used to generate the kernel\n");
	Printf("  neg_seqfile: negative sequence file used to generate the kernel\n");
	Printf("  out_prefix:  prefix of the output files {PREFIX}_svalpha.out, {PREFIX}_svseq.fa\n");
	Printf("               and {PREFIX}.gkmmodel\n");
	Printf("\n");
	Printf(" Options:\n");
	Printf("  -c C        regularisation parameter, default=1\n");
	Printf("  -w W        weight of the positive class (negative class weight is 1), default=1\n");
	Printf("  -e eps      stopping tolerance, default=0.001\n");
	Printf("  -v nfold    report nfold cross-validation accuracy before training\n");
	Printf("  -S          disable shrinking\n");
	Printf("  -q          quiet (no solver output)\n");
	Printf("  -n niter20  accepted and ignored (parameter of the pre-4b solver)\n");
	Printf("\n");
}

static void svm_print_via_Printf(const char *s) { Printf(s); }
static void svm_print_quiet(const char *) {}

static std::string trimmed(const char *s)
{
	std::string r(s);
	while (!r.empty() && (r.back() == '\n' || r.back() == '\r' || r.back() == ' ' || r.back() == '\t')) r.pop_back();
	return r;
}

// Reads the lower triangle (including the diagonal) of a text or binary kernel into a full N x N
// matrix. Returns 0 on success. N (if > 0) is the expected size; for binary files it is checked.
static int readKernel(const std::string &fn, int &N, std::vector<double> &K)
{
	FILE *f = fopen(fn.c_str(), "rb");
	if (f == NULL) { sprintf(globtmpstr, "\n ERROR: cannot open kernel file %s\n", fn.c_str()); Printf(globtmpstr); return 1; }
	unsigned char magic[4] = {0, 0, 0, 0};
	size_t got = fread(magic, 1, 4, f);
	if (got == 4 && memcmp(magic, "GKMK", 4) == 0) {
		GkmkReader rd;
		int rc = rd.read(f);
		fclose(f);
		if (rc != 0) { sprintf(globtmpstr, "\n ERROR: %s: %s\n", fn.c_str(), rd.error.c_str()); Printf(globtmpstr); return 1; }
		if (N > 0 && rd.hdr.n != N) { sprintf(globtmpstr, "\n ERROR: kernel has %d rows but the sequence files have %d sequences\n", rd.hdr.n, N); Printf(globtmpstr); return 1; }
		N = rd.hdr.n;
		K.assign((size_t)N * N, 0.0);
		size_t k = 0;
		for (int i = 0; i < N; i++) for (int j = 0; j <= i; j++) { K[(size_t)i * N + j] = K[(size_t)j * N + i] = rd.values[k++]; }
		return 0;
	}
	rewind(f);
	K.assign((size_t)N * N, 0.0);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			double v;
			if (fscanf(f, "%lf", &v) != 1) { sprintf(globtmpstr, "\n ERROR: error reading kernel %s at row %d (expected %d rows)\n", fn.c_str(), i, N); Printf(globtmpstr); fclose(f); return 1; }
			K[(size_t)i * N + j] = K[(size_t)j * N + i] = v;
		}
	}
	double extra;
	if (fscanf(f, "%lf", &extra) == 1) { sprintf(globtmpstr, "\n ERROR: kernel %s has more than %d rows\n", fn.c_str(), N); Printf(globtmpstr); fclose(f); return 1; }
	fclose(f);
	return 0;
}

static int readFasta(const std::string &fn, int maxseqlen, std::vector<std::string> &names, std::vector<std::string> &seqs)
{
	CSequence sgi(maxseqlen + 3);
	FILE *sfi = fopen(fn.c_str(), "r");
	if (sfi == NULL) { sprintf(globtmpstr, "\n ERROR: cannot open file %s\n", fn.c_str()); Printf(globtmpstr); return 1; }
	while (!feof(sfi)) {
		if (sgi.readFsa(sfi, true) < 0) { fclose(sfi); return 1; }
		if (sgi.getLength() > 0) { names.push_back(sgi.getName()); seqs.push_back(trimmed(sgi.getSeq())); }
	}
	fclose(sfi);
	return 0;
}

int gkmTrainRun(OptsGkmTrain &opt)
{
	if (opt.C <= 0 || opt.posWeight <= 0 || opt.eps <= 0) { Printf("\n ERROR: C, the class weight and eps must be positive.\n"); return 1; }

	std::vector<std::string> names, seqs;
	if (readFasta(opt.posfile, opt.maxseqlen, names, seqs) != 0) return 1;
	int npos = (int)names.size();
	if (readFasta(opt.negfile, opt.maxseqlen, names, seqs) != 0) return 1;
	int N = (int)names.size();
	int nneg = N - npos;
	sprintf(globtmpstr, "npos=%d, nneg=%d, N=%d\n", npos, nneg, N); Printf(globtmpstr);
	if (npos == 0 || nneg == 0) { Printf("\n ERROR: both classes must contain sequences.\n"); return 1; }

	std::vector<double> K;
	if (readKernel(opt.kernelfile, N, K) != 0) return 1;

	// the legacy model pair identifies support vectors by name: use the positional ids when names repeat
	std::vector<std::string> ids(N);
	{
		std::set<std::string> seen; bool dup = false;
		for (int i = 0; i < N; i++) { if (!seen.insert(names[i]).second) dup = true; ids[i] = seqRecordId(i, npos); }
		if (dup) { Printf("Note: sequence names are not unique; the model files use the row ids pos1.., neg1.. instead\n"); }
		else ids = names;
	}

	// LIBSVM problem with a precomputed kernel: x[i][0].value is the 1-based row index,
	// x[i][j].value = K(i, j-1) for j = 1..N, terminated by index -1.
	std::vector<svm_node> nodes((size_t)N * (N + 2));
	std::vector<svm_node *> x(N);
	std::vector<double> y(N);
	for (int i = 0; i < N; i++) {
		svm_node *row = &nodes[(size_t)i * (N + 2)];
		row[0].index = 0; row[0].value = i + 1;
		for (int j = 0; j < N; j++) { row[j + 1].index = j + 1; row[j + 1].value = K[(size_t)i * N + j]; }
		row[N + 1].index = -1; row[N + 1].value = 0;
		x[i] = row;
		y[i] = (i < npos) ? 1.0 : -1.0;
	}
	svm_problem prob; prob.l = N; prob.y = y.data(); prob.x = x.data();

	svm_parameter param;
	memset(&param, 0, sizeof param);
	param.svm_type = C_SVC; param.kernel_type = PRECOMPUTED;
	param.C = opt.C; param.eps = opt.eps; param.cache_size = opt.cacheMB;
	param.shrinking = opt.shrinking ? 1 : 0; param.probability = 0;
	int wlabel[1] = {1}; double wval[1] = {opt.posWeight};
	if (opt.posWeight != 1.0) { param.nr_weight = 1; param.weight_label = wlabel; param.weight = wval; }
	const char *perr = svm_check_parameter(&prob, &param);
	if (perr) { sprintf(globtmpstr, "\n ERROR: %s\n", perr); Printf(globtmpstr); return 1; }
	svm_set_print_string_function(opt.quiet ? &svm_print_quiet : &svm_print_via_Printf);

	if (opt.nfold > 1) {
		std::vector<double> target(N);
		svm_cross_validation(&prob, &param, opt.nfold, target.data());
		int correct = 0; for (int i = 0; i < N; i++) if (target[i] == y[i]) correct++;
		sprintf(globtmpstr, "Cross validation (%d folds) accuracy = %g%%\n", opt.nfold, 100.0 * correct / N); Printf(globtmpstr);
	}

	svm_model *model = svm_train(&prob, &param);
	int nsv = svm_get_nr_sv(model);
	std::vector<int> svidx(nsv);
	svm_get_sv_indices(model, svidx.data());
	// sv_coef is y_i*alpha_i with y = +1 for model->label[0]; make +1 the positive class regardless
	double sign = (model->label[0] == 1) ? 1.0 : -1.0;
	double rho = sign * model->rho[0];
	sprintf(globtmpstr, "nSV=%d, rho=%.10e\n", nsv, rho); Printf(globtmpstr);

	std::string alphaFN = opt.outprefix + "_svalpha.out", svFN = opt.outprefix + "_svseq.fa", modelFN = opt.outprefix + ".gkmmodel";
	FILE *fo_alpha = NULL, *fo_sv = NULL, *fo_model = NULL;
	if (opt.legacyPair) {
		fo_alpha = fopen(alphaFN.c_str(), "w"); fo_sv = fopen(svFN.c_str(), "w");
		if (fo_alpha == NULL || fo_sv == NULL) { sprintf(globtmpstr, "\n ERROR: cannot write %s / %s\n", alphaFN.c_str(), svFN.c_str()); Printf(globtmpstr); svm_free_and_destroy_model(&model); return 1; }
	}
	if (opt.modelFile) {
		fo_model = fopen(modelFN.c_str(), "w");
		if (fo_model == NULL) { sprintf(globtmpstr, "\n ERROR: cannot write %s\n", modelFN.c_str()); Printf(globtmpstr); svm_free_and_destroy_model(&model); return 1; }
		fprintf(fo_model, "#gkmmodel 1\n#rho %.10e\n#nsv %d\n#npos %d\n#nneg %d\n#C %g\n#solver libsvm-%d\n", rho, nsv, npos, nneg, opt.C, LIBSVM_VERSION);
	}
	for (int k = 0; k < nsv; k++) {
		int i = svidx[k] - 1;                      // 1-based training index -> row
		double coef = sign * model->sv_coef[0][k]; // signed alpha: > 0 for positives, < 0 for negatives
		if (fabs(coef) <= 1e-10) continue;
		if (fo_alpha) { fprintf(fo_alpha, "%s\t%e\n", ids[i].c_str(), coef); fprintf(fo_sv, ">%s\n%s\n", ids[i].c_str(), seqs[i].c_str()); }
		if (fo_model) fprintf(fo_model, ">%s\t%e\n%s\n", ids[i].c_str(), coef, seqs[i].c_str());
	}
	if (fo_alpha) { fclose(fo_alpha); fclose(fo_sv); }
	if (fo_model) fclose(fo_model);
	svm_free_and_destroy_model(&model);
	return 0;
}

int mainSVMtrain(int argc, char *argv[]) //mainSVMtrain
{
	OptsGkmTrain opt;
	int nfp = 0, perr = 0;
	for (int i = 1; i < argc && !perr; i++) {
		std::string a = argv[i];
		bool hasArg = (i + 1 < argc);
		if (a == "-c" && hasArg) opt.C = atof(argv[++i]);
		else if (a == "-w" && hasArg) opt.posWeight = atof(argv[++i]);
		else if (a == "-e" && hasArg) opt.eps = atof(argv[++i]);
		else if (a == "-v" && hasArg) opt.nfold = atoi(argv[++i]);
		else if (a == "-n" && hasArg) { i++; Printf("Note: -n niter20 belongs to the previous solver and is ignored; use -c C.\n"); }
		else if (a == "-S") opt.shrinking = false;
		else if (a == "-q") opt.quiet = true;
		else if (a.size() > 0 && a[0] == '-') { sprintf(globtmpstr, "\n parameter not recognized: %s \n", argv[i]); Printf(globtmpstr); perr = 1; }
		else {
			if (nfp == 0) opt.kernelfile = a; else if (nfp == 1) opt.posfile = a; else if (nfp == 2) opt.negfile = a; else if (nfp == 3) opt.outprefix = a;
			else { sprintf(globtmpstr, "\n parameter not recognized: %s \n", argv[i]); Printf(globtmpstr); perr = 1; }
			nfp++;
		}
	}
	if (nfp != 4) perr = 1;
	if (perr) { print_usage_train(); return 1; }
	return gkmTrainRun(opt);
}
