LIBSVM 3.37 (https://github.com/cjlin1/libsvm, tag v337), BSD-3-Clause — see COPYRIGHT.

Files: svm.h and COPYRIGHT are verbatim. svm.cpp carries one small, clearly delimited patch
("gkmSVM patch") so that the R package build (-DRPACKAGE) complies with CRAN policy:
  * rand() -> R_unif_index()          (used by svm_cross_validation's shuffling and the probability code)
  * fputs(stdout) -> Rprintf           (the default print function; gkmSVM installs its own anyway)
  * fprintf(stderr, ...) -> gkm_warn   (REprintf under RPACKAGE, stderr otherwise)
Without RPACKAGE the behaviour is exactly upstream's. To upgrade: copy the new svm.cpp/svm.h/COPYRIGHT
over these, re-apply the patch (dev/REFACTORING_PLAN.md Phase 4b), and update this note.

gkmSVM uses only kernel_type = PRECOMPUTED with svm_type = C_SVC (src/mainSVMtrain.cpp).
