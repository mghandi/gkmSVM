/* libsvm_svm.cpp : compiles the vendored LIBSVM (src/libsvm/svm.cpp) as part of the package.
 * R compiles only the top level of src/, and listing sub-directory sources in Makevars needs GNU
 * make extensions that CRAN does not allow; this one-line include is the portable alternative. */
#include "libsvm/svm.cpp"
