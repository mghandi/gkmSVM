/* main_kernel.cpp : standalone entry point for gkmsvm_kernel.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Thin shim over mainGkmKernel() so the CLI binary can be built from the same
 * source tree as the R package (see Makefile at the repository root). Not part
 * of the R package build: R compiles only the top level of src/, not src/cli/.
 */
int mainGkmKernel(int argc, char** argv);

int main(int argc, char** argv) { return mainGkmKernel(argc, argv); }
