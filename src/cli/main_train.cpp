/* main_train.cpp : standalone entry point for gkmsvm_train.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Thin shim over mainSVMtrain(); see main_kernel.cpp. Phase 4b of
 * dev/REFACTORING_PLAN.md replaces the solver behind it with LIBSVM.
 */
int mainSVMtrain(int argc, char* argv[]);

int main(int argc, char** argv) { return mainSVMtrain(argc, argv); }
