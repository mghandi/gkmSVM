/* main_classify.cpp : standalone entry point for gkmsvm_classify.
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Thin shim over mainSVMclassify(); see main_kernel.cpp.
 */
int mainSVMclassify(int argc, char** argv);

int main(int argc, char** argv) { return mainSVMclassify(argc, argv); }
