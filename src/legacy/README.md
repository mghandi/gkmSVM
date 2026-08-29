Quarantined code (Phase 2 of dev/REFACTORING_PLAN.md)

Nothing in this directory is compiled: R builds only the top level of src/, and the Makefile lists
src/*.cpp explicitly. These files were unreachable from both R entry points (gkmsvm_kernel,
gkmsvm_classify) and from gkmsvm_train:

  CountKLmers, CountKLmersGeneral, CountKLmersH, MLEstimKLmers, MLEstimKLmersLog, KLmer,
  EstimLogRatio, SequenceData (empty stubs), LKTree, GTree/GTreeLeafData (never constructed),
  GTree2/GTreeLeafData2 (only reachable through the literal `int UseGTree = 0;` branch that was
  deleted from mainGkmKernel.cpp together with CLTreeS::addToGTree).

They are kept for one release for reference and are deleted in Phase 6 (git keeps the history).
