# Builds the standalone gkmSVM command-line binaries from the same sources as the R package.
#
#   make                 -> build/gkmsvm_kernel build/gkmsvm_classify build/gkmsvm_train
#   make ASAN=1          -> AddressSanitizer + UndefinedBehaviorSanitizer build (into build-asan/)
#   make test            -> build, then run the golden tests (tests/golden/golden.py check)
#   make oracle          -> cross-check coefficient tables and small kernels against dev/oracle
#   make bench           -> DNA kernel wall-clock benchmark (dev/bench.sh)
#   make clean
#
# Variables: CXX, CXXFLAGS, BUILD (output dir), ASAN=1, PYTHON.
# This Makefile is NOT used by R CMD INSTALL (see .Rbuildignore); R compiles src/*.cpp itself.

CXX      ?= c++
PYTHON   ?= python3
CXXSTD   ?= -std=c++17
OPT      ?= -O2
WARN     ?= -Wall -Wno-deprecated-declarations -Wno-char-subscripts -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-warning-option -Wno-format-overflow -Wno-format-truncation
BUILD    ?= build

ifeq ($(ASAN),1)
  BUILD    := build-asan
  OPT      := -O1 -g -fno-omit-frame-pointer
  SANFLAGS := -fsanitize=address,undefined -fno-sanitize-recover=undefined
endif

CXXFLAGS ?= $(CXXSTD) $(OPT) $(WARN)
CXXFLAGS += $(SANFLAGS) -Isrc
LDFLAGS  += $(SANFLAGS) -pthread

# Everything in src/ except the Rcpp glue; each src/cli/*.cpp supplies one main().
LIB_SRCS := $(filter-out src/RcppExports.cpp src/gkmsvm_kernel.cpp src/gkmsvm_classify.cpp,$(wildcard src/*.cpp))
LIB_OBJS := $(patsubst src/%.cpp,$(BUILD)/obj/%.o,$(LIB_SRCS))

BINS := $(BUILD)/gkmsvm_kernel $(BUILD)/gkmsvm_classify $(BUILD)/gkmsvm_train

.PHONY: all test oracle bench clean
all: $(BINS)

$(BUILD)/obj/%.o: src/%.cpp | $(BUILD)/obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/obj/cli_%.o: src/cli/main_%.cpp | $(BUILD)/obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/gkmsvm_%: $(BUILD)/obj/cli_%.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD)/obj:
	mkdir -p $@

test: all
	$(PYTHON) tests/golden/golden.py check --bin $(BUILD)

oracle: all
	$(PYTHON) tests/oracle/oracle_check.py --bin $(BUILD)

bench: all
	dev/bench.sh $(BUILD)

clean:
	rm -rf build build-asan
