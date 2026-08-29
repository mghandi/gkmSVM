/* trie_b32.cpp : instantiation of the trie classes and the tree-path drivers for alphabets of up to
 * 32 symbols (namespace gkm_b32). See the note in global.h. Used for any alphabet with more than 4 symbols (up to GKM_MAX_ALPHABET).
 */
#define MAX_ALPHABET_SIZE 32
#define GKM_NS gkm_b32
#include "LTreeS_impl.h"
#include "LTreef_impl.h"
#include "kernel_tree_impl.h"
#include "classify_tree_impl.h"
