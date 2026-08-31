// gkm_metal.mm -- Route 2 ("G3") prototype of the gkm / general-B kernel on Apple Metal.
//
// The mismatch profile N_m[i][j] (window pairs of sequences i,j with per-block mismatch counts m)
// is computed by classifying every pair of DISTINCT l-mers on the GPU (bit-packed codes, XOR +
// per-field popcount per block) and scattering cnt_u (x) cnt_v into per-class n x n counters with
// integer atomics.  Only classes with |m| <= d are kept (d = -1: all).  The kernel
// K = sum_m c(m) N_m (filter or gkm coefficients) is formed on the CPU in double precision,
// normalised, and written in the .gkmk binary format of gkmsvm_kernel for direct comparison.
//
// Build:  make            (clang++ -std=c++17 -O2 -fobjc-arc -framework Metal -framework Foundation)
// Usage:  gkm_metal -l L -k K -d D [-t 0|2] -A dna,=01 [-R] [-o out.gkmk] [-L light] [-C chunk] in.mfa
// Limits (prototype): l * ceil(log2 max b) <= 64 bits per word, <= 8 blocks, kinds filter (0) / gkm (2).
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <chrono>
#include <algorithm>

static const char *kShader = R"MSL(
#include <metal_stdlib>
using namespace metal;

struct Params {
    uint U, n, fieldBits, nBlocks, light, rowLo, rowHi, ovCap;
    ulong lowMask;                // one bit per position (bit p*f)
    ulong blockMask[8];           // lowMask restricted to each block
    uint  blockStride[8];         // mixed-radix strides of the class index
    uint  nClassesTotal;
    uint  selfPairs;              // 1: this class group holds class 0 -> add the self pairs (u, u)
};

inline uint tri(uint i, uint j) {          // lower-triangle index, i >= j
    return i * (i + 1) / 2 + j;
}

inline uint classify(ulong a, ulong b, constant Params &P, thread bool &ok) {
    ulong x = a ^ b;
    ulong m = x;
    for (uint s = 1; s < P.fieldBits; s++) m |= (x >> s);
    m &= P.lowMask;
    uint cls = 0;
    for (uint bb = 0; bb < P.nBlocks; bb++) {
        ulong mb = m & P.blockMask[bb];
        uint c = popcount(uint(mb & 0xffffffffu)) + popcount(uint(mb >> 32));
        cls += c * P.blockStride[bb];
    }
    ok = true;
    return cls;
}

inline void scatter(device atomic_uint *N, uint triSize, uint cidx,
                    device const uint *postStart, device const uint *postSeq, device const uint *postCnt,
                    uint u, uint v, uint /*unused*/ factor) {
    device atomic_uint *Nc = N + (ulong)cidx * triSize;
    for (uint a = postStart[u]; a < postStart[u + 1]; a++) {
        uint i = postSeq[a], ci = postCnt[a];
        for (uint b = postStart[v]; b < postStart[v + 1]; b++) {
            uint j = postSeq[b];
            uint idx = (i >= j) ? tri(i, j) : tri(j, i);
            // ordered window pairs (w in W_i, w' in W_j): the two terms cnt_u[i]cnt_v[j] and cnt_v[i]cnt_u[j]
            // both arrive at the symmetric cell {i,j} from this double loop; on the diagonal they coincide
            atomic_fetch_add_explicit(&Nc[idx], ((i == j) ? 2u : 1u) * ci * postCnt[b], memory_order_relaxed);
        }
    }
}

// 2-D grid: thread (r, c) handles row u = rowLo + r and columns v in [c*VT, min((c+1)*VT, u))
kernel void classifyRows(device const ulong *codes        [[buffer(0)]],
                         device const uint  *postStart    [[buffer(1)]],
                         device const uint  *postSeq      [[buffer(2)]],
                         device const uint  *postCnt      [[buffer(3)]],
                         device const ushort *classTable  [[buffer(4)]],   // total class index -> compact index or 0xFFFF
                         device atomic_uint *N            [[buffer(5)]],
                         device uint *overflow            [[buffer(6)]],   // triples (u, v, cidx)
                         device atomic_uint *ovCount      [[buffer(7)]],
                         constant Params &P               [[buffer(8)]],
                         uint2 tid [[thread_position_in_grid]])
{
    const uint VT = 512;
    uint u = P.rowLo + tid.x;
    if (u >= P.rowHi) return;
    uint v0 = tid.y * VT;
    if (v0 >= u && tid.y > 0) return;
    uint triSize = P.n * (P.n + 1) / 2;
    ulong cu = codes[u];
    uint pu = postStart[u + 1] - postStart[u];
    if (tid.y == 0 && P.selfPairs) {          // self pair (class 0 is always reachable and always compact index 0 of group 0)
        device atomic_uint *N0 = N;
        for (uint a = postStart[u]; a < postStart[u + 1]; a++)
            for (uint b = a; b < postStart[u + 1]; b++) {           // unordered posting pairs: cell {i,j} once
                uint i = postSeq[a], j = postSeq[b];
                uint idx = (i >= j) ? tri(i, j) : tri(j, i);
                atomic_fetch_add_explicit(&N0[idx], postCnt[a] * postCnt[b], memory_order_relaxed);
            }
    }
    uint v1 = min(v0 + VT, u);
    for (uint v = v0; v < v1; v++) {
        bool ok;
        uint cls = classify(cu, codes[v], P, ok);
        uint cidx = classTable[cls];
        if (cidx == 0xFFFF) continue;
        uint pv = postStart[v + 1] - postStart[v];
        if (pu * pv <= P.light) {
            scatter(N, triSize, cidx, postStart, postSeq, postCnt, u, v, 1u);
        } else {
            uint slot = atomic_fetch_add_explicit(ovCount, 1u, memory_order_relaxed);
            if (slot < P.ovCap) {                       // longer posting list first: it becomes the thread dimension
                bool sw = pv > pu;
                overflow[3 * slot] = sw ? v : u; overflow[3 * slot + 1] = sw ? u : v; overflow[3 * slot + 2] = cidx;
            }
        }
    }
}

// heavy pairs: one thread per (overflow entry, posting of u)
kernel void heavyPairs(device const uint *postStart   [[buffer(1)]],
                       device const uint *postSeq     [[buffer(2)]],
                       device const uint *postCnt     [[buffer(3)]],
                       device atomic_uint *N          [[buffer(5)]],
                       device const uint *overflow    [[buffer(6)]],
                       device const uint *entryStart  [[buffer(7)]],   // prefix sum of |post(u)| over entries
                       constant Params &P             [[buffer(8)]],
                       uint t [[thread_position_in_grid]])
{
    // binary search the entry that owns thread t
    uint lo = 0, hi = P.ovCap;                 // ovCap reused as number of entries here
    while (hi - lo > 1) { uint mid = (lo + hi) / 2; if (entryStart[mid] <= t) lo = mid; else hi = mid; }
    uint e = lo, off = t - entryStart[e];
    uint u = overflow[3 * e], v = overflow[3 * e + 1], cidx = overflow[3 * e + 2];
    uint a = postStart[u] + off;
    if (a >= postStart[u + 1]) return;
    uint triSize = P.n * (P.n + 1) / 2;
    device atomic_uint *Nc = N + (ulong)cidx * triSize;
    uint i = postSeq[a], ci = postCnt[a];
    for (uint b = postStart[v]; b < postStart[v + 1]; b++) {
        uint j = postSeq[b];
        uint idx = (i >= j) ? tri(i, j) : tri(j, i);
        atomic_fetch_add_explicit(&Nc[idx], ((i == j) ? 2u : 1u) * ci * postCnt[b], memory_order_relaxed);
    }
}
)MSL";

struct Params {
    uint32_t U, n, fieldBits, nBlocks, light, rowLo, rowHi, ovCap;
    uint64_t lowMask;
    uint64_t blockMask[8];
    uint32_t blockStride[8];
    uint32_t nClassesTotal;
    uint32_t selfPairs;
};

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ---------------------------------------------------------------- alphabets
struct Alphabet { std::string syms; std::string comp; };   // comp empty = no complement
static Alphabet parseSpec(const std::string &s) {
    Alphabet a;
    if (s == "dna") { a.syms = "ACGT"; a.comp = "TGCA"; }
    else if (s == "rna") { a.syms = "ACGU"; a.comp = "UGCA"; }
    else if (s == "protein") { a.syms = "ACDEFGHIKLMNPQRSTVWY"; }
    else if (!s.empty() && s[0] == '=') { a.syms = s.substr(1); }
    else { fprintf(stderr, "unsupported alphabet spec %s (prototype: dna, rna, protein, =SYMBOLS)\n", s.c_str()); exit(1); }
    return a;
}

// ---------------------------------------------------------------- coefficients
static double binom(int n, int k) { if (k < 0 || k > n) return 0; double r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }
// filter coefficient H(m) = (1/prod b) sum_{j<=k} e_j(w), w = b-1 at matches, -1 at mismatches
static double Hcoef(const std::vector<int> &blockB, const std::vector<int> &blockLen, const std::vector<int> &m, int k) {
    std::vector<double> w;
    for (size_t b = 0; b < blockB.size(); b++) {
        for (int i = 0; i < blockLen[b] - m[b]; i++) w.push_back(blockB[b] - 1);
        for (int i = 0; i < m[b]; i++) w.push_back(-1);
    }
    std::vector<double> e(w.size() + 1, 0.0); e[0] = 1.0;
    for (double x : w) for (int j = (int)w.size(); j >= 1; j--) e[j] += x * e[j - 1];
    double s = 0; for (int j = 0; j <= k && j < (int)e.size(); j++) s += e[j];
    double denom = 1; for (size_t b = 0; b < blockB.size(); b++) denom *= std::pow((double)blockB[b], blockLen[b]);
    return s / denom;
}

// ---------------------------------------------------------------- crc32 + gkmk writer
static uint32_t crc32buf(const unsigned char *p, size_t n) {
    static uint32_t table[256]; static bool init = false;
    if (!init) { for (uint32_t i = 0; i < 256; i++) { uint32_t c = i; for (int k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1); table[i] = c; } init = true; }
    uint32_t crc = ~0u; for (size_t i = 0; i < n; i++) crc = table[(crc ^ p[i]) & 0xFF] ^ (crc >> 8); return ~crc;
}
static void put32(std::vector<unsigned char> &o, uint32_t v) { for (int i = 0; i < 4; i++) o.push_back((v >> (8 * i)) & 0xFF); }
static void writeGkmk(const std::string &fn, const std::vector<float> &tri, int n, int L, int K, int d, int kind, const std::string &alph) {
    std::vector<unsigned char> h = {'G', 'K', 'M', 'K', 1, 0, 0, 0};
    put32(h, n); put32(h, n);
    for (int v : {L, K, d, kind, 4, 1, 0}) put32(h, (uint32_t)v);
    put32(h, (uint32_t)alph.size()); h.insert(h.end(), alph.begin(), alph.end());
    std::vector<unsigned char> payload(tri.size() * 4);
    memcpy(payload.data(), tri.data(), payload.size());
    uint32_t crc = crc32buf(payload.data(), payload.size());
    FILE *f = fopen(fn.c_str(), "wb");
    fwrite(h.data(), 1, h.size(), f); fwrite(payload.data(), 1, payload.size(), f);
    unsigned char t[4] = {(unsigned char)(crc & 0xFF), (unsigned char)((crc >> 8) & 0xFF), (unsigned char)((crc >> 16) & 0xFF), (unsigned char)(crc >> 24)};
    fwrite(t, 1, 4, f); fclose(f);
}

int main(int argc, char **argv) {
    int L = 10, K = 6, D = -1, kind = 0, light = 16, chunk = 4096, groups = 0;
    double budgetGB = 0;          // -M: memory for the profile buffer; classes are processed in groups that fit (0 = 60 % of maxBufferLength)
    bool addRC = true;
    std::string specs = "dna", outFn = "out.gkmk", inFn;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "-l") L = atoi(argv[++i]); else if (a == "-k") K = atoi(argv[++i]); else if (a == "-d") D = atoi(argv[++i]);
        else if (a == "-t") kind = atoi(argv[++i]); else if (a == "-A") specs = argv[++i]; else if (a == "-R") addRC = false;
        else if (a == "-o") outFn = argv[++i]; else if (a == "-L") light = atoi(argv[++i]); else if (a == "-C") chunk = atoi(argv[++i]);
        else if (a == "-M") budgetGB = atof(argv[++i]); else if (a == "-g") groups = atoi(argv[++i]);
        else inFn = a;
    }
    if (inFn.empty()) { fprintf(stderr, "usage: gkm_metal -l L -k K -d D [-t 0|2] -A specs [-R] [-o out.gkmk] [-M GB] [-g groups] [-L light] [-C chunk] in.mfa\n"); return 1; }
    if (kind != 0 && kind != 2) { fprintf(stderr, "prototype supports -t 0 (filter) and -t 2 (gkm)\n"); return 1; }

    // alphabets and positions
    std::vector<Alphabet> alph;
    { size_t p = 0; while (p <= specs.size()) { size_t q = specs.find(',', p); if (q == std::string::npos) q = specs.size(); alph.push_back(parseSpec(specs.substr(p, q - p))); p = q + 1; } }
    int T = (int)alph.size(), ell = T * L;
    if (!addRC && alph[0].comp.empty()) addRC = false;
    if (addRC && alph[0].comp.empty()) addRC = false;
    std::vector<int> posB(ell); for (int t = 0; t < T; t++) for (int i = 0; i < L; i++) posB[t * L + i] = (int)alph[t].syms.size();
    int f = 1; for (int b : posB) { int need = 1; while ((1 << need) < b) need++; f = std::max(f, need); }
    if (ell * f > 64) { fprintf(stderr, "prototype limit: l*fieldBits = %d > 64\n", ell * f); return 1; }
    // blocks in order of first appearance of each alphabet size
    std::vector<int> blockB, blockLen, posBlock(ell);
    for (int p = 0; p < ell; p++) { int bi = -1; for (size_t b = 0; b < blockB.size(); b++) if (blockB[b] == posB[p]) bi = (int)b; if (bi < 0) { bi = (int)blockB.size(); blockB.push_back(posB[p]); blockLen.push_back(0); } blockLen[bi]++; posBlock[p] = bi; }
    int nB = (int)blockB.size(); if (nB > 8) { fprintf(stderr, "prototype limit: > 8 blocks\n"); return 1; }
    std::vector<uint32_t> stride(nB); uint32_t nClassesTotal = 1; for (int b = 0; b < nB; b++) { stride[b] = nClassesTotal; nClassesTotal *= (blockLen[b] + 1); }
    int dmax = (D < 0) ? ell : D;
    // compact class table and coefficients
    std::vector<uint16_t> classTable(nClassesTotal, 0xFFFF); std::vector<double> coef; std::vector<std::vector<int>> classM;
    for (uint32_t c = 0; c < nClassesTotal; c++) {
        std::vector<int> m(nB); uint32_t r = c; int tot = 0; for (int b = 0; b < nB; b++) { m[b] = r % (blockLen[b] + 1); r /= (blockLen[b] + 1); tot += m[b]; }
        if (tot > dmax) continue;
        classTable[c] = (uint16_t)coef.size(); classM.push_back(m);
        coef.push_back(kind == 0 ? Hcoef(blockB, blockLen, m, K) : binom(ell - tot, K));
    }
    int R = (int)coef.size();
    fprintf(stderr, "l=%d (T=%d x L=%d), k=%d, d=%d, %d blocks, %u classes, %d reachable, %d bits/position, RC=%d\n", ell, T, L, K, D, nB, nClassesTotal, R, f, (int)addRC);

    // ---- read .mfa, build distinct l-mers with postings
    double t0 = now();
    std::vector<uint64_t> codes; std::vector<uint32_t> postStart(1, 0), postSeq, postCnt;
    std::unordered_map<uint64_t, uint32_t> ids;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> post;
    std::vector<std::string> names;
    {
        FILE *fh = fopen(inFn.c_str(), "r"); if (!fh) { perror(inFn.c_str()); return 1; }
        char *line = NULL; size_t cap = 0; ssize_t len;
        std::vector<std::string> tracks; std::string name; int nseq = 0;
        auto flush = [&]() {
            if (name.empty()) return;
            if ((int)tracks.size() != T) { fprintf(stderr, "record %s: %zu tracks, expected %d\n", name.c_str(), tracks.size(), T); exit(1); }
            int n0 = (int)tracks[0].size();
            std::vector<std::vector<int>> code(T, std::vector<int>(n0, -1));
            for (int t = 0; t < T; t++) for (int i = 0; i < n0 && i < (int)tracks[t].size(); i++) { size_t q = alph[t].syms.find(tracks[t][i]); code[t][i] = (q == std::string::npos) ? -1 : (int)q; }
            auto addWord = [&](const std::vector<int> &w) {
                uint64_t c = 0; for (int p = 0; p < ell; p++) c |= (uint64_t)w[p] << (p * f);
                auto it = ids.find(c); uint32_t id;
                if (it == ids.end()) { id = (uint32_t)codes.size(); ids[c] = id; codes.push_back(c); post.emplace_back(); } else id = it->second;
                auto &pl = post[id];
                if (!pl.empty() && pl.back().first == (uint32_t)nseq) pl.back().second++; else pl.push_back({(uint32_t)nseq, 1});
            };
            for (int i = 0; i + L <= n0; i++) {
                std::vector<int> w(ell); bool ok = true;
                for (int t = 0; t < T && ok; t++) for (int q = 0; q < L; q++) { int v = code[t][i + q]; if (v < 0) { ok = false; break; } w[t * L + q] = v; }
                if (!ok) continue;
                addWord(w);
                if (addRC) {
                    std::vector<int> r(ell);
                    for (int q = 0; q < L; q++) { int v = w[q]; char cc = alph[0].comp[v]; r[L - 1 - q] = (int)alph[0].syms.find(cc); }
                    for (int t = 1; t < T; t++) for (int q = 0; q < L; q++) r[t * L + (L - 1 - q)] = w[t * L + q];
                    addWord(r);
                }
            }
            names.push_back(name); nseq++; name.clear(); tracks.clear();
        };
        while ((len = getline(&line, &cap, fh)) > 0) {
            while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
            if (len == 0) continue;
            if (line[0] == '>') { flush(); name = std::string(line + 1); size_t sp = name.find_first_of(" \t"); if (sp != std::string::npos) name = name.substr(0, sp); }
            else tracks.push_back(line);
        }
        flush(); free(line); fclose(fh);
    }
    uint32_t U = (uint32_t)codes.size(), n = (uint32_t)names.size();
    for (uint32_t u = 0; u < U; u++) { for (auto &pc : post[u]) { postSeq.push_back(pc.first); postCnt.push_back(pc.second); } postStart.push_back((uint32_t)postSeq.size()); }
    uint64_t triSize = (uint64_t)n * (n + 1) / 2;
    fprintf(stderr, "%u sequences, %u distinct l-mers, %zu postings; prep %.2f s\n", n, U, postSeq.size(), now() - t0);

    // ---- Metal
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev) { fprintf(stderr, "no Metal device\n"); return 1; }
    fprintf(stderr, "device: %s, maxBufferLength %.1f GB\n", dev.name.UTF8String, dev.maxBufferLength / 1e9);
    NSError *err = nil;
    MTLCompileOptions *opt = [MTLCompileOptions new];
    id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:kShader] options:opt error:&err];
    if (!lib) { fprintf(stderr, "shader compile failed: %s\n", err.localizedDescription.UTF8String); return 1; }
    id<MTLComputePipelineState> psoRows = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"classifyRows"] error:&err];
    id<MTLComputePipelineState> psoHeavy = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"heavyPairs"] error:&err];
    if (!psoRows || !psoHeavy) { fprintf(stderr, "pipeline failed: %s\n", err.localizedDescription.UTF8String); return 1; }
    id<MTLCommandQueue> queue = [dev newCommandQueue];
    auto mk = [&](const void *p, size_t bytes) { return [dev newBufferWithBytes:p length:std::max<size_t>(bytes, 16) options:MTLResourceStorageModeShared]; };
    id<MTLBuffer> bCodes = mk(codes.data(), codes.size() * 8), bStart = mk(postStart.data(), postStart.size() * 4),
                  bSeq = mk(postSeq.data(), postSeq.size() * 4), bCnt = mk(postCnt.data(), postCnt.size() * 4),
                  bTable = mk(classTable.data(), classTable.size() * 2);
    // ---- class groups: the profile buffer holds Rg <= R classes at a time; the pair space is classified once per
    // group (classification is cheap, the scatter is what costs), and each group's counts are folded into the
    // kernel on the CPU before the next group. Compact class indices are contiguous per group and folded in
    // increasing order, so the result is identical to the single-group computation.
    if (budgetGB <= 0) budgetGB = dev.maxBufferLength * 0.6 / 1e9;
    int Rg = (int)std::min<double>(R, std::floor(budgetGB * 1e9 / ((double)triSize * 4)));
    if (groups > 0) Rg = std::max(1, (R + groups - 1) / groups);
    if (Rg < 1) { fprintf(stderr, "even one class (%.2f GB) exceeds the memory budget %.1f GB\n", triSize * 4 / 1e9, budgetGB); return 1; }
    int G = (R + Rg - 1) / Rg;
    size_t Nbytes = (size_t)Rg * triSize * 4;
    fprintf(stderr, "profile buffer: %d of %d classes x %llu entries = %.2f GB, %d class group%s\n", Rg, R, (unsigned long long)triSize, Nbytes / 1e9, G, G == 1 ? "" : "s");
    if (Nbytes > dev.maxBufferLength) { fprintf(stderr, "profile buffer exceeds maxBufferLength\n"); return 1; }
    id<MTLBuffer> bN = [dev newBufferWithLength:Nbytes options:MTLResourceStorageModeShared];
    uint32_t ovCap = 1u << 27;
    id<MTLBuffer> bOv = [dev newBufferWithLength:(size_t)ovCap * 12 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bOvCount = [dev newBufferWithLength:16 options:MTLResourceStorageModeShared];
    memset(bOvCount.contents, 0, 16);
    Params P{}; P.U = U; P.n = n; P.fieldBits = f; P.nBlocks = nB; P.light = light; P.ovCap = ovCap; P.nClassesTotal = nClassesTotal;
    for (int p = 0; p < ell; p++) { P.lowMask |= 1ull << (p * f); P.blockMask[posBlock[p]] |= 1ull << (p * f); }
    for (int b = 0; b < nB; b++) P.blockStride[b] = stride[b];

    if ((uint64_t)chunk * U > ovCap) chunk = std::max<int>(64, (int)(ovCap / U));   // every pair of a chunk may be deferred
    fprintf(stderr, "row chunk %d, light threshold %d\n", chunk, light);
    std::vector<double> Kt(triSize, 0.0);
    unsigned long long tot = 0;
    double t1 = now(); double heavySec = 0, foldSec = 0; unsigned long long nOvTotal = 0;
    for (int g = 0; g < G; g++) {
    int c0 = g * Rg, c1 = std::min(R, c0 + Rg);
    std::vector<uint16_t> groupTable(nClassesTotal, 0xFFFF);
    for (uint32_t c = 0; c < nClassesTotal; c++) if (classTable[c] != 0xFFFF && (int)classTable[c] >= c0 && (int)classTable[c] < c1) groupTable[c] = (uint16_t)(classTable[c] - c0);
    id<MTLBuffer> bTable = mk(groupTable.data(), groupTable.size() * 2);
    memset(bN.contents, 0, Nbytes);
    P.selfPairs = (g == 0) ? 1u : 0u;
    for (uint32_t lo = 0; lo < U; lo += chunk) {
        uint32_t hi = std::min(U, lo + (uint32_t)chunk);
        P.rowLo = lo; P.rowHi = hi;
        id<MTLBuffer> bP = mk(&P, sizeof(P));
        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:psoRows];
        [enc setBuffer:bCodes offset:0 atIndex:0]; [enc setBuffer:bStart offset:0 atIndex:1]; [enc setBuffer:bSeq offset:0 atIndex:2];
        [enc setBuffer:bCnt offset:0 atIndex:3]; [enc setBuffer:bTable offset:0 atIndex:4]; [enc setBuffer:bN offset:0 atIndex:5];
        [enc setBuffer:bOv offset:0 atIndex:6]; [enc setBuffer:bOvCount offset:0 atIndex:7]; [enc setBuffer:bP offset:0 atIndex:8];
        uint32_t ctiles = (hi + 511) / 512;
        [enc dispatchThreads:MTLSizeMake(hi - lo, ctiles, 1) threadsPerThreadgroup:MTLSizeMake(32, 8, 1)];
        [enc endEncoding]; [cb commit]; [cb waitUntilCompleted];
        if (cb.error) { fprintf(stderr, "GPU error: %s\n", cb.error.localizedDescription.UTF8String); return 1; }
        uint32_t nOv = *(uint32_t *)bOvCount.contents;
        if (nOv > ovCap) { fprintf(stderr, "overflow list capacity exceeded in chunk %u (%u > %u): rerun with smaller -C or larger -L\n", lo, nOv, ovCap); return 1; }
        if (nOv > 0) {
            double th = now();
            std::vector<uint32_t> entryStart(nOv + 1, 0); const uint32_t *ov = (const uint32_t *)bOv.contents;
            for (uint32_t e = 0; e < nOv; e++) { uint32_t u = ov[3 * e]; entryStart[e + 1] = entryStart[e] + (postStart[u + 1] - postStart[u]); }
            uint32_t total = entryStart[nOv];
            id<MTLBuffer> bES = mk(entryStart.data(), entryStart.size() * 4);
            Params P2 = P; P2.ovCap = nOv; id<MTLBuffer> bP2 = mk(&P2, sizeof(P2));
            id<MTLCommandBuffer> cb2 = [queue commandBuffer]; id<MTLComputeCommandEncoder> enc2 = [cb2 computeCommandEncoder];
            [enc2 setComputePipelineState:psoHeavy];
            [enc2 setBuffer:bStart offset:0 atIndex:1]; [enc2 setBuffer:bSeq offset:0 atIndex:2]; [enc2 setBuffer:bCnt offset:0 atIndex:3];
            [enc2 setBuffer:bN offset:0 atIndex:5]; [enc2 setBuffer:bOv offset:0 atIndex:6]; [enc2 setBuffer:bES offset:0 atIndex:7]; [enc2 setBuffer:bP2 offset:0 atIndex:8];
            NSUInteger tg2 = std::min<NSUInteger>(psoHeavy.maxTotalThreadsPerThreadgroup, 256);
            [enc2 dispatchThreads:MTLSizeMake(total, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg2, 1, 1)];
            [enc2 endEncoding]; [cb2 commit]; [cb2 waitUntilCompleted];
            if (cb2.error) { fprintf(stderr, "GPU error (heavy): %s\n", cb2.error.localizedDescription.UTF8String); return 1; }
            heavySec += now() - th; nOvTotal += nOv;
            *(uint32_t *)bOvCount.contents = 0;
        }
    }
    // fold this group's counts into the kernel (double; classes in increasing compact order, as before)
    double tf = now();
    const uint32_t *Nbuf = (const uint32_t *)bN.contents;
    for (int c = c0; c < c1; c++) { const uint32_t *Nc = Nbuf + (uint64_t)(c - c0) * triSize; double w = coef[c]; unsigned long long sc = 0; for (uint64_t x = 0; x < triSize; x++) { Kt[x] += w * Nc[x]; sc += Nc[x]; } tot += sc; }
    foldSec += now() - tf;
    }
    double t2 = now();
    fprintf(stderr, "classify + scatter: %.2f s (%.3g pairs x %d group%s, %.3g pairs/s); heavy pairs: %llu in %.2f s; fold %.2f s\n", t2 - t1 - foldSec, (double)U * U / 2, G, G == 1 ? "" : "s", (double)U * U / 2 * G / (t2 - t1 - foldSec), nOvTotal, heavySec, foldSec);
    // ---- normalise, write
    double t3 = now();
    std::vector<double> diag(n); for (uint32_t i = 0; i < n; i++) diag[i] = Kt[(uint64_t)i * (i + 1) / 2 + i];
    std::vector<float> out(triSize);
    for (uint32_t i = 0; i < n; i++) for (uint32_t j = 0; j <= i; j++) { double dn = diag[i] * diag[j]; out[(uint64_t)i * (i + 1) / 2 + j] = (i == j) ? 1.0f : (dn < 1e-50 ? 0.0f : (float)(Kt[(uint64_t)i * (i + 1) / 2 + j] / std::sqrt(dn))); }
    writeGkmk(outFn, out, n, L, K, D, kind, specs);
    fprintf(stderr, "normalise + write: %.2f s; total %.2f s\n", now() - t3, now() - t0);
    fprintf(stderr, "profile entries counted (lower triangle, ordered pairs): %llu\n", tot);
    return 0;
}
