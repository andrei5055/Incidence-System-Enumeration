// =============================================================================
// k16gpu.cu -- CUDA engine for the K16 P1F full enumeration (K16P1F, UseKSolve=1).
//
// One WARP per root (canonical (r4,r5) pair). The kernel runs the S4-space
// iterative DFS validated on CPU as K16P1F::iter_solve_s4 (k16p1f.cpp): slot-MRV
// with cached counts, per-slot wipeout propagation, lazy edge-coverage prune with
// forced moves. The pool bitsets are LANE-SLICED: lane L owns words {L, L+32,
// L+64, ...} of every bitset, so global reads of adjacency rows are coalesced and
// all reductions are warp-synchronous (__reduce_*_sync / __shfl_*). No shared
// memory, no inter-warp communication; hits are appended with atomicAdd.
//
// Host API (extern "C", loaded at runtime via LoadLibrary from k16p1f.cpp, env
// K16_GPU=1): k16gpu_init / k16gpu_solve_group / k16gpu_close. Each host thread
// gets its own stream (cudaStreamPerThread), so the caller's existing
// parallel_for over r4-groups gives concurrent uploads+launches for free --
// ~12 groups x ~700 roots in flight saturates an RTX 5090.
//
// Build (no GPU needed to compile; PTX JITs on the target driver):
//   build_k16gpu.bat        (nvcc -arch=compute_90 -> K16GPU.dll)
// On the RTX 5090 box with CUDA 12.8+/13.x, prefer native code: -arch=sm_120.
// =============================================================================
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>

#define NSEARCH        12          // rows to place per root (r4,r5 occupy depths 0,1)
#define NEDGEBITS      120         // edges of K16
#define MRV_LIMIT      8           // K16_MRV_DEPTH_LIMIT
#define MRV_EARLY      1           // K16_MRV_EARLY_EXIT_THRESHOLD
#define PRUNE_START    9           // K16_EDGE_PRUNE_START
#define PRUNE_END      12          // K16_EDGE_PRUNE_END
#define PRUNE_THRESH   500         // K16_PRUNE_THRESHOLD_LOW
#define FORCED_THRESH  50
#define PERLANE        24          // max ceil(s4w4/32): supports s4w4 <= 768 (S4K <= ~49k)
#define WARPS_PER_BLOCK 4

struct K16GpuHit { int root; unsigned short v[NSEARCH]; };   // v[0..9] = chosen S4 indices, depths 2..11

struct K16GpuParams {
    const unsigned long long* s4_adj;     // S4K rows x s4w4 words
    const unsigned long long* presence;   // 120 x s4w4
    const unsigned long long* emasks;     // S4K x 2 (edge mask per candidate)
    const int*                ranges;     // 12 x 2 ints (start_word, end_word per slot)
    const unsigned long long* rootPools;  // nRoots x s4w4 (initial pool per root)
    const unsigned long long* rootUsed;   // nRoots x 2   (initial used_edges per root)
    K16GpuHit*                hits;
    int*                      hitCount;
    unsigned long long*       nodesTotal;
    int S4K, s4w4, nRoots, maxHits;
};

// warp-reduce helpers (full-warp, warp-uniform results)
__device__ __forceinline__ int warpSumI(int v) {
    for (int o = 16; o; o >>= 1) v += __shfl_xor_sync(0xFFFFFFFFu, v, o);
    return v;
}
__device__ __forceinline__ int warpMinI(int v) {
    for (int o = 16; o; o >>= 1) { int t = __shfl_xor_sync(0xFFFFFFFFu, v, o); v = t < v ? t : v; }
    return v;
}
__device__ __forceinline__ unsigned long long warpOrULL(unsigned long long v) {
    for (int o = 16; o; o >>= 1) v |= __shfl_xor_sync(0xFFFFFFFFu, v, o);
    return v;
}

__global__ void __launch_bounds__(WARPS_PER_BLOCK * 32, 8) k16Kernel(K16GpuParams p) {
    const int wid = blockIdx.x * WARPS_PER_BLOCK + (threadIdx.x >> 5);
    if (wid >= p.nRoots) return;
    const int lane = threadIdx.x & 31;
    const int W = p.s4w4;
    const int PL = (W + 31) >> 5;             // lane-owned words

    // ---- lane-sliced pool stack (word k*32+lane of depth d) ------------------
    unsigned long long pool[NSEARCH + 1][PERLANE];
    {
        const unsigned long long* rp = p.rootPools + (size_t)wid * W;
        for (int k = 0; k < PL; k++) { int w = (k << 5) + lane; pool[2][k] = (w < W) ? rp[w] : 0ULL; }
        // deeper levels: zero the padding words once — propagation only ever writes
        // slot windows, and the full-range prune scans must not see garbage
        for (int dd = 3; dd < NSEARCH; dd++)
            for (int k = 0; k < PL; k++) pool[dd][k] = 0ULL;
    }
    unsigned long long used0 = p.rootUsed[2 * (size_t)wid];
    unsigned long long used1 = p.rootUsed[2 * (size_t)wid + 1];

    // warp-uniform state, replicated in every lane
    unsigned char slots[NSEARCH];
    for (int s = 0; s < NSEARCH; s++) slots[s] = (unsigned char)s;
    unsigned short mrv[NSEARCH + 1][NSEARCH];
    bool cvalid[NSEARCH + 1];
    for (int i = 0; i <= NSEARCH; i++) cvalid[i] = false;

    int  f_bi[NSEARCH + 1];
    int  f_next[NSEARCH + 1];                 // next candidate bit index to try (within slot bit range)
    int  f_endbit[NSEARCH + 1];               // end of candidate bit range (exclusive)
    bool f_pushed[NSEARCH + 1];
    unsigned long long f_u0[NSEARCH + 1], f_u1[NSEARCH + 1];
    unsigned short vstack[NSEARCH + 1];
    unsigned long long nodes = 0;

    int d = 2; bool enter = true;
    while (d >= 2) {
        unsigned long long (&P)[PERLANE] = pool[d];

        if (enter) {
            enter = false;
            nodes++;
            f_pushed[d] = false;
            f_bi[d] = d;

            // ---- MRV slot selection ------------------------------------------
            if (d < MRV_LIMIT) {
                int best_idx = d, min_count = 1000000;
                if (cvalid[d]) {
                    for (int i2 = d; i2 < NSEARCH; i2++) {
                        int c = mrv[d][slots[i2]];
                        if (c < min_count) { min_count = c; best_idx = i2; if (min_count <= MRV_EARLY) break; }
                    }
                    cvalid[d] = false;
                } else {
                    for (int i2 = d; i2 < NSEARCH; i2++) {
                        int slot = slots[i2];
                        int ws = p.ranges[2 * slot], we = p.ranges[2 * slot + 1];
                        int c = 0;
                        for (int k = 0; k < PL; k++) {
                            int w = (k << 5) + lane;
                            if (w >= ws && w < we) c += __popcll(P[k]);
                        }
                        c = warpSumI(c);
                        if (c < min_count) { min_count = c; best_idx = i2; if (min_count <= MRV_EARLY) break; }
                    }
                }
                f_bi[d] = best_idx;
                unsigned char t = slots[d]; slots[d] = slots[best_idx]; slots[best_idx] = t;
                if (min_count == 0) { d--; continue; }
            }

            const int s = slots[d];
            f_next[d] = p.ranges[2 * s] * 64;
            f_endbit[d] = p.ranges[2 * s + 1] * 64;

            // ---- lazy edge-coverage prune + forced move ----------------------
            if (d >= PRUNE_START && d <= PRUNE_END) {
                int total_remaining = 0;
                for (int k = 0; k < PL && total_remaining <= PRUNE_THRESH; k++)
                    total_remaining += warpSumI(__popcll(P[k]));
                if (total_remaining <= PRUNE_THRESH) {
                    // coverage: OR of edge masks over all available candidates
                    unsigned long long u0 = 0, u1 = 0;
                    for (int k = 0; k < PL; k++) {
                        unsigned long long word = P[k];
                        while (word) {
                            int b = __ffsll((long long)word) - 1; word &= word - 1;
                            int v = ((k << 5) + lane) * 64 + b;
                            u0 |= p.emasks[2 * (size_t)v];
                            u1 |= p.emasks[2 * (size_t)v + 1];
                        }
                    }
                    u0 = warpOrULL(u0); u1 = warpOrULL(u1);
                    unsigned long long need0 = ~used0;
                    unsigned long long need1 = (~used1) & 0x00FFFFFFFFFFFFFFULL;
                    if ((need0 & ~u0) || (need1 & ~u1)) { d--; continue; }   // some edge unreachable

                    if (total_remaining <= FORCED_THRESH) {
                        int forced_v = -1; bool conflict = false;
                        for (int e = 0; e < NEDGEBITS; e++) {
                            bool is_needed = (e < 64) ? ((need0 >> e) & 1) : ((need1 >> (e - 64)) & 1);
                            if (!is_needed) continue;
                            const unsigned long long* ep = p.presence + (size_t)e * W;
                            int cnt = 0, myv = 0x7fffffff;
                            for (int k = 0; k < PL; k++) {
                                int w = (k << 5) + lane;
                                unsigned long long word = (w < W) ? (P[k] & ep[w]) : 0ULL;
                                if (word) {
                                    cnt += __popcll(word);
                                    int b = __ffsll((long long)word) - 1;
                                    int v = w * 64 + b; if (v < myv) myv = v;
                                }
                            }
                            int tot = warpSumI(cnt);
                            if (tot == 1) {
                                int v1 = warpMinI(myv);
                                if (v1 / 64 >= p.ranges[2 * s] && v1 / 64 < p.ranges[2 * s + 1]) {
                                    if (forced_v != -1 && forced_v != v1) { conflict = true; break; }
                                    forced_v = v1;
                                }
                            }
                        }
                        if (conflict) { d--; continue; }
                        if (forced_v != -1) { f_next[d] = forced_v; f_endbit[d] = forced_v + 1; }
                    }
                }
            }
        }

        // ---- undo previously committed child --------------------------------
        if (f_pushed[d]) { used0 = f_u0[d]; used1 = f_u1[d]; f_pushed[d] = false; }

        // ---- next candidate: smallest set bit in pool[d] within [f_next, f_endbit)
        int v;
        {
            const int from = f_next[d], endb = f_endbit[d];
            int myv = 0x7fffffff;
            for (int k = 0; k < PL; k++) {
                int w = (k << 5) + lane;
                int wb = w * 64;
                if (wb + 64 <= from || wb >= endb) continue;
                unsigned long long word = P[k];
                if (wb < from) word &= ~0ULL << (from - wb);
                if (wb + 64 > endb) word &= (endb - wb == 64) ? ~0ULL : ((1ULL << (endb - wb)) - 1);
                if (word) { int b = __ffsll((long long)word) - 1; int cand = wb + b; if (cand < myv) myv = cand; }
            }
            v = warpMinI(myv);
            if (v == 0x7fffffff) {
                unsigned char t = slots[d]; slots[d] = slots[f_bi[d]]; slots[f_bi[d]] = t;
                d--;
                continue;
            }
            f_next[d] = v + 1;
        }

        const unsigned long long fm0 = p.emasks[2 * (size_t)v];
        const unsigned long long fm1 = p.emasks[2 * (size_t)v + 1];
        if ((used0 & fm0) | (used1 & fm1)) continue;    // safety (mirrors CPU)

        if (d + 1 < NSEARCH) {
            // propagate pool & per-slot counts for the child
            const unsigned long long* rel = p.s4_adj + (size_t)v * W;
            unsigned long long (&NP)[PERLANE] = pool[d + 1];
            bool wipe = false;
            for (int n = d + 1; n < NSEARCH; n++) {
                int ns = slots[n];
                int ws = p.ranges[2 * ns], we = p.ranges[2 * ns + 1];
                int c = 0;
                for (int k = 0; k < PL; k++) {
                    int w = (k << 5) + lane;
                    if (w >= ws && w < we) {
                        unsigned long long r = P[k] & rel[w];
                        NP[k] = r;
                        c += __popcll(r);
                    }
                }
                c = warpSumI(c);
                mrv[d + 1][ns] = (unsigned short)c;
                if (c == 0) { wipe = true; break; }
            }
            if (wipe) continue;
            cvalid[d + 1] = true;

            f_u0[d] = used0; f_u1[d] = used1;
            used0 |= fm0; used1 |= fm1;
            vstack[d] = (unsigned short)v;
            f_pushed[d] = true;
            d++; enter = true;
        } else {
            // complete: report hit (host re-checks canonicity)
            vstack[d] = (unsigned short)v;
            nodes++;
            if (lane == 0) {
                int idx = atomicAdd(p.hitCount, 1);
                if (idx < p.maxHits) {
                    K16GpuHit& h = p.hits[idx];
                    h.root = wid;
                    for (int j = 0; j < NSEARCH - 2; j++) h.v[j] = vstack[j + 2];
                    h.v[NSEARCH - 2] = 0; h.v[NSEARCH - 1] = 0;
                }
            }
        }
    }
    if (lane == 0) atomicAdd(p.nodesTotal, nodes);
}

// ============================== host API =====================================
#define CUCHECK(x) do { cudaError_t _e = (x); if (_e != cudaSuccess) { \
    fprintf(stderr, "[K16GPU] %s failed: %s (%s:%d)\n", #x, cudaGetErrorString(_e), __FILE__, __LINE__); return -1; } } while (0)

extern "C" __declspec(dllexport) int k16gpu_init(int device) {
    CUCHECK(cudaSetDevice(device));
    cudaDeviceProp prop;
    CUCHECK(cudaGetDeviceProperties(&prop, device));
    printf("[K16GPU] device %d: %s, %d SMs, %.1f GB\n", device, prop.name,
           prop.multiProcessorCount, prop.totalGlobalMem / 1073741824.0);
    fflush(stdout);
    return 0;
}

// Solve every root of one r4-group. Synchronous w.r.t. the calling host thread
// (cudaStreamPerThread), so concurrent host threads pipeline naturally.
extern "C" __declspec(dllexport) int k16gpu_solve_group(
    const unsigned long long* h_s4_adj, size_t adjWords,
    const unsigned long long* h_presence,
    const unsigned long long* h_emasks,
    const int* h_ranges,
    const unsigned long long* h_rootPools,
    const unsigned long long* h_rootUsed,
    int S4K, int s4w4, int nRoots,
    K16GpuHit* h_hits, int maxHits, int* out_nHits,
    unsigned long long* out_nodes)
{
    if (s4w4 > PERLANE * 32) { fprintf(stderr, "[K16GPU] s4w4=%d exceeds kernel cap %d\n", s4w4, PERLANE * 32); return -2; }
    cudaStream_t st = cudaStreamPerThread;

    unsigned long long* d_adj = nullptr; unsigned long long* d_pres = nullptr;
    unsigned long long* d_em = nullptr; int* d_rng = nullptr;
    unsigned long long* d_rp = nullptr; unsigned long long* d_ru = nullptr;
    K16GpuHit* d_hits = nullptr; int* d_hc = nullptr; unsigned long long* d_nodes = nullptr;

    CUCHECK(cudaMallocAsync(&d_adj, adjWords * 8, st));
    CUCHECK(cudaMallocAsync(&d_pres, (size_t)120 * s4w4 * 8, st));
    CUCHECK(cudaMallocAsync(&d_em, (size_t)S4K * 16, st));
    CUCHECK(cudaMallocAsync(&d_rng, 12 * 2 * sizeof(int), st));
    CUCHECK(cudaMallocAsync(&d_rp, (size_t)nRoots * s4w4 * 8, st));
    CUCHECK(cudaMallocAsync(&d_ru, (size_t)nRoots * 16, st));
    CUCHECK(cudaMallocAsync(&d_hits, (size_t)maxHits * sizeof(K16GpuHit), st));
    CUCHECK(cudaMallocAsync(&d_hc, sizeof(int), st));
    CUCHECK(cudaMallocAsync(&d_nodes, 8, st));

    CUCHECK(cudaMemcpyAsync(d_adj, h_s4_adj, adjWords * 8, cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemcpyAsync(d_pres, h_presence, (size_t)120 * s4w4 * 8, cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemcpyAsync(d_em, h_emasks, (size_t)S4K * 16, cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemcpyAsync(d_rng, h_ranges, 12 * 2 * sizeof(int), cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemcpyAsync(d_rp, h_rootPools, (size_t)nRoots * s4w4 * 8, cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemcpyAsync(d_ru, h_rootUsed, (size_t)nRoots * 16, cudaMemcpyHostToDevice, st));
    CUCHECK(cudaMemsetAsync(d_hc, 0, sizeof(int), st));
    CUCHECK(cudaMemsetAsync(d_nodes, 0, 8, st));

    K16GpuParams prm;
    prm.s4_adj = d_adj; prm.presence = d_pres; prm.emasks = d_em; prm.ranges = d_rng;
    prm.rootPools = d_rp; prm.rootUsed = d_ru;
    prm.hits = d_hits; prm.hitCount = d_hc; prm.nodesTotal = d_nodes;
    prm.S4K = S4K; prm.s4w4 = s4w4; prm.nRoots = nRoots; prm.maxHits = maxHits;

    int blocks = (nRoots + WARPS_PER_BLOCK - 1) / WARPS_PER_BLOCK;
    k16Kernel<<<blocks, WARPS_PER_BLOCK * 32, 0, st>>>(prm);
    CUCHECK(cudaGetLastError());

    int hc = 0; unsigned long long nn = 0;
    CUCHECK(cudaMemcpyAsync(&hc, d_hc, sizeof(int), cudaMemcpyDeviceToHost, st));
    CUCHECK(cudaMemcpyAsync(&nn, d_nodes, 8, cudaMemcpyDeviceToHost, st));
    CUCHECK(cudaStreamSynchronize(st));
    if (hc > maxHits) { fprintf(stderr, "[K16GPU] hit buffer overflow (%d > %d)\n", hc, maxHits); return -3; }
    if (hc > 0)
        CUCHECK(cudaMemcpy(h_hits, d_hits, (size_t)hc * sizeof(K16GpuHit), cudaMemcpyDeviceToHost));
    *out_nHits = hc;
    *out_nodes = nn;

    cudaFreeAsync(d_adj, st); cudaFreeAsync(d_pres, st); cudaFreeAsync(d_em, st);
    cudaFreeAsync(d_rng, st); cudaFreeAsync(d_rp, st); cudaFreeAsync(d_ru, st);
    cudaFreeAsync(d_hits, st); cudaFreeAsync(d_hc, st); cudaFreeAsync(d_nodes, st);
    return 0;
}

extern "C" __declspec(dllexport) int k16gpu_close() {
    cudaDeviceSynchronize();
    return 0;
}
