#ifndef K18A2_H
#define K18A2_H

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include "kBase.h"
#include <unordered_map>

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#define VECTOR_CALL __vectorcall
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define VECTOR_CALL
#endif
#define K18_Use_rdtsc 0
#define K18_N 18
#define K18_MATCH 17
#define K18_FIXED 3
#define K18_SEARCH 14
#define K18_M_MAX 900000
#define K18_SORT_INPUT 0
#define K18_USE_ROOT_SORT        1

#define K18_DISABLE_IS_CANONICAL_R45_CHECK 0
#define K18_BENCHMARK_EXIT     0

#define K18_EDGE_PRUNE_ENABLED 1
#define K18_EDGE_PRUNE_START   9 
#define K18_EDGE_PRUNE_END     12 
#define K18_PRUNE_THRESHOLD_HIGH 20000
#define K18_PRUNE_THRESHOLD_LOW  500 

#define K18_USE_MRV              1
#define K18_MRV_DEPTH_LIMIT      8 
#define K18_MRV_EARLY_EXIT       1
#define K18_MRV_EARLY_EXIT_THRESHOLD 1

struct Mask18_C { uint64_t m[4] = { 0, 0, 0, 0 }; };

class K18A2 : public KBase<Mask18_C> {
public:

    struct alignas(64) PackedAdj { 
        Mask18_C edge_mask;
        uint8_t adj[K18_N];
    };

    struct FastSortedFactor : FastSortedFactorBase<32> {};
    struct FastRowTriplet { FastSortedFactor r[3]; };

    struct Factor {
        uint8_t adj[K18_N];
        uint8_t src[K18_N];
        Mask18_C edge_mask;
        FastSortedFactor fs;
    };

    // Bitset state for candidates in the search pool
    #define K18_WORDS (((K18_M_MAX + 255) / 256 * 4) + 4)
    struct State { alignas(64) uint64_t bits[K18_WORDS]; };

    struct SearchContext {
        alignas(64) State pool[K18_SEARCH + 1];
        Mask18_C used_edges;
        uint8_t slots[K18_SEARCH];
        int r4_idx = 0;
        int r5_idx = 0;
        int root_idx = 0;
        uint16_t mrv_counts[K18_SEARCH + 1][K18_SEARCH]; // [depth][slot_idx]
        bool counts_valid[K18_SEARCH + 1];
    };

    struct CycleMeter {
        uint64_t total_cycles = 0;
        uint64_t prop_cycles = 0;
        uint64_t mask_cycles = 0;
        uint64_t calls = 0;
        uint64_t masks_checked = 0;
        uint64_t props_done = 0;
    };

    struct alignas(64) ThreadLocalBuffers : LocalBufers<Mask18_C> {

        State local_edge_presence[153];
        // Search Phase Memory
        SearchContext local_ctx;
        uint64_t r5_row_mask_in_s4[K18_SEARCH][K18_WORDS];
        int dirty_s4_words[K18_SEARCH][K18_WORDS];
        int dirty_s4_count[K18_SEARCH];
        SubSamplePlanItem sub_sampling_plans[K18_SEARCH][K18_WORDS];
        int sub_sampling_plan_sizes[K18_SEARCH];
        int s4_bit_cursor_table[K18_SEARCH];

        // Constant pointers for the search (to keep SearchContext small)
        const RowRange* ranges = nullptr;
        const uint64_t* adj = nullptr;
        const size_t* offsets = nullptr;
        const Mask18_C* edge_masks = nullptr;
        const int* s_to_f = nullptr;
        double local_compression = 0.0;

        CycleMeter meter;

        ThreadLocalBuffers() {
            memset(dirty_s4_count, 0, sizeof(dirty_s4_count));
            memset(sub_sampling_plan_sizes, 0, sizeof(sub_sampling_plan_sizes));
            memset(dirty_s4_words, 0, sizeof(dirty_s4_words));
            memset(r5_row_mask_in_s4, 0, sizeof(r5_row_mask_in_s4));
        }
    };

    struct TransInfo { Permutation perms[512]; int count = 0; };
    struct CycleUnion { int cycles[18][18]; int lens[18]; int count = 0; };

    K18A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr = NULL, bool bPrint = true);
    ~K18A2() {}

    bool addRow(int iRow, const unsigned char* source);
    void solve(int mode = 0);

private:
    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr);

    void* cbClass = NULL;
    int call_counter = 0;
    int kThreads;
    ResultCallback resultCallback;
    Factor fixedRows[K18_FIXED];
    PackedAdj fixed_packed[K18_FIXED];
    Mask18_C fixedEdgesMask;
    int edge_id_table[K18_N][K18_N];
    uint8_t edge_to_u[256];
    uint8_t edge_to_v[256];

    std::vector<Factor> global_pool;
    std::vector<PackedAdj> packed_pool;

    struct VectorHash {
        size_t operator()(const std::vector<uint8_t>& v) const {
            size_t hash = 17;
            for (uint8_t x : v) {
                hash = hash * 31 + x;
            }
            return hash;
        }
    };
    std::unordered_map<std::vector<uint8_t>, int, VectorHash> f_map_unordered;

    int restart_index;

    std::vector<std::unique_ptr<ThreadLocalBuffers>> thread_buffers;
    int thread_row4[256];
    int thread_row5[256];
    int thread_root_idx[256];

    bool bTimeSet = false;

    // Internal Help Methods
    void VECTOR_CALL internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx, ThreadLocalBuffers* buf);
    void diagnostic_printout(double current_compr);
    PackedAdj pack_factor_adj(const uint8_t* adj);

    FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if ((p1.edge_mask.m[0] & p2.edge_mask.m[0]) || 
            (p1.edge_mask.m[1] & p2.edge_mask.m[1]) ||
            (p1.edge_mask.m[2] & p2.edge_mask.m[2])) 
            return false;
        return is_perfect_scalar(p1.adj, p2.adj);
    }
    bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2);

    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    void get_transformations(const Factor& fi, const Factor& fj, TransInfo& info);
    void get_transformations_general(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl, TransInfo& info);
    void apply_perm_18(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj);
    FastSortedFactor get_fast_sorted(const uint8_t* adj);
    
    CycleUnion target_cu;
    TransInfo fixed_trans[3];
    FastSortedFactor r1_can, r2_can;

    // Statistics and Reporting
    bool compare_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)   { return memcmp(a.pairs, b.pairs, 18) < 0; }
    bool equal_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)     { return memcmp(a.pairs, b.pairs, 18) == 0; }
    bool compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b);

    std::mutex result_mutex;
};

#endif
