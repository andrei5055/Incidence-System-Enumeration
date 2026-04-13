#ifndef K16P1F_H
#define K16P1F_H
// Decoupled from TripleSys.h to prevent full rebuilds.

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include "kBase.h"

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#define VECTOR_CALL __vectorcall
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define VECTOR_CALL
#endif
#define K16_Use_rdtsc 0
#define K16_N 16
#define K16_MATCH 15
#define K16_FIXED 3
#define K16_SEARCH 12
#define K16_M_MAX 60000
#define K16_SORT_INPUT 0
#define K16_USE_ROOT_SORT        1

#define K16_DISABLE_IS_CANONICAL_R45_CHECK 0
#define K16_BENCHMARK_EXIT     0

#define K16_EDGE_PRUNE_ENABLED 1
#define K16_EDGE_PRUNE_START   9 // 9
#define K16_EDGE_PRUNE_END     12 // 12
#define K16_PRUNE_THRESHOLD_HIGH 20000
#define K16_PRUNE_THRESHOLD_LOW  500 // 500

#define K16_USE_MRV              1
#define K16_MRV_DEPTH_LIMIT      8 // Only reorder for depths < 8 
#define K16_MRV_EARLY_EXIT       1
#define K16_MRV_EARLY_EXIT_THRESHOLD 1

struct Mask256_C { uint64_t m[2] = { 0, 0 }; };

class K16P1F : public KBase<Mask256_C> {
public:

    struct alignas(64) PackedAdj { 
        Mask256_C edge_mask;
        uint8_t adj[K16_N];
        __m128i v_e2o;
        __m128i v_o2e;
    };

    struct FastSortedFactor : FastSortedFactorBase<32> {};
    struct FastRowTriplet { FastSortedFactor r[3]; };

    struct Factor {
        uint8_t adj[K16_N];
        uint8_t src[K16_N];
        Mask256_C edge_mask;
        FastSortedFactor fs;
        __m128i v_e2o;
        __m128i v_o2e;
    };

    // Bitset state for candidates in the search pool
    #define K16_WORDS (((K16_M_MAX + 255) / 256 * 4) + 4)
    struct State { alignas(64) uint64_t bits[K16_WORDS]; };

    struct SearchContext {
        alignas(64) State pool[K16_SEARCH + 1];
        Mask256_C used_edges;
        uint8_t slots[K16_SEARCH];
        int r4_idx = 0;
        int r5_idx = 0;
        int root_idx = 0;
        uint16_t mrv_counts[K16_SEARCH + 1][K16_SEARCH]; // [depth][slot_idx]
        bool counts_valid[K16_SEARCH + 1];
    };

    struct CycleMeter {
        uint64_t total_cycles = 0;
        uint64_t prop_cycles = 0;
        uint64_t mask_cycles = 0;
        uint64_t calls = 0;
        uint64_t masks_checked = 0;
        uint64_t props_done = 0;
    };

    struct alignas(64) ThreadLocalBuffers : LocalBufers<Mask256_C> {

        State local_edge_presence[120];
        // Search Phase Memory
        SearchContext local_ctx;
        uint64_t r5_row_mask_in_s4[K16_SEARCH][K16_WORDS];
        int dirty_s4_words[K16_SEARCH][K16_WORDS];
        int dirty_s4_count[K16_SEARCH];
        SubSamplePlanItem sub_sampling_plans[K16_SEARCH][K16_WORDS];
        int sub_sampling_plan_sizes[K16_SEARCH];

        // Constant pointers for the search (to keep SearchContext small)
        const RowRange* ranges = nullptr;
        const uint64_t* adj = nullptr;
        const size_t* offsets = nullptr;
        const Mask256_C* edge_masks = nullptr;
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
    struct CycleUnion { int cycles[16][16]; int lens[16]; int count = 0; };

    K16P1F(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr = NULL, bool bPrint = true);
    ~K16P1F() {}

    bool addRow(int iRow, const unsigned char* source);
    void solve(int mode = 0);

private:
    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr);

    void* cbClass = NULL;
    int call_counter = 0;
    int kThreads;
    ResultCallback resultCallback;
    Factor fixedRows[K16_FIXED];
    PackedAdj fixed_packed[K16_FIXED];
    Mask256_C fixedEdgesMask;
    int edge_id_table[K16_N][K16_N];
    uint8_t edge_to_u[128];
    uint8_t edge_to_v[128];

    std::vector<Factor> global_pool;
    std::vector<PackedAdj> packed_pool;

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
    static FORCE_INLINE bool check_cycle_simd(__m128i v_e2o1, __m128i v_o2e2, __m128i v_id) {
        __m128i v_sigma = _mm_shuffle_epi8(v_o2e2, v_e2o1);
        // We want a single cycle of length 8 in the edge permutation for K16 Hamiltonian.
        uint32_t mask1 = _mm_movemask_epi8(_mm_cmpeq_epi8(v_sigma, v_id)) & 0x00FF;
        if (mask1) return false;
        
        __m128i v_s2 = _mm_shuffle_epi8(v_sigma, v_sigma);
        uint32_t mask2 = _mm_movemask_epi8(_mm_cmpeq_epi8(v_s2, v_id)) & 0x00FF;
        if (mask2) return false;
        
        __m128i v_s4 = _mm_shuffle_epi8(v_s2, v_s2);
        uint32_t mask4 = _mm_movemask_epi8(_mm_cmpeq_epi8(v_s4, v_id)) & 0x00FF;
        if (mask4) return false;
        
        __m128i v_s8 = _mm_shuffle_epi8(v_s4, v_s4);
        uint32_t mask8 = _mm_movemask_epi8(_mm_cmpeq_epi8(v_s8, v_id)) & 0x00FF;
        return mask8 == 0x00FF;
    }

    FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if ((p1.edge_mask.m[0] & p2.edge_mask.m[0]) || (p1.edge_mask.m[1] & p2.edge_mask.m[1])) 
            return false;
        return is_perfect_scalar(p1.adj, p2.adj);
    }
    bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2);

    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    void get_transformations(const Factor& fi, const Factor& fj, TransInfo& info);
    void apply_perm_16(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj);
    FastSortedFactor get_fast_sorted(const uint8_t* adj);
    
    bool is_canonical(int r4_fid, int r5_fid, const TransInfo* r4_dependent_trans);
    bool is_canonical_stab(int r5_fid, const Permutation* stab, int stab_count);

    CycleUnion target_cu;
    TransInfo fixed_trans[3];
    FastSortedFactor r1_can, r2_can;

    // Statistics and Reporting
    static bool compare_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)   { return memcmp(a.pairs, b.pairs, 16) < 0; }
    static bool equal_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)     { return memcmp(a.pairs, b.pairs, 16) == 0; }
    static bool compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b);

    // Threading Helpers
    void parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size = 1);
    std::mutex result_mutex;
};

#endif
