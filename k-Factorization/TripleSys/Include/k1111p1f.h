#ifndef K1111P1F_H
#define K1111P1F_H
// Decoupled from TripleSys.h to prevent full rebuilds.

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include <thread>
#include <mutex>
#include <functional>
#include <map>
#include <set>
#include <memory>
#include <cstring>

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif

#define K1111_DISABLE_IS_CANONICAL_R45_CHECK 0
#define K1111_N 22
#define K1111_MATCH 11
#define K1111_FIXED 3
#define K1111_SEARCH (K1111_MATCH - K1111_FIXED)
#define K1111_M_MAX 250000
#define K1111_SORT_INPUT 0
#define K1111_USE_SIMD_CYCLE_CHECK 1
#define K1111_USE_ROOT_SORT        1
#define K1111_USE_BRANCH_FREE_EDGE 0
#define K1111_USE_STATIC_V_ID      1

typedef bool (*ResultCallback)(const void* cbClass, const unsigned char* results, int r4, int r5, int mode);

class K1111P1F {
public:
    struct Mask256 { uint64_t m[3] = { 0, 0, 0 }; };

    struct alignas(64) PackedAdj { 
        uint64_t m[4] = { 0, 0, 0, 0 }; 
        __m128i v_e2o;
        __m128i v_o2e;
    };

    struct alignas(64) FastSortedFactor { uint8_t pairs[64]; };
    struct FastRowTriplet { FastSortedFactor r[3]; };

    struct Factor {
        uint8_t adj[K1111_N];
        uint8_t src[K1111_N];
        Mask256 edge_mask;
        __m128i v_e2o;
        __m128i v_o2e;
        FastSortedFactor fs;
    };

    // Bitset state for candidates in the search pool
    #define K1111_WORDS (((K1111_M_MAX + 255) / 256 * 4) + 4)
    struct State { alignas(64) uint64_t bits[K1111_WORDS]; };

    struct RowRange { int start_word; int end_word; };

    struct SubSamplePlanItem {
        uint16_t s4_word_idx;
        uint16_t l_word_off;
        uint64_t mask;
        uint8_t l_bit_shift;
        uint8_t bits_pushed;
        uint8_t split; 
    };

    struct SearchContext {
        alignas(64) State pool[K1111_SEARCH + 1];
        Mask256 used_edges;

        const RowRange* ranges;
        const uint64_t* adj;
        size_t adj_size;
        const size_t* offsets;
        const Mask256* edge_masks;
        const int* s_to_f;

        int r4_idx = 0;
        int r5_idx = 0;
        double local_compression = 1.0;
    };

    struct alignas(64) ThreadLocalBuffers {
        std::vector<int> g_to_l;
        std::vector<int> local_to_global;
        std::vector<int> cl;
        std::vector<RowRange> local_ranges;
        std::vector<uint64_t> local_adj;
        std::vector<size_t> local_offsets;
        std::vector<Mask256> local_edge_masks;
        std::vector<int> local_s_to_f;
        struct ActiveWord { int global_word_idx; uint64_t mask; int local_bit_start; };
        std::vector<ActiveWord> active_words;

        // Phase 2: Semi-local cache for a fixed r4
        std::vector<int> s4_to_global;
        std::vector<int> global_to_s4;
        std::vector<uint64_t> s4_adj;
        std::vector<size_t> s4_offsets;
        std::vector<uint8_t> s4_canonical; 
        
        // Search Phase Memory
        SearchContext local_ctx;
        int sub_sampling_plan_sizes[K1111_SEARCH];
        uint64_t r5_row_mask_in_s4[K1111_SEARCH][K1111_WORDS];
        int dirty_s4_words[K1111_SEARCH][K1111_WORDS];
        int dirty_s4_count[K1111_SEARCH];
        SubSamplePlanItem sub_sampling_plans[K1111_SEARCH][K1111_WORDS];

        ThreadLocalBuffers() {
            memset(dirty_s4_count, 0, sizeof(dirty_s4_count));
            memset(sub_sampling_plan_sizes, 0, sizeof(sub_sampling_plan_sizes));
            memset(dirty_s4_words, 0, sizeof(dirty_s4_words));
            memset(r5_row_mask_in_s4, 0, sizeof(r5_row_mask_in_s4));
        }
    };

    // Canonicity structures
    struct alignas(32) Permutation22 { 
        uint8_t p[32]; 
        uint8_t p_inv[32]; 
        __m256i pv;
        __m256i pinvv;
        int swaps = 0;
        int padding[7]; // Ensure sizeof=160 (multiple of 32)
    };
    struct TransInfo { Permutation22 perms[512]; int count = 0; };
    struct CycleUnion { int cycles[24][24]; int lens[24]; int count = 0; };

    K1111P1F();
    ~K1111P1F();

    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr = NULL, bool bPrint = true);
    bool addRow(int iRow, const unsigned char* source);
    void solve();

private:
    void* cbClass = NULL;
    int kThreads;
    int fixed3RowsIndex = 0;
    ResultCallback resultCallback;
    bool bPrint = true;
    Factor fixedRows[K1111_FIXED];
    PackedAdj fixed_packed[K1111_FIXED];
    Mask256 fixedEdgesMask;
    int edge_id_table[K1111_N][K1111_N];

    std::vector<Factor> global_pool;
    std::vector<PackedAdj> packed_pool;
    std::vector<int> temp_slot_ids[K1111_SEARCH];
    std::map<std::vector<uint8_t>, int> f_map;
    std::set<int> row_factor_ids[K1111_SEARCH];

    // Search Structures
    std::vector<int> search_to_factor;
    std::vector<uint64_t> adj_matrix;
    std::vector<Mask256> factor_edge_masks;
    std::vector<RowRange> row_ranges;
    std::vector<size_t> factor_offsets;
    int roots_done[256] = { 0 };
    std::atomic<int> num_results{ 0 };
    std::atomic<int> num_notCanon{ 0 };
    int total_roots = 0;
    int first_root = -1;
    std::chrono::steady_clock::time_point solve_start_time;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point last_print_time;

    // Thread Local State
    std::vector<std::unique_ptr<ThreadLocalBuffers>> thread_buffers;
    int thread_row4[256];
    int thread_row5[256];
    int thread_root_idx[256];
    int diag_cnt = 0;
    std::vector<std::vector<unsigned char>> results_to_sort;

    // Internal Help Methods
    void internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx);
    void diagnostic_printout(double current_compr);
    PackedAdj pack_factor_adj(const uint8_t* adj);
    static FORCE_INLINE bool check_cycle_simd(__m128i v_e2o1, __m128i v_o2e2, __m128i v_id) {
        // v_e2o and v_o2e contain (v-offset)/2 values to map vertices to indices 0..10
        __m128i v_sigma = _mm_shuffle_epi8(v_o2e2, v_e2o1); // σ acts on 11 elements
        
        // P1F for n=11 requires a single cycle of length 11.
        // This is true if σ^k has no fixed points for k=1..5.
        __m128i v_s2 = _mm_shuffle_epi8(v_sigma, v_sigma);
        __m128i v_s3 = _mm_shuffle_epi8(v_s2, v_sigma);
        __m128i v_s4 = _mm_shuffle_epi8(v_s2, v_s2);
        __m128i v_s5 = _mm_shuffle_epi8(v_s4, v_sigma);

        __m128i v_fix = _mm_cmpeq_epi8(v_sigma, v_id);
        v_fix = _mm_or_si128(v_fix, _mm_cmpeq_epi8(v_s2, v_id));
        v_fix = _mm_or_si128(v_fix, _mm_cmpeq_epi8(v_s3, v_id));
        v_fix = _mm_or_si128(v_fix, _mm_cmpeq_epi8(v_s4, v_id));
        v_fix = _mm_or_si128(v_fix, _mm_cmpeq_epi8(v_s5, v_id));

        // Use mask for indices 0..10
        uint32_t mask = _mm_movemask_epi8(v_fix) & 0x07FF;
        return mask == 0;
    }

    FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if (!_mm256_testz_si256(_mm256_load_si256((const __m256i*)&p1), _mm256_load_si256((const __m256i*)&p2))) return false;
#if K1111_USE_SIMD_CYCLE_CHECK
        const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        return check_cycle_simd(_mm_load_si128(&p1.v_e2o), _mm_load_si128(&p2.v_o2e), v_id);
#else
        return is_perfect_scalar(p1.adj, p2.adj);
#endif
    }
    bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2);

    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    void get_transformations(const Factor& fi, const Factor& fj, TransInfo& info);
    void apply_perm_22(const uint8_t* src_adj, const Permutation22& perm, uint8_t* dst_adj);
    FastSortedFactor get_fast_sorted(const uint8_t* adj);
    
    // Canonicity Modes
#define K1111_USE_STABILIZER_CANON 1 // 0=Brute Force Mode, 1=New Automorphism Mode

    bool is_canonical(int r4_fid, int r5_fid, const TransInfo* r4_dependent_trans);
    bool is_canonical_stab(int r5_fid, const Permutation22* stab, int stab_count);

    CycleUnion target_cu;
    TransInfo fixed_trans[3];
    FastSortedFactor r1_can, r2_can;
    __m128i gfs_table_adj[256];
    __m128i gfs_table_idx[256];

    // Statistics and Reporting
    std::atomic<int> g_total_pairs_checked{ 0 };
    std::atomic<int> g_rejected_edges{ 0 };
    std::atomic<int> g_rejected_cycles{ 0 };
    std::atomic<int> g_rejected_canon{ 0 };
    std::atomic<int> g_total_rejected{ 0 };
    std::atomic<int> g_total_kept{ 0 };
    
    std::atomic<int> min_stab_size{ 1000000 };
    std::atomic<int> max_stab_size{ 0 };

    static bool compare_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b);
    static bool equal_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b);
    static bool compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b);

    // Threading Helpers
    void parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size = 1);
    std::mutex result_mutex;
};

#endif
