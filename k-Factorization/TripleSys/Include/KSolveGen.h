#ifndef KSOLVEGEN_H
#define KSOLVEGEN_H

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include <cstring>
#include <mutex>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include <thread>
#include <memory>
#include "kBase.h"

extern int g_printCount;
#define SAFE_PRINTF(...) { \
    if (g_printCount++ >= 100) { \
        printf("\nERROR: Print limit (100) reached. Exiting to prevent console flood.\n"); \
        exit(1); \
    } \
    printf(__VA_ARGS__); \
}

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#define VECTOR_CALL __vectorcall
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define VECTOR_CALL
#endif

struct alignas(64) Mask512 {
    uint64_t m[8] = { 0 }; 
    FORCE_INLINE void clear() { memset(m, 0, sizeof(m)); }
    FORCE_INLINE void set_bit(int bit) { if (bit >= 0 && bit < 512) m[bit >> 6] |= (1ULL << (bit & 63)); }
    
    FORCE_INLINE void union_with(const Mask512& other) {
        __m256i* dest = (__m256i*)m;
        const __m256i* src = (const __m256i*)other.m;
        dest[0] = _mm256_or_si256(dest[0], _mm256_loadu_si256(&src[0]));
        dest[1] = _mm256_or_si256(dest[1], _mm256_loadu_si256(&src[1]));
    }
    FORCE_INLINE void and_with(const Mask512& other) {
        __m256i* dest = (__m256i*)m;
        const __m256i* src = (const __m256i*)other.m;
        dest[0] = _mm256_and_si256(dest[0], _mm256_loadu_si256(&src[0]));
        dest[1] = _mm256_and_si256(dest[1], _mm256_loadu_si256(&src[1]));
    }
    FORCE_INLINE void and_not(const Mask512& other) {
        __m256i* dest = (__m256i*)m;
        const __m256i* src = (const __m256i*)other.m;
        dest[0] = _mm256_andnot_si256(_mm256_loadu_si256(&src[0]), dest[0]);
        dest[1] = _mm256_andnot_si256(_mm256_loadu_si256(&src[1]), dest[1]);
    }
    FORCE_INLINE void remove_mask(const Mask512& other) { and_not(other); }

    FORCE_INLINE bool is_empty() const {
        for (int i = 0; i < 8; i++) if (m[i]) return false;
        return true;
    }
    FORCE_INLINE bool has_overlap(const Mask512& other) const {
        __m256i t0 = _mm256_loadu_si256((const __m256i*)&m[0]);
        __m256i t1 = _mm256_loadu_si256((const __m256i*)&m[4]);
        __m256i o0 = _mm256_loadu_si256((const __m256i*)&other.m[0]);
        __m256i o1 = _mm256_loadu_si256((const __m256i*)&other.m[4]);
        return !_mm256_testz_si256(t0, o0) || !_mm256_testz_si256(t1, o1);
    }
    static FORCE_INLINE bool any_missing(const Mask512& target, const Mask512& available) {
        __m256i t0 = _mm256_loadu_si256((const __m256i*)&target.m[0]);
        __m256i a0 = _mm256_loadu_si256((const __m256i*)&available.m[0]);
        __m256i res0 = _mm256_andnot_si256(a0, t0);
        if (!_mm256_testz_si256(res0, res0)) return true;
        __m256i t1 = _mm256_loadu_si256((const __m256i*)&target.m[4]);
        __m256i a1 = _mm256_loadu_si256((const __m256i*)&available.m[4]);
        __m256i res1 = _mm256_andnot_si256(a1, t1);
        return !_mm256_testz_si256(res1, res1);
    }
};

#define KGEN_M_MAX 250000
#define KGEN_WORDS (((KGEN_M_MAX + 255) / 256 * 4) + 4)

class KSolveGen : public KBase<Mask512> {
public:
    struct State { alignas(64) uint64_t bits[KGEN_WORDS]; };
    struct Permutation { alignas(32) uint8_t p[32]; uint8_t p_inv[32]; __m256i pv; __m256i pinvv; };
    struct TransInfo { int count; Permutation perms[512]; };
    struct RootPair { int r4_v, r5_v; };
    struct CompRange { int start_word, end_word; };
    struct CycleUnion { int cycles[32][32]; int lens[32]; int count = 0; };
    
    struct SubSamplePlanItem {
        uint16_t global_word_idx;
        uint64_t mask;
        uint16_t local_word_off;
        uint8_t local_bit_shift;
        uint8_t bits_pushed;
        uint8_t split;
    };

    struct RootGroup {
        int start_idx;
        int end_idx;
        int r4_v;
    };
    
    struct Factor {
        uint8_t adj[32];
        uint8_t src[32];
        Mask512 edge_mask;
        FastSortedFactorBase<32> fs;
    };

    struct PackedAdj {
        Mask512 edge_mask;
        uint8_t adj[32];
        __m128i v_e2o, v_o2e;
    };

    struct SearchContext {
        alignas(64) State pool[32];

        // MRV State
        uint8_t slots[32];
        uint16_t mrv_counts[32][32]; // [depth][slot]
        bool counts_valid[32];

        int r4_idx = 0;
        int r5_idx = 0;
        Mask512 used_edges[32];
    };

    struct ThreadLocalBuffers {
        SearchContext local_ctx;
        TransInfo r5_trans;
        TransInfo r4_trans[3];
        
        // Local Matrix State
        uint64_t* local_adj_aligned = nullptr;
        size_t local_adj_size = 0;
        std::vector<size_t> local_offsets;
        std::vector<Mask512> local_edge_masks;
        std::vector<int> s_to_f;
        std::vector<int> local_to_global;
        std::vector<int> g_to_l; 

        SubSamplePlanItem sub_sampling_plans[32][4096];
        int sub_sampling_plan_sizes[32];

        // Specialized Transposed Index for K16-style SIMD Pruning
        // Maps edge_id -> bitset of valid candidates covering that edge
        alignas(64) State local_edge_presence[512]; 

        CompRange local_ranges[32];
        uint64_t* adj_ptr = nullptr;
        const CompRange* ranges_ptr = nullptr;
        const size_t* offsets_ptr = nullptr;

        ThreadLocalBuffers() { 
            memset(&local_ctx, 0, sizeof(local_ctx)); 
            memset(sub_sampling_plan_sizes, 0, sizeof(sub_sampling_plan_sizes));
        }
        ~ThreadLocalBuffers() {
            if (local_adj_aligned) _aligned_free(local_adj_aligned);
        }

        void reserve_aligned_adj(size_t size) {
            if (size > local_adj_size) {
                if (local_adj_aligned) _aligned_free(local_adj_aligned);
                local_adj_aligned = (uint64_t*)_aligned_malloc(size * sizeof(uint64_t), 64);
                local_adj_size = size;
            }
            if (local_adj_aligned) memset(local_adj_aligned, 0, size * sizeof(uint64_t));
        }
    };

    struct SortItem { int global_id, row_id, degree; };

    KSolveGen(const FactorParams& factParam, bool bipartite, const unsigned char* targetCycles, int cycleCount,
              int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, 
              ResultCallback callback, void* cbClassPtr, bool bPrint,
              const uint8_t* v4Row = nullptr, const uint8_t* v4Val = nullptr);
    
    virtual ~KSolveGen() {}

    bool addRow(int iRow, const unsigned char* source) override;

    /**
     * @brief Performs internal candidate generation and initial filtering.
     * Clears existing pools and generates all valid rows for all slots using optimized 
     * on-the-fly routines, then applies the Iterative Global Filter (Arc Consistency).
     */
    void preSolve() override;

    /**
     * @brief Returns the number of initial candidates in a specific slot.
     * @param slot The slot index (0 = Row 4, 1 = Row 5, etc.).
     */
    int getCandidateCount(int slot) override;

    /**
     * @brief Retrieves a candidate result in the standard pairs format.
     * @param slot The slot index.
     * @param index The candidate index (0 to count-1).
     * @param out_src Output buffer (m_n bytes) to receive the pair-format data (e.g. 0,7, 1,4...).
     * @return True if candidate exists and is active; False if out of bounds or marked deleted.
     */
    bool getCandidate(int slot, int index, unsigned char* out_src) override;

    /**
     * @brief Marks a candidate for removal from the final search pool.
     * Note: Does not shift indices immediately; removal occurs inside solve().
     * @param slot The slot index.
     * @param index The candidate index to delete.
     */
    void deleteCandidate(int slot, int index) override;

    /**
     * @brief Starts the final clique search for solutions.
     * @param mode 0 = Manual Capture Mode: Filters captured rows (addRow) using AC before search.
     *             1 = Presolve Mode: Compresses slots (removes deleteCandidate entries) before search.
     */
    void solve(int mode = 0) override;

private:
    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr);
    void runArcConsistency();
    
    bool m_bipartite;
    bool m_useOptimized;
    int m_cycleCount;
    uint8_t m_targetCycles[32];
    int m_fixedCount = 3;
    Mask512 total_edges_mask;
    bool m_enableCoveragePruning = false;
    bool m_enableUnitProp = false;
    bool m_generateOnTheFly = false;

    void* cbClass = NULL;
    int kThreads;
    ResultCallback resultCallback;
    Factor fixedRows[8];
    PackedAdj fixed_packed[8];
    Mask512 fixedEdgesMask;
    int edge_id_table[32][32];
    
    TransInfo fixed_trans[3];
    bool m_presolveDone = false;
    FastSortedFactorBase<32> r1_can, r2_can;

    std::vector<Factor> global_pool;
    std::vector<PackedAdj> packed_pool;
    std::map<std::vector<uint8_t>, int> f_map;

    std::vector<std::unique_ptr<ThreadLocalBuffers>> thread_buffers;
    std::vector<CompRange> row_ranges;
    std::vector<size_t> factor_offsets;
    std::vector<uint64_t> adj_matrix;
    std::vector<int> search_to_factor;
    std::vector<std::vector<unsigned char>> results_to_sort;

    // Parity & Reference Tracking
    std::vector<std::set<std::vector<uint8_t>>> m_referenceSets;
    bool m_parityReportDone = false;

    void VECTOR_CALL internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx, ThreadLocalBuffers* buf);
    bool checkCycleStructure(const uint8_t* adj1, const uint8_t* adj2);
    bool check_cycle_bipartite_simd(const __m128i& v_e2o_i, const __m128i& v_o2e_j);
    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    void generateSlotCandidates(int slotIndex, int startIdx, int endIdx);
    void generateRecursive(int depth, uint64_t usedMask, uint8_t* currentAdj, int slotIndex, int startPairIdx, int endPairIdx, int rowK);
    PackedAdj pack_factor_adj(const uint8_t* adj);
    
    int m_bipartiteFixedPoints[32];
    FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if (p1.edge_mask.has_overlap(p2.edge_mask)) return false;
        if (m_bipartite && m_nPlayers <= 32) return check_cycle_bipartite_simd(p1.v_e2o, p2.v_o2e);
        return checkCycleStructure(p1.adj, p2.adj);
    }

    void get_transformations(const Factor& fi, const Factor& fj, TransInfo& info);
    FastSortedFactorBase<32> get_fast_sorted(const uint8_t* adj);
    void apply_perm_generic(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj);
    bool is_canonical(int r4_fid, int r5_fid, const TransInfo* r4_dependent_trans, ThreadLocalBuffers* buf);
    bool is_canonical_stab(int r5_fid, const Permutation* stab, int stab_count);

    std::mutex result_mutex;
    std::mutex pool_mutex;
    std::mutex stdout_mutex;
    std::atomic<int> m_resultsFound{ 0 };
    std::atomic<int> m_rootsDone{ 0 };
    std::chrono::steady_clock::time_point last_report_time;
    std::chrono::steady_clock::time_point start_time;

    uint8_t m_v4Row[4];
    uint8_t m_v4Val[4];

    void parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size = 1);
    void diagnostic_printout(int total_roots);
};

#endif
