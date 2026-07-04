#ifndef K16A2_H
#define K16A2_H

#include <vector>
#include <array>
#include <cstdint>
#include <functional>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include "kBase.h"
#include <unordered_map>
#include <new>
#include <malloc.h>
#include <crtdbg.h>

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#define VECTOR_CALL __vectorcall
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define VECTOR_CALL
#endif
#define K16_Use_rdtsc 0
#define K16_SORT_INPUT 0
#define K16_USE_ROOT_SORT        1

#define K16_DISABLE_IS_CANONICAL_R45_CHECK 0
#define K16_BENCHMARK_EXIT     0

#define K16_EDGE_PRUNE_ENABLED 1
#define K16_EDGE_PRUNE_START   9 
#define K16_EDGE_PRUNE_END     12 
#define K16_PRUNE_THRESHOLD_HIGH 20000
#define K16_PRUNE_THRESHOLD_LOW  500 

#define K16_USE_MRV              1
#define K16_MRV_DEPTH_LIMIT      8 
#define K16_MRV_EARLY_EXIT       1
#define K16_MRV_EARLY_EXIT_THRESHOLD 1

// Representative-method classifier (k16A2rep.cpp; N=16 sibling of k18a2rep.cpp).
// When K16_USE_REP_METHOD == 1, K16A2::runExhaustiveSearch runs the unseeded
// orderly-generation classifier instead of the cyclic search. Default order is
// K16_REP_ORDER, overridable at runtime via env REP_ORDER (single) or REP_ORDERS
// (comma list). K16_REP_SYMCAP bounds centralizer enumeration. Set to 0 to
// restore the original cyclic K16A2 behavior.
#define K16_USE_REP_METHOD 1          // 1 = run representative method; 0 = cyclic search
#define K16_REP_ORDER      3          // default automorphism order to classify (env REP_ORDER overrides)
#define K16_REP_SYMCAP     2000000LL  // enumerate C(alpha) only when |C(alpha)| <= this

struct Mask16_C { uint64_t m[4] = { 0, 0, 0, 0 }; };

#pragma push_macro("new")
#undef new

class K16A2 : public KBase<Mask16_C> {
public:
    // ---- N-derived constants (single knob: NP). K16 port differs only in NP. ----
    static constexpr int NP = 16;                      // number of players
    static constexpr int NM = NP - 1;                  // factors / match rows (was NM=17)
    static constexpr int NFIXED = 3;                   // fixed starter rows (was NFIXED)
    static constexpr int NSEARCH = NP - NFIXED - 1;    // search depth (was NSEARCH=14)
    static constexpr int M_MAX = 60000;               // pool cap (size tuning; differs for K16)
    static constexpr int NHALF = NP / 2;               // half (was literal 9)
    static constexpr int NEDGES = NP * (NP - 1) / 2;   // complement edges (was literal NEDGES)
    static constexpr int NWORDS = (((M_MAX + 255) / 256 * 4) + 4); // bitset words (was NWORDS)
    void* operator new(size_t size) {
        void* p = _aligned_malloc(size, 64);
        if (!p) throw std::bad_alloc();
        return p;
    }
    void operator delete(void* p) {
        _aligned_free(p);
    }
#if defined(_DEBUG) && defined(_MSC_VER)
    void* operator new(size_t size, int blockType, const char* filename, int linenumber) {
        void* p = _aligned_malloc_dbg(size, 64, filename, linenumber);
        if (!p) throw std::bad_alloc();
        return p;
    }
    void operator delete(void* p, int blockType, const char* filename, int linenumber) {
        _aligned_free(p);
    }
#endif

    struct alignas(64) PackedAdj { 
        Mask16_C edge_mask;
        uint8_t adj[NP];
    };

    struct Factor {
        uint8_t adj[NP];
        uint8_t src[NP];
        Mask16_C edge_mask;
    };

    // Bitset state for candidates in the search pool
    struct State { alignas(64) uint64_t bits[NWORDS]; };

    struct SearchContext {
        alignas(64) State pool[NSEARCH + 1];
        Mask16_C used_edges;
        uint8_t slots[NSEARCH];
        int r4_idx = 0;
        int r5_idx = 0;
        int root_idx = 0;
        uint16_t mrv_counts[NSEARCH + 1][NSEARCH]; // [depth][slot_idx]
        bool counts_valid[NSEARCH + 1];
    };

    struct CycleMeter {
        uint64_t total_cycles = 0;
        uint64_t prop_cycles = 0;
        uint64_t mask_cycles = 0;
        uint64_t calls = 0;
        uint64_t masks_checked = 0;
        uint64_t props_done = 0;
    };

    struct alignas(64) ThreadLocalBuffers : LocalBufers<Mask16_C> {

        State local_edge_presence[NEDGES];
        // Search Phase Memory
        SearchContext local_ctx;
        uint64_t r5_row_mask_in_s4[NSEARCH][NWORDS];
        int dirty_s4_words[NSEARCH][NWORDS];
        int dirty_s4_count[NSEARCH];
        SubSamplePlanItem sub_sampling_plans[NSEARCH][NWORDS];
        int sub_sampling_plan_sizes[NSEARCH];
        int s4_bit_cursor_table[NSEARCH];

        // Constant pointers for the search (to keep SearchContext small)
        const RowRange* ranges = nullptr;
        const uint64_t* adj = nullptr;
        const size_t* offsets = nullptr;
        const Mask16_C* edge_masks = nullptr;
        const int* s_to_f = nullptr;
        double local_compression = 0.0;

        CycleMeter meter;

        ThreadLocalBuffers() {
            memset(dirty_s4_count, 0, sizeof(dirty_s4_count));
            memset(sub_sampling_plan_sizes, 0, sizeof(sub_sampling_plan_sizes));
            memset(dirty_s4_words, 0, sizeof(dirty_s4_words));
            memset(r5_row_mask_in_s4, 0, sizeof(r5_row_mask_in_s4));
        }
        void* operator new(size_t size) {
            void* p = _aligned_malloc(size, 64);
            if (!p) throw std::bad_alloc();
            return p;
        }
        void operator delete(void* p) {
            _aligned_free(p);
        }
#if defined(_DEBUG) && defined(_MSC_VER)
        void* operator new(size_t size, int blockType, const char* filename, int linenumber) {
            void* p = _aligned_malloc_dbg(size, 64, filename, linenumber);
            if (!p) throw std::bad_alloc();
            return p;
        }
        void operator delete(void* p, int blockType, const char* filename, int linenumber) {
            _aligned_free(p);
        }
#endif
    };

    struct CycleUnion { int cycles[NP][NP]; int lens[NP]; int count = 0; };

    K16A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr = NULL, bool bPrint = true);
    ~K16A2() {}

    bool addRow(int iRow, const unsigned char* source);
    void solve(int mode = 0);
    
private:

    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr);

    void* cbClass = NULL;
    int call_counter = 0;
    int kThreads;
    ResultCallback resultCallback;
    Factor fixedRows[NFIXED];
    std::vector<std::unique_ptr<ThreadLocalBuffers>> thread_buffers;
    int thread_row4[256];
    int thread_row5[256];
    int thread_root_idx[256];

    bool bTimeSet = false;

    // Internal Help Methods
    void VECTOR_CALL internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx, ThreadLocalBuffers* buf);
    void diagnostic_printout(double current_compr);

    FORCE_INLINE bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
        uint8_t curr = 0;
        for (int i = 0; i < NHALF; i++) {
            curr = adj1[curr];
            curr = adj2[curr];
            if (curr == 0 && i < NHALF - 1) return false;
        }
        return curr == 0;
    }

    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    
    FORCE_INLINE void apply_perm_18(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
        for (int i = 0; i < NP; i++)
            dst_adj[perm.p[i]] = perm.p[src_adj[i]];
    }

    


    struct CycleBacktrackState {
        K16A2* self;
        uint64_t total_generated = 0;
        uint64_t total_passed_p1f = 0;
        std::vector<std::array<uint8_t, NP>> valid_alphas;
        uint8_t F[9][NP];
        uint8_t pair_elements[8][2];
        int vertex_to_pair[NP];
        int vertex_to_pos[NP];
        uint8_t v0;

        struct SearchNode {
            uint8_t c[NP];
            bool used[NP];
            int pairs_visited;
        };
        bool is_collecting = false;
        int target_depth = 0;
        std::vector<SearchNode>* p_nodes = nullptr;

        // Parallel remaining-cycle generation. Used for small L, where the
        // main-cycle search produces too few work items to fill all threads:
        // we run the few main cycles sequentially but fan out each candidate's
        // generate_remaining_cycles across the thread pool. Explores the exact
        // same branches as the sequential path (results are identical).
        struct RemTask {
            uint8_t alpha_p[NP];
            uint8_t rem_used[NP];
            int start_idx;
        };
        bool parallel_remaining = false;
        bool collecting_rem = false;
        int rem_target_depth = 0;
        std::vector<RemTask>* p_rem_tasks = nullptr;

        CycleBacktrackState() = default;
        CycleBacktrackState(const CycleBacktrackState& other) {
            self = other.self;
            total_generated = 0;
            total_passed_p1f = 0;
            memcpy(F, other.F, sizeof(F));
            memcpy(pair_elements, other.pair_elements, sizeof(pair_elements));
            memcpy(vertex_to_pair, other.vertex_to_pair, sizeof(vertex_to_pair));
            memcpy(vertex_to_pos, other.vertex_to_pos, sizeof(vertex_to_pos));
            v0 = other.v0;
            is_collecting = false;
            target_depth = 0;
            p_nodes = nullptr;
            parallel_remaining = false;
            collecting_rem = false;
            rem_target_depth = 0;
            p_rem_tasks = nullptr;
        }
        CycleBacktrackState& operator=(const CycleBacktrackState& other) {
            if (this != &other) {
                self = other.self;
                total_generated = 0;
                total_passed_p1f = 0;
                valid_alphas.clear();
                memcpy(F, other.F, sizeof(F));
                memcpy(pair_elements, other.pair_elements, sizeof(pair_elements));
                memcpy(vertex_to_pair, other.vertex_to_pair, sizeof(vertex_to_pair));
                memcpy(vertex_to_pos, other.vertex_to_pos, sizeof(vertex_to_pos));
                v0 = other.v0;
                is_collecting = false;
                target_depth = 0;
                p_nodes = nullptr;
                parallel_remaining = false;
                collecting_rem = false;
                rem_target_depth = 0;
                p_rem_tasks = nullptr;
            }
            return *this;
        }

        void apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj);
        bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2);
        void backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L);
        void backtrackRecurse(int depth, int pairs_visited, uint8_t* c, bool* used, int L);
        void tryVertexForCycle(int v, int depth, int pairs_visited, uint8_t* c, bool* used, int L);
        void recurseWithVertex(int v, int depth, int next_pairs, uint8_t* c, bool* used, int L);
        void processCandidate(uint8_t* c, bool* used, int L);
        void buildPermutation(uint8_t* c, uint8_t* alpha_p, int L);
        bool checkPermutationPassed(const uint8_t* alpha_p, bool* used, int L);
        bool validateCandidateL(const uint8_t* alpha_p, bool* used, int L);
        void saveAlpha(const uint8_t* alpha_p);
        void generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size, bool* rem_used, uint8_t* alpha_p, int L, int depth);
        void generateRemainingParallel(const uint8_t* rem, int rem_size, uint8_t* alpha_p, int L);
        void checkTimeoutAndReport(int L, long long checked_reps);
        bool validateDefinedOrbits(const uint8_t* alpha_p, int L, const bool* defined, const uint8_t F0[][NP]);
    };

    struct CycleLengthStats {
        int L = 0;
        long long inputs = 0;
        long long passed = 0;
        size_t unique_classes = 0;
    };

    void runExhaustiveSearch();
    // Representative-method classifier (defined in k16A2rep.cpp). Private; called only
    // from runExhaustiveSearch when K16_USE_REP_METHOD==1. Uses kThreads workers and
    // gates all output on m_bPrint.
    void runRepresentativeMethod(int order, int target = 0);  // target>0 = harvest mode: stop after `target` distinct classes
    void searchCycleLength(int L, std::set<std::vector<uint8_t>>& unique_results, CycleLengthStats& stats);
    bool checkCyclesCompatibility(const uint8_t G[][NP], int L);
    bool decomposeMissingEdges(const uint8_t G[][NP], int L, const uint8_t* alpha_p, uint8_t H[][NP],
                               const std::function<void(const uint8_t[][NP])>& onCover = {});
    bool backtrackColor(int edge_idx, int num_edges, int num_colors,
                        const std::pair<uint8_t, uint8_t>* edges,
                        uint8_t matchings[][NP], const uint8_t* G0);
    void recordIsomorphicResults(const uint8_t H[][NP],
                                 std::set<std::vector<uint8_t>>& unique_results);
    void setupBacktrackState(CycleBacktrackState& state, int search_type);
    void constructFullH(const uint8_t G[][NP], int L, const uint8_t* alpha, uint8_t H[][NP]);
    void processAutomorphism(const std::array<uint8_t, NP>& alpha_arr, int L, std::set<std::vector<uint8_t>>& unique_results);
    void adj_to_src(const uint8_t* adj, unsigned char* src);
    void reportTotalResults(const std::set<std::vector<uint8_t>>& unique_results, std::chrono::high_resolution_clock::time_point start);
    void sendResultsToCallback(const std::set<std::vector<uint8_t>>& unique_results);
    void setupPairsTable(CycleBacktrackState& state, int search_type);
    void fillPairsTable(CycleBacktrackState& state, int search_type);
    void storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v);
    void constructFullHFromAut(const uint8_t* alpha, int L, uint8_t H[][NP]);
    void findMissingEdges(const uint8_t G[][NP], int L, std::pair<uint8_t, uint8_t>* edges);
    bool isEdgeMissing(const uint8_t G[][NP], int L, int u, int v);
    bool checkMatchingsCompatibility(uint8_t matchings[][NP], int num_colors, const uint8_t* G0);
    bool tryColoringEdge(int edge_idx, int num_edges, int num_colors, const std::pair<uint8_t, uint8_t>* edges, uint8_t matchings[][NP], const uint8_t* G0);
    int getSymmetryBreakingLimit(uint8_t matchings[][NP], int num_colors);
    bool is_color_compatible(int c, int num_colors, uint8_t matchings[][NP], const uint8_t* G0);
    bool isColorUsed(const uint8_t* matching);
    void copyMatchingsToH(uint8_t matchings[][NP], int num_colors, const uint8_t G[][NP], int L, uint8_t H[][NP]);
    void tryIsomorphicMapping(const uint8_t H[][NP], const CycleUnion& cu_H, int v, std::set<std::vector<uint8_t>>& unique_results);
    void buildStarterCycle(uint8_t* cyc_R);
    void buildHCycle(const uint8_t H[][NP], int v, uint8_t* cyc_H);
    void buildMappingPermutation(const uint8_t* cyc_H, const uint8_t* cyc_R, uint8_t* p);
    void checkAndRecordPermutedH(const uint8_t H[][NP], const uint8_t* p, std::set<std::vector<uint8_t>>& unique_results);
    void applyPermToH(const uint8_t H[][NP], const uint8_t* p, uint8_t S[][NP]);
    bool doesSMatchFixedRows(const uint8_t S[][NP]);
    void recordS(const uint8_t S[][NP], std::set<std::vector<uint8_t>>& unique_results);
    void verifyL16Pairing(const std::set<std::vector<uint8_t>>& local_unique);
    std::chrono::steady_clock::time_point case_start_time;
    std::chrono::steady_clock::time_point last_print_time;
    double case_timeout_seconds = 1e9; // Very large timeout value (effectively infinite)
    int min_cycle_length = 3; // K16 aut>2 validation: L=3..15 (skip slow L=2 for now)
    bool case_timed_out = false;
    long long current_checked_reps = 0;
    int current_top_branch_idx = 0;
    int total_top_branches = 0;

    void printEstimatedTime(int L, long long checked_reps, int search_type);

    std::mutex result_mutex;
};

#pragma pop_macro("new")

#endif
