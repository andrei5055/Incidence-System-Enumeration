#ifndef K16A2_H
#define K16A2_H

#include <vector>
#include <array>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <immintrin.h>
#include "kBase.h"
#include <unordered_map>
#include <set>
#include <mutex>
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
#define K16_N 16
#define K16_MATCH 15
#define K16_FIXED 3
#define K16_SEARCH 12
#define K16_M_MAX 60000

#ifndef K16P1F_H
struct Mask256_C { uint64_t m[2] = { 0, 0 }; };
#endif

FORCE_INLINE bool is_disjoint(const Mask256_C& a, const Mask256_C& b) {
    return ((a.m[0] & b.m[0]) | (a.m[1] & b.m[1])) == 0;
}

#pragma push_macro("new")
#undef new

class K16A2 : public KBase<Mask256_C> {
public:
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
        Mask256_C edge_mask;
        uint8_t adj[K16_N];
    };

    struct FastSortedFactor : FastSortedFactorBase<32> {};
    struct FastRowTriplet { FastSortedFactor r[3]; };

    struct Factor {
        uint8_t adj[K16_N];
        uint8_t src[K16_N];
        Mask256_C edge_mask;
        FastSortedFactor fs;
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
        int s4_bit_cursor_table[K16_SEARCH];

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

    struct TransInfo { Permutation perms[512]; int count = 0; };
    struct CycleUnion { int cycles[16][16]; int lens[16]; int count = 0; };

    K16A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr = NULL, bool bPrint = true);
    ~K16A2() {}

    bool addRow(int iRow, const unsigned char* source);
    void solve(int mode = 0);
    
    void setTransitionConfig(const int config[5]) {
        for (int i = 0; i < 5; i++) trans_config[i] = config[i];
    }

private:
    int trans_config[5] = { 1, 2, 3, 4, 3 };

    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr);

    struct ValidChain {
        Permutation alpha;
        uint8_t factors_adj[16][16];
        PackedAdj factors_packed[16];
        bool step_present[16];
    };
    std::vector<ValidChain> valid_chains;
    std::vector<Permutation> saved_transitions;
    int max_row_seen;
    uint64_t rejected_candidates[16];

    void* cbClass = NULL;
    int call_counter = 0;
    int kThreads;
    ResultCallback resultCallback;
    Factor fixedRows[K16_FIXED];
    PackedAdj fixed_packed[K16_FIXED];
    Mask256_C fixedEdgesMask;
    int edge_id_table[K16_N][K16_N];
    uint8_t edge_to_u[256];
    uint8_t edge_to_v[256];

    std::vector<Factor> global_pool;
    std::vector<PackedAdj> packed_pool;

    struct ArrayHash {
        size_t operator()(const std::array<uint8_t, 16>& a) const {
            size_t hash = 17;
            for (uint8_t x : a) {
                hash = hash * 31 + x;
            }
            return hash;
        }
    };
    std::unordered_map<std::array<uint8_t, 16>, int, ArrayHash> f_map_unordered;

    // Lazy Candidate Caching States
    bool slot_generated[K16_SEARCH];
    std::mutex slot_mutex[K16_SEARCH];
    std::mutex pool_mutex;

    // Lazy generation methods
    void ensureSlotGenerated(int slotIndex);
    void generateSlotCandidates(int slotIndex);
    void generateRecursive(int depth, uint64_t usedMask, uint8_t* factor, int slotIndex, int P0);

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

    FORCE_INLINE bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
        uint8_t curr = 0;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr]; if (curr == 0) return false;
        curr = adj1[curr]; curr = adj2[curr];
        return curr == 0;
    }

    FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if ((p1.edge_mask.m[0] & p2.edge_mask.m[0]) || 
            (p1.edge_mask.m[1] & p2.edge_mask.m[1])) 
            return false;
        return is_perfect_scalar(p1.adj, p2.adj);
    }

    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    void get_transformations(const Factor& fi, const Factor& fj, TransInfo& info);
    void get_transformations_general(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl, TransInfo& info);
    
    FORCE_INLINE void apply_perm_16(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
        dst_adj[perm.p[0]] = perm.p[src_adj[0]];
        dst_adj[perm.p[1]] = perm.p[src_adj[1]];
        dst_adj[perm.p[2]] = perm.p[src_adj[2]];
        dst_adj[perm.p[3]] = perm.p[src_adj[3]];
        dst_adj[perm.p[4]] = perm.p[src_adj[4]];
        dst_adj[perm.p[5]] = perm.p[src_adj[5]];
        dst_adj[perm.p[6]] = perm.p[src_adj[6]];
        dst_adj[perm.p[7]] = perm.p[src_adj[7]];
        dst_adj[perm.p[8]] = perm.p[src_adj[8]];
        dst_adj[perm.p[9]] = perm.p[src_adj[9]];
        dst_adj[perm.p[10]] = perm.p[src_adj[10]];
        dst_adj[perm.p[11]] = perm.p[src_adj[11]];
        dst_adj[perm.p[12]] = perm.p[src_adj[12]];
        dst_adj[perm.p[13]] = perm.p[src_adj[13]];
        dst_adj[perm.p[14]] = perm.p[src_adj[14]];
        dst_adj[perm.p[15]] = perm.p[src_adj[15]];
    }

    FastSortedFactor get_fast_sorted(const uint8_t* adj);
    
    CycleUnion target_cu;
    TransInfo fixed_trans[3];
    FastSortedFactor r1_can, r2_can;

    // Statistics and Reporting
    bool compare_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)   { return memcmp(a.pairs, b.pairs, 16) < 0; }
    bool equal_fast_sorted(const FastSortedFactor& a, const FastSortedFactor& b)     { return memcmp(a.pairs, b.pairs, 16) == 0; }
    bool compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b);

    struct CycleBacktrackState {
        K16A2* self;
        uint64_t total_generated = 0;
        uint64_t total_passed_p1f = 0;
        std::vector<std::array<uint8_t, 16>> valid_alphas;
        uint8_t F[9][16];
        uint8_t pair_elements[8][2];
        int vertex_to_pair[16];
        int vertex_to_pos[16];
        uint8_t v0;
        uint8_t fixed_point = 0;
        int transposed_point = 1;
        int search_type = 1;

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
        bool is_partial_permutation_ok(const uint8_t* alpha_p, const bool* defined);
        void generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size, bool* rem_used, bool* defined, uint8_t* alpha_p, int L);
    };

    struct CycleLengthStats {
        int L = 0;
        long long inputs = 0;
        long long passed = 0;
        size_t unique_classes = 0;
    };

    void runExhaustiveSearch();
    void searchCycleLength(int L, std::set<std::vector<uint8_t>>& unique_results, CycleLengthStats& stats);
    bool checkCyclesCompatibility(const uint8_t G[][16], int L);
    bool decomposeMissingEdges(const uint8_t G[][16], int L, uint8_t H[][16]);
    bool backtrackColor(int edge_idx, int num_edges, int num_colors,
                        const std::pair<uint8_t, uint8_t>* edges,
                        uint8_t matchings[][16], const uint8_t* G0);
    void recordIsomorphicResults(const uint8_t H[][16],
                                 std::set<std::vector<uint8_t>>& unique_results);
    void setupBacktrackState(CycleBacktrackState& state, int search_type, int fixed_point = 0, int transposed_point = 1, int override_v0 = -1);
    void constructFullH(const uint8_t G[][16], int L, const uint8_t* alpha, uint8_t H[][16]);
    void processAutomorphism(const std::array<uint8_t, 16>& alpha_arr, int L, std::set<std::vector<uint8_t>>& unique_results);
    void adj_to_src(const uint8_t* adj, unsigned char* src);
    void reportTotalResults(const std::set<std::vector<uint8_t>>& unique_results, std::chrono::high_resolution_clock::time_point start);
    void sendResultsToCallback(const std::set<std::vector<uint8_t>>& unique_results);
    void setupPairsTable(CycleBacktrackState& state, int search_type);
    void fillPairsTable(CycleBacktrackState& state, int search_type);
    void storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v);
    bool constructFullHFromAut(const uint8_t* alpha, int L, uint8_t H[][16]);
    void findMissingEdges(const uint8_t G[][16], int L, std::pair<uint8_t, uint8_t>* edges);
    bool isEdgeMissing(const uint8_t G[][16], int L, int u, int v);
    bool checkMatchingsCompatibility(uint8_t matchings[][16], int num_colors, const uint8_t* G0);
    bool tryColoringEdge(int edge_idx, int num_edges, int num_colors, const std::pair<uint8_t, uint8_t>* edges, uint8_t matchings[][16], const uint8_t* G0);
    int getSymmetryBreakingLimit(uint8_t matchings[][16], int num_colors);
    bool is_color_compatible(int c, int num_colors, uint8_t matchings[][16], const uint8_t* G0);
    bool isColorUsed(const uint8_t* matching);
    void copyMatchingsToH(uint8_t matchings[][16], int num_colors, const uint8_t G[][16], int L, uint8_t H[][16]);
    void tryIsomorphicMapping(const uint8_t H[][16], int a, int b, const CycleUnion& cu_H, int v, std::set<std::vector<uint8_t>>& unique_results);
    void buildStarterCycle(uint8_t* cyc_R);
    void buildHCycle(const uint8_t H[][16], int a, int b, int v, uint8_t* cyc_H);
    void buildMappingPermutation(const uint8_t* cyc_H, const uint8_t* cyc_R, uint8_t* p);
    void checkAndRecordPermutedH(const uint8_t H[][16], const uint8_t* p, std::set<std::vector<uint8_t>>& unique_results);
    void applyPermToH(const uint8_t H[][16], const uint8_t* p, uint8_t S[][16]);
    bool doesSMatchFixedRows(const uint8_t S[][16]);
    void recordS(const uint8_t S[][16], std::set<std::vector<uint8_t>>& unique_results);
    void verifyL14Pairing(const std::set<std::vector<uint8_t>>& local_unique);

    std::chrono::steady_clock::time_point case_start_time;
    std::chrono::steady_clock::time_point last_print_time;
    double case_timeout_seconds = 1800.0;
    bool case_timed_out = false;
    long long current_checked_reps = 0;
    int current_top_branch_idx = 0;
    int total_top_branches = 0;

    void printEstimatedTime(int L, long long checked_reps);

    std::mutex result_mutex;
};

#pragma pop_macro("new")

#endif
