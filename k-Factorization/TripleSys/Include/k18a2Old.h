#ifndef K18A2OLD_H
#define K18A2OLD_H

// =============================================================================
// K18A2Old — starter-based (pool + arc-consistency + orbit-clique) solver for K18.
//
// Purpose (see MD/k18a2old_spec.md §2): a FAST completion engine for finding K18
// Perfect-1-Factorizations that admit an order-4 automorphism (L=4), where the
// cyclic K18A2's per-alpha `decomposeMissingEdges` is intractable (30 cores x 50h+).
//
// Design (benchmark-driven):
//   * Build ONCE an alpha-INDEPENDENT candidate pool = all factor rows pairwise
//     Hamiltonian with the 3 fixed starter rows (amortized across all alphas).
//   * Arc-consistency prune the pool.
//   * For each candidate order-4 automorphism alpha: group the pool into alpha-orbits
//     and run an orbit clique exact-cover to emit every alpha-invariant P1F.
//
// This is a NP=18 scalar adaptation of K16A2Old. Key differences from K16A2Old:
//   * 153 edges (vs 120) -> 256-bit Mask18O_C (m[0..2] used; K16 used 128-bit m[0..1]).
//   * NO __m128i/__m256i kernels (18 vertices don't fit a 16-byte shuffle and 4 does
//     not divide 18 so no AVX-512); reuse scalar find_cycles/is_perfect_scalar/apply_perm.
//   * alpha source is NOT get_transformations (an 18-cycle has no order-4 rotation,
//     4 does not divide 18); order-4 alphas come from a dedicated generator (see .cpp).
//
// Zero risk to the committed cyclic K18A2 / decompose path: this is a new class in
// new files; only TripleSys.cpp dispatch + the vcxproj get small additive edits.
// =============================================================================

#include <vector>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <array>
#include "kBase.h"

#ifdef _MSC_VER
#define K18O_FORCE_INLINE __forceinline
#else
#define K18O_FORCE_INLINE inline __attribute__((always_inline))
#endif

#define K18O_N        18                 // players
#define K18O_MATCH    17                 // factors / match rows (NP-1)
#define K18O_FIXED    3                  // fixed starter rows
#define K18O_SEARCH   14                 // rows to place = MATCH - FIXED
#define K18O_EDGES    153                // NP*(NP-1)/2
#define K18O_HALF     9                  // NP/2 (edges per perfect matching)
#define K18O_M_MAX    2097152            // candidate-pool cap (K18 single-starter pool > 131072; sized large to discover/hold true size)
#define K18O_AUT_ORDER 4                 // target automorphism order (production: L=4)

// 256-bit edge mask. 153 edges need words m[0],m[1],m[2]; m[3] is headroom.
// Distinct type/name from K18A2's Mask18_C to avoid any ODR coupling.
struct Mask18O_C { uint64_t m[4] = { 0, 0, 0, 0 }; };

class K18A2Old : public KBase<Mask18O_C> {
public:
    struct alignas(64) PackedAdj {
        Mask18O_C edge_mask;
        uint8_t adj[K18O_N];
    };

    struct Factor {
        uint8_t adj[K18O_N];
        uint8_t src[K18O_N];
        Mask18O_C edge_mask;
    };

    struct CycleUnion { int cycles[K18O_N][K18O_N]; int lens[K18O_N]; int count = 0; };

    K18A2Old(const FactorParams& factParam, int fixed3RowsIndex, int kThreads,
             const unsigned char* first3Rows, ResultCallback callback,
             void* cbClassPtr = NULL, bool bPrint = true);
    ~K18A2Old() {}

    bool addRow(int iRow, const unsigned char* source);   // append a candidate row to the pool
    void solve(int mode = 0);

private:
    void init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows,
              ResultCallback callback, void* cbClassPtr);

    // ---- scalar primitives (NP=18) ----
    K18O_FORCE_INLINE bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
        // two factors are "perfect" (compatible) iff their union is a single Hamiltonian cycle
        uint8_t curr = 0;
        for (int i = 0; i < K18O_HALF; i++) {
            curr = adj1[curr];
            curr = adj2[curr];
            if (curr == 0 && i < K18O_HALF - 1) return false;
        }
        return curr == 0;
    }
    K18O_FORCE_INLINE bool is_perfect_packed(const PackedAdj& p1, const PackedAdj& p2) {
        if ((p1.edge_mask.m[0] & p2.edge_mask.m[0]) ||
            (p1.edge_mask.m[1] & p2.edge_mask.m[1]) ||
            (p1.edge_mask.m[2] & p2.edge_mask.m[2]))
            return false;                              // share an edge -> not edge-disjoint
        return is_perfect_scalar(p1.adj, p2.adj);
    }
    void apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj) {
        // permute a factor's adjacency by vertex permutation `perm`
        for (int i = 0; i < K18O_N; i++)
            dst_adj[perm[i]] = perm[src_adj[i]];
    }
    CycleUnion find_cycles(const uint8_t* adj1, const uint8_t* adj2);
    PackedAdj pack_factor_adj(const uint8_t* adj);
    void set_edge_mask(Mask18O_C& mask, int eid) {     // 153 edges -> 3 words
        if (eid < 64)       mask.m[0] |= (1ULL << eid);
        else if (eid < 128) mask.m[1] |= (1ULL << (eid - 64));
        else                mask.m[2] |= (1ULL << (eid - 128));
    }

    // ---- engine phases (defined in k18a2Old.cpp) ----
    void arcConsistencyPrune(std::vector<uint8_t>& active);                 // shrink pool
    // group the (pruned) pool into alpha-orbits and clique-cover -> emit P1Fs for this alpha
    void completeForAlpha(const uint8_t* alpha, const std::vector<uint8_t>& active,
                          const std::vector<int>& active_slots);
    // enumerate candidate order-4 automorphisms consistent with the 3 fixed rows and
    // drive completeForAlpha for each. (porting target: K18A2 CycleBacktrackState, L=4)
    void generateOrder4Alphas(const std::vector<uint8_t>& active,
                              const std::vector<int>& active_slots);
    // Factor-alignment alpha generation (scalar NP=18 port of K16A2Old::get_transformations_general):
    // if fi-union-fj and fk-union-fl are each a single Hamiltonian cycle, every rotation/reflection
    // alignment of the source cycle onto the target cycle yields a candidate permutation mapping
    // fi->fk, fj->fl. Appends them to `out`.
    void alignmentAlphas(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl,
                         std::vector<std::array<uint8_t, K18O_N>>& out);
    bool has_order(const uint8_t* p, int target_order);  // perm order == target_order ?

    // ---- order-4 automorphism generator (ported from K18A2::CycleBacktrackState, L=4) ----
    // 4 does NOT divide 18, so order-4 alphas cannot be Hamiltonian-cycle rotations
    // (K16A2Old's get_transformations). Instead we enumerate, by backtracking, every vertex
    // permutation alpha whose main cycle (through v0) has length L=K18O_AUT_ORDER and whose
    // remaining cycles have lengths dividing L (so order(alpha) == L exactly), consistent with
    // the first fixed row F[0]. This is the proven cyclic-K18A2 generator with its decompose
    // existence-check removed: every structurally valid alpha is saved and later completed by
    // completeForAlpha (orbit clique) instead of decomposeMissingEdges.
    struct CycleBacktrackState {
        K18A2Old* self = nullptr;
        uint64_t total_generated = 0;     // main-cycle candidates examined
        uint64_t total_passed_p1f = 0;    // alphas saved (powers of F[0] all perfect)
        std::vector<std::array<uint8_t, K18O_N>> valid_alphas;
        uint8_t F[K18O_HALF][K18O_N];     // only F[0] is used (= fixedRows[0].adj)
        uint8_t pair_elements[K18O_HALF][2];
        int vertex_to_pair[K18O_N];
        int vertex_to_pos[K18O_N];
        uint8_t v0 = 0;
        int search_type = 1;              // 1 = fixed out-of-cycle pts, 2 = transposed (0 1)
        long long n_reach_leaf = 0;       // [DIAG] alphas fully built (reached checkPermutationPassed)
        long long n_pass_step1 = 0;       // [DIAG] passed F[0]-power perfect check
        long long n_pass_step2 = 0;       // [DIAG] passed 3-fixed-row pool check (== saved)
        long long enum_steps = 0;         // cycle-extension steps (for the safety cap)
        bool enum_aborted = false;        // set if enum_steps exceeds the cap (enumeration too large)

        void apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj);
        bool is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2);
        void backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L);
        void tryVertexForCycle(int v, int depth, int pairs_visited, uint8_t* c, bool* used, int L);
        void recurseWithVertex(int v, int depth, int next_pairs, uint8_t* c, bool* used, int L);
        void processCandidate(uint8_t* c, bool* used, int L);
        void buildPermutation(uint8_t* c, uint8_t* alpha_p, int L);
        bool checkPermutationPassed(const uint8_t* alpha_p, bool* used, int L);
        void saveAlpha(const uint8_t* alpha_p);
        void generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size,
                                       bool* rem_used, uint8_t* alpha_p, int L, int depth);
    };
    void setupBacktrackState(CycleBacktrackState& state, int search_type);
    void setupPairsTable(CycleBacktrackState& state, int search_type);
    void fillPairsTable(CycleBacktrackState& state, int search_type);
    void storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v);

    // ---- state ----
    void* cbClass = NULL;
    int kThreads = 1;
    ResultCallback resultCallback = nullptr;

    Factor fixedRows[K18O_FIXED];
    PackedAdj fixed_packed[K18O_FIXED];
    Mask18O_C fixedEdgesMask;
    int edge_id_table[K18O_N][K18O_N];

    std::vector<Factor>    global_pool;       // candidate rows (index 0..2 = fixed rows)
    std::vector<PackedAdj> packed_pool;
    std::map<std::vector<uint8_t>, int> f_map;       // adj -> pool index (orbit lookup)
    std::vector<std::vector<int>> temp_slot_ids;     // per-search-slot candidate ids

    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point last_print_time;
    std::mutex result_mutex;
    bool bTimeSet = false;
    std::atomic<uint64_t> solutions{0};

    // [PERF-MEASURE] split candidate-generation (main SW) vs addRow body cost
    double m_addRowSeconds = 0.0;
    long long m_addRowCalls = 0;
};

#endif // K18A2OLD_H
