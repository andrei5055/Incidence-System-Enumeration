#pragma once

#include <thread>
#include <mutex>
#include <functional>
#include <map>
#include <set>
#include <memory>

typedef bool (*ResultCallback)(const void* cbClass, const unsigned char* results, int r4, int r5, int mode);

struct RowRange { int start_word; int end_word; };

struct SubSamplePlanItem {
    uint16_t s4_word_idx;
    uint16_t l_word_off;
    uint64_t mask;
    uint8_t l_bit_shift;
    uint8_t bits_pushed;
    uint8_t split;
};

// Canonicity structures
struct alignas(32) Permutation {
    uint8_t p[32];
    uint8_t p_inv[32];
    __m256i pv;
    __m256i pinvv;
    int swaps = 0;
    int padding[7]; // Ensure sizeof=160 (multiple of 32)
};

template<typename T>
struct LocalBufers {
    std::vector<int> g_to_l;
    std::vector<int> local_to_global;
    std::vector<int> cl;
    std::vector<RowRange> local_ranges;
    std::vector<uint64_t> local_adj;
    std::vector<size_t> local_offsets;
    std::vector<T> local_edge_masks;
    std::vector<int> local_s_to_f;
    struct ActiveWord { int global_word_idx; uint64_t mask; int local_bit_start; };
    std::vector<ActiveWord> active_words;

    // Phase 2: Semi-local cache for a fixed r4
    std::vector<int> s4_to_global;
    std::vector<int> global_to_s4;
    std::vector<uint64_t> s4_adj;
    std::vector<size_t> s4_offsets;
    std::vector<uint8_t> s4_canonical;
};

class KSolver {
public:
    virtual ~KSolver() = default;
    virtual void solve(int mode = 0) = 0;
    virtual void preSolve() {}
    virtual int getCandidateCount(int slot) { return 0; }
    virtual bool getCandidate(int slot, int index, unsigned char* out_src) { return false; }
    virtual void deleteCandidate(int slot, int index) {}
    virtual bool addRow(int iRow, const unsigned char* source) = 0;
};

template<std::size_t N>
struct alignas(64) FastSortedFactorBase { uint8_t pairs[N]; };

struct FactorParams {
    const int m_nPlayers;
	const int m_nMatched;
    const int m_nFixedRows;
    const int m_nMax;
    const int m_nWords;
    const int m_nSearch;
    FactorParams(int numPlayers, int numMatched, int numFixedRows, int nMax) : 
        m_nPlayers(numPlayers), m_nMatched(numMatched), m_nFixedRows(numFixedRows), m_nMax(nMax), 
        m_nWords((((nMax + 255) / 256 * 4) + 4)), m_nSearch(numMatched - numFixedRows) {}
};

template<typename T>
class KBase : public KSolver, public FactorParams {
protected:
    KBase(const FactorParams& factParam, bool bPrint) : m_bPrint(bPrint), FactorParams(factParam) {
        temp_slot_ids = new std::vector<int>[m_nSearch];
        row_factor_ids = new std::set<int>[m_nSearch];
        m_pResults = new unsigned char[m_nMatched * m_nPlayers];
	}
    ~KBase() {
		delete[] temp_slot_ids;
		delete[] row_factor_ids;
        delete[] m_pResults;
	}

    void init() {
        num_results = 0;
        num_notCanon = 0;
        g_total_pairs_checked = 0;
        g_rejected_edges = 0;
        g_rejected_cycles = 0;
        g_rejected_canon = 0;
        g_total_rejected = 0;
        g_total_kept = 0;
        min_stab_size = 1000000;
        max_stab_size = 0;
        results_to_sort.clear();
        for (int i = m_nSearch; i--;) {
            temp_slot_ids[i].clear();
            row_factor_ids[i].clear();
        }
	}

    void addFactorToSlot(int slot, int factor_id) {
        if (row_factor_ids[slot].find(factor_id) == row_factor_ids[slot].end()) {
            temp_slot_ids[slot].push_back(factor_id);
            row_factor_ids[slot].insert(factor_id);
        }
    }

    // Search Structures
    const bool m_bPrint;
    std::vector<int> search_to_factor;
    std::vector<uint64_t> adj_matrix;
    std::vector<T> factor_edge_masks;
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

    // Statistics and Reporting
    std::atomic<int> g_total_pairs_checked{ 0 };
    std::atomic<int> g_rejected_edges{ 0 };
    std::atomic<int> g_rejected_cycles{ 0 };
    std::atomic<int> g_rejected_canon{ 0 };
    std::atomic<int> g_total_rejected{ 0 };
    std::atomic<int> g_total_kept{ 0 };

    std::atomic<int> min_stab_size{ 1000000 };
    std::atomic<int> max_stab_size{ 0 };
    std::vector<std::vector<unsigned char>> results_to_sort;
    std::map<std::vector<uint8_t>, int> f_map;

    std::vector<int>* temp_slot_ids;
    std::set<int>* row_factor_ids;

    int fixed3RowsIndex = 0;

//    FastSortedFactor<N> r1_can, r2_can;
    __m128i gfs_table_adj[256];
    __m128i gfs_table_idx[256];

    unsigned char *m_pResults = NULL;
};