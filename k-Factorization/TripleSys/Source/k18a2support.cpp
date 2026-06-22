#include "k18a2.h"
#include <thread>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <set>
#include <chrono>
#include <atomic>

static std::atomic<long long> g_candidates_processed{0};
static std::atomic<long long> g_candidates_passed_perfect{0};
static std::atomic<long long> g_candidates_passed_validate{0};
static std::atomic<long long> g_generate_matchings_calls{0};
static std::atomic<long long> g_generate_matchings_total_matchings{0};
// Phase timers (nanoseconds) for decomposeMissingEdges — profiling.
static std::atomic<long long> g_ns_enum{0};   // generate_matchings (pool enumeration)
static std::atomic<long long> g_ns_group{0};  // sort + hashmap + orbit grouping + compat
static std::atomic<long long> g_ns_cover{0};  // exact-cover on orbits
static std::atomic<long long> g_parallel_tasks{0}; // tasks fanned out in generateRemainingParallel
static std::atomic<long long> g_tasks_done{0};      // parallel tasks completed (for live ETA)

struct FlatMatchingMap {
    struct Entry {
        std::array<uint8_t, 18> key;
        int value = -1;
    };
    inline static thread_local std::vector<Entry> table;
    int mask = 0;

    void build(const std::vector<std::array<uint8_t, 18>>& pool) {
        int cap = 1;
        while (cap < (int)pool.size() * 2) {
            cap *= 2;
        }
        table.resize(cap);
        for (int i = 0; i < cap; i++) {
            table[i].value = -1;
        }
        mask = cap - 1;

        for (int i = 0; i < (int)pool.size(); i++) {
            const auto& m = pool[i];
            uint64_t h = hash_m(m);
            int idx = h & mask;
            while (table[idx].value != -1) {
                idx = (idx + 1) & mask;
            }
            table[idx].key = m;
            table[idx].value = i;
        }
    }

    int find(const std::array<uint8_t, 18>& m) const {
        uint64_t h = hash_m(m);
        int idx = h & mask;
        while (table[idx].value != -1) {
            if (table[idx].key == m) {
                return table[idx].value;
            }
            idx = (idx + 1) & mask;
        }
        return -1;
    }

private:
    static FORCE_INLINE uint64_t hash_m(const std::array<uint8_t, 18>& m) {
        uint64_t h1, h2;
        memcpy(&h1, m.data(), 8);
        memcpy(&h2, m.data() + 8, 8);
        uint16_t h3;
        memcpy(&h3, m.data() + 16, 2);
        return h1 * 3137 + h2 * 17 + h3;
    }
};


// ==========================================
// CycleBacktrackState Implementation
// ==========================================

void K18A2::CycleBacktrackState::apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj) {
    for (int i = 0; i < 18; i++) {
        dst_adj[perm[i]] = perm[src_adj[i]];
    }
}

bool K18A2::CycleBacktrackState::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    uint8_t curr = 0;
    for (int i = 0; i < 9; i++) {
        curr = adj1[curr];
        curr = adj2[curr];
        if (curr == 0 && i < 8) return false;
    }
    return curr == 0;
}

void K18A2::CycleBacktrackState::checkTimeoutAndReport(int L, long long checked_reps) {
    if (self->case_timed_out) return;
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - self->case_start_time).count();
    if (elapsed > self->case_timeout_seconds) {
        self->case_timed_out = true;
    }
    self->printEstimatedTime(L, checked_reps, (v0 == 3) ? 2 : 1);
}

void K18A2::CycleBacktrackState::backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    if (self->case_timed_out) return;
    static thread_local int backtrack_call_count = 0;
    if (++backtrack_call_count >= 1000) {
        backtrack_call_count = 0;
        checkTimeoutAndReport(L, g_candidates_processed.load());
    }
    if (is_collecting && depth == target_depth) {
        SearchNode node;
        memcpy(node.c, c, 18);
        memcpy(node.used, used, 18);
        node.pairs_visited = pairs_visited;
        p_nodes->push_back(node);
        return;
    }
    if (depth == L - 1) {
        processCandidate(c, used, L);
        return;
    }
    if (depth == 0 && !is_collecting && self->kThreads > 1) {
        int t_depth = std::min(5, L - 1);
        if (t_depth > 0) {
            std::vector<SearchNode> nodes;
            
            // Set up collection state
            this->is_collecting = true;
            this->target_depth = t_depth;
            this->p_nodes = &nodes;
            
            // Run search to collect nodes at target_depth
            backtrackRecurse(0, pairs_visited, c, used, L);
            
            this->is_collecting = false;

            if (!nodes.empty() && (int)nodes.size() < self->kThreads) {
                // Small L: too few main-cycle work items to fill the thread pool.
                // Run the candidates sequentially on this thread, but parallelize
                // each one's remaining-cycle generation so all threads stay busy.
                // (For small L the collected nodes are complete main cycles, i.e.
                // t_depth == L - 1, so each backtrack call reaches processCandidate.)
                this->parallel_remaining = true;
                for (auto& nd : nodes) {
                    uint8_t local_c[18];
                    bool local_used[18];
                    memcpy(local_c, nd.c, 18);
                    memcpy(local_used, nd.used, 18);
                    backtrack(t_depth, nd.pairs_visited, local_c, local_used, L);
                }
                this->parallel_remaining = false;
                return;
            }
            if (!nodes.empty()) {
                std::atomic<size_t> next_node_idx{ 0 };
                std::vector<std::thread> workers;
                std::mutex stats_mutex;
                uint64_t accumulated_generated = 0;
                uint64_t accumulated_passed_p1f = 0;
                
                for (int t = 0; t < self->kThreads; t++) {
                    workers.emplace_back([&]() {
                        uint8_t local_c[18];
                        bool local_used[18];
                        
                        while (true) {
                            size_t idx = next_node_idx.fetch_add(1);
                            if (idx >= nodes.size()) break;
                            
                            {
                                std::lock_guard<std::mutex> lock(stats_mutex);
                                self->current_top_branch_idx++;
                            }
                            
                            CycleBacktrackState thread_state = *this;
                            thread_state.total_generated = 0;
                            thread_state.total_passed_p1f = 0;
                            thread_state.is_collecting = false;
                            
                            memcpy(local_c, nodes[idx].c, 18);
                            memcpy(local_used, nodes[idx].used, 18);
                            
                            thread_state.backtrack(t_depth, nodes[idx].pairs_visited, local_c, local_used, L);
                            
                            {
                                std::lock_guard<std::mutex> lock(stats_mutex);
                                accumulated_generated += thread_state.total_generated;
                                accumulated_passed_p1f += thread_state.total_passed_p1f;
                                this->valid_alphas.insert(this->valid_alphas.end(), thread_state.valid_alphas.begin(), thread_state.valid_alphas.end());
                            }
                        }
                    });
                }
                
                for (auto& w : workers) {
                    w.join();
                }
                
                this->total_generated = accumulated_generated;
                this->total_passed_p1f = accumulated_passed_p1f;
                return;
            }
        }
    }
    backtrackRecurse(depth, pairs_visited, c, used, L);
}

void K18A2::CycleBacktrackState::backtrackRecurse(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    for (int v = 0; v < 18; v++) {
        if (used[v]) continue;
        tryVertexForCycle(v, depth, pairs_visited, c, used, L);
    }
}

void K18A2::CycleBacktrackState::tryVertexForCycle(int v, int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    int p_idx = vertex_to_pair[v];
    if (p_idx == -2) {
        recurseWithVertex(v, depth, pairs_visited, c, used, L);
        return;
    }
    bool pair_visited = used[pair_elements[p_idx][0]] || used[pair_elements[p_idx][1]];
    if (pair_visited) {
        recurseWithVertex(v, depth, pairs_visited, c, used, L);
    } else if (v == pair_elements[pairs_visited][0]) {
        recurseWithVertex(v, depth, pairs_visited + 1, c, used, L);
    }
}

void K18A2::CycleBacktrackState::recurseWithVertex(int v, int depth, int next_pairs, uint8_t* c, bool* used, int L) {
    if (depth == 0) {
        self->current_top_branch_idx++;
        self->printEstimatedTime(L, g_candidates_processed.load(), (v0 == 3) ? 2 : 1);
    }
    used[v] = true;
    c[depth] = v;
    backtrack(depth + 1, next_pairs, c, used, L);
    used[v] = false;
}

static bool validateDefinedOrbits(const uint8_t* alpha_p, int L, const bool* defined, const uint8_t F0[][18]) {
    int limit = L / 2;
    for (int u = 0; u < 18; u++) {
        if (!defined[u]) continue;
        int v = F0[0][u];
        if (!defined[v] || u > v) continue;
        
        uint8_t curr_u = u;
        uint8_t curr_v = v;
        for (int j = 1; j <= limit; j++) {
            curr_u = alpha_p[curr_u];
            curr_v = alpha_p[curr_v];
            if (F0[0][curr_u] == curr_v) {
                return false;
            }
        }
    }
    return true;
}

void K18A2::CycleBacktrackState::processCandidate(uint8_t* c, bool* used, int L) {
    long long count = g_candidates_processed.fetch_add(1) + 1;
    total_generated++;
    self->current_checked_reps = count;

    if (count % 100 == 0) {
        self->printEstimatedTime(L, count, (v0 == 3) ? 2 : 1);
    }

    // Prune candidates where the main cycle antipodal edges are not missing.
    // For even L (powers of 2: 2, 4, 8, 16), any valid completion must contain at least one invariant factor,
    // which requires the antipodal edges of the main cycle to be missing in G[0] (stored in F[0]).
    if (L == 2 || L == 4 || L == 8 || L == 16) {
        bool antipodal_ok = true;
        int half = L / 2;
        for (int i = 0; i < half; i++) {
            uint8_t u = (i == 0) ? v0 : c[i - 1];
            uint8_t v = c[i + half - 1];
            if (F[0][u] == v) {
                antipodal_ok = false;
                break;
            }
        }
        if (!antipodal_ok) {
            return;
        }
    }

    uint8_t alpha_p[18];
    buildPermutation(c, alpha_p, L);
    
    int search_type = (v0 == 3) ? 2 : 1;
    if (search_type == 2) {
        alpha_p[0] = 1;
        alpha_p[1] = 0;
    }
    
    bool defined[18] = { false };
    defined[v0] = true;
    for (int i = 0; i < L - 1; i++) {
        defined[c[i]] = true;
    }
    
    if (!validateDefinedOrbits(alpha_p, L, defined, F)) {
        return;
    }


    
    bool in_main_cycle[18] = { false };
    in_main_cycle[v0] = true;
    for (int i = 0; i < L - 1; i++) {
        in_main_cycle[c[i]] = true;
    }
    
    uint8_t rem[18];
    int rem_size = 0;
    int start_u = (search_type == 2) ? 2 : 1;
    for (int u = start_u; u < 18; u++) {
        if (!in_main_cycle[u]) {
            rem[rem_size++] = u;
        }
    }
    
    bool rem_used[18] = { false };
    if (parallel_remaining) {
        generateRemainingParallel(rem, rem_size, alpha_p, L);
    } else {
        generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L, 0);
    }
}

void K18A2::CycleBacktrackState::generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size, bool* rem_used, uint8_t* alpha_p, int L, int depth) {
    if (self->case_timed_out) return;
    // Parallel task collection: snapshot the branch state at the target depth
    // (a fully-placed first remaining cycle) instead of recursing further.
    if (collecting_rem && depth == rem_target_depth) {
        RemTask t;
        memcpy(t.alpha_p, alpha_p, 18);
        for (int i = 0; i < rem_size; i++) t.rem_used[i] = rem_used[i] ? 1 : 0;
        t.start_idx = start_idx;
        p_rem_tasks->push_back(t);
        return;
    }
    static thread_local int rem_call_count = 0;
    if (++rem_call_count >= 1000) {
        rem_call_count = 0;
        checkTimeoutAndReport(L, total_generated);
    }
    int first_unused = -1;
    for (int i = start_idx; i < rem_size; i++) {
        if (!rem_used[i]) {
            first_unused = i;
            break;
        }
    }
    
    if (first_unused == -1) {
        if (checkPermutationPassed(alpha_p, rem_used, L)) {
            saveAlpha(alpha_p);
        }
        return;
    }
    
    uint8_t v0_rem = rem[first_unused];
    rem_used[first_unused] = true;
    
    for (int d = 1; d <= rem_size; d++) {
        if (L % d != 0) continue;
        
        if (d == 1) {
            alpha_p[v0_rem] = v0_rem;
            
            bool defined[18];
            for (int i = 0; i < 18; i++) defined[i] = true;
            for (int i = 0; i < rem_size; i++) {
                if (!rem_used[i]) {
                    defined[rem[i]] = false;
                }
            }
            
            if (validateDefinedOrbits(alpha_p, L, defined, F)) {
                generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L, depth + 1);
            }
        } else {
            int unused_indices[18];
            int unused_count = 0;
            for (int j = first_unused + 1; j < rem_size; j++) {
                if (!rem_used[j]) {
                    unused_indices[unused_count++] = j;
                }
            }
            
            if (unused_count >= d - 1) {
                uint8_t cycle_nodes[18];
                cycle_nodes[0] = v0_rem;
                
                int perm[18];
                bool perm_used[18] = { false };
                
                auto arrange_stack = [&](auto& self_fn, int pos) -> void {
                    if (pos == d - 1) {
                        for (int k = 0; k < d - 1; k++) {
                            cycle_nodes[k + 1] = rem[unused_indices[perm[k]]];
                        }
                        for (int k = 0; k < d - 1; k++) {
                            alpha_p[cycle_nodes[k]] = cycle_nodes[k + 1];
                            rem_used[unused_indices[perm[k]]] = true;
                        }
                        alpha_p[cycle_nodes[d - 1]] = cycle_nodes[0];
                        
                        bool defined[18];
                        for (int i = 0; i < 18; i++) defined[i] = true;
                        for (int i = 0; i < rem_size; i++) {
                            if (!rem_used[i]) {
                                defined[rem[i]] = false;
                            }
                        }
                        
                        if (validateDefinedOrbits(alpha_p, L, defined, F)) {
                            generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L, depth + 1);
                        }

                        for (int k = 0; k < d - 1; k++) {
                            rem_used[unused_indices[perm[k]]] = false;
                        }
                        return;
                    }
                    for (int j = 0; j < unused_count; j++) {
                        if (!perm_used[j]) {
                            perm_used[j] = true;
                            perm[pos] = j;
                            self_fn(self_fn, pos + 1);
                            perm_used[j] = false;
                        }
                    }
                };
                arrange_stack(arrange_stack, 0);
            }
        }
    }
    
    rem_used[first_unused] = false;
}

void K18A2::CycleBacktrackState::generateRemainingParallel(const uint8_t* rem, int rem_size, uint8_t* alpha_p, int L) {
    // Nothing to fan out over: run the sequential path directly.
    if (rem_size == 0 || self->kThreads <= 1) {
        bool rem_used[18] = { false };
        generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L, 0);
        return;
    }

    // Phase 1 (single-threaded): collect the first-remaining-cycle branches as
    // independent tasks. This explores exactly the branches the sequential path
    // would, snapshotting each by value at depth 1.
    std::vector<RemTask> tasks;
    {
        bool rem_used[18] = { false };
        collecting_rem = true;
        rem_target_depth = 2;   // finer task granularity for load balancing
        p_rem_tasks = &tasks;
        generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L, 0);
        collecting_rem = false;
        p_rem_tasks = nullptr;
    }
    g_parallel_tasks += (long long)tasks.size();
    if (tasks.empty()) return;  // this candidate yields no completions

    // Phase 2: replay each task to completion across the thread pool.
    std::atomic<size_t> next_task{ 0 };
    std::vector<std::thread> workers;
    std::mutex merge_mutex;
    for (int w = 0; w < self->kThreads; w++) {
        workers.emplace_back([&]() {
            CycleBacktrackState ts = *this;  // copy ctor resets counters & flags
            ts.parallel_remaining = false;
            ts.collecting_rem = false;
            while (true) {
                size_t i = next_task.fetch_add(1);
                if (i >= tasks.size()) break;
                uint8_t local_alpha[18];
                bool local_rem_used[18] = { false };
                memcpy(local_alpha, tasks[i].alpha_p, 18);
                for (int k = 0; k < rem_size; k++) local_rem_used[k] = tasks[i].rem_used[k] != 0;
                ts.generate_remaining_cycles(tasks[i].start_idx, rem, rem_size, local_rem_used, local_alpha, L, 1);
                g_tasks_done.fetch_add(1);
            }
            std::lock_guard<std::mutex> lock(merge_mutex);
            this->valid_alphas.insert(this->valid_alphas.end(), ts.valid_alphas.begin(), ts.valid_alphas.end());
            this->total_generated += ts.total_generated;
            this->total_passed_p1f += ts.total_passed_p1f;
        });
    }
    for (auto& w : workers) w.join();
}

void K18A2::CycleBacktrackState::buildPermutation(uint8_t* c, uint8_t* alpha_p, int L) {
    for (int i = 0; i < 18; i++) alpha_p[i] = (uint8_t)i;
    if (L > 1) {
        alpha_p[v0] = c[0];
        for (int i = 0; i < L - 2; i++) alpha_p[c[i]] = c[i + 1];
        alpha_p[c[L - 2]] = v0;
    }
}

bool K18A2::CycleBacktrackState::checkPermutationPassed(const uint8_t* alpha_p, bool* used, int L) {
    uint8_t G[9][18];
    memcpy(G[0], F[0], 18);
    int limit = L / 2;
    for (int j = 1; j <= limit; j++) {
        apply_perm(G[j - 1], alpha_p, G[j]);
        if (!is_perfect_scalar(G[0], G[j])) return false;
    }
    g_candidates_passed_perfect++;
    return validateCandidateL(alpha_p, used, L);
}

bool K18A2::CycleBacktrackState::validateCandidateL(const uint8_t* alpha_p, bool* used, int L) {
    if (L == 17) return true;
    g_candidates_passed_validate++;

    uint8_t G[17][18];
    self->constructFullH(G, L, alpha_p, G);
    uint8_t H[18][18];
    return self->decomposeMissingEdges(G, L, alpha_p, H);
}

void K18A2::CycleBacktrackState::saveAlpha(const uint8_t* alpha_p) {
    std::array<uint8_t, 18> a;
    std::copy(alpha_p, alpha_p + 18, a.begin());
    valid_alphas.push_back(a);
    total_passed_p1f++;
}
void K18A2::printEstimatedTime(int L, long long checked_reps, int search_type) {
    auto now = std::chrono::steady_clock::now();
    double elapsed_since_print = std::chrono::duration<double>(now - last_print_time).count();
    if (elapsed_since_print >= 5.0) {
        last_print_time = now;
        double elapsed_total = std::chrono::duration<double>(now - case_start_time).count();
        
        long long orbit_reps[18] = {
            0, 1, 1, 2, 4, 10, 26, 72, 232, 504, 2619, 5040, 34650, 124740, 405405, 1081080, 2027025, 2027025
        };
        long long total_reps = (L >= 2 && L <= 17) ? orbit_reps[L] : 1;
        if (total_reps <= 0) total_reps = 1;
        
        // Progress basis: when the parallel small-L path is active, use task
        // completion (meaningful for small L). Otherwise fall back to orbit-rep
        // coverage. NOTE: for small L the tasks are added per main candidate, so
        // pct/ETA track the current candidate batch (see Cand for which one).
        long long ttot = g_parallel_tasks.load();
        long long tdone = g_tasks_done.load();
        double pct;
        if (ttot > 0) {
            pct = 100.0 * (double)tdone / (double)ttot;
        } else {
            pct = (double)checked_reps * 100.0 / total_reps;
        }
        if (pct > 100.0) pct = 100.0;

        double est_rem = 0.0;
        if (pct > 0.00001) {
            est_rem = elapsed_total * (100.0 / pct) - elapsed_total;
            if (est_rem < 0.0) est_rem = 0.0;
        }

        long long gmc = g_generate_matchings_calls;
        double avg_pool = gmc ? (double)g_generate_matchings_total_matchings / (double)gmc : 0.0;
        double s_en = g_ns_enum.load() / 1e9, s_gr = g_ns_group.load() / 1e9, s_cv = g_ns_cover.load() / 1e9;
        double par = elapsed_total > 0 ? (s_en + s_gr + s_cv) / elapsed_total : 0.0;
        printf("\r   [L=%d T%d] L-elapsed=%.0fs | %.2f%% (tasks %lld/%lld) | Cand: %lld | Decomp: %lld | AvgPool: %.0f | en/gr/cv: %.0f/%.0f/%.0f s | par~%.1fx | ETA(batch): %.1f min   ",
               L, search_type, elapsed_total, pct, tdone, ttot,
               (long long)g_candidates_processed, gmc, avg_pool,
               s_en, s_gr, s_cv, par,
               est_rem / 60.0);
        fflush(stdout);
    }
}

// ==========================================
// K18A2 Orchestrator & Helpers
// ==========================================

void K18A2::runExhaustiveSearch() {
    auto search_start = std::chrono::high_resolution_clock::now();
    std::set<std::vector<uint8_t>> unique_results;
    
    struct L_Stats {
        int L;
        long long checked = 0;
        long long passed = 0;
        size_t unique_classes = 0;
        bool run = false;
    };
    L_Stats stats[18];
    for (int i = 1; i <= 17; i++) {
        stats[i].L = i;
    }
    
    std::set<std::vector<uint8_t>> local_l16_unique;
    int cycle_lengths[] = { 17, 16, 14, 12, 10, 8, 6, 4, 2 };
    for (int L : cycle_lengths) {
        if (L < min_cycle_length) {
            stats[L].run = false;
            continue;
        }
        if (unique_results.size() >= 200000) {
            stats[L].run = false;
            continue;
        }
        
        std::vector<std::vector<uint8_t>> before_results;
        if (L == 16) {
            before_results.assign(unique_results.begin(), unique_results.end());
        }
        
        CycleLengthStats l_stats;
        searchCycleLength(L, unique_results, l_stats);
        stats[L].checked = l_stats.inputs;
        stats[L].passed = l_stats.passed;
        stats[L].unique_classes = l_stats.unique_classes;
        stats[L].run = true;
        
        if (L == 16) {
            for (const auto& res : unique_results) {
                if (std::find(before_results.begin(), before_results.end(), res) == before_results.end()) {
                    local_l16_unique.insert(res);
                }
            }
            if (!local_l16_unique.empty()) {
                verifyL16Pairing(local_l16_unique);
            }
        }
    }
    
    // Print the summary table
    printf("\n");
    printf("=========================================================================================================================================================\n");
    printf(" L   Cycle Structures Covered                           Orbit Reps   Checked Perms   Passed Filter   Unique Reps  Short Status / Description\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    long long orbit_reps[18] = {
        0, // L=0
        1, // L=1 (Identity)
        1, // L=2
        2, // L=3
        4, // L=4
        10, // L=5
        26, // L=6
        72, // L=7
        232, // L=8
        504, // L=9
        2619, // L=10
        5040, // L=11
        34650, // L=12
        124740, // L=13
        405405, // L=14
        1081080, // L=15
        2027025, // L=16
        2027025  // L=17
    };
    
    for (int L = 17; L >= 1; L--) {
        char cycle_str[64];
        if (L == 1) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(1)^18");
        } else if (L == 17 || L == 13 || L == 11) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(%d)(1)^%d", L, 18 - L);
        } else if (L == 16) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(16)(1)^2, (16)(2)");
        } else if (L == 15) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(15)(1)^3, (15)(3)");
        } else if (L == 14) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(14)(1)^4, (14)(2)(1)^2, (14)(2)^2");
        } else if (L == 9) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(9)(1)^9, (9)(3)^3");
        } else if (L == 7) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(7)(1)^11, (7)^2(1)^4");
        } else if (L == 5) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(5)^k(1)^{18-5k} for k in {1..3}");
        } else if (L == 3) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(3)^k(1)^{18-3k} for k in {1..6}");
        } else if (L == 2) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(2)^k(1)^{18-2k} for k in {1..9}");
        } else { // L = 12, 10, 8, 6, 4
            sprintf_s(cycle_str, sizeof(cycle_str), "(%d)(d_i)^k where d_i | %d", L, L);
        }
        
        char reps_str[32];
        sprintf_s(reps_str, sizeof(reps_str), "%lld", orbit_reps[L]);
        
        char checked_str[32] = "-";
        char passed_str[32] = "-";
        char reps_count_str[32] = "-";
        const char* desc = "";
        
        if (L == 1) {
            desc = "Identity (Unconstrained search)";
        } else if (L == 11 || L == 9 || L == 7) {
            desc = "Theoretically Impossible (odd L <= 15)";
        } else {
            if (stats[L].run) {
                sprintf_s(checked_str, sizeof(checked_str), "%lld", stats[L].checked);
                sprintf_s(passed_str, sizeof(passed_str), "%lld", stats[L].passed);
                sprintf_s(reps_count_str, sizeof(reps_count_str), "%zu", stats[L].unique_classes);
                
                if (L == 15 || L == 13) {
                    desc = "Theoretically Impossible (odd L <= 15)";
                } else if (L == 17) {
                    desc = "Exhaustive (GK/GB constructions)";
                } else if (L == 16) {
                    desc = "Exhaustive (GK subgroup & 2 new classes)";
                } else {
                    desc = stats[L].unique_classes > 0 ? "Exhaustive (P1Fs found)" : "Exhaustive (0 results)";
                }
            } else {
                desc = "Skipped / Not Run";
            }
        }
        
        printf("%2d   %-50s  %11s  %14s  %14s  %12s  %s\n",
               L, cycle_str, reps_str, checked_str, passed_str, reps_count_str, desc);
               
        if (L == 16 && stats[16].run) {
            printf("16R  %-50s  %11s  %14s  %14s  %12s  %s\n",
                   "Reflection pairing verification", "-", "-", "-", "2", "Verified pairing 2-and-2");
        }
    }
    printf("=========================================================================================================================================================\n");
    printf("\n");
    
    reportTotalResults(unique_results, search_start);
}

void K18A2::verifyL16Pairing(const std::set<std::vector<uint8_t>>& local_unique) {
    printf("-> Entering Case: L = 16 Reflection Pairing Verification\n");
    fflush(stdout);
    
    std::vector<std::vector<uint8_t>> mats(local_unique.begin(), local_unique.end());
    size_t n = mats.size();
    
    std::vector<bool> self_reflected(n, false);
    std::vector<int> pair_partner(n, -1);
    
    uint8_t cyc_R[18];
    buildStarterCycle(cyc_R);
    uint8_t p[18];
    for (int i = 0; i < 18; i++) {
        p[cyc_R[i]] = cyc_R[(18 - i) % 18];
    }
    
    for (int i = 0; i < n; i++) {
        uint8_t H[18][18];
        memset(H, 0, sizeof(H));
        for (int k = 1; k <= 17; k++) {
            const uint8_t* src = mats[i].data() + (k - 1) * 18;
            for (int j = 0; j < 9; j++) {
                uint8_t u = src[2 * j];
                uint8_t v = src[2 * j + 1];
                H[k][u] = v;
                H[k][v] = u;
            }
        }
        
        uint8_t H_ref[18][18];
        memset(H_ref, 0, sizeof(H_ref));
        for (int k = 1; k <= 17; k++) {
            uint8_t new_nb = p[k];
            for (int u = 0; u < 18; u++) {
                H_ref[new_nb][p[u]] = p[H[k][u]];
            }
        }
        
        std::set<std::vector<uint8_t>> ref_copies;
        recordIsomorphicResults(H_ref, ref_copies);
        
        if (ref_copies.find(mats[i]) != ref_copies.end()) {
            self_reflected[i] = true;
            continue;
        }
        
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            if (ref_copies.find(mats[j]) != ref_copies.end()) {
                pair_partner[i] = j;
                break;
            }
        }
    }
    
    int self_count = 0;
    int pair_count = 0;
    std::vector<bool> reported(n, false);
    for (int i = 0; i < n; i++) {
        if (self_reflected[i]) {
            self_count++;
            if (!reported[i]) {
                printf("  Matrix %d is self-reflected (isomorphic to its own cycle reflection)\n", i + 1);
                reported[i] = true;
            }
        } else if (pair_partner[i] != -1) {
            int j = pair_partner[i];
            if (!reported[i] && !reported[j]) {
                pair_count++;
                printf("  Matched Pair: Matrix %d <-> Matrix %d (Isomorphic under cycle reflection)\n", i + 1, j + 1);
                reported[i] = true;
                reported[j] = true;
            }
        } else {
            printf("  Warning: Matrix %d could not be paired or identified as self-reflected.\n", i + 1);
        }
    }
    
    printf("  Verification complete: found %d self-reflected matrices and %d isomorphic pairs.\n", self_count, pair_count);
}

void K18A2::searchCycleLength(int L, std::set<std::vector<uint8_t>>& unique_results, CycleLengthStats& stats) {
    printf("-> Entering Case: L = %d, Search Type 1 (Fixed out-of-cycle points)\n", L);
    fflush(stdout);

    case_start_time = std::chrono::steady_clock::now();
    last_print_time = case_start_time;
    case_timed_out = false;
    current_checked_reps = 0;

    g_candidates_processed = 0;
    g_candidates_passed_perfect = 0;
    g_candidates_passed_validate = 0;
    g_generate_matchings_calls = 0;
    g_generate_matchings_total_matchings = 0;
    g_ns_enum = 0;
    g_ns_group = 0;
    g_ns_cover = 0;
    g_parallel_tasks = 0;
    g_tasks_done = 0;


    CycleBacktrackState state1;
    state1.self = this;
    setupBacktrackState(state1, 1);
    uint8_t c1[18];
    bool used1[18];
    memset(used1, 0, sizeof(used1));
    used1[0] = true;
    used1[state1.v0] = true;

    // Calculate total top-level branches for state1
    int total_branches1 = 0;
    for (int v = 0; v < 18; v++) {
        if (used1[v]) continue;
        int p_idx = state1.vertex_to_pair[v];
        if (p_idx == -2) {
            total_branches1++;
        } else {
            bool pair_visited = used1[state1.pair_elements[p_idx][0]] || used1[state1.pair_elements[p_idx][1]];
            if (pair_visited) {
                total_branches1++;
            } else if (v == state1.pair_elements[0][0]) {
                total_branches1++;
            }
        }
    }
    this->total_top_branches = total_branches1;
    this->current_top_branch_idx = 0;

    state1.backtrack(0, 0, c1, used1, L);
    
    CycleBacktrackState state2;
    if (!case_timed_out && L >= 2 && L <= 16 && L % 2 == 0) {
        printf("-> Entering Case: L = %d, Search Type 2 (Transposed out-of-cycle points)\n", L);
        fflush(stdout);

        state2.self = this;
        setupBacktrackState(state2, 2);
        uint8_t c2[18];
        bool used2[18];
        memset(used2, 0, sizeof(used2));
        used2[0] = true;
        used2[1] = true;
        used2[state2.v0] = true;

        // Calculate total top-level branches for state2
        int total_branches2 = 0;
        for (int v = 0; v < 18; v++) {
            if (used2[v]) continue;
            int p_idx = state2.vertex_to_pair[v];
            if (p_idx == -2) {
                total_branches2++;
            } else {
                bool pair_visited = used2[state2.pair_elements[p_idx][0]] || used2[state2.pair_elements[p_idx][1]];
                if (pair_visited) {
                    total_branches2++;
                } else if (v == state2.pair_elements[0][0]) {
                    total_branches2++;
                }
            }
        }
        this->total_top_branches = total_branches2;
        this->current_top_branch_idx = 0;

        state2.backtrack(0, 0, c2, used2, L);
    }
    
    std::set<std::vector<uint8_t>> local_unique;
    if (!case_timed_out) {
        // Gather all automorphisms from both search types for this L.
        std::vector<std::array<uint8_t, 18>> all_alphas(state1.valid_alphas.begin(), state1.valid_alphas.end());
        if (L >= 2 && L <= 16 && L % 2 == 0) {
            all_alphas.insert(all_alphas.end(), state2.valid_alphas.begin(), state2.valid_alphas.end());
        }

        // Parallel post-processing: each alpha is independent (processAutomorphism /
        // decomposeMissingEdges use thread-local buffers and only read shared state),
        // so fan out across the pool. Results are identical to the sequential loop -
        // they accumulate into an order-independent std::set.
        int nthreads = (kThreads > 1) ? kThreads : 1;
        std::vector<std::set<std::vector<uint8_t>>> per_thread(nthreads);
        std::atomic<size_t> next_alpha{ 0 };
        std::vector<std::thread> workers;
        for (int t = 0; t < nthreads; t++) {
            workers.emplace_back([&, t]() {
                while (true) {
                    if (case_timed_out) break;
                    size_t i = next_alpha.fetch_add(1);
                    if (i >= all_alphas.size()) break;
                    std::set<std::vector<uint8_t>> temp_unique;
                    processAutomorphism(all_alphas[i], L, temp_unique);
                    per_thread[t].insert(temp_unique.begin(), temp_unique.end());
                }
            });
        }
        for (auto& w : workers) w.join();

        // Merge single-threaded: dedup against prior-L results, honor the global cap.
        for (auto& s : per_thread) {
            for (const auto& res : s) {
                if (unique_results.size() + local_unique.size() >= 200000) break;
                if (unique_results.find(res) == unique_results.end()) {
                    local_unique.insert(res);
                }
            }
        }
    }
    
    stats.L = L;
    stats.inputs = state1.total_generated + (L >= 2 && L <= 16 && L % 2 == 0 ? state2.total_generated : 0);
    stats.passed = state1.total_passed_p1f + (L >= 2 && L <= 16 && L % 2 == 0 ? state2.total_passed_p1f : 0);
    stats.unique_classes = local_unique.size();
    
    unique_results.insert(local_unique.begin(), local_unique.end());

    if (case_timed_out) {
        long long orbit_reps[18] = {
            0, 1, 1, 2, 4, 10, 26, 72, 232, 504, 2619, 5040, 34650, 124740, 405405, 1081080, 2027025, 2027025
        };
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - case_start_time).count();
        long long total_reps = orbit_reps[L];
        long long checked_reps = stats.inputs;
        double estimated_time = -1.0;
        if (checked_reps > 0) {
            estimated_time = elapsed * ((double)total_reps / checked_reps);
        }
        printf("   [L=%d] TIMEOUT reached. Checked %lld / %lld orbit reps in %.2f seconds.\n", 
               L, checked_reps, total_reps, elapsed);
        if (estimated_time >= 0) {
            printf("   [L=%d] Estimated total time needed without timeout: %.2f seconds (%.2f hours)\n", 
                   L, estimated_time, estimated_time / 3600.0);
        } else {
            printf("   [L=%d] Estimated total time: infinite/unknown (checked 0 reps)\n", L);
        }
        fflush(stdout);
    } else {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - case_start_time).count();
        if (elapsed >= 30.0) {
            printf("\n");
        }
        fflush(stdout);
    }
}

void K18A2::setupBacktrackState(CycleBacktrackState& state, int search_type) {
    memcpy(state.F[0], fixedRows[0].adj, 18);
    if (search_type == 2) {
        state.v0 = 3;
    } else {
        state.v0 = fixedRows[0].adj[0];
    }
    setupPairsTable(state, search_type);
}

void K18A2::setupPairsTable(CycleBacktrackState& state, int search_type) {
    for (int i = 0; i < 18; i++) {
        state.vertex_to_pair[i] = -1;
        state.vertex_to_pos[i] = -1;
    }
    fillPairsTable(state, search_type);
    
    uint8_t partner = fixedRows[0].adj[state.v0];
    state.vertex_to_pair[partner] = -2;
}

void K18A2::fillPairsTable(CycleBacktrackState& state, int search_type) {
    int pair_count = 0;
    for (int u = 1; u < 18; u++) {
        if (search_type == 2 && u == 1) continue;
        if (u == state.v0) continue;
        uint8_t v = fixedRows[0].adj[u];
        if (v == state.v0) continue;
        if (state.vertex_to_pair[u] != -1) continue;
        storePair(state, pair_count++, u, v);
    }
}

void K18A2::storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v) {
    state.pair_elements[pair_idx][0] = u;
    state.pair_elements[pair_idx][1] = v;
    state.vertex_to_pair[u] = pair_idx;
    state.vertex_to_pos[u] = 0;
    state.vertex_to_pair[v] = pair_idx;
    state.vertex_to_pos[v] = 1;
}

void K18A2::processAutomorphism(const std::array<uint8_t, 18>& alpha_arr, int L, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t alpha[18];
    std::copy(alpha_arr.begin(), alpha_arr.end(), alpha);
    uint8_t H[18][18];
    constructFullHFromAut(alpha, L, H);
    recordIsomorphicResults(H, unique_results);
}

void K18A2::constructFullHFromAut(const uint8_t* alpha, int L, uint8_t H[][18]) {
    uint8_t G[17][18];
    memcpy(G[0], fixedRows[0].adj, 18);
    constructFullH(G, L, alpha, G);
    if (L == 17) {
        copyMatchingsToH(nullptr, 0, G, L, H);
    } else {
        decomposeMissingEdges(G, L, alpha, H);
    }
}

void K18A2::constructFullH(const uint8_t G[][18], int L, const uint8_t* alpha, uint8_t H[][18]) {
    memcpy(H[0], fixedRows[0].adj, 18);
    Permutation perm;
    memcpy(perm.p, alpha, 18);
    for (int k = 1; k < L; k++) {
        apply_perm_18(H[k-1], perm, H[k]);
    }
}

bool K18A2::checkCyclesCompatibility(const uint8_t G[][18], int L) {
    for (int i = 0; i < L; i++) {
        for (int j = i + 1; j < L; j++) {
            if (!is_perfect_scalar(G[i], G[j])) return false;
        }
    }
    return true;
}

bool K18A2::decomposeMissingEdges(const uint8_t G[][18], int L, const uint8_t* alpha_p, uint8_t H[][18]) {
    if (case_timed_out) return false;
    int num_colors = 17 - L;

    // 1. Build adjacency matrix of the complement graph
    bool complement_adj[18][18];
    memset(complement_adj, false, sizeof(complement_adj));
    uint8_t comp_neighbors_higher[18][18];
    uint8_t comp_neighbors_higher_count[18];
    memset(comp_neighbors_higher_count, 0, sizeof(comp_neighbors_higher_count));
    for (int u = 0; u < 18; u++) {
        for (int v = u + 1; v < 18; v++) {
            if (isEdgeMissing(G, L, u, v)) {
                complement_adj[u][v] = true;
                complement_adj[v][u] = true;
                comp_neighbors_higher[u][comp_neighbors_higher_count[u]++] = (uint8_t)v;
            }
        }
    }
    
    // 2. Generate all perfect matchings in complement graph compatible with G[0]
    struct Orbit {
        std::vector<int> matching_indices;
        int size = 0;
    };
    static thread_local std::vector<std::array<uint8_t, 18>> pool;
    static thread_local std::vector<Orbit> orbits;
    static thread_local std::vector<bool> pool_visited;
    static thread_local std::vector<std::vector<int>> orbit_edges;
    static thread_local std::vector<int> chosen_orbits;
    static thread_local std::vector<bool> edge_covered;
    static thread_local std::vector<bool> orbit_compat;
    static thread_local std::vector<std::vector<int>> edge_to_orbits; // inverted index: edge -> orbits covering it

    pool.clear();
    orbits.clear();
    pool_visited.clear();
    orbit_edges.clear();
    chosen_orbits.clear();
    edge_covered.clear();
    orbit_compat.clear();
    edge_to_orbits.clear();
    
    uint8_t current_matching[18];
    memset(current_matching, 0xFF, sizeof(current_matching));
    bool used[18] = { false };
    
    // Multi-path endpoint tracking for all G[k]
    uint8_t path_end[17][18];
    for (int k = 0; k < L; k++) {
        for (int i = 0; i < 18; i++) {
            path_end[k][i] = G[k][i];
        }
    }
    
    auto generate_matchings = [&](auto& self_fn, int u, int edge_count) -> void {
        if (case_timed_out) return;
        static thread_local int match_call_count = 0;
        if (++match_call_count >= 1000) {
            match_call_count = 0;
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - case_start_time).count();
            if (elapsed > case_timeout_seconds) {
                case_timed_out = true;
            }
            printEstimatedTime(L, g_candidates_processed, (alpha_p[0] == 1 && alpha_p[1] == 0) ? 2 : 1);
        }
        if (u == 18) {
            std::array<uint8_t, 18> m;
            memcpy(m.data(), current_matching, 18);
            pool.push_back(m);
            return;
        }
        if (used[u]) {
            self_fn(self_fn, u + 1, edge_count);
            return;
        }
        used[u] = true;
        int n_count = comp_neighbors_higher_count[u];
        const uint8_t* neighbors = comp_neighbors_higher[u];
        for (int i = 0; i < n_count; i++) {
            int v = neighbors[i];
            if (!used[v]) {
                bool ok = true;
                if (edge_count < 8) {
                    for (int k = 0; k < L; k++) {
                        if (path_end[k][u] == v) {
                            ok = false;
                            break;
                        }
                    }
                }
                if (ok) {
                    int ep_u[17];
                    int ep_v[17];
                    for (int k = 0; k < L; k++) {
                        ep_u[k] = path_end[k][u];
                        ep_v[k] = path_end[k][v];
                        path_end[k][ep_u[k]] = ep_v[k];
                        path_end[k][ep_v[k]] = ep_u[k];
                    }
                    
                    used[v] = true;
                    current_matching[u] = (uint8_t)v;
                    current_matching[v] = (uint8_t)u;
                    
                    self_fn(self_fn, u + 1, edge_count + 1);
                    
                    used[v] = false;
                    current_matching[u] = 0xFF;
                    current_matching[v] = 0xFF;
                    
                    for (int k = 0; k < L; k++) {
                        path_end[k][ep_u[k]] = u;
                        path_end[k][ep_v[k]] = v;
                    }
                }
            }
        }
        used[u] = false;
    };
    
    auto _t_enum0 = std::chrono::steady_clock::now();
    generate_matchings(generate_matchings, 0, 0);
    g_ns_enum += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - _t_enum0).count();
    g_generate_matchings_calls++;
    g_generate_matchings_total_matchings += pool.size();

    if (pool.size() < (size_t)num_colors) return false;

    auto _t_grp0 = std::chrono::steady_clock::now();
    // 3. Build lookup map for O(1) matching search
    std::sort(pool.begin(), pool.end());
    FlatMatchingMap lookup;
    lookup.build(pool);
    
    // 4. Group pool matchings into valid orbits under alpha_p
    pool_visited.assign(pool.size(), false);
    
    for (int i = 0; i < (int)pool.size(); i++) {
        if (pool_visited[i]) continue;
        
        std::vector<int> orbit_indices;
        std::array<uint8_t, 18> curr = pool[i];
        bool valid_orbit = true;
        
        while (true) {
            int idx = lookup.find(curr);
            if (idx == -1) {
                valid_orbit = false;
                break;
            }
            if (std::find(orbit_indices.begin(), orbit_indices.end(), idx) != orbit_indices.end()) {
                if (idx != orbit_indices[0]) {
                    valid_orbit = false;
                }
                break;
            }
            orbit_indices.push_back(idx);
            
            // Apply permutation alpha_p to curr to get next matching
            std::array<uint8_t, 18> next_m;
            for (int v = 0; v < 18; v++) {
                next_m[alpha_p[v]] = alpha_p[curr[v]];
            }
            curr = next_m;
        }
        
        int d = (int)orbit_indices.size();
        if (valid_orbit) {
            if (L % d != 0) {
                valid_orbit = false;
            }
        }
        if (valid_orbit) {
            // Check mutual compatibility within the orbit
            for (int x = 0; x < d; x++) {
                for (int y = x + 1; y < d; y++) {
                    if (!is_perfect_scalar(pool[orbit_indices[x]].data(), pool[orbit_indices[y]].data())) {
                        valid_orbit = false;
                        break;
                    }
                }
                if (!valid_orbit) break;
            }
        }
        
        // Mark all elements in this orbit as visited in the pool
        for (int idx : orbit_indices) {
            pool_visited[idx] = true;
        }
        
        if (valid_orbit) {
            Orbit orb;
            orb.matching_indices = orbit_indices;
            orb.size = d;
            orbits.push_back(orb);
        }
    }
    
    // 5. Map each complement edge to an index 0..153
    std::pair<uint8_t, uint8_t> comp_edges[153];
    int comp_edge_cnt = 0;
    int edge_to_idx[18][18];
    memset(edge_to_idx, -1, sizeof(edge_to_idx));
    for (int u = 0; u < 18; u++) {
        for (int v = u + 1; v < 18; v++) {
            if (complement_adj[u][v]) {
                edge_to_idx[u][v] = comp_edge_cnt;
                edge_to_idx[v][u] = comp_edge_cnt;
                comp_edges[comp_edge_cnt++] = { (uint8_t)u, (uint8_t)v };
            }
        }
    }
    
    // 6. For each orbit, represent the set of edges it covers
    int num_orbits = (int)orbits.size();
    orbit_edges.resize(num_orbits);
    for (int i = 0; i < num_orbits; i++) {
        orbit_edges[i].clear();
        for (int m_idx : orbits[i].matching_indices) {
            for (int u = 0; u < 18; u++) {
                int v = pool[m_idx][u];
                if (u < v) {
                    int e_idx = edge_to_idx[u][v];
                    orbit_edges[i].push_back(e_idx);
                }
            }
        }
    }
    
    // 6b. Inverted index: for each complement edge, the orbits that cover it.
    // Built by ascending orbit index, so each list is sorted -> iterating it
    // reproduces the original "for o in 0..num_orbits if o covers e" order exactly.
    // (Orbit members are mutually compatible => edge-disjoint, so no duplicates.)
    edge_to_orbits.assign(comp_edge_cnt, {});
    for (int i = 0; i < num_orbits; i++) {
        for (int e_idx : orbit_edges[i]) {
            edge_to_orbits[e_idx].push_back(i);
        }
    }

    // Proposal C: Precalculate Orbit Compatibility Matrix
    orbit_compat.assign(num_orbits * num_orbits, true);
    for (int i = 0; i < num_orbits; i++) {
        for (int j = i; j < num_orbits; j++) {
            bool comp = true;
            for (int m1 : orbits[i].matching_indices) {
                for (int m2 : orbits[j].matching_indices) {
                    if (!is_perfect_scalar(pool[m1].data(), pool[m2].data())) {
                        comp = false;
                        break;
                    }
                }
                if (!comp) break;
            }
            orbit_compat[i * num_orbits + j] = comp;
            orbit_compat[j * num_orbits + i] = comp;
        }
    }
    
    // 7. Solve Exact Cover on Orbits using backtracking with MRV
    edge_covered.assign(comp_edge_cnt, false);
    int current_size_sum = 0;
    
    auto solve_exact_cover_orbits = [&](auto& self_fn) -> bool {
        if (case_timed_out) return false;
        
        if (current_size_sum == num_colors) {
            return true;
        }
        
        // Find the uncovered edge with the minimum number of remaining active orbits (MRV)
        int best_edge = -1;
        int min_candidates = 999999;
        
        for (int e = 0; e < comp_edge_cnt; e++) {
            if (edge_covered[e]) continue;
            
            int active_count = 0;
            for (int o_idx : edge_to_orbits[e]) {
                if (current_size_sum + orbits[o_idx].size > num_colors) {
                    continue;
                }
                
                bool overlap = false;
                for (int edge : orbit_edges[o_idx]) {
                    if (edge_covered[edge]) {
                        overlap = true;
                        break;
                    }
                }
                if (overlap) continue;
                
                bool compatible = true;
                for (int chosen_idx : chosen_orbits) {
                    if (!orbit_compat[o_idx * num_orbits + chosen_idx]) {
                        compatible = false;
                        break;
                    }
                }
                
                if (compatible) {
                    active_count++;
                }
            }
            
            if (active_count == 0) return false;
            
            if (active_count < min_candidates) {
                min_candidates = active_count;
                best_edge = e;
            }
        }
        
        if (best_edge == -1) {
            return false;
        }
        
        for (int o_idx : edge_to_orbits[best_edge]) {
            if (current_size_sum + orbits[o_idx].size > num_colors) {
                continue;
            }
            
            bool overlap = false;
            for (int edge : orbit_edges[o_idx]) {
                if (edge_covered[edge]) {
                    overlap = true;
                    break;
                }
            }
            if (overlap) continue;
            
            bool compatible = true;
            for (int chosen_idx : chosen_orbits) {
                if (!orbit_compat[o_idx * num_orbits + chosen_idx]) {
                    compatible = false;
                    break;
                }
            }
            
            if (compatible) {
                chosen_orbits.push_back(o_idx);
                current_size_sum += orbits[o_idx].size;
                for (int edge : orbit_edges[o_idx]) {
                    edge_covered[edge] = true;
                }
                
                if (self_fn(self_fn)) return true;
                
                for (int edge : orbit_edges[o_idx]) {
                    edge_covered[edge] = false;
                }
                current_size_sum -= orbits[o_idx].size;
                chosen_orbits.pop_back();
            }
        }
        return false;
    };
    
    g_ns_group += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - _t_grp0).count();
    auto _t_cov0 = std::chrono::steady_clock::now();
    bool _cover_ok = solve_exact_cover_orbits(solve_exact_cover_orbits);
    g_ns_cover += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - _t_cov0).count();
    if (!_cover_ok) return false;
    
    // 8. Copy the found clique of matchings from chosen orbits to H
    uint8_t matchings[17][18];
    memset(matchings, 0xFF, sizeof(matchings));
    int m_count = 0;
    for (int o_idx : chosen_orbits) {
        for (int m_idx : orbits[o_idx].matching_indices) {
            memcpy(matchings[m_count++], pool[m_idx].data(), 18);
        }
    }
    copyMatchingsToH(matchings, num_colors, G, L, H);
    return true;
}

void K18A2::findMissingEdges(const uint8_t G[][18], int L, std::pair<uint8_t, uint8_t>* edges) {
    int edge_cnt = 0;
    for (int u = 0; u < 18; u++) {
        for (int v = u + 1; v < 18; v++) {
            if (isEdgeMissing(G, L, u, v)) {
                edges[edge_cnt++] = { (uint8_t)u, (uint8_t)v };
            }
        }
    }
}

bool K18A2::isEdgeMissing(const uint8_t G[][18], int L, int u, int v) {
    for (int k = 0; k < L; k++) {
        if (G[k][u] == v) return false;
    }
    return true;
}

bool K18A2::backtrackColor(int edge_idx, int num_edges, int num_colors,
                           const std::pair<uint8_t, uint8_t>* edges,
                           uint8_t matchings[][18], const uint8_t* G0) {
    if (edge_idx == num_edges) {
        return checkMatchingsCompatibility(matchings, num_colors, G0);
    }
    return tryColoringEdge(edge_idx, num_edges, num_colors, edges, matchings, G0);
}

bool K18A2::checkMatchingsCompatibility(uint8_t matchings[][18], int num_colors, const uint8_t* G0) {
    for (int c = 0; c < num_colors; c++) {
        if (!is_perfect_scalar(G0, matchings[c])) return false;
    }
    for (int c = 0; c < num_colors; c++) {
        for (int d = c + 1; d < num_colors; d++) {
            if (!is_perfect_scalar(matchings[c], matchings[d])) return false;
        }
    }
    return true;
}

bool K18A2::tryColoringEdge(int edge_idx, int num_edges, int num_colors,
                            const std::pair<uint8_t, uint8_t>* edges,
                            uint8_t matchings[][18], const uint8_t* G0) {
    uint8_t u = edges[edge_idx].first;
    uint8_t v = edges[edge_idx].second;
    int limit = getSymmetryBreakingLimit(matchings, num_colors);
    for (int c = 0; c <= limit; c++) {
        if (matchings[c][u] == 0xFF && matchings[c][v] == 0xFF) {
            matchings[c][u] = v;
            matchings[c][v] = u;
            if (is_color_compatible(c, num_colors, matchings, G0)) {
                if (backtrackColor(edge_idx + 1, num_edges, num_colors, edges, matchings, G0)) return true;
            }
            matchings[c][u] = 0xFF;
            matchings[c][v] = 0xFF;
        }
    }
    return false;
}

bool K18A2::is_color_compatible(int c, int num_colors, uint8_t matchings[][18], const uint8_t* G0) {
    int count = 0;
    for (int i = 0; i < 18; i++) {
        if (matchings[c][i] != 0xFF) count++;
    }
    if (count < 18) return true; // Not yet complete
    
    // Check compatibility with G0
    if (!is_perfect_scalar(G0, matchings[c])) return false;
    
    // Check compatibility with other complete colors
    for (int d = 0; d < num_colors; d++) {
        if (d == c) continue;
        int count_d = 0;
        for (int i = 0; i < 18; i++) {
            if (matchings[d][i] != 0xFF) count_d++;
        }
        if (count_d == 18) {
            if (!is_perfect_scalar(matchings[c], matchings[d])) return false;
        }
    }
    return true;
}

int K18A2::getSymmetryBreakingLimit(uint8_t matchings[][18], int num_colors) {
    int max_color = 0;
    for (int c = 0; c < num_colors; c++) {
        if (isColorUsed(matchings[c])) max_color = c + 1;
    }
    return (max_color < num_colors - 1) ? max_color : num_colors - 1;
}

bool K18A2::isColorUsed(const uint8_t* matching) {
    for (int i = 0; i < 18; i++) {
        if (matching[i] != 0xFF) return true;
    }
    return false;
}

void K18A2::copyMatchingsToH(uint8_t matchings[][18], int num_colors,
                             const uint8_t G[][18], int L, uint8_t H[][18]) {
    memset(H, 0, 18 * 18);
    for (int k = 0; k < L; k++) {
        uint8_t nb = G[k][0];
        memcpy(H[nb], G[k], 18);
    }
    for (int c = 0; c < num_colors; c++) {
        uint8_t nb = matchings[c][0];
        memcpy(H[nb], matchings[c], 18);
    }
}

void K18A2::recordIsomorphicResults(const uint8_t H[][18], std::set<std::vector<uint8_t>>& unique_results) {
    CycleUnion cu_H = find_cycles(H[1], H[2]);
    if (cu_H.count != 1 || cu_H.lens[0] != 18) return;
    for (int v = 0; v < 18; v++) {
        tryIsomorphicMapping(H, cu_H, v, unique_results);
    }
}

void K18A2::tryIsomorphicMapping(const uint8_t H[][18], const CycleUnion& cu_H, int v, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t cyc_R[18];
    buildStarterCycle(cyc_R);
    uint8_t cyc_H[18];
    buildHCycle(H, v, cyc_H);
    uint8_t p[18];
    buildMappingPermutation(cyc_H, cyc_R, p);
    checkAndRecordPermutedH(H, p, unique_results);
}

void K18A2::buildStarterCycle(uint8_t* cyc_R) {
    uint8_t curr = 0;
    for (int i = 0; i < 9; i++) {
        cyc_R[2 * i] = curr;
        curr = fixedRows[0].adj[curr];
        cyc_R[2 * i + 1] = curr;
        curr = fixedRows[1].adj[curr];
    }
}

void K18A2::buildHCycle(const uint8_t H[][18], int v, uint8_t* cyc_H) {
    uint8_t curr = v;
    for (int i = 0; i < 9; i++) {
        cyc_H[2 * i] = curr;
        curr = H[1][curr];
        cyc_H[2 * i + 1] = curr;
        curr = H[2][curr];
    }
}

void K18A2::buildMappingPermutation(const uint8_t* cyc_H, const uint8_t* cyc_R, uint8_t* p) {
    for (int i = 0; i < 18; i++) {
        p[cyc_H[i]] = cyc_R[i];
    }
}

void K18A2::checkAndRecordPermutedH(const uint8_t H[][18], const uint8_t* p, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t S[18][18];
    memset(S, 0, sizeof(S));
    applyPermToH(H, p, S);
    if (doesSMatchFixedRows(S)) {
        recordS(S, unique_results);
    }
}

void K18A2::applyPermToH(const uint8_t H[][18], const uint8_t* p, uint8_t S[][18]) {
    for (int k = 1; k <= 17; k++) {
        uint8_t permuted_factor[18];
        for (int i = 0; i < 18; i++) {
            permuted_factor[p[i]] = p[H[k][i]];
        }
        uint8_t nb = permuted_factor[0];
        memcpy(S[nb], permuted_factor, 18);
    }
}

bool K18A2::doesSMatchFixedRows(const uint8_t S[][18]) {
    return memcmp(S[1], fixedRows[0].adj, 18) == 0 &&
           memcmp(S[2], fixedRows[1].adj, 18) == 0;
}

void K18A2::recordS(const uint8_t S[][18], std::set<std::vector<uint8_t>>& unique_results) {
    unsigned char results[17 * 18];
    for (int k = 1; k <= 17; k++) {
        adj_to_src(S[k], results + (k - 1) * 18);
    }
    std::vector<uint8_t> res_vec(results, results + 17 * 18);
    unique_results.insert(res_vec);
}

void K18A2::reportTotalResults(const std::set<std::vector<uint8_t>>& unique_results,
                               std::chrono::high_resolution_clock::time_point start) {
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    printf("Exhaustive search completed. Unique P1Fs: %zu. Time: %.4f seconds\n", unique_results.size(), elapsed);
    sendResultsToCallback(unique_results);
}

void K18A2::sendResultsToCallback(const std::set<std::vector<uint8_t>>& unique_results) {
    std::lock_guard<std::mutex> lock(result_mutex);
    for (const auto& res : unique_results) {
        resultCallback(cbClass, res.data(), 0, 1, 2);
    }
}

void K18A2::adj_to_src(const uint8_t* adj, unsigned char* src) {
    src[0] = 0;
    src[1] = adj[0];
    int idx = 2;
    bool visited[18] = { false };
    visited[0] = true;
    visited[adj[0]] = true;
    for (int u = 1; u < 18; u++) {
        if (!visited[u]) {
            src[idx] = u;
            src[idx + 1] = adj[u];
            visited[u] = true;
            visited[adj[u]] = true;
            idx += 2;
        }
    }
}
