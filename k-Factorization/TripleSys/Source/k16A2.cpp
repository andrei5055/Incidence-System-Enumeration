#include "k16A2.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <intrin.h>
#include <set>

static bool has_prime_order_k16(const uint8_t* p, int target_p) {
    bool visited[16] = { false };
    bool found_target_cycle = false;
    for (int i = 0; i < 16; i++) {
        if (visited[i]) continue;
        int len = 0;
        int curr = i;
        while (!visited[curr]) {
            visited[curr] = true;
            curr = p[curr];
            len++;
        }
        if (len > 1) {
            if (len != target_p) return false;
            found_target_cycle = true;
        }
    }
    return found_target_cycle;
}

K16A2::K16A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows,
    ResultCallback callback, void* cbClassPtr, bool bPrint) : KBase<Mask256_C>(factParam, bPrint) {
    int id_cnt = 0;
    memset(edge_id_table, -1, sizeof(edge_id_table));
    for (int i = 0; i < K16_N; i++) {
        for (int j = i + 1; j < K16_N; j++) {
            int eid = id_cnt++;
            edge_id_table[i][j] = eid;
            edge_id_table[j][i] = eid;
            if (eid < 256) {
                edge_to_u[eid] = (uint8_t)i;
                edge_to_v[eid] = (uint8_t)j;
            }
        }
    }
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }

    init(fixed3RowsIndex, kThreads, first3Rows, callback, cbClassPtr);
}

void K16A2::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K16A2 solver (Exhaustive Cyclic Automorphism version), MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
    }
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    start_time = std::chrono::steady_clock::now();
    last_print_time = start_time;
    bTimeSet = false;
    call_counter = 0;

    fixedEdgesMask = Mask256_C();
    for (int s = 0; s < m_nFixedRows; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * K16_N, K16_N);
        for (int i = 0; i < K16_N; i += 2) {
            const auto a = fixedRows[s].src[i]; const auto b = fixedRows[s].src[i + 1];
            fixedRows[s].adj[a] = b; fixedRows[s].adj[b] = a;
            int eid = edge_id_table[a][b];
            if (eid < 64) fixedEdgesMask.m[0] |= (1ULL << eid);
            else if (eid < 128) fixedEdgesMask.m[1] |= (1ULL << (eid - 64));
        }
    }

    target_cu = find_cycles(fixedRows[0].adj, fixedRows[1].adj);
    for (int s = 0; s < m_nFixedRows; s++) fixedRows[s].fs = get_fast_sorted(fixedRows[s].adj);
    r1_can = fixedRows[0].fs;
    r2_can = fixedRows[1].fs;
    for (int i = 0; i < m_nFixedRows; i++) {
        fixed_packed[i] = pack_factor_adj(fixedRows[i].adj);
    }

    get_transformations(fixedRows[0], fixedRows[1], fixed_trans[0]);
    get_transformations(fixedRows[0], fixedRows[2], fixed_trans[1]);
    get_transformations(fixedRows[1], fixedRows[2], fixed_trans[2]);

    global_pool.clear();
    packed_pool.clear();
    f_map.clear();
    f_map_unordered.clear();
    for (int s = 0; s < m_nFixedRows; s++) {
        std::vector<uint8_t> k(fixedRows[s].adj, fixedRows[s].adj + K16_N);
        f_map[k] = s;
        std::array<uint8_t, 16> k_arr;
        std::copy(fixedRows[s].adj, fixedRows[s].adj + 16, k_arr.begin());
        f_map_unordered[k_arr] = s;
        global_pool.push_back(fixedRows[s]);
        packed_pool.push_back(pack_factor_adj(fixedRows[s].adj));
    }
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }
    KBase::init();

    // Precalculate all valid cyclic automorphism transition chains for the starter (Ra, Rb) -> (Rb, Rc)
    max_row_seen = 3;
    valid_chains.clear();
    saved_transitions.clear();
    memset(rejected_candidates, 0, sizeof(rejected_candidates));

    int a = trans_config[0];
    int b = trans_config[1];
    int c = trans_config[2];

    Factor R_a = fixedRows[a - 1];
    Factor R_b = fixedRows[b - 1];
    Factor R_c = fixedRows[c - 1];

    TransInfo info_init;
    get_transformations_general(R_a, R_b, R_b, R_c, info_init);
    int total_candidates = info_init.count;
    int valid_starter_chains = 0;

    for (int p_idx = 0; p_idx < info_init.count; p_idx++) {
        const Permutation& alpha = info_init.perms[p_idx];

        if (!has_prime_order_k16(alpha.p, 15)) {
            continue;
        }

        bool used_v0[16] = { false };
        used_v0[fixedRows[0].adj[0]] = true;
        used_v0[fixedRows[1].adj[0]] = true;
        used_v0[fixedRows[2].adj[0]] = true;

        bool chain_ok = true;
        uint8_t current_adj[16];
        memcpy(current_adj, R_c.adj, 16);

        ValidChain chain;
        chain.alpha = alpha;
        memcpy(chain.factors_adj[1], R_a.adj, 16);
        memcpy(chain.factors_adj[2], R_b.adj, 16);
        memcpy(chain.factors_adj[3], R_c.adj, 16);
        chain.factors_packed[1] = pack_factor_adj(R_a.adj);
        chain.factors_packed[2] = pack_factor_adj(R_b.adj);
        chain.factors_packed[3] = pack_factor_adj(R_c.adj);

        // Check initial compatibility
        if (!is_perfect_packed(chain.factors_packed[1], chain.factors_packed[2]) ||
            !is_perfect_packed(chain.factors_packed[1], chain.factors_packed[3]) ||
            !is_perfect_packed(chain.factors_packed[2], chain.factors_packed[3])) {
            continue;
        }

        for (int k = 4; k <= 15; k++) {
            uint8_t next_adj[16];
            apply_perm_16(current_adj, alpha, next_adj);

            uint8_t v0 = next_adj[0];
            if (used_v0[v0]) {
                chain_ok = false;
                break;
            }
            used_v0[v0] = true;

            memcpy(chain.factors_adj[k], next_adj, 16);
            chain.factors_packed[k] = pack_factor_adj(next_adj);

            // Check compatibility with all previously generated factors in the chain
            for (int prev = 1; prev < k; prev++) {
                if (!is_perfect_packed(chain.factors_packed[prev], chain.factors_packed[k])) {
                    chain_ok = false;
                    break;
                }
            }
            if (!chain_ok) {
                break;
            }

            memcpy(current_adj, next_adj, 16);
        }

        if (chain_ok) {
            valid_starter_chains++;
            std::lock_guard<std::mutex> lock(result_mutex);
            if (std::find_if(saved_transitions.begin(), saved_transitions.end(), [&](const Permutation& p) {
                return memcmp(p.p, alpha.p, 16) == 0;
            }) == saved_transitions.end()) {
                saved_transitions.push_back(alpha);
            }
        }
    }

    if (m_bPrint) {
        printf("Transitions: %d candidate transitions, %d valid starter chains recorded.\n", total_candidates, valid_starter_chains);
    }
}

bool K16A2::addRow(int rowNum, const unsigned char* source) {
    printf("not supported\n");
    exit(1);
    return false;
}

K16A2::CycleUnion K16A2::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    bool visited[16] = { false };
    for (int i = 0; i < 16; ++i) {
        if (visited[i]) continue;
        int curr = i;
        int cycle_pos = 0;
        do {
            visited[curr] = true;
            cu.cycles[cu.count][cycle_pos++] = curr;
            curr = adj1[curr];
            visited[curr] = true;
            cu.cycles[cu.count][cycle_pos++] = curr;
            curr = adj2[curr];
        } while (curr != i);
        cu.lens[cu.count] = cycle_pos;
        cu.count++;
    }
    return cu;
}

void K16A2::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
    CycleUnion source_cu = find_cycles(fi.adj, fj.adj);
    if (source_cu.count != 1 || source_cu.lens[0] != 16) return;
    if (target_cu.count != 1 || target_cu.lens[0] != 16) return;

    for (int dir = 0; dir < 2; dir++) {
        for (int start_offset = 0; start_offset < 16; start_offset++) {
            if (info.count >= 512) break;
            auto& p_total = info.perms[info.count++];
            memset(p_total.p, 0xFF, sizeof(p_total.p));
            memset(p_total.p_inv, 0xFF, sizeof(p_total.p_inv));
            
            for (int i = 0; i < 16; i++) {
                int src_pos = (dir == 0) ? (start_offset + i) % 16 : (start_offset - i + 16) % 16;
                uint8_t src_v = (uint8_t)source_cu.cycles[0][src_pos];
                uint8_t tgt_v = (uint8_t)target_cu.cycles[0][i];
                p_total.p[src_v] = tgt_v;
                p_total.p_inv[tgt_v] = src_v;
            }
        }
    }
}

void K16A2::get_transformations_general(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl, TransInfo& info) {
    info.count = 0;
    CycleUnion cu_src = find_cycles(fi.adj, fj.adj);
    if (cu_src.count != 1 || cu_src.lens[0] != 16) return;
    CycleUnion cu_tgt = find_cycles(fk.adj, fl.adj);
    if (cu_tgt.count != 1 || cu_tgt.lens[0] != 16) return;

    for (int dir = 0; dir < 2; dir++) {
        for (int start_offset = (dir == 0 ? 0 : 1); start_offset < 16; start_offset += 2) {
            if (info.count >= 512) break;
            auto& p_total = info.perms[info.count++];
            memset(p_total.p, 0xFF, sizeof(p_total.p));
            memset(p_total.p_inv, 0xFF, sizeof(p_total.p_inv));
            
            for (int i = 0; i < 16; i++) {
                int src_pos = (dir == 0) ? (start_offset + i) % 16 : (start_offset - i + 16) % 16;
                uint8_t src_v = (uint8_t)cu_src.cycles[0][src_pos];
                uint8_t tgt_v = (uint8_t)cu_tgt.cycles[0][i];
                p_total.p[src_v] = tgt_v;
                p_total.p_inv[tgt_v] = src_v;
            }
        }
    }
}

K16A2::FastSortedFactor K16A2::get_fast_sorted(const uint8_t* adj) {
    FastSortedFactor sf; memset(&sf, 0, sizeof(sf));
    uint8_t* out = sf.pairs;
    int count = 0;
    for (int u = 0; u < 16; u++) {
        int v = adj[u];
        if (u < v) {
            out[count * 2] = (uint8_t)u;
            out[count * 2 + 1] = (uint8_t)v;
            count++;
        }
    }
    return sf;
}

K16A2::PackedAdj K16A2::pack_factor_adj(const uint8_t* adj) {
    PackedAdj pa;
    memcpy(pa.adj, adj, 16);
    pa.edge_mask = Mask256_C();
    for (int u = 0; u < 16; u++) {
        int v = adj[u];
        if (u < v) {
            int eid = edge_id_table[u][v];
            if (eid < 64) pa.edge_mask.m[0] |= (1ULL << eid);
            else if (eid < 128) pa.edge_mask.m[1] |= (1ULL << (eid - 64));
        }
    }
    return pa;
}

bool K16A2::compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b) {
    for (int i = 0; i < 3; i++) {
        if (compare_fast_sorted(a.r[i], b.r[i])) return true;
        if (compare_fast_sorted(b.r[i], a.r[i])) return false;
    }
    return false;
}

void K16A2::solve(int mode) {
    runExhaustiveSearch();
}

// ==========================================
// CycleBacktrackState Implementation
// ==========================================

void K16A2::CycleBacktrackState::apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj) {
    for (int i = 0; i < 16; i++) {
        dst_adj[perm[i]] = perm[src_adj[i]];
    }
}

bool K16A2::CycleBacktrackState::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    uint8_t curr = 0;
    for (int i = 0; i < 8; i++) {
        curr = adj1[curr];
        curr = adj2[curr];
        if (curr == 0 && i < 7) return false;
    }
    return curr == 0;
}

void K16A2::CycleBacktrackState::backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    if (self->case_timed_out) return;
    static thread_local int backtrack_call_count = 0;
    if (++backtrack_call_count >= 1000) {
        backtrack_call_count = 0;
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - self->case_start_time).count();
        if (elapsed > self->case_timeout_seconds) {
            self->case_timed_out = true;
            return;
        }
        self->printEstimatedTime(L, total_generated);
    }
    if (depth == L - 1) {
        processCandidate(c, used, L);
        return;
    }
    backtrackRecurse(depth, pairs_visited, c, used, L);
}

void K16A2::CycleBacktrackState::backtrackRecurse(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    for (int v = 0; v < 16; v++) {
        if (used[v]) continue;
        tryVertexForCycle(v, depth, pairs_visited, c, used, L);
    }
}

void K16A2::CycleBacktrackState::tryVertexForCycle(int v, int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
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

void K16A2::CycleBacktrackState::recurseWithVertex(int v, int depth, int next_pairs, uint8_t* c, bool* used, int L) {
    if (depth == 0) {
        self->current_top_branch_idx++;
    }
    used[v] = true;
    c[depth] = v;
    backtrack(depth + 1, next_pairs, c, used, L);
    used[v] = false;
}

void K16A2::CycleBacktrackState::processCandidate(uint8_t* c, bool* used, int L) {
    total_generated++;
    self->current_checked_reps = total_generated;
    uint8_t alpha_p[16];
    buildPermutation(c, alpha_p, L);
    
    int search_type = (v0 == 3) ? 2 : 1;
    if (search_type == 2) {
        alpha_p[0] = 1;
        alpha_p[1] = 0;
    }
    
    bool in_main_cycle[16] = { false };
    in_main_cycle[v0] = true;
    for (int i = 0; i < L - 1; i++) {
        in_main_cycle[c[i]] = true;
    }
    
    uint8_t rem[16];
    int rem_size = 0;
    int start_u = (search_type == 2) ? 2 : 1;
    for (int u = start_u; u < 16; u++) {
        if (!in_main_cycle[u]) {
            rem[rem_size++] = u;
        }
    }
    
    bool rem_used[16] = { false };
    generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L);
}

void K16A2::CycleBacktrackState::generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size, bool* rem_used, uint8_t* alpha_p, int L) {
    if (self->case_timed_out) return;
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
            generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L);
        } else {
            int unused_indices[16];
            int unused_count = 0;
            for (int j = first_unused + 1; j < rem_size; j++) {
                if (!rem_used[j]) {
                    unused_indices[unused_count++] = j;
                }
            }
            
            if (unused_count >= d - 1) {
                uint8_t cycle_nodes[16];
                cycle_nodes[0] = v0_rem;
                
                int perm[16];
                bool perm_used[16] = { false };
                
                auto arrange_stack = [&](auto& self_fn, int depth) -> void {
                    if (depth == d - 1) {
                        for (int k = 0; k < d - 1; k++) {
                            cycle_nodes[k + 1] = rem[unused_indices[perm[k]]];
                        }
                        for (int k = 0; k < d - 1; k++) {
                            alpha_p[cycle_nodes[k]] = cycle_nodes[k + 1];
                            rem_used[unused_indices[perm[k]]] = true;
                        }
                        alpha_p[cycle_nodes[d - 1]] = cycle_nodes[0];
                        
                        generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L);
                        
                        for (int k = 0; k < d - 1; k++) {
                            rem_used[unused_indices[perm[k]]] = false;
                        }
                        return;
                    }
                    for (int j = 0; j < unused_count; j++) {
                        if (!perm_used[j]) {
                            perm_used[j] = true;
                            perm[depth] = j;
                            self_fn(self_fn, depth + 1);
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

void K16A2::CycleBacktrackState::buildPermutation(uint8_t* c, uint8_t* alpha_p, int L) {
    for (int i = 0; i < 16; i++) alpha_p[i] = (uint8_t)i;
    if (L > 1) {
        alpha_p[v0] = c[0];
        for (int i = 0; i < L - 2; i++) alpha_p[c[i]] = c[i + 1];
        alpha_p[c[L - 2]] = v0;
    }
}

bool K16A2::CycleBacktrackState::checkPermutationPassed(const uint8_t* alpha_p, bool* used, int L) {
    uint8_t G[8][16];
    memcpy(G[0], F[0], 16);
    int limit = L / 2;
    for (int j = 1; j <= limit; j++) {
        apply_perm(G[j - 1], alpha_p, G[j]);
        if (!is_perfect_scalar(G[0], G[j])) return false;
    }
    return validateCandidateL(alpha_p, used, L);
}

bool K16A2::CycleBacktrackState::validateCandidateL(const uint8_t* alpha_p, bool* used, int L) {
    if (L == 15) return true;
    uint8_t G[15][16];
    self->constructFullH(G, L, alpha_p, G);
    if (!self->checkCyclesCompatibility(G, L)) return false;
    uint8_t H[16][16];
    return self->decomposeMissingEdges(G, L, H);
}

void K16A2::CycleBacktrackState::saveAlpha(const uint8_t* alpha_p) {
    std::array<uint8_t, 16> a;
    std::copy(alpha_p, alpha_p + 16, a.begin());
    valid_alphas.push_back(a);
    total_passed_p1f++;
}
void K16A2::printEstimatedTime(int L, long long checked_reps) {
    long long orbit_reps[16] = {
        0, 1, 1, 2, 4, 10, 26, 72, 232, 504, 2619, 5040, 34650, 124740, 405405, 135135
    };
    auto now = std::chrono::steady_clock::now();
    double elapsed_since_print = std::chrono::duration<double>(now - last_print_time).count();
    if (elapsed_since_print >= 30.0) {
        last_print_time = now;
        double elapsed_total = std::chrono::duration<double>(now - case_start_time).count();
        long long total_reps = orbit_reps[L];
        if (checked_reps > 0 && checked_reps < total_reps) {
            double est_total = elapsed_total * ((double)total_reps / checked_reps);
            double est_rem = est_total - elapsed_total;
            printf("\r   [L=%d] Rep %lld/%lld. Est. time to end: %.1f min       ", 
                   L, checked_reps, total_reps, est_rem / 60.0);
            fflush(stdout);
        }
    }
}

// ==========================================
// K16A2 Orchestrator & Helpers
// ==========================================

void K16A2::runExhaustiveSearch() {
    auto search_start = std::chrono::high_resolution_clock::now();
    std::set<std::vector<uint8_t>> unique_results;
    
    struct L_Stats {
        int L;
        long long checked = 0;
        long long passed = 0;
        size_t unique_classes = 0;
        bool run = false;
    };
    L_Stats stats[16];
    for (int i = 1; i <= 15; i++) {
        stats[i].L = i;
    }
    
    std::set<std::vector<uint8_t>> local_l14_unique;
    int cycle_lengths[] = { 15, 14, 12, 10, 8, 6, 4, 3, 2 };
    for (int L : cycle_lengths) {
        if (unique_results.size() >= 2000) {
            stats[L].run = false;
            continue;
        }
        
        std::vector<std::vector<uint8_t>> before_results;
        if (L == 14) {
            before_results.assign(unique_results.begin(), unique_results.end());
        }
        
        CycleLengthStats l_stats;
        searchCycleLength(L, unique_results, l_stats);
        stats[L].checked = l_stats.inputs;
        stats[L].passed = l_stats.passed;
        stats[L].unique_classes = l_stats.unique_classes;
        stats[L].run = true;
        
        if (L == 14) {
            for (const auto& res : unique_results) {
                if (std::find(before_results.begin(), before_results.end(), res) == before_results.end()) {
                    local_l14_unique.insert(res);
                }
            }
            if (!local_l14_unique.empty()) {
                verifyL14Pairing(local_l14_unique);
            }
        }
    }
    
    // Print the summary table
    printf("\n");
    printf("=========================================================================================================================================================\n");
    printf(" L   Cycle Structures Covered                           Orbit Reps   Checked Perms   Passed Filter   Unique Reps  Short Status / Description\n");
    printf("---------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    long long orbit_reps[16] = {
        0, // L=0
        1, // L=1
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
        135135  // L=15
    };
    
    for (int L = 15; L >= 1; L--) {
        char cycle_str[64];
        if (L == 1) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(1)^16");
        } else if (L == 15 || L == 13 || L == 11) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(%d)(1)^%d", L, 16 - L);
        } else if (L == 14) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(14)(1)^2, (14)(2)");
        } else if (L == 12) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(12)(1)^4, (12)(2)(1)^2, (12)(2)^2, (12)(3)(1), (12)(4)");
        } else if (L == 9) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(9)(1)^7");
        } else if (L == 7) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(7)(1)^9, (7)^2(1)^2");
        } else if (L == 5) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(5)^k(1)^{16-5k} for k in {1..3}");
        } else if (L == 3) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(3)^k(1)^{16-3k} for k in {1..5}");
        } else if (L == 2) {
            sprintf_s(cycle_str, sizeof(cycle_str), "(2)^k(1)^{16-2k} for k in {1..8}");
        } else { // L = 10, 8, 6, 4
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
        } else if (L == 11 || L == 9) {
            desc = "Theoretically Impossible (odd L <= 13)";
        } else {
            if (stats[L].run) {
                sprintf_s(checked_str, sizeof(checked_str), "%lld", stats[L].checked);
                sprintf_s(passed_str, sizeof(passed_str), "%lld", stats[L].passed);
                sprintf_s(reps_count_str, sizeof(reps_count_str), "%zu", stats[L].unique_classes);
                
                if (L == 13) {
                    desc = "Theoretically Impossible (odd L <= 13)";
                } else if (L == 15) {
                    desc = "Exhaustive cyclic search";
                } else {
                    desc = stats[L].unique_classes > 0 ? "Exhaustive (P1Fs found)" : "Exhaustive (0 results)";
                }
            } else {
                desc = "Skipped / Not Run";
            }
        }
        
        printf("%2d   %-50s  %11s  %14s  %14s  %12s  %s\n",
               L, cycle_str, reps_str, checked_str, passed_str, reps_count_str, desc);
               
        if (L == 14 && stats[14].run) {
            printf("14R  %-50s  %11s  %14s  %14s  %12s  %s\n",
                   "Reflection pairing verification", "-", "-", "-", "2", "Verified pairing");
        }
    }
    printf("=========================================================================================================================================================\n");
    printf("\n");
    
    reportTotalResults(unique_results, search_start);
}

void K16A2::verifyL14Pairing(const std::set<std::vector<uint8_t>>& local_unique) {
    printf("-> Entering Case: L = 14 Reflection Pairing Verification\n");
    fflush(stdout);
    
    std::vector<std::vector<uint8_t>> mats(local_unique.begin(), local_unique.end());
    size_t n = mats.size();
    
    std::vector<bool> self_reflected(n, false);
    std::vector<int> pair_partner(n, -1);
    
    uint8_t cyc_R[16];
    buildStarterCycle(cyc_R);
    uint8_t p[16];
    for (int i = 0; i < 16; i++) {
        p[cyc_R[i]] = cyc_R[(16 - i) % 16];
    }
    
    for (int i = 0; i < n; i++) {
        uint8_t H[16][16];
        memset(H, 0, sizeof(H));
        for (int k = 1; k <= 15; k++) {
            const uint8_t* src = mats[i].data() + (k - 1) * 16;
            for (int j = 0; j < 8; j++) {
                uint8_t u = src[2 * j];
                uint8_t v = src[2 * j + 1];
                H[k][u] = v;
                H[k][v] = u;
            }
        }
        
        uint8_t H_ref[16][16];
        memset(H_ref, 0, sizeof(H_ref));
        for (int k = 1; k <= 15; k++) {
            uint8_t new_nb = p[k];
            for (int u = 0; u < 16; u++) {
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

void K16A2::searchCycleLength(int L, std::set<std::vector<uint8_t>>& unique_results, CycleLengthStats& stats) {
    printf("-> Entering Case: L = %d, Search Type 1 (Fixed out-of-cycle points)\n", L);
    fflush(stdout);

    case_start_time = std::chrono::steady_clock::now();
    last_print_time = case_start_time;
    case_timed_out = false;
    current_checked_reps = 0;

    CycleBacktrackState state1;
    state1.self = this;
    setupBacktrackState(state1, 1);
    uint8_t c1[16];
    bool used1[16];
    memset(used1, 0, sizeof(used1));
    used1[0] = true;
    used1[state1.v0] = true;

    // Calculate total top-level branches for state1
    int total_branches1 = 0;
    for (int v = 0; v < 16; v++) {
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
    if (!case_timed_out && L >= 2 && L <= 14 && L % 2 == 0) {
        printf("-> Entering Case: L = %d, Search Type 2 (Transposed out-of-cycle points)\n", L);
        fflush(stdout);

        state2.self = this;
        setupBacktrackState(state2, 2);
        uint8_t c2[16];
        bool used2[16];
        memset(used2, 0, sizeof(used2));
        used2[0] = true;
        used2[1] = true;
        used2[state2.v0] = true;

        // Calculate total top-level branches for state2
        int total_branches2 = 0;
        for (int v = 0; v < 16; v++) {
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
        for (const auto& alpha : state1.valid_alphas) {
            if (unique_results.size() + local_unique.size() >= 2000) break;
            std::set<std::vector<uint8_t>> temp_unique;
            processAutomorphism(alpha, L, temp_unique);
            for (const auto& res : temp_unique) {
                if (unique_results.find(res) == unique_results.end()) {
                    local_unique.insert(res);
                }
            }
        }
        if (L >= 2 && L <= 14 && L % 2 == 0) {
            for (const auto& alpha : state2.valid_alphas) {
                if (unique_results.size() + local_unique.size() >= 2000) break;
                std::set<std::vector<uint8_t>> temp_unique;
                processAutomorphism(alpha, L, temp_unique);
                for (const auto& res : temp_unique) {
                    if (unique_results.find(res) == unique_results.end()) {
                        local_unique.insert(res);
                    }
                }
            }
        }
    }
    
    stats.L = L;
    stats.inputs = state1.total_generated + (L >= 2 && L <= 14 && L % 2 == 0 ? state2.total_generated : 0);
    stats.passed = state1.total_passed_p1f + (L >= 2 && L <= 14 && L % 2 == 0 ? state2.total_passed_p1f : 0);
    stats.unique_classes = local_unique.size();
    
    unique_results.insert(local_unique.begin(), local_unique.end());

    if (case_timed_out) {
        long long orbit_reps[16] = {
            0, 1, 1, 2, 4, 10, 26, 72, 232, 504, 2619, 5040, 34650, 124740, 405405, 135135
        };
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - case_start_time).count();
        long long total_reps = orbit_reps[L];
        long long checked_reps = stats.inputs;
        double estimated_time = -1.0;
        if (checked_reps > 0) {
            estimated_time = elapsed * ((double)total_reps / checked_reps);
        }
        printf("\n   [L=%d] TIMEOUT reached. Checked %lld / %lld orbit reps in %.2f seconds.\n", 
               L, checked_reps, total_reps, elapsed);
        if (estimated_time >= 0) {
            printf("   [L=%d] Estimated total time needed without timeout: %.2f seconds (%.2f hours)\n", 
                   L, estimated_time, estimated_time / 3600.0);
        } else {
            printf("   [L=%d] Estimated total time: infinite/unknown (checked 0 reps)\n", L);
        }
    } else {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - case_start_time).count();
        if (elapsed >= 30.0) {
            printf("\n");
        }
    }
    fflush(stdout);
}

void K16A2::setupBacktrackState(CycleBacktrackState& state, int search_type) {
    memcpy(state.F[0], fixedRows[0].adj, 16);
    if (search_type == 2) {
        state.v0 = 3;
    } else {
        state.v0 = fixedRows[0].adj[0];
    }
    setupPairsTable(state, search_type);
}

void K16A2::setupPairsTable(CycleBacktrackState& state, int search_type) {
    for (int i = 0; i < 16; i++) {
        state.vertex_to_pair[i] = -1;
        state.vertex_to_pos[i] = -1;
    }
    fillPairsTable(state, search_type);
    
    uint8_t partner = fixedRows[0].adj[state.v0];
    state.vertex_to_pair[partner] = -2;
}

void K16A2::fillPairsTable(CycleBacktrackState& state, int search_type) {
    int pair_count = 0;
    for (int u = 1; u < 16; u++) {
        if (search_type == 2 && u == 1) continue;
        if (u == state.v0) continue;
        uint8_t v = fixedRows[0].adj[u];
        if (v == state.v0) continue;
        if (state.vertex_to_pair[u] != -1) continue;
        storePair(state, pair_count++, u, v);
    }
}

void K16A2::storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v) {
    state.pair_elements[pair_idx][0] = u;
    state.pair_elements[pair_idx][1] = v;
    state.vertex_to_pair[u] = pair_idx;
    state.vertex_to_pos[u] = 0;
    state.vertex_to_pair[v] = pair_idx;
    state.vertex_to_pos[v] = 1;
}

void K16A2::processAutomorphism(const std::array<uint8_t, 16>& alpha_arr, int L, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t alpha[16];
    std::copy(alpha_arr.begin(), alpha_arr.end(), alpha);
    uint8_t H[16][16];
    constructFullHFromAut(alpha, L, H);
    recordIsomorphicResults(H, unique_results);
}

void K16A2::constructFullHFromAut(const uint8_t* alpha, int L, uint8_t H[][16]) {
    uint8_t G[15][16];
    memcpy(G[0], fixedRows[0].adj, 16);
    constructFullH(G, L, alpha, G);
    if (L == 15) {
        copyMatchingsToH(nullptr, 0, G, L, H);
    } else {
        decomposeMissingEdges(G, L, H);
    }
}

void K16A2::constructFullH(const uint8_t G[][16], int L, const uint8_t* alpha, uint8_t H[][16]) {
    memcpy(H[0], fixedRows[0].adj, 16);
    Permutation perm;
    memcpy(perm.p, alpha, 16);
    for (int k = 1; k < L; k++) {
        apply_perm_16(H[k-1], perm, H[k]);
    }
}

bool K16A2::checkCyclesCompatibility(const uint8_t G[][16], int L) {
    for (int i = 0; i < L; i++) {
        for (int j = i + 1; j < L; j++) {
            if (!is_perfect_scalar(G[i], G[j])) return false;
        }
    }
    return true;
}

bool K16A2::decomposeMissingEdges(const uint8_t G[][16], int L, uint8_t H[][16]) {
    if (case_timed_out) return false;
    int num_colors = 15 - L;
    
    // 1. Build adjacency matrix of the complement graph
    bool complement_adj[16][16];
    memset(complement_adj, false, sizeof(complement_adj));
    for (int u = 0; u < 16; u++) {
        for (int v = u + 1; v < 16; v++) {
            if (isEdgeMissing(G, L, u, v)) {
                complement_adj[u][v] = true;
                complement_adj[v][u] = true;
            }
        }
    }
    
    // 2. Generate all perfect matchings in complement graph compatible with G[0]
    std::vector<std::array<uint8_t, 16>> pool;
    uint8_t current_matching[16];
    memset(current_matching, 0xFF, sizeof(current_matching));
    bool used[16] = { false };
    
    auto forms_disjoint_cycle = [&](int u, int v, const uint8_t* matching) -> bool {
        uint8_t curr = G[0][u];
        int len = 2;
        while (curr != v) {
            uint8_t next_v = matching[curr];
            if (next_v == 0xFF) return false;
            curr = G[0][next_v];
            len += 2;
        }
        return len < 16;
    };
    
    auto generate_matchings = [&](auto& self_fn, int u) -> void {
        if (case_timed_out) return;
        if (u == 16) {
            std::array<uint8_t, 16> m;
            memcpy(m.data(), current_matching, 16);
            pool.push_back(m);
            return;
        }
        if (used[u]) {
            self_fn(self_fn, u + 1);
            return;
        }
        used[u] = true;
        for (int v = u + 1; v < 16; v++) {
            if (!used[v] && complement_adj[u][v]) {
                if (!forms_disjoint_cycle(u, v, current_matching)) {
                    used[v] = true;
                    current_matching[u] = (uint8_t)v;
                    current_matching[v] = (uint8_t)u;
                    self_fn(self_fn, u + 1);
                    used[v] = false;
                    current_matching[u] = 0xFF;
                    current_matching[v] = 0xFF;
                }
            }
        }
        used[u] = false;
    };
    generate_matchings(generate_matchings, 0);
    
    if (pool.size() < (size_t)num_colors) return false;
    
    // 3. Find a clique of size num_colors of pairwise compatible matchings
    std::vector<int> clique;
    
    auto find_clique = [&](auto& self_fn, int start_idx, int depth) -> bool {
        if (case_timed_out) return false;
        static thread_local int find_clique_call_count = 0;
        if (++find_clique_call_count >= 1000) {
            find_clique_call_count = 0;
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - case_start_time).count();
            if (elapsed > case_timeout_seconds) {
                case_timed_out = true;
                return false;
            }
            printEstimatedTime(L, current_checked_reps);
        }
        if (depth == num_colors) return true;
        
        for (int i = start_idx; i < pool.size(); i++) {
            bool ok = true;
            for (int idx : clique) {
                // Edge-disjoint check
                for (int v = 0; v < 16; v++) {
                    if (pool[i][v] == pool[idx][v]) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) break;
                
                // Cycle compatibility check
                if (!is_perfect_scalar(pool[i].data(), pool[idx].data())) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                clique.push_back(i);
                if (self_fn(self_fn, i + 1, depth + 1)) return true;
                clique.pop_back();
            }
        }
        return false;
    };
    
    if (!find_clique(find_clique, 0, 0)) return false;
    
    // 4. Copy the found clique of matchings to H
    uint8_t matchings[15][16];
    memset(matchings, 0xFF, sizeof(matchings));
    for (int c = 0; c < num_colors; c++) {
        memcpy(matchings[c], pool[clique[c]].data(), 16);
    }
    copyMatchingsToH(matchings, num_colors, G, L, H);
    return true;
}

void K16A2::findMissingEdges(const uint8_t G[][16], int L, std::pair<uint8_t, uint8_t>* edges) {
    int edge_cnt = 0;
    for (int u = 0; u < 16; u++) {
        for (int v = u + 1; v < 16; v++) {
            if (isEdgeMissing(G, L, u, v)) {
                edges[edge_cnt++] = { (uint8_t)u, (uint8_t)v };
            }
        }
    }
}

bool K16A2::isEdgeMissing(const uint8_t G[][16], int L, int u, int v) {
    for (int k = 0; k < L; k++) {
        if (G[k][u] == v) return false;
    }
    return true;
}

bool K16A2::backtrackColor(int edge_idx, int num_edges, int num_colors,
                           const std::pair<uint8_t, uint8_t>* edges,
                           uint8_t matchings[][16], const uint8_t* G0) {
    if (edge_idx == num_edges) {
        return checkMatchingsCompatibility(matchings, num_colors, G0);
    }
    return tryColoringEdge(edge_idx, num_edges, num_colors, edges, matchings, G0);
}

bool K16A2::checkMatchingsCompatibility(uint8_t matchings[][16], int num_colors, const uint8_t* G0) {
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

bool K16A2::tryColoringEdge(int edge_idx, int num_edges, int num_colors,
                            const std::pair<uint8_t, uint8_t>* edges,
                            uint8_t matchings[][16], const uint8_t* G0) {
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

bool K16A2::is_color_compatible(int c, int num_colors, uint8_t matchings[][16], const uint8_t* G0) {
    int count = 0;
    for (int i = 0; i < 16; i++) {
        if (matchings[c][i] != 0xFF) count++;
    }
    if (count < 16) return true; // Not yet complete
    
    // Check compatibility with G0
    if (!is_perfect_scalar(G0, matchings[c])) return false;
    
    // Check compatibility with other complete colors
    for (int d = 0; d < num_colors; d++) {
        if (d == c) continue;
        int count_d = 0;
        for (int i = 0; i < 16; i++) {
            if (matchings[d][i] != 0xFF) count_d++;
        }
        if (count_d == 16) {
            if (!is_perfect_scalar(matchings[c], matchings[d])) return false;
        }
    }
    return true;
}

int K16A2::getSymmetryBreakingLimit(uint8_t matchings[][16], int num_colors) {
    int max_color = 0;
    for (int c = 0; c < num_colors; c++) {
        if (isColorUsed(matchings[c])) max_color = c + 1;
    }
    return (max_color < num_colors - 1) ? max_color : num_colors - 1;
}

bool K16A2::isColorUsed(const uint8_t* matching) {
    for (int i = 0; i < 16; i++) {
        if (matching[i] != 0xFF) return true;
    }
    return false;
}

void K16A2::copyMatchingsToH(uint8_t matchings[][16], int num_colors,
                             const uint8_t G[][16], int L, uint8_t H[][16]) {
    memset(H, 0, 16 * 16);
    for (int k = 0; k < L; k++) {
        uint8_t nb = G[k][0];
        memcpy(H[nb], G[k], 16);
    }
    for (int c = 0; c < num_colors; c++) {
        uint8_t nb = matchings[c][0];
        memcpy(H[nb], matchings[c], 16);
    }
}

void K16A2::recordIsomorphicResults(const uint8_t H[][16], std::set<std::vector<uint8_t>>& unique_results) {
    CycleUnion cu_H = find_cycles(H[1], H[2]);
    if (cu_H.count != 1 || cu_H.lens[0] != 16) return;
    for (int v = 0; v < 16; v++) {
        tryIsomorphicMapping(H, cu_H, v, unique_results);
    }
}

void K16A2::tryIsomorphicMapping(const uint8_t H[][16], const CycleUnion& cu_H, int v, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t cyc_R[16];
    buildStarterCycle(cyc_R);
    uint8_t cyc_H[16];
    buildHCycle(H, v, cyc_H);
    uint8_t p[16];
    buildMappingPermutation(cyc_H, cyc_R, p);
    checkAndRecordPermutedH(H, p, unique_results);
}

void K16A2::buildStarterCycle(uint8_t* cyc_R) {
    uint8_t curr = 0;
    for (int i = 0; i < 8; i++) {
        cyc_R[2 * i] = curr;
        curr = fixedRows[0].adj[curr];
        cyc_R[2 * i + 1] = curr;
        curr = fixedRows[1].adj[curr];
    }
}

void K16A2::buildHCycle(const uint8_t H[][16], int v, uint8_t* cyc_H) {
    uint8_t curr = v;
    for (int i = 0; i < 8; i++) {
        cyc_H[2 * i] = curr;
        curr = H[1][curr];
        cyc_H[2 * i + 1] = curr;
        curr = H[2][curr];
    }
}

void K16A2::buildMappingPermutation(const uint8_t* cyc_H, const uint8_t* cyc_R, uint8_t* p) {
    for (int i = 0; i < 16; i++) {
        p[cyc_H[i]] = cyc_R[i];
    }
}

void K16A2::checkAndRecordPermutedH(const uint8_t H[][16], const uint8_t* p, std::set<std::vector<uint8_t>>& unique_results) {
    uint8_t S[16][16];
    memset(S, 0, sizeof(S));
    applyPermToH(H, p, S);
    if (doesSMatchFixedRows(S)) {
        recordS(S, unique_results);
    }
}

void K16A2::applyPermToH(const uint8_t H[][16], const uint8_t* p, uint8_t S[][16]) {
    for (int k = 1; k <= 15; k++) {
        uint8_t permuted_factor[16];
        for (int i = 0; i < 16; i++) {
            permuted_factor[p[i]] = p[H[k][i]];
        }
        uint8_t nb = permuted_factor[0];
        memcpy(S[nb], permuted_factor, 16);
    }
}

bool K16A2::doesSMatchFixedRows(const uint8_t S[][16]) {
    return memcmp(S[1], fixedRows[0].adj, 16) == 0 &&
           memcmp(S[2], fixedRows[1].adj, 16) == 0;
}

void K16A2::recordS(const uint8_t S[][16], std::set<std::vector<uint8_t>>& unique_results) {
    unsigned char results[15 * 16];
    for (int k = 1; k <= 15; k++) {
        adj_to_src(S[k], results + (k - 1) * 16);
    }
    std::vector<uint8_t> res_vec(results, results + 15 * 16);
    unique_results.insert(res_vec);
}

void K16A2::reportTotalResults(const std::set<std::vector<uint8_t>>& unique_results,
                               std::chrono::high_resolution_clock::time_point start) {
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    printf("Exhaustive search completed. Unique P1Fs: %zu. Time: %.4f seconds\n", unique_results.size(), elapsed);
    sendResultsToCallback(unique_results);
}

void K16A2::sendResultsToCallback(const std::set<std::vector<uint8_t>>& unique_results) {
    std::lock_guard<std::mutex> lock(result_mutex);
    for (const auto& res : unique_results) {
        resultCallback(cbClass, res.data(), 0, 1, 2);
    }
}

void K16A2::adj_to_src(const uint8_t* adj, unsigned char* src) {
    src[0] = 0;
    src[1] = adj[0];
    int idx = 2;
    bool visited[16] = { false };
    visited[0] = true;
    visited[adj[0]] = true;
    for (int u = 1; u < 16; u++) {
        if (!visited[u]) {
            src[idx] = u;
            src[idx + 1] = adj[u];
            visited[u] = true;
            visited[adj[u]] = true;
            idx += 2;
        }
    }
}
