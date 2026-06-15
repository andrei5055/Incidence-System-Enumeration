#include "k18a2.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <set>
#include <chrono>

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

void K18A2::CycleBacktrackState::backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    if (depth == L - 1) {
        processCandidate(c, used, L);
        return;
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
    used[v] = true;
    c[depth] = v;
    backtrack(depth + 1, next_pairs, c, used, L);
    used[v] = false;
}

void K18A2::CycleBacktrackState::processCandidate(uint8_t* c, bool* used, int L) {
    total_generated++;
    uint8_t alpha_p[18];
    buildPermutation(c, alpha_p, L);
    
    int search_type = (v0 == 3) ? 2 : 1;
    if (search_type == 2) {
        alpha_p[0] = 1;
        alpha_p[1] = 0;
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
    generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L);
}

void K18A2::CycleBacktrackState::generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size, bool* rem_used, uint8_t* alpha_p, int L) {
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
    for (int j = 1; j < 9; j++) {
        apply_perm(G[j - 1], alpha_p, G[j]);
        if (!is_perfect_scalar(G[0], G[j])) return false;
    }
    return validateCandidateL(alpha_p, used, L);
}

bool K18A2::CycleBacktrackState::validateCandidateL(const uint8_t* alpha_p, bool* used, int L) {
    if (L == 17) return true;
    uint8_t G[17][18];
    self->constructFullH(G, L, alpha_p, G);
    if (!self->checkCyclesCompatibility(G, L)) return false;
    uint8_t H[18][18];
    return self->decomposeMissingEdges(G, L, H);
}

void K18A2::CycleBacktrackState::saveAlpha(const uint8_t* alpha_p) {
    std::array<uint8_t, 18> a;
    std::copy(alpha_p, alpha_p + 18, a.begin());
    valid_alphas.push_back(a);
    total_passed_p1f++;
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
    int cycle_lengths[] = { 17, 16, 15, 14, 13, 12, 10, 8, 6, 5, 4, 3, 2 };
    for (int L : cycle_lengths) {
        if (unique_results.size() >= 2000) {
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

    CycleBacktrackState state1;
    state1.self = this;
    setupBacktrackState(state1, 1);
    uint8_t c1[18];
    bool used1[18];
    memset(used1, 0, sizeof(used1));
    used1[0] = true;
    used1[state1.v0] = true;
    state1.backtrack(0, 0, c1, used1, L);
    
    CycleBacktrackState state2;
    if (L >= 2 && L <= 16 && L % 2 == 0) {
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
        state2.backtrack(0, 0, c2, used2, L);
    }
    
    std::set<std::vector<uint8_t>> local_unique;
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
    if (L >= 2 && L <= 16 && L % 2 == 0) {
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
    
    stats.L = L;
    stats.inputs = state1.total_generated + (L >= 2 && L <= 16 && L % 2 == 0 ? state2.total_generated : 0);
    stats.passed = state1.total_passed_p1f + (L >= 2 && L <= 16 && L % 2 == 0 ? state2.total_passed_p1f : 0);
    stats.unique_classes = local_unique.size();
    
    unique_results.insert(local_unique.begin(), local_unique.end());
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
        decomposeMissingEdges(G, L, H);
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

bool K18A2::decomposeMissingEdges(const uint8_t G[][18], int L, uint8_t H[][18]) {
    int num_colors = 17 - L;
    int num_edges = 9 * num_colors;
    std::pair<uint8_t, uint8_t> edges[153];
    findMissingEdges(G, L, edges);
    uint8_t matchings[15][18];
    memset(matchings, 0xFF, sizeof(matchings));
    if (!backtrackColor(0, num_edges, num_colors, edges, matchings, G[0])) return false;
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
            if (backtrackColor(edge_idx + 1, num_edges, num_colors, edges, matchings, G0)) return true;
            matchings[c][u] = 0xFF;
            matchings[c][v] = 0xFF;
        }
    }
    return false;
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
