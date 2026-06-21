#include "k18a2.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <set>
#include <chrono>
#include <vector>
#include <array>

// ============================================================================
// Lazy candidate generation methods
// ============================================================================

void K18A2::ensureSlotGenerated(int slotIndex) {
    if (slot_generated[slotIndex]) return;
    std::lock_guard<std::mutex> lock(slot_mutex[slotIndex]);
    if (slot_generated[slotIndex]) return;

    generateSlotCandidates(slotIndex);
    slot_generated[slotIndex] = true;
}

void K18A2::generateSlotCandidates(int slotIndex) {
    int P0 = slotIndex + 4;
    uint8_t factor[18];
    factor[0] = 0; factor[1] = (uint8_t)P0;
    uint64_t usedMask = (1ULL << 0) | (1ULL << P0);
    generateRecursive(0, usedMask, factor, slotIndex, P0);
}

void K18A2::generateRecursive(int depth, uint64_t usedMask, uint8_t* factor, int slotIndex, int P0) {
    if (depth == K18_N / 2 - 1) {
        addRow(P0, factor);
        return;
    }

    unsigned long u;
    if (!_BitScanForward64(&u, ~usedMask)) return;

    for (int v = 0; v < K18_N; ++v) {
        if (usedMask & (1ULL << v)) continue;
        if (v <= (int)u) continue;

        int eid = edge_id_table[u][v];
        if (eid != -1 && fixedEdgesMask.m[eid >> 6] & (1ULL << (eid & 63))) continue;

        factor[(depth + 1) * 2] = (uint8_t)u;
        factor[(depth + 1) * 2 + 1] = (uint8_t)v;
        generateRecursive(depth + 1, usedMask | (1ULL << u) | (1ULL << v), factor, slotIndex, P0);
    }
}

// ============================================================================
// Transition-based Automorphism Search (K16A2Old algorithm scaled for K18)
// ============================================================================

static bool has_prime_order_k18(const uint8_t* p, int target_p) {
    bool visited[18] = { false };
    bool found_target_cycle = false;
    for (int i = 0; i < 18; i++) {
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

void K18A2::runTransitionSearch() {
    // 1. Generate all candidate pools internally
    for (int s = 0; s < K18_SEARCH; s++) {
        ensureSlotGenerated(s);
    }
    
    int Global_M_total = (int)global_pool.size();
    if (m_bPrint) {
        printf("K18A2 solver starting (Transition-based search with internal candidate generation)...\n");
        printf("Row slots size: ");
        for (int s = 0; s < K18_SEARCH; s++) printf("%d:%zd ", s + 4, temp_slot_ids[s].size());
        printf("\n");
        printf("Total candidates: %d\n", Global_M_total);
    }
    
    std::vector<uint8_t> active(Global_M_total, 1);
    std::vector<uint32_t> f_rows_mask(Global_M_total, 0);
    for (int s = 0; s < K18_SEARCH; s++) {
        for (int id : temp_slot_ids[s]) {
            f_rows_mask[id] |= (1 << s);
        }
    }
    
    // Arc consistency check / Filter
    bool changed = true;
    int iter = 0;
    while (changed) {
        changed = false;
        for (int i = 0; i < Global_M_total; i++) {
            if (!active[i]) continue;
            if (i < 3) continue;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < K18_SEARCH; r++) {
                if (f_rows_mask[i] & (1 << r)) continue; 
                bool row_found = false;
                for (int id : temp_slot_ids[r]) {
                    if (active[id] && is_perfect_packed(pi, packed_pool[id])) { 
                        row_found = true; 
            break;
        }
    }
                if (!row_found) {
                    active[i] = 0;
                    changed = true;
                    break;
        }
            }
        }
        int rem = 0; for (int i = 0; i < Global_M_total; i++) if (active[i]) rem++;
        iter++;
        if (m_bPrint) { printf(" Remaining(%d):%d ", iter, rem); }
        if (rem == 0) {
            printf("\nAll candidates filtered by Arc consistency check\n");
        return;
    }
    }
    if (m_bPrint) {
        printf("\n");
        printf("DEBUG: Rebuilding active slots...\n");
    }
    
    // Rebuild active slots map for O(1) lookup
    std::vector<int> active_slots(Global_M_total, -1);
    std::vector<std::vector<int>> new_slot_ids(K18_SEARCH);
    for (int s = 0; s < K18_SEARCH; s++) {
        for (int id : temp_slot_ids[s]) {
            if (active[id]) {
                new_slot_ids[s].push_back(id);
                active_slots[id] = s;
            }
                }
            }
            
    std::mutex results_mutex;
    std::vector<int> target_primes = {17, 13, 11, 7, 5, 3, 2};
    int pairs_processed = 0;
    uint64_t total_pairs_to_process = 6;
    auto search_start = std::chrono::steady_clock::now();
                
    auto run_search_for_pair = [&](int type1, int val1, int type2, int val2) {
        Factor T1 = global_pool[val1];
        Factor T2 = global_pool[val2];
                
        if (m_bPrint) {
            printf("\nDEBUG: Starting search for pair transition %d <-> %d...\n", val1, val2);
            fflush(stdout);
                        }
                        
        uint64_t total_clique_calls = 0;
        uint64_t total_lookahead_slots = 0;
        uint64_t total_lookahead_orbits = 0;
        uint64_t total_lookahead_failures = 0;

        auto process_permutation = [&](const Permutation& alpha) {
            for (int p : target_primes) {
                if (!has_prime_order_k18(alpha.p, p)) continue;
                        
                if (m_bPrint) {
                    printf("  DEBUG: Permutation order %d: building orbits...\n", p);
                    fflush(stdout);
                        }

                // Build orbits under alpha for all factors
                struct Orbit {
                    std::vector<int> fids;
                    std::vector<int> slots;
                    uint16_t slots_mask = 0;
                    Mask18_C edge_mask;
                };

                std::vector<Orbit> fixed_orbits;
                std::vector<Orbit> valid_orbits;
                std::vector<uint8_t> visited_factor(Global_M_total, 0);
                bool alpha_valid = true;

                for (int fid = 0; fid < Global_M_total; fid++) {
                    if (fid >= 3 && !active[fid]) continue;
                    if (visited_factor[fid]) continue;

                    Orbit orb;
                    int curr_fid = fid;
                    bool orb_ok = true;

                    while (true) {
                        visited_factor[curr_fid] = 1;
                        orb.fids.push_back(curr_fid);
                        orb.slots.push_back(active_slots[curr_fid]);

                        uint8_t tr_f[18];
                        apply_perm_18(global_pool[curr_fid].adj, alpha, tr_f);

                        std::array<uint8_t, 18> key_arr;
                        std::copy(tr_f, tr_f + 18, key_arr.begin());
                        auto it_f = f_map_unordered.find(key_arr);
                        if (it_f == f_map_unordered.end()) {
                            orb_ok = false;
                            break;
                    }
                        int next_fid = it_f->second;
                        if (next_fid >= 3 && !active[next_fid]) {
                            orb_ok = false;
                            break;
                        }
                        if (next_fid == fid) break;
                        if (visited_factor[next_fid]) {
                            orb_ok = false;
                            break;
                    }
                        curr_fid = next_fid;
            }

                    // Check pairwise compatibility of factors within the orbit using automorphism symmetry
                    bool cycle_ok = true;
                    int f1 = orb.fids[0];
                    for (size_t j = 1; j < orb.fids.size(); j++) {
                        int f2 = orb.fids[j];
                        if (!is_perfect_packed(packed_pool[f1], packed_pool[f2])) {
                            cycle_ok = false;
                            break;
        }
    }
    
                    // Check slot conflicts and build slots_mask
                    orb.slots_mask = 0;
                    bool slot_conflict = false;
                    for (int sl : orb.slots) {
                        if (sl >= 0) {
                            if (orb.slots_mask & (1 << sl)) {
                                slot_conflict = true;
                                break;
    }
                            orb.slots_mask |= (1 << sl);
    }
    }
    
                    bool contains_fixed = false;
                    for (int f_id : orb.fids) {
                        if (f_id < 3) contains_fixed = true;
        }
        
                    if (contains_fixed) {
                        if (!orb_ok || (orb.fids.size() != 1 && orb.fids.size() != p) || !cycle_ok || slot_conflict) {
                            alpha_valid = false;
                            break;
        }
        
                        // Check internal edge conflict
                        orb.edge_mask = Mask18_C();
                        bool conflict = false;
                        for (int f_id : orb.fids) {
                            if ((orb.edge_mask.m[0] & packed_pool[f_id].edge_mask.m[0]) ||
                                (orb.edge_mask.m[1] & packed_pool[f_id].edge_mask.m[1]) ||
                                (orb.edge_mask.m[2] & packed_pool[f_id].edge_mask.m[2])) {
                                conflict = true;
                                break;
            }
                            orb.edge_mask.m[0] |= packed_pool[f_id].edge_mask.m[0];
                            orb.edge_mask.m[1] |= packed_pool[f_id].edge_mask.m[1];
                            orb.edge_mask.m[2] |= packed_pool[f_id].edge_mask.m[2];
        }
                        if (conflict) {
                            alpha_valid = false;
                            break;
    }
                        fixed_orbits.push_back(orb);
                } else {
                        if (orb_ok && (orb.fids.size() == 1 || orb.fids.size() == p) && cycle_ok && !slot_conflict) {
                            // Check internal edge conflict
                            orb.edge_mask = Mask18_C();
                            bool conflict = false;
                            for (int f_id : orb.fids) {
                                if ((orb.edge_mask.m[0] & packed_pool[f_id].edge_mask.m[0]) ||
                                    (orb.edge_mask.m[1] & packed_pool[f_id].edge_mask.m[1]) ||
                                    (orb.edge_mask.m[2] & packed_pool[f_id].edge_mask.m[2])) {
                                    conflict = true;
                                    break;
                }
                                orb.edge_mask.m[0] |= packed_pool[f_id].edge_mask.m[0];
                                orb.edge_mask.m[1] |= packed_pool[f_id].edge_mask.m[1];
                                orb.edge_mask.m[2] |= packed_pool[f_id].edge_mask.m[2];
                            }
                            if (!conflict) {
                                valid_orbits.push_back(orb);
            }
        }
        }
    }
    
                if (!alpha_valid) continue;
    
                // Reset visited markers
                for (int f_id = 0; f_id < Global_M_total; f_id++) visited_factor[f_id] = 0;
    
                // Check compatibility among fixed orbits and initialize search state
                Mask18_C initial_edges;
                uint16_t initial_slots_mask = 0;
                bool fixed_compatible = true;
    
                std::vector<int> fixed_fids;
                for (const auto& orb : fixed_orbits) {
                    // Check edge conflicts
                    if ((initial_edges.m[0] & orb.edge_mask.m[0]) ||
                        (initial_edges.m[1] & orb.edge_mask.m[1]) ||
                        (initial_edges.m[2] & orb.edge_mask.m[2])) {
                        fixed_compatible = false;
                        break;
            }
                    initial_edges.m[0] |= orb.edge_mask.m[0];
                    initial_edges.m[1] |= orb.edge_mask.m[1];
                    initial_edges.m[2] |= orb.edge_mask.m[2];

                    // Check slot conflicts and set mask
                    if (initial_slots_mask & orb.slots_mask) {
                        fixed_compatible = false;
                        break;
        }
                    initial_slots_mask |= orb.slots_mask;
        
                    for (int fid : orb.fids) {
                        fixed_fids.push_back(fid);
            }
        }
        
                if (!fixed_compatible) continue;
        
                // Check pairwise compatibility between all factors in fixed_orbits
                for (size_t i = 0; i < fixed_fids.size(); i++) {
                    int f1 = fixed_fids[i];
                    for (size_t j = i + 1; j < fixed_fids.size(); j++) {
                        int f2 = fixed_fids[j];
                        if (!is_perfect_packed(packed_pool[f1], packed_pool[f2])) {
                            fixed_compatible = false;
                            break;
                        }
        }
                    if (!fixed_compatible) break;
                }

                if (!fixed_compatible) continue;
        
                // Filter valid_orbits to only those compatible with fixed_orbits
                std::vector<Orbit> filtered_valid;
                for (const auto& orb : valid_orbits) {
                    bool compatible = true;
                    // Check edge conflict with fixed orbits
                    if ((initial_edges.m[0] & orb.edge_mask.m[0]) ||
                        (initial_edges.m[1] & orb.edge_mask.m[1]) ||
                        (initial_edges.m[2] & orb.edge_mask.m[2])) {
                        compatible = false;
                    } else {
                        // Check slot overlap with fixed orbits
                        if (initial_slots_mask & orb.slots_mask) {
                            compatible = false;
                        } else {
                            // Check pairwise compatibility with fixed_fids
                            for (int fid_fixed : fixed_fids) {
                                for (int fid_new : orb.fids) {
                                    if (!is_perfect_packed(packed_pool[fid_fixed], packed_pool[fid_new])) {
                                        compatible = false;
                break;
            }
        }
                                if (!compatible) break;
    }
            }
            }
                    if (compatible) {
                        filtered_valid.push_back(orb);
        }
    }
    
                int num_orbits = (int)filtered_valid.size();
                if (num_orbits == 0) continue;

                // Local Arc Consistency on factors belonging to valid orbits
                std::vector<uint8_t> local_factor_active(Global_M_total, 0);
                for (int fid : fixed_fids) {
                    local_factor_active[fid] = 1;
    }
                for (const auto& orb : filtered_valid) {
                    for (int fid : orb.fids) {
                        local_factor_active[fid] = 1;
            }
        }

                // Map each factor to its orbit index in filtered_valid
                std::vector<int> factor_to_orbit(Global_M_total, -1);
                for (int o_idx = 0; o_idx < num_orbits; o_idx++) {
                    for (int fid : filtered_valid[o_idx].fids) {
                        factor_to_orbit[fid] = o_idx;
    }
                }

                bool local_changed = true;
                int local_iter = 0;
                while (local_changed) {
                    local_changed = false;
                    
                    // Build list of active factors per slot for fast lookup
                    std::vector<int> active_factors_in_slot[K18_SEARCH];
                    for (int s = 0; s < K18_SEARCH; s++) {
                        for (int id : new_slot_ids[s]) {
                            if (local_factor_active[id]) {
                                active_factors_in_slot[s].push_back(id);
            }
        }
    }
    
                    for (int i = 3; i < Global_M_total; i++) {
                        if (!local_factor_active[i]) continue;
    
                        const PackedAdj& pi = packed_pool[i];
                        int my_slot = active_slots[i];

                        for (int r = 0; r < K18_SEARCH; r++) {
                            if (r == my_slot) continue;

                            bool row_found = false;
                            for (int id : active_factors_in_slot[r]) {
                                if (is_perfect_packed(pi, packed_pool[id])) {
                                    row_found = true;
                                    break;
    }
    }

                            if (!row_found) {
                                int o_idx = factor_to_orbit[i];
                                if (o_idx != -1) {
                                    for (int fid : filtered_valid[o_idx].fids) {
                                        local_factor_active[fid] = 0;
                                    }
    } else {
                                    local_factor_active[i] = 0;
    }
                                local_changed = true;
                                break;
    }
        }
    }
                    local_iter++;
                    if (m_bPrint && local_iter % 2 == 0) {
                        int rem = 0; 
                        for (int i = 3; i < Global_M_total; i++) if (local_factor_active[i]) rem++;
                        printf("      Local AC iteration %d: remaining active factors = %d\n", local_iter, rem);
                    }
                }

                // Filter valid orbits to only those whose factors are still active
                std::vector<Orbit> local_filtered;
                for (const auto& orb : filtered_valid) {
                    bool all_active = true;
                    for (int fid : orb.fids) {
                        if (!local_factor_active[fid]) {
                            all_active = false;
                            break;
            }
        }
                    if (all_active) {
                        local_filtered.push_back(orb);
    }
    }
                filtered_valid = std::move(local_filtered);
                num_orbits = (int)filtered_valid.size();

                if (m_bPrint) {
                    printf("    DEBUG: Orbits built. fixed=%zd, valid=%zd, filtered=%zd. Starting clique search...\n",
                        fixed_orbits.size(), valid_orbits.size(), filtered_valid.size());
                    fflush(stdout);
    }

                // Pre-index the valid orbits by slot to enable exact cover backtracking search
                std::vector<int> orbits_by_slot[K18_SEARCH];
                for (int o_idx = 0; o_idx < num_orbits; o_idx++) {
                    for (int s : filtered_valid[o_idx].slots) {
                        if (s >= 0) {
                            orbits_by_slot[s].push_back(o_idx);
    }
        }
    }

                // Sort slot indices by number of candidate orbits (static MRV)
                std::vector<int> slot_order(K18_SEARCH);
                for (int i = 0; i < K18_SEARCH; i++) slot_order[i] = i;
                std::sort(slot_order.begin(), slot_order.end(), [&](int a, int b) {
                    return orbits_by_slot[a].size() < orbits_by_slot[b].size();
                });

                // Clique of orbits search using dynamic MRV (Minimum Remaining Values)
                std::vector<int> chosen_orbits;
                auto find_clique = [&](auto& self, uint16_t slots_mask, Mask18_C edges_mask) -> void {
                    uint64_t calls = ++total_clique_calls;
                    int d = (int)chosen_orbits.size();

                    if (m_bPrint && calls % 10000 == 0) {
                        auto current_time = std::chrono::steady_clock::now();
                        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - search_start).count();
                        double elapsed_min = (double)elapsed_ms / 60000.0;
                        printf("    CLIQUE_PROGRESS: depth=%d, slots_mask=0x%03X, calls=%llu. Elapsed: %.2f min.\n",
                            d, slots_mask, calls, elapsed_min);
                        fflush(stdout);
    }

                    if (slots_mask == 0x3FFF) {
                        std::vector<int> solution_fids;
                        for (int fid : fixed_fids) {
                            solution_fids.push_back(fid);
                        }
                        for (int o_idx : chosen_orbits) {
                            for (int fid : filtered_valid[o_idx].fids) {
                                solution_fids.push_back(fid);
    }
    }

                        bool sol_perfect = true;
                        for (size_t i = 0; i < solution_fids.size(); i++) {
                            const uint8_t* adj_i = global_pool[solution_fids[i]].adj;
                            for (size_t j = i + 1; j < solution_fids.size(); j++) {
                                const uint8_t* adj_j = global_pool[solution_fids[j]].adj;
                                if (!is_perfect_scalar(adj_i, adj_j)) {
                                    sol_perfect = false;
                                    break;
    }
    }
                            if (!sol_perfect) break;
    }
                        if (!sol_perfect) return;

                        unsigned char results[K18_MATCH * K18_N];
                        std::vector<Factor> sol_factors;
                        for (int fid : solution_fids) {
                            sol_factors.push_back(global_pool[fid]);
    }
                        std::sort(sol_factors.begin(), sol_factors.end(), [](const Factor& a, const Factor& b) { return a.src[1] < b.src[1]; });
                        for (int i = 0; i < K18_MATCH; i++) memcpy(results + i * K18_N, sol_factors[i].src, K18_N);

                        std::lock_guard<std::mutex> lock(results_mutex);
                        resultCallback(cbClass, results, val1, val2, 2);
                        return;
    }

                    // Dynamic MRV: Find the uncovered slot with the absolute minimum number of compatible candidate orbits
                    int s = -1;
                    int min_compatible = 9999999;
                    std::vector<int> compatible_orbits_for_s;

                    for (int i = 0; i < K18_SEARCH; i++) {
                        int slot = slot_order[i];
                        if (slots_mask & (1 << slot)) continue;

                        int compatible_count = 0;
                        std::vector<int> temp_compatible;

                        for (int o_idx : orbits_by_slot[slot]) {
                            const auto& orb = filtered_valid[o_idx];

                            // 1. Slot conflict
                            if (slots_mask & orb.slots_mask) continue;

                            // 2. Check edge conflict
                            if (!is_disjoint(edges_mask, orb.edge_mask)) continue;

                            // 3. Check cycle compatibility with chosen orbits
                            bool compatible = true;
                            int f1 = orb.fids[0];
                            for (int cho_o_idx : chosen_orbits) {
                                for (int f2 : filtered_valid[cho_o_idx].fids) {
                                    if (!is_perfect_scalar(packed_pool[f1].adj, packed_pool[f2].adj)) {
                                        compatible = false;
                                        break;
                                    }
                                }
                                if (!compatible) break;
                            }
                            if (compatible) {
                                compatible_count++;
                                temp_compatible.push_back(o_idx);
                                if (compatible_count >= min_compatible) break; // Break early if we exceed min
                            }
    }

                        if (compatible_count < min_compatible) {
                            min_compatible = compatible_count;
                            s = slot;
                            compatible_orbits_for_s = std::move(temp_compatible);
                            if (min_compatible == 0) break; // Lookahead prune!
        }
    }

                    if (s == -1) return;
                    if (min_compatible == 0) return; // Prune branch

                    // Branch on the compatible orbits covering slot s
                    for (int o_idx : compatible_orbits_for_s) {
                        const auto& orb = filtered_valid[o_idx];

                        // Recurse
                        chosen_orbits.push_back(o_idx);
                        Mask18_C new_edges;
                        new_edges.m[0] = edges_mask.m[0] | orb.edge_mask.m[0];
                        new_edges.m[1] = edges_mask.m[1] | orb.edge_mask.m[1];
                        new_edges.m[2] = edges_mask.m[2] | orb.edge_mask.m[2];
                        self(self, slots_mask | orb.slots_mask, new_edges);
                        chosen_orbits.pop_back();
                    }
                };

                find_clique(find_clique, initial_slots_mask, initial_edges);
    }
        };

        TransInfo info;
        get_transformations_general(global_pool[0], global_pool[1], T1, T2, info);
        for (int p_idx = 0; p_idx < info.count; p_idx++) {
            const Permutation& alpha = info.perms[p_idx];
            process_permutation(alpha);
        }

        int processed = ++pairs_processed;
        if (m_bPrint) {
            auto current_time = std::chrono::steady_clock::now();
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - search_start).count();
            double elapsed_min = (double)elapsed_ms / 60000.0;
            double pct = (double)processed * 100.0 / total_pairs_to_process;
            double est_rem_min = (processed > 0) ? (elapsed_min * (total_pairs_to_process - processed) / processed) : 0.0;
            
            char t1_str[64] = "";
            char t2_str[64] = "";
            int p1 = 0, p2 = 0;
            for (int i = 0; i < 18; i += 2) {
                p1 += sprintf_s(t1_str + p1, 64 - p1, "%d-%d ", T1.src[i], T1.src[i + 1]);
                p2 += sprintf_s(t2_str + p2, 64 - p2, "%d-%d ", T2.src[i], T2.src[i + 1]);
            }
            printf("Processed pair transition %d/%lld (%.1f%% executed). Current T1: %s | T2: %s. Elapsed: %.1f min. Est. remaining: %.1f min\n",
                processed, total_pairs_to_process, pct, t1_str, t2_str, elapsed_min, est_rem_min);
            printf("  STATS [pair %d/%d]: clique_calls=%llu, lookahead_slots=%llu, lookahead_orbits=%llu, lookahead_failures=%llu\n",
                val1, val2, total_clique_calls, total_lookahead_slots, total_lookahead_orbits, total_lookahead_failures);
        }
    };

    if (m_bPrint)
        printf("DEBUG: Starting simplified transition search (6 starter transitions)...\n");

    // Process starter row target pairs only (6 transitions)
    run_search_for_pair(0, 0, 0, 1); 
    run_search_for_pair(0, 0, 0, 2); 
    run_search_for_pair(0, 1, 0, 0); 
    run_search_for_pair(0, 1, 0, 2); 
    run_search_for_pair(0, 2, 0, 0); 
    run_search_for_pair(0, 2, 0, 1); 

    if (m_bPrint) {
        auto end_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(end_time - search_start).count();
        printf("Transition search completed in %.4f seconds.\n", elapsed);
    }
}
