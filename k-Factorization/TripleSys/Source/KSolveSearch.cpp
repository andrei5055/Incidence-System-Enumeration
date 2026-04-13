#include "KSolveGen.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <immintrin.h>
#include <nmmintrin.h>
#include <set>

#define KGEN_VERIFY_PARITY 1
#define KGEN_USE_ROOT_CANONICITY 0 // Set to 1 to enable full root canonicity filtering (Slow)

void KSolveGen::runArcConsistency() {
    int Global_M_total = (int)global_pool.size();
    if (Global_M_total == 0) return;

    if (m_bPrint) {
        printf("KSolveGen: Running Iterative Global Filtering (Arc Consistency)...\n");
        printf("Total candidates: %d ", Global_M_total);
    }

    std::vector<uint8_t> active(Global_M_total, 1);
    std::vector<uint32_t> f_rows_mask(Global_M_total, 0);
    for (int s = 0; s < m_nSearch; s++) {
        for (int id : temp_slot_ids[s]) {
            if (id >= 0 && id < Global_M_total) f_rows_mask[id] |= (1 << s);
        }
    }

    bool changed = true;
    int iter = 0;
    while (changed) {
        changed = false;
        parallel_for(0, Global_M_total, [this, &active, &changed, &f_rows_mask](int i, int tid) {
            if (!active[i]) return;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < m_nSearch; r++) {
                if (f_rows_mask[i] & (1 << r)) continue;
                bool found = false;
                for (int id : temp_slot_ids[r]) {
                    if (active[id] && is_perfect_packed(pi, packed_pool[id])) {
                        found = true; break;
                    }
                }
                if (!found) { active[i] = 0; changed = true; break; }
            }
        }, 128);

        if (m_bPrint) {
            int rem = 0;
            for (int i = 0; i < Global_M_total; i++) if (active[i]) rem++;
            printf(" Remaining(%d):%d ", ++iter, rem);
        }
    }
    if (m_bPrint) printf("\n");

    // Re-index remaining active candidates to make pool compact
    std::vector<int> old_to_new(Global_M_total, -1);
    std::vector<Factor> f_new; std::vector<PackedAdj> p_new;
    for (int i = 0; i < Global_M_total; i++) {
        if (active[i]) {
            old_to_new[i] = (int)f_new.size();
            f_new.push_back(global_pool[i]);
            p_new.push_back(packed_pool[i]);
        }
    }
    global_pool.swap(f_new);
    packed_pool.swap(p_new);

    for (int s = 0; s < m_nSearch; s++) {
        std::vector<int> filtered;
        for (int id : temp_slot_ids[s]) {
            if (id >= 0 && id < (int)active.size() && active[id]) filtered.push_back(old_to_new[id]);
        }
        temp_slot_ids[s] = filtered;
    }
}

void KSolveGen::preSolve() {
    auto gen_start_time = std::chrono::steady_clock::now();
    if (m_bPrint) printf("KSolveGen: Presolving (Internal Generation + Arc Consistency)...\n");

    // 1. Capture Legacy Candidates if pool isn't empty (for parity verification)
    std::vector<std::set<std::vector<uint8_t>>> legacy_sets(m_nSearch);
    bool has_legacy = !global_pool.empty();
    if (has_legacy) {
        if (m_bPrint) printf("KSolveGen: Capturing %zd legacy candidates for parity check...\n", global_pool.size());
        for (int s = 0; s < m_nSearch; s++) {
            for (int id : temp_slot_ids[s]) {
                if (id >= 0 && id < (int)global_pool.size()) {
                    legacy_sets[s].insert(std::vector<uint8_t>(global_pool[id].adj, global_pool[id].adj + m_nPlayers));
                }
            }
        }
    }

    m_generateOnTheFly = true; 
    global_pool.clear(); 
    packed_pool.clear(); 
    f_map.clear();
    for (int s = 0; s < m_nSearch; s++) {
        temp_slot_ids[s].clear();
        row_factor_ids[s].clear();
    }

    parallel_for(0, m_nSearch, [this](int s, int tid) {
        generateSlotCandidates(s, -1, -1); 
    }, 1);

    if (m_bPrint) {
        printf("KSolveGen: Candidates generated (slot:count):");
        size_t total = 0;
        for (int s = 0; s < m_nSearch; s++) {
            size_t cnt = temp_slot_ids[s].size();
            printf(" %d:%zd", s, cnt);
            total += cnt;
        }
        printf("  Total: %zd\n", total);
    }

    // 2. Perform Parity Check immediately if we had legacy data
    if (has_legacy) {
        printf("\nKSolveGen: --- Parity Check (Legacy vs Generator) ---\n");
        int total_missing = 0;
        std::vector<std::set<std::vector<uint8_t>>> new_slot_sets(m_nSearch);
        for (int s = 0; s < m_nSearch; s++) {
            for (int new_id : temp_slot_ids[s]) {
                new_slot_sets[s].insert(std::vector<uint8_t>(global_pool[new_id].adj, global_pool[new_id].adj + m_nPlayers));
            }
        }

        for (int s = 0; s < m_nSearch; s++) {
            int s_missing = 0;
            for (const auto& leg_adj : legacy_sets[s]) {
                if (new_slot_sets[s].find(leg_adj) == new_slot_sets[s].end()) s_missing++;
            }
            printf("Slot %2d: Legacy=%6zd, Gen=%6zd, Missing=%d%s\n", 
                s, legacy_sets[s].size(), new_slot_sets[s].size(), s_missing, (s_missing > 0) ? " [FAIL]" : "");
            total_missing += s_missing;
        }
        if (total_missing == 0) printf("SUCCESS: All legacy candidates accounted for in generator.\n");
        printf("--------------------------------------------------\n\n");
        m_parityReportDone = true;
    }

    // 3. Store reference sets for solve(0) difference reporting
    m_referenceSets.clear();
    m_referenceSets.resize(m_nSearch);
    for (int s = 0; s < m_nSearch; s++) {
        for (int id : temp_slot_ids[s]) {
            m_referenceSets[s].insert(std::vector<uint8_t>(global_pool[id].adj, global_pool[id].adj + m_nPlayers));
        }
    }

    runArcConsistency();

    m_presolveDone = true;
    if (m_bPrint) {
        auto gen_end_time = std::chrono::steady_clock::now();
        auto gen_ms = std::chrono::duration_cast<std::chrono::milliseconds>(gen_end_time - gen_start_time).count();
        printf("KSolveGen: Presolve complete in %lld ms. Active candidates in pool: %zd\n", gen_ms, global_pool.size());
    }
}

int KSolveGen::getCandidateCount(int slot) {
    if (slot < 0 || slot >= m_nSearch) return 0;
    return (int)temp_slot_ids[slot].size();
}

bool KSolveGen::getCandidate(int slot, int index, unsigned char* out_src) {
    if (slot < 0 || slot >= m_nSearch || index < 0 || index >= (int)temp_slot_ids[slot].size()) return false;
    int id = temp_slot_ids[slot][index];
    if (id == -1) return false;
    memcpy(out_src, global_pool[id].src, m_nPlayers);
    return true;
}

void KSolveGen::deleteCandidate(int slot, int index) {
    if (slot >= 0 && slot < m_nSearch && index >= 0 && index < (int)temp_slot_ids[slot].size()) {
        temp_slot_ids[slot][index] = -1;
    }
}

void KSolveGen::solve(int mode) {
    auto solve_start = std::chrono::steady_clock::now();

    if (mode == 1) { 
        // Workflow B: Manual Slot Filtering (Canonization)
        if (m_bPrint) printf("KSolveGen: solve(1) - Compressing slots after external canonization...\n");
        for (int s = 0; s < m_nSearch; s++) {
            std::vector<int> compressed;
            for (int id : temp_slot_ids[s]) if (id != -1) compressed.push_back(id);
            temp_slot_ids[s] = compressed;
        }
    } else { 
        // Workflow A: Manual Capture & Verification
#if KGEN_VERIFY_PARITY
        if (m_bPrint) printf("KSolveGen: solve(0) - Running Parity Verification (Capture vs Reference)...\n");

        if (m_presolveDone && !m_referenceSets.empty()) {
            // Optimization: If preSolve() already happened, we can verify the CURRENT status 
            // against the PRE-SOLVE reference without re-generating everything.
            printf("\nKSolveGen: --- Candidate Difference Report (Reference vs Current) ---\n");
            int total_filtered = 0;
            for (int s = 0; s < m_nSearch; s++) {
                std::set<std::vector<uint8_t>> current_set;
                for (int id : temp_slot_ids[s]) {
                    if (id != -1 && id < (int)global_pool.size()) {
                        current_set.insert(std::vector<uint8_t>(global_pool[id].adj, global_pool[id].adj + m_nPlayers));
                    }
                }

                int s_missing = 0;
                for (const auto& ref_adj : m_referenceSets[s]) {
                    if (current_set.find(ref_adj) == current_set.end()) s_missing++;
                }
                printf("Slot %2d: Baseline=%6zd, Current=%6zd, Filtered Out=%d\n", 
                    s, m_referenceSets[s].size(), current_set.size(), s_missing);
                total_filtered += s_missing;
            }
            printf("----------------------------------------------------------\n\n");
        } 
        else {
            // Fallback: Full Capture & Re-generation (for cases where solve(0) is called without preSolve)
            std::vector<std::set<std::vector<uint8_t>>> legacy_sets(m_nSearch);
            size_t total_legacy_count = global_pool.size();
            if (total_legacy_count > 0) {
                for (int s = 0; s < m_nSearch; s++) {
                    for (int id : temp_slot_ids[s]) {
                        if (id >= 0 && id < (int)global_pool.size()) {
                            std::vector<uint8_t> adj(global_pool[id].adj, global_pool[id].adj + m_nPlayers);
                            legacy_sets[s].insert(adj);
                        }
                    }
                }
            }

            m_generateOnTheFly = true; 
            global_pool.clear(); 
            packed_pool.clear(); 
            f_map.clear();
            for (int s = 0; s < m_nSearch; s++) {
                temp_slot_ids[s].clear();
                row_factor_ids[s].clear();
            }

            auto gen_start_time = std::chrono::steady_clock::now();
            parallel_for(0, m_nSearch, [this](int s, int tid) {
                generateSlotCandidates(s, -1, -1); 
            }, 1);

            if (m_bPrint) {
                printf("KSolveGen: Candidates generated (slot:count):");
                size_t total = 0;
                for (int s = 0; s < m_nSearch; s++) {
                    size_t cnt = temp_slot_ids[s].size();
                    printf(" %d:%zd", s, cnt);
                    total += cnt;
                }
                printf("  Total: %zd\n", total);
            }

            if (total_legacy_count > 0) {
                printf("\nKSolveGen: --- Candidate Comparison Report ---\n");
                int total_missing = 0;
                std::vector<std::set<std::vector<uint8_t>>> new_slot_sets(m_nSearch);
                for (int s = 0; s < m_nSearch; s++) {
                    for (int new_id : temp_slot_ids[s]) {
                        new_slot_sets[s].insert(std::vector<uint8_t>(global_pool[new_id].adj, global_pool[new_id].adj + m_nPlayers));
                    }
                }

                for (int s = 0; s < m_nSearch; s++) {
                    int s_missing = 0;
                    for (const auto& leg_adj : legacy_sets[s]) {
                        if (new_slot_sets[s].find(leg_adj) == new_slot_sets[s].end()) s_missing++;
                    }
                    size_t new_count = temp_slot_ids[s].size();
                    size_t old_count = legacy_sets[s].size();
                    printf("Slot %2d: Legacy=%6zd, New=%6zd, Missing=%d%s\n", 
                        s, old_count, new_count, s_missing, (s_missing > 0) ? " [FAIL]" : "");
                    total_missing += s_missing;
                }
                if (total_missing == 0) {
                    auto gen_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - gen_start_time).count();
                    printf("SUCCESS: All legacy candidates are accounted for. (Gen: %lld ms)\n", gen_ms);
                }
                printf("-------------------------------------------\n\n");
            }
        }
#else
        if (m_bPrint) printf("KSolveGen: solve(0) - Running Arc Consistency on captured candidates...\n");
#endif
        runArcConsistency();
    }

    m_enableCoveragePruning = true; 
    m_enableUnitProp = true; 

    // Finalize items for search
    std::vector<SortItem> items;
    if (m_bPrint) printf("Slots (Active Candidates): (");
    for (int s = 0; s < m_nSearch; s++) {
        int s_active = 0;
        for (int id : temp_slot_ids[s]) {
            items.push_back({ id, s, 0 });
            s_active++;
        }
        if (m_bPrint) printf("%d%s", s_active, (s == m_nSearch - 1) ? "" : ", ");
        while (items.size() % 256 != 0) items.push_back({ -1, s, 0 });
    }
    if (m_bPrint) printf(")\n");
    
    int MS = (int)items.size();
    if (MS == 0) {
        if (m_bPrint) printf("KSolveGen: No candidates remaining after filtering. Search aborted.\n");
        return;
    }
    
    int Global_M_total = (int)global_pool.size();
    search_to_factor.assign(MS, -1);
    std::vector<int> g_to_item_idx(Global_M_total, -1);
    std::vector<int> fid_to_row(Global_M_total, -1);
    for (int i = 0; i < MS; ++i) {
        if (items[i].global_id != -1) {
            search_to_factor[i] = items[i].global_id;
            g_to_item_idx[items[i].global_id] = i;
            fid_to_row[items[i].global_id] = items[i].row_id;
        }
    }

    row_ranges.assign(m_nSearch + 1, { 0, 0 });
    int cp = 0;
    for (int s = 0; s < m_nSearch; ++s) {
        row_ranges[s].start_word = cp / 64;
        int count = 0; for (const auto& it : items) if (it.row_id == s) count++;
        cp += count; row_ranges[s].end_word = (cp + 63) / 64;
    }

    factor_offsets.assign(MS, 0); size_t coff = 0;
    int total_words = (row_ranges[m_nSearch - 1].end_word + 3) / 4 * 4;
    for (int i = 0; i < MS; i++) {
        factor_offsets[i] = coff; coff += total_words;
    }
    adj_matrix.assign(coff, 0);

    // 3. Build Adjacency Matrix
    if (m_bPrint) {
        printf("KSolveGen: Populating adj_matrix (%%): ");
    }
    
    std::atomic<int> matrix_done{ 0 };
    auto last_matrix_report = std::chrono::steady_clock::now();
    int report_step = MS / 20; if (report_step < 1) report_step = 1;

    parallel_for(0, MS, [&](int i, int tid) {
        if (items[i].global_id == -1) return;
        const PackedAdj& pi = packed_pool[items[i].global_id];
        int slot_i = items[i].row_id;
        uint64_t* row_i = &adj_matrix[factor_offsets[i]];
        
        for (int j = 0; j < MS; ++j) {
            if (i == j) continue;
            if (items[j].global_id == -1) continue;
            if (items[j].row_id == slot_i) continue;

            if (is_perfect_packed(pi, packed_pool[items[j].global_id])) {
                row_i[j >> 6] |= (1ULL << (j & 63));
            }
        }
        
        int done = ++matrix_done;
        if (done % report_step == 0 || done == MS) {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_matrix_report).count() >= 5 || done == MS) {
                last_matrix_report = now;
                if (m_bPrint) {
                    printf("%d ", (int)(done * 100LL / MS));
                }
            }
        }
    }, 64);
    
    if (m_bPrint) {
        printf("100 Done.\n");
    }


    // 4. Root Filtering with Canonicity
    std::vector<RootPair> potential_roots;
    std::mutex pr_mutex;
    if (m_bPrint) {
        printf("KSolveGen: Root Filtering (%%): ");
    }
    
    std::atomic<int> root_done{ 0 };
    int root_total = row_ranges[0].end_word * 64 - row_ranges[0].start_word * 64;
    int r_report_step = root_total / 20; if (r_report_step < 1) r_report_step = 1;
    std::atomic<int> root_before_canonicity{ 0 };
    auto last_root_report = std::chrono::steady_clock::now();

    parallel_for(0, root_total, [&](int i_raw, int tid) {
        int i = i_raw + row_ranges[0].start_word * 64;
        if (i >= MS || items[i].row_id != 0 || items[i].global_id == -1) {
            int done = ++root_done;
            return;
        }
        
        ThreadLocalBuffers& b = *thread_buffers[tid];
        
        std::vector<RootPair> local_roots;
        for (int j = row_ranges[1].start_word * 64; j < row_ranges[1].end_word * 64; j++) {
            if (j >= MS || items[j].row_id != 1 || items[j].global_id == -1) continue;
            if (memcmp(global_pool[items[i].global_id].fs.pairs, 
                       global_pool[items[j].global_id].fs.pairs, m_nPlayers) >= 0) continue;
            
            if (adj_matrix[factor_offsets[i] + (j >> 6)] & (1ULL << (j & 63))) {
                root_before_canonicity++;
                
                bool ok = true;
                /**
                 * AI AGENT PERSISTENCE NOTE: 
                 * DO NOT REMOVE THIS BLOCK. Root canonicity filtering is essential for 
                 * non-bipartite search space reduction. It is currently disabled via 
                 * macro for performance, but must remain available for complex test cases.
                 */
                if (KGEN_USE_ROOT_CANONICITY || !m_bipartite) {
                    ok = is_canonical(items[i].global_id, items[j].global_id, nullptr, &b);
                }

                if (ok) {
                    local_roots.push_back({ i, j });
                }
            }
        }
        if (!local_roots.empty()) {
            std::lock_guard<std::mutex> lock(pr_mutex);
            potential_roots.insert(potential_roots.end(), local_roots.begin(), local_roots.end());
        }

        int done = ++root_done;
        if (done % r_report_step == 0 || done == root_total) {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_root_report).count() >= 5 || done == root_total) {
                last_root_report = now;
                if (m_bPrint) {
                    printf("%d ", (int)(done * 100LL / root_total));
                }
            }
        }
    }, 1);

    int total_roots = (int)potential_roots.size();
    if (m_bPrint) {
        printf("100 Done.\n");
        printf("KSolveGen: %d/%d root pairs before/after canonicity filtering.\n", (int)root_before_canonicity, total_roots);
    }
    
    // 6. Parallel Search with Matrix Compression
    std::vector<RootGroup> r4_groups;
    if (!potential_roots.empty()) {
        std::sort(potential_roots.begin(), potential_roots.end(), [](const RootPair& a, const RootPair& b) {
            return a.r4_v < b.r4_v;
        });
        int start = 0;
        for (int i = 1; i < (int)potential_roots.size(); i++) {
            if (potential_roots[i].r4_v != potential_roots[start].r4_v) {
                r4_groups.push_back({ start, i, potential_roots[start].r4_v });
                start = i;
            }
        }
        r4_groups.push_back({ start, (int)potential_roots.size(), potential_roots[start].r4_v });
    }

    start_time = std::chrono::steady_clock::now();
    const uint64_t* p_adj_matrix = adj_matrix.data();
    const size_t* p_factor_offsets = factor_offsets.data();
    const int* p_g_to_item_idx = g_to_item_idx.data();
    const int* p_fid_to_row = fid_to_row.data();

    // 6. Parallel Search with Matrix Compression
    parallel_for(0, (int)r4_groups.size(), [&](int g_idx, int tid) {
        ThreadLocalBuffers& b = *thread_buffers[tid];
        const RootGroup& group = r4_groups[g_idx];
        int r4_idx = group.r4_v;
        int r4_fid = items[r4_idx].global_id;
        const uint64_t* r4_bits_global = &adj_matrix[factor_offsets[r4_idx]];

        // 6.1 Identify Compatible Subset for this R4 group
        if (b.g_to_l.empty()) b.g_to_l.assign(Global_M_total, -1);
        b.local_to_global.clear();
        int l_r_start[32], l_r_end[32];
        for (int r = 1; r < m_nSearch; r++) { // R5 to Rn
            while (b.local_to_global.size() % 256 != 0) b.local_to_global.push_back(-1);
            l_r_start[r] = (int)b.local_to_global.size();
            for (int w = row_ranges[r].start_word; w < row_ranges[r].end_word; w++) {
                uint64_t mask = r4_bits_global[w];
                while (mask) {
                    int bit = (int)_tzcnt_u64(mask);
                    int g_idx_f = w * 64 + bit;
                    if (g_idx_f < MS && items[g_idx_f].global_id != -1) {
                        int fid = items[g_idx_f].global_id;
                        b.g_to_l[fid] = (int)b.local_to_global.size();
                        b.local_to_global.push_back(fid);
                    }
                    mask &= (mask - 1);
                }
            }
            l_r_end[r] = (int)b.local_to_global.size();
        }

        int LK = (int)b.local_to_global.size();
        if (LK == 0) { m_rootsDone += (group.end_idx - group.start_idx); return; }

        // 6.2 Build Sub-Sampling Plans for this Root group
        for (int r = 1; r < m_nSearch; r++) {
            b.sub_sampling_plan_sizes[r] = 0;
            int bit_cursor = l_r_start[r];
            for (int w = row_ranges[r].start_word; w < row_ranges[r].end_word; w++) {
                uint64_t mask = r4_bits_global[w];
                if (mask) {
                    if (b.sub_sampling_plan_sizes[r] >= 4096) {
                        printf("ERROR: Sub-sampling plan overflow in slot %d (tid %d). Matrix too sparse/large.\n", r, tid);
                        exit(1);
                    }
                    SubSamplePlanItem& item = b.sub_sampling_plans[r][b.sub_sampling_plan_sizes[r]++];
                    item.global_word_idx = (uint16_t)w;
                    item.mask = mask;
                    item.local_word_off = (uint16_t)(bit_cursor / 64);
                    item.local_bit_shift = (uint8_t)(bit_cursor % 64);
                    item.bits_pushed = (uint8_t)_mm_popcnt_u64(mask);
                    item.split = (item.local_bit_shift + item.bits_pushed > 64) ? 1 : 0;
                    bit_cursor += item.bits_pushed;
                }
            }
        }

        // 6.3 Build Local Matrix & Edge Presence (Padded for SIMD safety)
        int l_words = (LK + 63) / 64;
        l_words = (l_words + 3) / 4 * 4; // Pad to 256-bit boundary
        b.reserve_aligned_adj((size_t)LK * l_words + 4); 
        
        b.local_offsets.assign(LK, 0);
        b.local_edge_masks.assign(LK, Mask512());
        b.s_to_f.assign(LK, -1);
        
        // Clear local edge presence
        for (int e = 0; e < 512; e++) memset(b.local_edge_presence[e].bits, 0, l_words * 8);

        for (int i = 0; i < LK; i++) {
            int fid = b.local_to_global[i];
            if (fid == -1) continue;
            b.local_offsets[i] = (size_t)i * l_words;
            const Mask512& m = packed_pool[fid].edge_mask;
            b.local_edge_masks[i] = m;
            b.s_to_f[i] = fid;
            
            // Populate transposed index (Optimized for sparse edge masks)
            for (int w = 0; w < 8; w++) {
                uint64_t mask = m.m[w];
                while (mask) {
                    unsigned long e_bit;
                    _BitScanForward64(&e_bit, mask);
                    int e = (w << 6) | e_bit;
                    b.local_edge_presence[e].bits[i >> 6] |= (1ULL << (i & 63));
                    mask &= (mask - 1);
                }
            }
        }
        
        uint64_t* l_adj_raw = b.local_adj_aligned;
        for (int i = 0; i < LK; i++) {
            int fid = b.local_to_global[i];
            if (fid == -1) continue;
            uint64_t* l_row = &l_adj_raw[b.local_offsets[i]];
            const uint64_t* g_row_base = &p_adj_matrix[p_factor_offsets[p_g_to_item_idx[fid]]];
            for (int r = 1; r < m_nSearch; r++) {
                for (int p_idx = 0; p_idx < b.sub_sampling_plan_sizes[r]; p_idx++) {
                    const auto& p = b.sub_sampling_plans[r][p_idx];
                    uint64_t extracted = _pext_u64(g_row_base[p.global_word_idx], p.mask);
                    l_row[p.local_word_off] |= (extracted << p.local_bit_shift);
                    if (p.split) l_row[p.local_word_off + 1] |= (extracted >> (64 - p.local_bit_shift));
                }
            }
        }

        for (int r = 1; r < m_nSearch; r++) {
            b.local_ranges[r].start_word = l_r_start[r] / 64;
            b.local_ranges[r].end_word = (l_r_end[r] + 63) / 64;
        }

        b.ranges_ptr = b.local_ranges;
        b.adj_ptr = l_adj_raw;
        b.offsets_ptr = b.local_offsets.data();
        
        // 6.4 Process Root Pairs in this group
        for (int root_idx = group.start_idx; root_idx < group.end_idx; root_idx++) {
            int r5_idx = potential_roots[root_idx].r5_v;
            int r5_fid = items[r5_idx].global_id;
            int l5_idx = b.g_to_l[r5_fid];
            if (l5_idx == -1) continue;
            
            SearchContext& ctx = b.local_ctx;
            ctx.r4_idx = group.r4_v;
            ctx.r5_idx = r5_idx;
            ctx.slots[0] = 0; ctx.slots[1] = 1; 
            for (int s = 0; s < m_nSearch - 2; s++) ctx.slots[s + 2] = (uint8_t)(s + 2); // Variable rows
            memset(ctx.counts_valid, 0, sizeof(ctx.counts_valid));
            memset(ctx.used_edges, 0, sizeof(ctx.used_edges));
            
            ctx.used_edges[0] = fixedEdgesMask;
            ctx.used_edges[0].union_with(global_pool[r4_fid].edge_mask);
            ctx.used_edges[1] = ctx.used_edges[0];
            ctx.used_edges[1].union_with(global_pool[r5_fid].edge_mask);
            
            const uint64_t* l5_adj = &b.local_adj_aligned[b.local_offsets[l5_idx]];
            for (int v = 0; v < l_words; v++) ctx.pool[2].bits[v] = l5_adj[v];
            
            std::vector<int> clique = { r4_fid, r5_fid };
            clique.reserve(m_nSearch);
            
            internal_solve(2, clique, ctx, &b);
            m_rootsDone++;
            if(m_bPrint)
                diagnostic_printout(total_roots);
        }
        for (int fid : b.local_to_global) if (fid != -1) b.g_to_l[fid] = -1;
    }, 1);


    if (m_bPrint) { printf("Reporting %zd sorted results...\n", results_to_sort.size()); }
    // Sort canonical results lexicographically
    std::sort(results_to_sort.begin(), results_to_sort.end(), [](const std::vector<unsigned char>& a, const std::vector<unsigned char>& b) {
        return memcmp(a.data(), b.data(), a.size()) < 0;
        });

    // Final callback loop with mode=0
    for (const auto& res : results_to_sort) {
        resultCallback(cbClass, res.data(), 0, 0, 0);
    }
    std::vector<std::vector<unsigned char>>().swap(results_to_sort);
}

void VECTOR_CALL KSolveGen::internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx, ThreadLocalBuffers* buf) {
    const uint64_t* l_adj = buf->adj_ptr;
    const CompRange* l_ranges = buf->ranges_ptr;
    const size_t* l_offsets = buf->offsets_ptr;
    const int* l_to_f = buf->s_to_f.data();

    if (depth == m_nSearch) { // Found a full factorization
        unsigned char res[1024];
        for (int i = 0; i < 3; i++) memcpy(res + i * m_nPlayers, fixedRows[i].src, m_nPlayers);
        
        // clique[0] is R4 (slot 0), clique[1] is R5 (slot 1)
        // for d >= 2, clique[d] is for slot ctx.slots[d]
        memcpy(res + (0 + 3) * m_nPlayers, global_pool[clique[0]].src, m_nPlayers);
        memcpy(res + (1 + 3) * m_nPlayers, global_pool[clique[1]].src, m_nPlayers);
        for (int i = 2; i < m_nSearch; i++) {
            memcpy(res + (ctx.slots[i] + 3) * m_nPlayers, global_pool[clique[i]].src, m_nPlayers);
        }
        
        std::lock_guard<std::mutex> lock(result_mutex);
        if (resultCallback(cbClass, (const unsigned char*)res, ctx.r4_idx, ctx.r5_idx, 1)) {
            m_resultsFound++;
            results_to_sort.push_back(std::vector<unsigned char>(res, res + (m_nSearch + 3) * m_nPlayers));
        }
        return;
    }

    int best_depth_idx = depth;
    int min_count = 1000000;

    // 1. MRV (Minimum Remaining Values) Slot Selection
    const int early_exit_threshold = m_useOptimized ? 2 : 1; 
    if (ctx.counts_valid[depth]) {
        for (int i = depth; i < m_nSearch; i++) {
            int slot = ctx.slots[i];
            int count = ctx.mrv_counts[depth][slot];
            if (count < min_count) {
                min_count = count;
                best_depth_idx = i;
                if (min_count <= early_exit_threshold) break;
            }
        }
    } else {
        const uint64_t* p_bits = ctx.pool[depth].bits;
        for (int i = depth; i < m_nSearch; i++) {
            int slot = ctx.slots[i];
            int count = 0;
            const int sw = l_ranges[slot].start_word;
            const int ew = l_ranges[slot].end_word;
            for (int w = sw; w < ew; w++) count += (int)__popcnt64(p_bits[w]);
            ctx.mrv_counts[depth][slot] = (uint16_t)count;
            if (count < min_count) {
                min_count = count;
                best_depth_idx = i;
                if (min_count <= early_exit_threshold) break;
            }
        }
        ctx.counts_valid[depth] = true;
    }

    if (min_count == 0) return;

    int target_slot = ctx.slots[best_depth_idx];
    int forced_candidate = -1;

    // 2. Advanced Pruning & Unit Propagation (Non-Bipartite Only)
    // Thresholds tuned to k16p1f strategy: Prune early, prune hard
    const int pruning_activation = m_useOptimized ? 500 : 150;
    const int pruning_start_depth = m_useOptimized ? 9 : 6;
    if (!m_bipartite && depth >= pruning_start_depth && depth < m_nSearch - 1 && m_enableCoveragePruning && min_count <= 60) {
        int total_remaining = 0;
        if (m_useOptimized) {
            // Faster total_remaining check
            const uint64_t* p_bits = ctx.pool[depth].bits;
            for (int i = depth; i < m_nSearch; i++) {
                total_remaining += ctx.mrv_counts[depth][ctx.slots[i]];
                if (total_remaining > pruning_activation) break;
            }
        } else {
            for (int i = depth; i < m_nSearch; i++) total_remaining += ctx.mrv_counts[depth][ctx.slots[i]];
        }

        if (total_remaining <= pruning_activation) {
            Mask512 needed;
            Mask512 current_used = ctx.used_edges[depth - 1];
            for (int i = 0; i < 8; i++) needed.m[i] = total_edges_mask.m[i] & ~current_used.m[i];

            int forced_v = -1;
            bool conflict = false;

            // SIMD-Accelerated Edge Coverage Check
            const uint64_t* p_bits = ctx.pool[depth].bits;
            const int LK = (int)buf->local_to_global.size();
            const int l_words = (LK + 63) / 64;
            const int l_v4 = (l_words + 3) / 4;

            const int max_edge_check = m_useOptimized ? 256 : 512;

            for (int e = 0; e < max_edge_check; e++) {
                if (!(needed.m[e >> 6] & (1ULL << (e & 63)))) continue;

                int e_count = 0;
                int last_v = -1;
                const uint64_t* e_presence = buf->local_edge_presence[e].bits;

                for (int v = 0; v < l_v4; v++) {
                    __m256i p_vec = _mm256_loadu_si256((const __m256i*)&p_bits[v * 4]);
                    __m256i e_vec = _mm256_loadu_si256((const __m256i*)&e_presence[v * 4]);
                    __m256i res_vec = _mm256_and_si256(p_vec, e_vec);
                    if (!_mm256_testz_si256(res_vec, res_vec)) {
                        uint64_t m0 = _mm256_extract_epi64(res_vec, 0);
                        uint64_t m1 = _mm256_extract_epi64(res_vec, 1);
                        uint64_t m2 = _mm256_extract_epi64(res_vec, 2);
                        uint64_t m3 = _mm256_extract_epi64(res_vec, 3);
                        e_count += (int)(__popcnt64(m0) + __popcnt64(m1) + __popcnt64(m2) + __popcnt64(m3));
                        if (e_count == 1) {
                            if (m0) last_v = (v * 4 + 0) * 64 + (int)_tzcnt_u64(m0);
                            else if (m1) last_v = (v * 4 + 1) * 64 + (int)_tzcnt_u64(m1);
                            else if (m2) last_v = (v * 4 + 2) * 64 + (int)_tzcnt_u64(m2);
                            else if (m3) last_v = (v * 4 + 3) * 64 + (int)_tzcnt_u64(m3);
                        }
                        if (e_count > 1) break;
                    }
                }

                if (e_count == 0) { conflict = true; break; }
                if (e_count == 1) { forced_v = last_v; break; }
            }

            if (conflict) return;

            if (forced_v != -1) {
                forced_candidate = forced_v;
                for (int i = depth; i < m_nSearch; i++) {
                    int s_idx = ctx.slots[i];
                    if (forced_v >= buf->local_ranges[s_idx].start_word * 64 && 
                        forced_v < (int)(buf->local_ranges[s_idx].end_word * 64)) {
                        best_depth_idx = i;
                        break;
                    }
                }
                goto found_forced;
            }
        }
    }
    found_forced:;

    // 3. Search Loop
    std::swap(ctx.slots[depth], ctx.slots[best_depth_idx]);
    int s = ctx.slots[depth];
    int sw = l_ranges[s].start_word;
    int ew = l_ranges[s].end_word;
    uint64_t* p_pool = ctx.pool[depth].bits;

    for (int w = sw; w < ew; w++) {
        uint64_t word = p_pool[w];
        while (word) {
            unsigned long bit;
            _BitScanForward64(&bit, word);
            int l_idx = w * 64 + (int)bit;
            if (forced_candidate != -1 && l_idx != forced_candidate) { word &= (word - 1); continue; }

            int fid = buf->s_to_f[l_idx];
            clique.push_back(fid);
            ctx.used_edges[depth] = ctx.used_edges[depth - 1];
            ctx.used_edges[depth].union_with(buf->local_edge_masks[l_idx]);

            const uint64_t* mat_row = &buf->adj_ptr[buf->offsets_ptr[l_idx]];
            const int next_depth = depth + 1;
            bool domain_wipeout = false;

            for (int i = next_depth; i < m_nSearch; i++) {
                int ns = ctx.slots[i];
                int n_sw = l_ranges[ns].start_word;
                int n_ew = l_ranges[ns].end_word;
                int words_v = (n_ew - n_sw + 3) / 4;
                __m256i row_acc = _mm256_setzero_si256();
                int slot_count = 0;

                for (int v = 0; v < words_v; v++) {
                    int off = n_sw + v * 4;
                    __m256i v1 = _mm256_loadu_si256((const __m256i*)&p_pool[off]);
                    __m256i v2 = _mm256_loadu_si256((const __m256i*)&mat_row[off]);
                    __m256i res = _mm256_and_si256(v1, v2);
                    _mm256_storeu_si256((__m256i*)&ctx.pool[next_depth].bits[off], res);
                    row_acc = _mm256_or_si256(row_acc, res);
                    
                    if (m_useOptimized && depth < 8) { // MRV limit from k16p1f
                        // Inline popcount for MRV
                        slot_count += (int)(__popcnt64(_mm256_extract_epi64(res, 0)) + 
                                            __popcnt64(_mm256_extract_epi64(res, 1)) + 
                                            __popcnt64(_mm256_extract_epi64(res, 2)) + 
                                            __popcnt64(_mm256_extract_epi64(res, 3)));
                    }
                }
                if (_mm256_testz_si256(row_acc, row_acc)) { domain_wipeout = true; break; }
                if (m_useOptimized && depth < 8) ctx.mrv_counts[next_depth][ns] = (uint16_t)slot_count;
            }

            if (!domain_wipeout) {
                // Ensure MRV cache is marked invalid for lower depths to force recalculation if needed
                // though usually we recalculate for depth < 8 anyway.
                internal_solve(next_depth, clique, ctx, buf);
            }
            clique.pop_back();
            if (forced_candidate != -1) break; 
            word &= (word - 1);
        }
    }
    std::swap(ctx.slots[depth], ctx.slots[best_depth_idx]);
}

void KSolveGen::diagnostic_printout(int total_roots) {
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::seconds>(now - last_report_time).count() < 30) return;
    std::lock_guard<std::mutex> lock(stdout_mutex);
    if (std::chrono::duration_cast<std::chrono::seconds>(now - last_report_time).count() < 30) return;
    last_report_time = now;
    auto elap = (int)std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
    int done = m_rootsDone.load();
    double pct = total_roots ? (done * 100.0 / total_roots) : 100.0;
    int expected = (pct > 0.001) ? (int)(elap * 100.0 / pct + 0.5) : 0;
    printf("[RUN] %04d | %%:%6.3f | T/E:%4d/%4d | Found:%d | Threads:%d\n", 
           fixed3RowsIndex, pct, elap, expected, m_resultsFound.load(), kThreads);
}
