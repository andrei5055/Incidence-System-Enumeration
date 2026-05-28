#include "k16p1f.h"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <format>
#include <iterator>
#include <utility>

K16P1F::K16P1F(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr, bool bPrint) : KBase<Mask256_C>(factParam, bPrint) {
    int id_cnt = 0;
    memset(edge_id_table, -1, sizeof(edge_id_table));
    for (int i = 0; i < K16_N; i++) {
        for (int j = i + 1; j < K16_N; j++) {
            int eid = id_cnt++;
            edge_id_table[i][j] = eid;
            edge_id_table[j][i] = eid;
            if (eid < 128) {
                edge_to_u[eid] = (uint8_t)i;
                edge_to_v[eid] = (uint8_t)j;
            }
        }
    }
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }

    for (int m = 0; m < 256; m++) {
        uint8_t adj_mask[16], idx_mask[16];
        memset(adj_mask, 0xFF, 16);
        memset(idx_mask, 0, 16);
        uint8_t count = 0;
        for (int i = 0; i < 8; i++) {
            if (m & (1 << i)) {
                adj_mask[count] = i;
                idx_mask[count] = i;
                count++;
            }
        }
        gfs_table_adj[m] = _mm_loadu_si128((const __m128i*)adj_mask);
        gfs_table_idx[m] = _mm_loadu_si128((const __m128i*)idx_mask);
    }

    init(fixed3RowsIndex, kThreads, first3Rows, callback, cbClassPtr);
}

void K16P1F::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K16 p1f, MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
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

    // Reset transient state for reuse
    global_pool.clear();
    packed_pool.clear();
    f_map.clear();
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }
    KBase::init();
}

bool K16P1F::addRow(int iRow, const unsigned char* source) {
    iRow--;
    if (iRow < m_nFixedRows || iRow >= K16_MATCH) return false;
    int slot = iRow - m_nFixedRows;

    Factor f; memset(&f, 0, sizeof(f));
    memcpy(f.src, source, K16_N);
    for (int i = 0; i < K16_N; i += 2) { 
        const auto a = f.src[i]; const auto b = f.src[i + 1]; 
        if (a >= K16_N || b >= K16_N) return false; // Safety check
        f.adj[a] = b; f.adj[b] = a; 
    }

    f.edge_mask = Mask256_C();
    for (int i = 0; i < K16_N; i += 2) {
        int eid = edge_id_table[f.src[i]][f.src[i+1]];
        if (eid < 64) f.edge_mask.m[0] |= (1ULL << eid);
        else if (eid < 128) f.edge_mask.m[1] |= (1ULL << (eid - 64));
    }

    if ((f.edge_mask.m[0] & fixedEdgesMask.m[0]) || (f.edge_mask.m[1] & fixedEdgesMask.m[1])) {
        printf("\n*** Error: Row has common edges with fixed rows!\n");
        exit(1);
    }

    // Check cycle structure with fixed rows (r1, r2, r3)
    PackedAdj pf = pack_factor_adj(f.adj);
    for (int i = 0; i < m_nFixedRows; i++) {
        if (!is_perfect_packed(pf, fixed_packed[i])) {
            printf("\n*** Error: Row has incorrect cycle structure with fixed row %d!\n", i + 1);
            exit(1);
        }
    }

    f.fs = get_fast_sorted(f.adj);

    std::vector<uint8_t> k(f.adj, f.adj + K16_N);
    int factor_id;
    auto it = f_map.find(k);
    if (it == f_map.end()) {
        factor_id = (int)global_pool.size();
        if (factor_id >= K16_M_MAX) return false;
        f_map[k] = factor_id;
        global_pool.push_back(f);
        packed_pool.push_back(pack_factor_adj(f.adj));
    } else {
        factor_id = it->second;
    }

    addFactorToSlot(slot, factor_id);
    return true;
}

void K16P1F::parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size) {
    if (kThreads <= 1 || end <= start) {
        for (int i = start; i < end; ++i) task(i, 0);
        return;
    }
    std::atomic<int> next_idx{ start };
    std::vector<std::thread> workers;
    for (int t = 0; t < kThreads; ++t) {
        workers.emplace_back([&, t]() {
            while (true) {
                int i = next_idx.fetch_add(grain_size);
                if (i >= end) break;
                int limit = std::min(i + grain_size, end);
                for (int j = i; j < limit; ++j) {
                    task(j, t);
                }
            }
        });
    }
    for (auto& w : workers) w.join();
}

void K16P1F::solve(int mode) {
    int Global_M_total = (int)global_pool.size();
    if (m_bPrint) {
        printf("Row slots size: ");
        for (int s = 0; s < K16_SEARCH; s++) printf("%d:%zd ", s * 2 + 7, temp_slot_ids[s].size());
        printf("\n");
        printf("Iterative Global Filtering (Arc Consistency):\n");
        printf("Total candidates: %d ", Global_M_total);
    }
    std::vector<uint8_t> active(Global_M_total, 1);
    std::vector<uint32_t> f_rows_mask(Global_M_total, 0);
    for (int s = 0; s < K16_SEARCH; s++) for (int id : temp_slot_ids[s]) f_rows_mask[id] |= (1 << s);

    std::atomic<bool> changed{ true };
    int iter = 0;
    while (changed) {
        changed = false;
        parallel_for(0, Global_M_total, [&](int i, int tid) {
            if (!active[i]) return;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < K16_SEARCH; r++) {
                if (f_rows_mask[i] & (1 << r)) continue; 
                bool row_found = false;
                for (int id : temp_slot_ids[r]) {
                    if (active[id] && is_perfect_packed(pi, packed_pool[id])) { row_found = true; break; }
                }
                if (!row_found) {
                    active[i] = 0;
                    changed = true;
                    break;
                }
            }
        }, 128);
        int rem = 0; for (int i = 0; i < Global_M_total; i++) if (active[i]) rem++;
        iter++;
        if (m_bPrint) { printf(" Remaining(%d):%d ", iter, rem); }
        if (rem == 0) {
            printf("\nAll candidates filtered by Arc consistency check\n");
        }
    }
    if (m_bPrint) { printf("\n"); }

    std::vector<int> old_to_new(Global_M_total, -1);
    std::vector<Factor> f_pool_new; std::vector<PackedAdj> p_pool_new;
    for (int i = 0; i < Global_M_total; i++) if (active[i]) { old_to_new[i] = (int)f_pool_new.size(); f_pool_new.push_back(global_pool[i]); p_pool_new.push_back(packed_pool[i]); }
    global_pool = std::move(f_pool_new); packed_pool = std::move(p_pool_new);
    
    std::vector<std::vector<int>> new_slot_ids(K16_SEARCH);
    if (m_bPrint) printf("Row: Candidates ");
    for (int s = 0; s < K16_SEARCH; s++) {
        for (int id : temp_slot_ids[s]) 
            if (active[id]) new_slot_ids[s].push_back(old_to_new[id]);
        if (m_bPrint) { printf("%d:%d ", s * 2 + 7, (int)new_slot_ids[s].size()); }
    }
    if (m_bPrint) { printf("\n"); }

    const int M = (int)global_pool.size();
    std::vector<int> deg(M, 0);
    if (m_bPrint) { printf("Calculating candidate degrees(%%): 0 "); }
    auto last_deg_report = std::chrono::steady_clock::now();
    parallel_for(0, M, [&](int i, int tid) {
        PackedAdj pi = packed_pool[i];
        for (int r = 0; r < K16_SEARCH; r++) for (int id : new_slot_ids[r]) if (is_perfect_packed(pi, packed_pool[id])) deg[i]++;
        
        if (tid == 0) {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_deg_report).count() >= 5) {
                last_deg_report = now;
                if (m_bPrint) { printf("%d ", (int)(i * 100LL / M)); }
            }
        }
    }, 256);
    if (m_bPrint) {
        printf("100%% Done.\n");
    }

    struct SortItem { int global_id, row_id, degree; };
    std::vector<SortItem> items_sorted;
    for (int s = 0; s < K16_SEARCH; ++s) {
        std::vector<SortItem> row_items;
        for (int id : new_slot_ids[s]) row_items.push_back({ id, s, deg[id] });
        /*
        std::sort(row_items.begin(), row_items.end(), [](const SortItem& a, const SortItem& b) {
            if (K16_SORT_INPUT == 1) return a.degree < b.degree;
            if (K16_SORT_INPUT == 2) return a.degree > b.degree;
            return false; 
            });*/
        while (items_sorted.size() % 256 != 0) items_sorted.push_back({ -1, s, 0 }); 
        items_sorted.insert(items_sorted.end(), row_items.begin(), row_items.end());
    }
    const std::vector<SortItem>& items = items_sorted;

    int MS = (int)items.size(); 
    search_to_factor.assign(MS, -1); 
    for (int i = 0; i < MS; ++i) if (items[i].global_id != -1) search_to_factor[i] = items[i].global_id;
    
    factor_edge_masks.assign(MS, Mask256_C());
    for (int i = 0; i < MS; i++) {
        if (items[i].global_id == -1) continue;
        const Factor& f = global_pool[items[i].global_id];
        for (int u = 0; u < K16_N; u++) { 
            int v = f.adj[u];
            if (u < v) {
                int id = edge_id_table[u][v];
                if (id != -1) {
                    if (id < 64) factor_edge_masks[i].m[0] |= (1ULL << id);
                    else if (id < 128) factor_edge_masks[i].m[1] |= (1ULL << (id - 64));
                }
            }
        }
    }

    row_ranges.assign(K16_SEARCH + 1, { 0, 0 });
    int cp = 0;
    for (int s = 0; s < K16_SEARCH; ++s) {
        row_ranges[s].start_word = cp / 64;
        int count = 0; for (const auto& it : items) if (it.row_id == s) count++;
        cp += count; row_ranges[s].end_word = (cp + 63) / 64;
    }
    row_ranges[K16_SEARCH].start_word = cp / 64;
    row_ranges[K16_SEARCH].end_word = (cp + 63) / 64;

    factor_offsets.assign(MS, 0); size_t coff = 0;
    const int total_words = (row_ranges[K16_SEARCH - 1].end_word + 3) / 4 * 4;
    for (int i = 0; i < MS; i++) {
        factor_offsets[i] = coff; coff += total_words;
    }
    adj_matrix.assign(coff, 0);
    if (m_bPrint) { printf("Populating adj_matrix (%%): "); }
    auto last_matrix_report = std::chrono::steady_clock::now();

    std::vector<PackedAdj> sorted_packed(MS);
    for (int i = 0; i < MS; i++) {
        if (items[i].global_id != -1) sorted_packed[i] = packed_pool[items[i].global_id];
        else memset(&sorted_packed[i], 0, sizeof(PackedAdj));
    }

    std::atomic<int> matrix_done{ 0 };
    int report_step = MS / 20; if (report_step < 1) report_step = 1;

    parallel_for(0, MS, [&](int i, int tid) {
        if (items[i].global_id == -1) return;
        const PackedAdj& pi = sorted_packed[i];
        const uint8_t* adj_i = global_pool[items[i].global_id].adj;
        __m256i mi = _mm256_load_si256((const __m256i*)&pi);
        
        uint64_t* row_i = &adj_matrix[factor_offsets[i]];
        int sw = 0;
        int j_start = 0;

        const Mask256_C& mi_e = factor_edge_masks[i];
        for (int j = j_start; j < MS; ++j) {
            if (items[j].global_id == -1) continue;
            const PackedAdj& pj = sorted_packed[j];
            if (is_perfect_packed(pi, pj)) {
                const Mask256_C& mj_e = factor_edge_masks[j];
                if (!((mi_e.m[0] & mj_e.m[0]) | (mi_e.m[1] & mj_e.m[1]))) {
                    row_i[(j / 64) - sw] |= (1ULL << (j % 64));
                }
            }
        }
        int done = ++matrix_done;
        if (done % report_step == 0 || done == MS) {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_matrix_report).count() >= 5 || done == MS) {
                last_matrix_report = now;
                if (m_bPrint) { printf("%d ", (int)(done * 100LL / MS)); }
        }
    }
    }, 64);
    if (m_bPrint) { printf("100%% Done.\n"); }


    std::vector<int> r4ids;
    for (int idx = row_ranges[0].start_word * 64; idx < row_ranges[0].end_word * 64; idx++) {
        if (idx < MS && items[idx].row_id == 0 && items[idx].global_id != -1) r4ids.push_back(idx);
    }

    struct RootPair { int r4_v; int r5_v; };
    std::vector<RootPair> potential_roots;
    if (m_bPrint) { printf("Gathering potential Root Pairs and filtering (%%): "); }
    auto last_canon_report = std::chrono::steady_clock::now();
    auto canon_start = last_canon_report;
    std::atomic<int> canon_done{ 0 };
    std::atomic<int> r4_branches_rejected{ 0 };
    int np = (int)r4ids.size();
    int c_report_step = np / 20; if (c_report_step < 1) c_report_step = 1;

    std::mutex pr_mutex;
    parallel_for(0, np, [&](int i, int tid) {
        int r4_v = r4ids[i];
        int r4_fid = search_to_factor[r4_v];
        
        auto r4_trans = std::make_unique<TransInfo[]>(3);
        for (int j = 0; j < 3; j++) get_transformations(fixedRows[j], global_pool[r4_fid], r4_trans[j]);

        PackedAdj p4 = packed_pool[r4_fid];
        const Mask256_C& m4 = factor_edge_masks[r4_v];

        std::vector<RootPair> local_roots;
        int local_checked = 0, local_kept = 0;
        int local_rej_edge = 0, local_rej_cycle = 0, local_rej_canon = 0;

        const __m128i v_e2o_4 = _mm_load_si128(&p4.v_e2o);
        const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        for (int r5_v = row_ranges[1].start_word * 64; r5_v < row_ranges[1].end_word * 64; r5_v++) {
            if (r5_v >= MS || items[r5_v].row_id != 1 || items[r5_v].global_id == -1) continue;
            if (r5_v <= r4_v) continue; // Keep triangle search

            int r5_fid = search_to_factor[r5_v];
            local_checked++;
            
            // 1. Edge compatibility check (Fastest)
            const Mask256_C& m5 = factor_edge_masks[r5_v];
            if ((m4.m[0] & m5.m[0]) || (m4.m[1] & m5.m[1])) {
                local_rej_edge++;
                continue;
            }

            // 2. Cycle structure check (SIMD)
            if (!is_perfect_scalar(global_pool[r4_fid].adj, global_pool[r5_fid].adj)) {
                local_rej_cycle++;
                continue;
            }

            // 3. Canonicity check
#if !K16_DISABLE_IS_CANONICAL_R45_CHECK
            if (!is_canonical(r4_fid, r5_fid, r4_trans.get())) {
                local_rej_canon++;
                continue;
            }
#endif
            local_kept++; 
            local_roots.push_back({ r4_v, r5_v });
        }
        
        g_total_pairs_checked += local_checked;
        g_rejected_edges += local_rej_edge;
        g_rejected_cycles += local_rej_cycle;
        g_rejected_canon += local_rej_canon;
        g_total_rejected += (local_rej_edge + local_rej_cycle + local_rej_canon);
        g_total_kept += local_kept;

        if (!local_roots.empty()) {
            std::lock_guard<std::mutex> lock(pr_mutex);
            potential_roots.insert(potential_roots.end(), local_roots.begin(), local_roots.end());
        }
        int done = ++canon_done;
        if (done % c_report_step == 0 || done == np) {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_canon_report).count() >= 5 || done == np) {
                last_canon_report = now;
                if (m_bPrint) { printf("%d ", (int)(done * 100LL / np)); }
            }
        }
    }, 128);
#if K16_USE_ROOT_SORT
    std::sort(potential_roots.begin(), potential_roots.end(), [](const RootPair& a, const RootPair& b) {
        if (a.r4_v != b.r4_v) return a.r4_v < b.r4_v;
        return a.r5_v < b.r5_v;
        });
#endif
    total_roots = (int)potential_roots.size();
    auto canon_end = std::chrono::steady_clock::now();
    auto canon_elap = std::chrono::duration_cast<std::chrono::seconds>(canon_end - canon_start).count();
    if (m_bPrint) {
        printf("100%% Done in %lld seconds.\n", canon_elap);
        printf("Total: %d, Kept: %d, Rejected: %d (Edges: %d, Cycles: %d, Canon: %d)\n", 
            (int)g_total_pairs_checked, (int)g_total_kept, (int)g_total_rejected,
            (int)g_rejected_edges, (int)g_rejected_cycles, (int)g_rejected_canon);
        printf("%d Root Pairs will be processed. (r4 branches rejected: %d) | StabMin/Max:%d/%d\n", total_roots, (int)r4_branches_rejected, (int)min_stab_size, (int)max_stab_size);
    }
    solve_start_time = std::chrono::steady_clock::now();
    last_print_time = solve_start_time;
    std::vector<std::pair<int, int>> r4_groups;
    if (!potential_roots.empty()) {
        int start = 0;
        for (int i = 1; i < (int)potential_roots.size(); i++) {
            if (potential_roots[i].r4_v != potential_roots[start].r4_v) {
                r4_groups.push_back({ start, i });
                start = i;
            }
        }
        r4_groups.push_back({ start, (int)potential_roots.size() });
    }

    solve_start_time = std::chrono::steady_clock::now();
    last_print_time = solve_start_time;
    if (m_bPrint) { printf("Solving with Group-by-r4 optimization. Groups: %zd\n", r4_groups.size()); }

    restart_index = 0;//3000000;//800;//2111111;
    if (restart_index)
        printf("[RESTART] Starting with restart_index = %d\n", restart_index);

    if (total_roots <= restart_index) {
        printf("*** number of roots(%d) <= restart index(%d), exit(1)\n", total_roots, restart_index);
        exit(1);
    }
    //total_roots -= restart_index;

    parallel_for(0, (int)r4_groups.size(), [&](int group_idx, int tid) {
        if (restart_index < r4_groups[group_idx].second) {
            const int start_i = r4_groups[group_idx].first < restart_index ? restart_index : r4_groups[group_idx].first;
        const int end_i = r4_groups[group_idx].second;
        const int r4_v = potential_roots[start_i].r4_v;
        const int r4_fid = search_to_factor[r4_v];
        const uint64_t* r4_bits = &adj_matrix[factor_offsets[r4_v]];
        const int r4_col_sw = row_ranges[1].start_word;

        ThreadLocalBuffers& buf = *thread_buffers[tid];
        if (buf.global_to_s4.empty()) buf.global_to_s4.assign(K16_M_MAX, -1);
        if (buf.g_to_l.empty()) buf.g_to_l.assign(K16_M_MAX, -1);
        
            // Defensive: ensure clean state when starting potentially in middle of a group
        buf.s4_to_global.clear();
        buf.active_words.clear();
        
        // 1. S4 = {f | f compatible with r4, row(f) >= 2}
        Mask256_C mask_r4;
        mask_r4.m[0] = fixedEdgesMask.m[0] | factor_edge_masks[r4_v].m[0];
        mask_r4.m[1] = fixedEdgesMask.m[1] | factor_edge_masks[r4_v].m[1];

        for (int r = 2; r < K16_SEARCH; r++) {
            while (buf.s4_to_global.size() % 64 != 0) buf.s4_to_global.push_back(-1); // Align row starts for bitwise logic
            const int sw = row_ranges[r].start_word;
            const int ew = row_ranges[r].end_word;
            for (int w = sw; w < ew; w++) {
                if (w < r4_col_sw) continue;
                    uint64_t word = r4_bits[w];
                if (word) {
                    uint64_t filtered_word = 0;
                    uint64_t temp_word = word;
                    while (temp_word) {
                        int bit = (int)_tzcnt_u64(temp_word);
                        int g_idx = w * 64 + bit;
                            if (g_idx < MS && items[g_idx].row_id == r && items[g_idx].global_id != -1) {
                            const Mask256_C& fm = factor_edge_masks[g_idx];
                            if (!((mask_r4.m[0] & fm.m[0]) || (mask_r4.m[1] & fm.m[1]))) {
                                filtered_word |= (1ULL << bit);
                            }
                        }
                        temp_word &= (temp_word - 1);
                    }
                    if (filtered_word) {
                        buf.active_words.push_back({ w, filtered_word, (int)buf.s4_to_global.size() });
                            uint64_t word_to_process = filtered_word;
                            while (word_to_process) {
                                int bit = (int)_tzcnt_u64(word_to_process);
                            int g_idx = w * 64 + bit;
                                buf.global_to_s4[g_idx] = (int)buf.s4_to_global.size();
                            buf.s4_to_global.push_back(g_idx);
                                word_to_process &= (word_to_process - 1);
                        }
                    }
                }
            }
        }

        const int S4K = (int)buf.s4_to_global.size();
        if (S4K == 0) {
                roots_done[tid] += (end_i - start_i);
                thread_root_idx[tid] = end_i - 1;
            return;
        }

        // 2. Build Semi-Local Matrix (S4 x S4) with fixed-width rows

        int s4_words = (S4K + 63) / 64 + 1; // PAD with 1 extra word for spill safety
        try {
            buf.s4_adj.assign((size_t)S4K * s4_words, 0);
            buf.s4_offsets.assign(S4K, 0);
            }
            catch (...) {
            printf("Error: Memory allocation failed for S4 matrix (S4K=%d, words=%d)\n", S4K, s4_words);
            exit(1);
        }
        for (int j = 0; j < S4K; j++) buf.s4_offsets[j] = (size_t)j * s4_words;

        int s4_l_cur = 0;
        for (int r = 2; r < K16_SEARCH; r++) {
            s4_l_cur = (s4_l_cur + 63) / 64 * 64;
            while (s4_l_cur < S4K && buf.s4_to_global[s4_l_cur] != -1 && items[buf.s4_to_global[s4_l_cur]].row_id == r) 
                s4_l_cur++;
        }

        for (int j = 0; j < S4K; j++) {
            int gj = buf.s4_to_global[j];
            if (gj == -1) continue; // Skip alignment padding
            int r_j = items[gj].row_id;
            uint64_t* s4_row = &buf.s4_adj[buf.s4_offsets[j]];
            const uint64_t* g_row = &adj_matrix[factor_offsets[gj]];
            const int g_sw = row_ranges[r_j + 1].start_word;
            
            for (const auto& aw : buf.active_words) {
                    uint64_t g_word = g_row[aw.global_word_idx];
                if (g_word & aw.mask) {
                    uint64_t extracted = _pext_u64(g_word, aw.mask);
                    int s4_bit_start = aw.local_bit_start;
                    s4_row[s4_bit_start / 64] |= (extracted << (s4_bit_start % 64));
                    if ((s4_bit_start % 64) + (int)__popcnt64(aw.mask) > 64) {
                        s4_row[s4_bit_start / 64 + 1] |= (extracted >> (64 - (s4_bit_start % 64)));
                    }
                }
            }
        }
        // 3. Process r5 mates
        for (int i = start_i; i < end_i; i++) {
                if (i < restart_index) { printf("internal error 3, exit(3)\n"); exit(3); }
            int r5_v = potential_roots[i].r5_v;
            int r5_fid = search_to_factor[r5_v];
            const uint64_t* r5_bits = &adj_matrix[factor_offsets[r5_v]];
            const int r5_col_sw = row_ranges[2].start_word; // Row 5 (Slot 1) has columns starting from Slot 2
            
            thread_row4[tid] = (r4_v - row_ranges[0].start_word * 64) + 1;
            thread_row5[tid] = (r5_v - row_ranges[1].start_word * 64) + 1;
            thread_root_idx[tid] = i;

            buf.local_ranges.assign(K16_SEARCH + 1, { 0, 0 });
            buf.local_to_global.clear();
            buf.local_s_to_f.clear();
            // Use heap-allocated buffers from ThreadLocalBuffers to avoid stack overflow
            SearchContext& local_ctx = buf.local_ctx;
            for (int r = 2; r < K16_SEARCH; r++) {
                for (int j = 0; j < buf.dirty_s4_count[r]; j++) buf.r5_row_mask_in_s4[r][buf.dirty_s4_words[r][j]] = 0;
                buf.dirty_s4_count[r] = 0;
                buf.sub_sampling_plan_sizes[r] = 0;
            }
            memset(buf.local_edge_presence, 0, sizeof(buf.local_edge_presence));
            int l_idx = 0;
            l_idx = 0; // Reset l_idx for the actual processing loop
            for (int gi : buf.local_to_global) if (gi != -1) buf.g_to_l[gi] = -1;
            buf.local_to_global.clear();
            buf.local_s_to_f.clear();

            const Mask256_C& m5 = factor_edge_masks[r5_v];

            for (int r = 2; r < K16_SEARCH; r++) {
                l_idx = (l_idx + 255) / 256 * 256;
                buf.local_ranges[r].start_word = l_idx / 64;
                const int g_sw = row_ranges[r].start_word;
                const int g_ew = row_ranges[r].end_word;

                for (int w = g_sw; w < g_ew; w++) {
                        uint64_t g_word = r5_bits[w];
                    while (g_word) {
                        int bit = (int)_tzcnt_u64(g_word);
                        int g_idx = w * 64 + bit;
                        if (g_idx < MS && items[g_idx].row_id == r && items[g_idx].global_id != -1) {
                            const Mask256_C& fm = factor_edge_masks[g_idx];
                            if (!((m5.m[0] & fm.m[0]) | (m5.m[1] & fm.m[1]))) {
                                int s4_idx = buf.global_to_s4[g_idx];
                                if (s4_idx != -1) {
                                    int w_idx = s4_idx / 64;
                                    if (buf.r5_row_mask_in_s4[r][w_idx] == 0) {
                                        if (buf.dirty_s4_count[r] < K16_WORDS) buf.dirty_s4_words[r][buf.dirty_s4_count[r]++] = w_idx;
                                    }
                                    buf.r5_row_mask_in_s4[r][w_idx] |= (1ULL << (s4_idx % 64));
                                    buf.g_to_l[g_idx] = l_idx;
                                    buf.local_to_global.push_back(g_idx);

                                    // Populate local edge presence
                                    uint64_t e0 = fm.m[0];
                                    uint64_t e1 = fm.m[1] & 0x00FFFFFFFFFFFFFFULL;
                                    while (e0) {
                                        int b = (int)_tzcnt_u64(e0);
                                            buf.local_edge_presence[b].bits[l_idx >> 6] |= (1ULL << (l_idx & 63));
                                        e0 &= (e0 - 1);
                                    }
                                    while (e1) {
                                        int b = (int)_tzcnt_u64(e1);
                                            buf.local_edge_presence[b + 64].bits[l_idx >> 6] |= (1ULL << (l_idx & 63));
                                        e1 &= (e1 - 1);
                                    }
                                        l_idx++;
                                }
                            }
                        }
                        g_word &= (g_word - 1);
                    }
                }
                buf.local_ranges[r].end_word = (l_idx + 63) / 64;
            }

            if (l_idx == 0) { roots_done[tid]++; continue; }

            buf.local_offsets.assign(l_idx, 0);
            buf.local_edge_masks.assign(l_idx, { 0, 0 });
            buf.local_s_to_f.assign(l_idx, -1);
            size_t l_coff = 0;
            int K = (int)buf.local_to_global.size();
            buf.local_compression = (double)Global_M_total / K;

            for (int j = 0; j < K; j++) {
                int gi = buf.local_to_global[j];
                int li = buf.g_to_l[gi];
                buf.local_edge_masks[li] = factor_edge_masks[gi];
                buf.local_s_to_f[li] = search_to_factor[gi];
                buf.local_offsets[li] = l_coff;
                int d = items[gi].row_id;
                    int wn = (buf.local_ranges[K16_SEARCH - 1].end_word + 3) / 4 * 4;
                    l_coff += (size_t)wn;

            }
            buf.local_adj.assign(l_coff + 16, 0);
            uint64_t* aligned_base = (uint64_t*)(((uintptr_t)buf.local_adj.data() + 63) & ~63);

            // Sub-sample S4_adj into local_adj
            // Generate sub-sampling plan once for this r5 mate
            for (int r = 2; r < K16_SEARCH; r++) {
                if (buf.dirty_s4_count[r] == 0) continue;
                std::sort(buf.dirty_s4_words[r], buf.dirty_s4_words[r] + buf.dirty_s4_count[r]);
                
                int l_bit_cursor = buf.local_ranges[r].start_word * 64;
                for (int j = 0; j < buf.dirty_s4_count[r]; j++) {
                    int w = buf.dirty_s4_words[r][j];
                    uint64_t mask = buf.r5_row_mask_in_s4[r][w];
                    SubSamplePlanItem& item = buf.sub_sampling_plans[r][buf.sub_sampling_plan_sizes[r]++];
                    item.s4_word_idx = (uint16_t)w;
                    item.mask = mask;
                    item.l_word_off = (uint16_t)(l_bit_cursor / 64);
                    item.l_bit_shift = (uint8_t)(l_bit_cursor % 64);
                    item.bits_pushed = (uint8_t)__popcnt64(mask);
                    item.split = (item.l_bit_shift + item.bits_pushed > 64) ? 1 : 0;
                    l_bit_cursor += item.bits_pushed;
                }
            }

            // Sub-sample S4_adj into local_adj using the plan
            for (int k_idx = 0; k_idx < K; k_idx++) {
                int gj = buf.local_to_global[k_idx];
                int li = buf.g_to_l[gj];
                int s4_j = buf.global_to_s4[gj];
                int d = items[gj].row_id;
                
                uint64_t* l_row = &aligned_base[buf.local_offsets[li]];
                const uint64_t* s4_row_ptr = &buf.s4_adj[buf.s4_offsets[s4_j]];
                const int l_row_sw = 0;
                const int l_row_ew = buf.local_ranges[K16_SEARCH - 1].end_word;
                const int l_row_wn = l_row_ew - l_row_sw;

                    for (int r = 2; r < K16_SEARCH; r++) {
                    int p_size = buf.sub_sampling_plan_sizes[r];
                    const auto* plan = buf.sub_sampling_plans[r];
                    for (int j = 0; j < p_size; j++) {
                        const auto& item = plan[j];
                        uint64_t extracted = _pext_u64(s4_row_ptr[item.s4_word_idx], item.mask);
                        int l_off = item.l_word_off - l_row_sw;
                        if (l_off >= 0 && l_off < l_row_wn) {
                            l_row[l_off] |= (extracted << item.l_bit_shift);
                            if (item.split && l_off + 1 < l_row_wn) {
                                l_row[l_off + 1] |= (extracted >> (64 - item.l_bit_shift));
                            }
                        }
                    }
                }
            }

            buf.ranges = buf.local_ranges.data();
            buf.adj = aligned_base;
            buf.offsets = buf.local_offsets.data();
            buf.edge_masks = buf.local_edge_masks.data();
            buf.s_to_f = buf.local_s_to_f.data();
            buf.local_compression = (double)Global_M_total / K;

                local_ctx.r4_idx = r4_v; local_ctx.r5_idx = r5_v; local_ctx.root_idx = i;
                memset(local_ctx.counts_valid, 0, sizeof(local_ctx.counts_valid));
                memset(local_ctx.mrv_counts, 0, sizeof(local_ctx.mrv_counts));

            memset(local_ctx.pool, 0, sizeof(local_ctx.pool));
                for (int j = 0; j < (int)buf.local_to_global.size(); j++) {
                    int gi = buf.local_to_global[j];
                    int li = buf.g_to_l[gi];
                    local_ctx.pool[2].bits[li / 64] |= (1ULL << (li % 64));
            }
            buf.cl.clear(); buf.cl.push_back(r4_fid); buf.cl.push_back(r5_fid);
            const Mask256_C& m4 = factor_edge_masks[r4_v];
            local_ctx.used_edges.m[0] = fixedEdgesMask.m[0] | m4.m[0] | m5.m[0];
            local_ctx.used_edges.m[1] = fixedEdgesMask.m[1] | m4.m[1] | m5.m[1];
                for (int s = 0; s < K16_SEARCH; s++) local_ctx.slots[s] = (uint8_t)s;

                local_ctx.r4_idx = r4_v; local_ctx.r5_idx = r5_v; local_ctx.root_idx = i;
                if (thread_root_idx[tid] < restart_index) { printf("internal error 4, exit(4)\n"); exit(4); }
                internal_solve(2, buf.cl, local_ctx, &buf);
            roots_done[tid]++;

            thread_root_idx[tid] = 0x7fffffff;
        }
        for (int gi : buf.s4_to_global) if (gi != -1) buf.global_to_s4[gi] = -1;
            // Also clear g_to_l mappings to be safe for next group/restart
            for (int gi : buf.local_to_global) if (gi != -1) buf.g_to_l[gi] = -1;
        }
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

void VECTOR_CALL K16P1F::internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx, ThreadLocalBuffers* pBuf) {
    ThreadLocalBuffers& buf = *pBuf;
    buf.meter.calls++;
#if K16_Use_rdtsc
    uint64_t t0 = __rdtsc();
#endif
    if (depth == 2 && m_bPrint) diagnostic_printout(buf.local_compression);
    
    if (depth == K16_SEARCH) {
        unsigned char results[K16_MATCH * K16_N];
        for (int i = 0; i < m_nFixedRows; i++) memcpy(results + i * K16_N, this->fixedRows[i].src, K16_N);
        std::vector<Factor> sol_factors;
        for (int fid : clique) sol_factors.push_back(this->global_pool[fid]);
        std::sort(sol_factors.begin(), sol_factors.end(), [](const Factor& a, const Factor& b) { return a.src[1] < b.src[1]; });
        for (int i = 0; i < K16_SEARCH; i++) memcpy(results + (i + m_nFixedRows) * K16_N, sol_factors[i].src, K16_N);
        
        {
            std::lock_guard<std::mutex> lock(this->result_mutex);
            // Call with mode=1 to check canonicity
            if (this->resultCallback(this->cbClass, results, ctx.r4_idx, ctx.r5_idx, 1)) {
                results_to_sort.push_back(std::vector<unsigned char>(results, results + K16_MATCH * K16_N));
                num_results++;
                if (0 && m_bPrint) {
                    printf("[CAN] Found canonical result at Root Index %d (r4:%d, r5:%d) | Total: %d\n",
                        ctx.root_idx, (int)ctx.r4_idx, (int)ctx.r5_idx, (int)num_results);
                }
            }
            else {
                num_notCanon++;
                if (0 && m_bPrint) {
                    printf("[NC] Found non-canonical result at Root Index %d (r4:%d, r5:%d) | Total NC: %d\n",
                        ctx.root_idx, (int)ctx.r4_idx, (int)ctx.r5_idx, (int)num_notCanon);
                }
            }
        }
        return;
    }

#if K16_USE_MRV
    if (depth >= 2 && depth < K16_MRV_DEPTH_LIMIT) {
        int best_idx = depth;
        int min_count = 1000000;
        
        if (ctx.counts_valid[depth]) {
            for (int i = depth; i < K16_SEARCH; i++) {
                int slot = ctx.slots[i];
                int count = ctx.mrv_counts[depth][slot];
                if (count < min_count) {
                    min_count = count;
                    best_idx = i;
                    if (min_count <= K16_MRV_EARLY_EXIT_THRESHOLD) break;
                }
            }
            ctx.counts_valid[depth] = false; // Consumption
        } else {
            for (int i = depth; i < K16_SEARCH; i++) {
                int slot = ctx.slots[i];
                int count = 0;
                const int wStart = buf.ranges[slot].start_word;
                const int wEnd = buf.ranges[slot].end_word;
                for (int w = wStart; w < wEnd; w++) {
                    count += (int)_mm_popcnt_u64(ctx.pool[depth].bits[w]);
                }
                if (count < min_count) {
                    min_count = count;
                    best_idx = i;
#if K16_MRV_EARLY_EXIT
                    if (min_count <= K16_MRV_EARLY_EXIT_THRESHOLD) break;
#endif
                }
            }
        }
        
        // MRV: Swap the most constrained slot to the front
        uint8_t tmp = ctx.slots[depth];
        ctx.slots[depth] = ctx.slots[best_idx];
        ctx.slots[best_idx] = tmp;
        
        if (min_count == 0) return;
    }
#endif

    int s = ctx.slots[depth];
    State& P = ctx.pool[depth];

    // --- High-Speed Bit-Sliced Edge Coverage Optimization ---
#if K16_EDGE_PRUNE_ENABLED
    if (depth >= K16_EDGE_PRUNE_START && depth <= K16_EDGE_PRUNE_END) {
        const int uStart = 0;
        const int uEnd = buf.ranges[K16_SEARCH - 1].end_word;

        // Ultra-Lazy Threshold: Zero overhead for main search, only prunes final branches
        const int threshold = K16_PRUNE_THRESHOLD_LOW;

        int total_remaining = 0;
        for (int w = uStart; w < uEnd; w++) {
            total_remaining += (int)__popcnt64(P.bits[w]);
            if (total_remaining > threshold) break; 
        }

        if (total_remaining <= threshold) {
            uint64_t u0 = 0, u1 = 0;
            const Mask256_C* masks = buf.edge_masks;
            uint64_t m0_needed = ~ctx.used_edges.m[0];
            uint64_t m1_needed = (~ctx.used_edges.m[1]) & 0x00FFFFFFFFFFFFFFULL;

            for (int w = uStart; w < uEnd; w++) {
                uint64_t word = P.bits[w];
                if (!word) continue;
                int base = w << 6;
                do {
                    int bit = (int)_tzcnt_u64(word);
                    const Mask256_C& fm = masks[base + bit];
                    u0 |= fm.m[0];
                    u1 |= fm.m[1];
                    // Early Exit: if all needed edges are now covered, we can stop
                    if (((u0 & m0_needed) == m0_needed) && ((u1 & m1_needed) == m1_needed)) break;
                    word &= (word - 1);
                } while (word);
                if (((u0 & m0_needed) == m0_needed) && ((u1 & m1_needed) == m1_needed)) break;
            }

            if ((m0_needed & ~u0) || (m1_needed & ~u1)) {
                return; 
            }

            // The Acceleration: Only if pool is very small
            if (total_remaining <= 50) {
                int forced_v = -1;
                bool conflict = false;

                for (int e = 0; e < 120; e++) {
                    bool is_needed = (e < 64) ? (m0_needed & (1ULL << e)) : (m1_needed & (1ULL << (e - 64)));
                    if (!is_needed) continue;

                    int total_e_count = 0;
                    int last_e_v = -1;
                    for (int w = uStart; w < uEnd; w += 4) {
                        __m256i p_vec = _mm256_load_si256((__m256i*) & P.bits[w]);
                        __m256i e_vec = _mm256_load_si256((__m256i*) & buf.local_edge_presence[e].bits[w]);
                        __m256i res_vec = _mm256_and_si256(p_vec, e_vec);
                        if (!_mm256_testz_si256(res_vec, res_vec)) {
                            uint64_t m0 = _mm256_extract_epi64(res_vec, 0);
                            uint64_t m1 = _mm256_extract_epi64(res_vec, 1);
                            uint64_t m2 = _mm256_extract_epi64(res_vec, 2);
                            uint64_t m3 = _mm256_extract_epi64(res_vec, 3);
                            total_e_count += (int)(__popcnt64(m0) + __popcnt64(m1) + __popcnt64(m2) + __popcnt64(m3));
                            if (total_e_count == 1) {
                                if (m0) last_e_v = w * 64 + (int)_tzcnt_u64(m0);
                                else if (m1) last_e_v = (w + 1) * 64 + (int)_tzcnt_u64(m1);
                                else if (m2) last_e_v = (w + 2) * 64 + (int)_tzcnt_u64(m2);
                                else if (m3) last_e_v = (w + 3) * 64 + (int)_tzcnt_u64(m3);
                            }
                            if (total_e_count > 1) break; 
                        }
                    }

                    if (total_e_count == 1) {
                        int r_s = ctx.slots[depth];
                        int r_start = buf.ranges[r_s].start_word;
                        int r_end = buf.ranges[r_s].end_word;
                        if (last_e_v / 64 >= r_start && last_e_v / 64 < r_end) {
                            if (forced_v != -1 && forced_v != last_e_v) { conflict = true; break; }
                            forced_v = last_e_v;
                        }
                    }
                }

                if (conflict) return;
                if (forced_v != -1) {
                    int v = forced_v;
                    const Mask256_C& fm = buf.edge_masks[v];
                    State& next_P = ctx.pool[depth + 1];
                    const uint64_t* rel_adj = &buf.adj[buf.offsets[v]];
                    
                    bool domain_wipeout = false;
                    for (int n = depth + 1; n < K16_SEARCH; n++) {
                        int next_s = ctx.slots[n];
                        __m256i row_acc = _mm256_setzero_si256();
                        const int wStart = buf.ranges[next_s].start_word;
                        const int wEnd = buf.ranges[next_s].end_word;
                        for (int w = wStart; w < wEnd; w += 4) {
                            __m256i p_vec = _mm256_load_si256((__m256i*) & P.bits[w]);
                            __m256i a_vec = _mm256_load_si256((__m256i*) & rel_adj[w]);
                            __m256i res_vec = _mm256_and_si256(p_vec, a_vec);
                            _mm256_store_si256((__m256i*) & next_P.bits[w], res_vec);
                            row_acc = _mm256_or_si256(row_acc, res_vec);
                        }
                        if (_mm256_testz_si256(row_acc, row_acc)) { domain_wipeout = true; break; }
                    }
                    if (!domain_wipeout) {
                        Mask256_C old = ctx.used_edges;
                        ctx.used_edges.m[0] |= fm.m[0]; ctx.used_edges.m[1] |= fm.m[1];
                        clique.push_back(buf.s_to_f[v]);
                        internal_solve(depth + 1, clique, ctx, &buf);
                        clique.pop_back();
                        ctx.used_edges = old;
                    }
                    return;
                }
            }
        }
    }
#endif

    int start = buf.ranges[s].start_word;
    int end = buf.ranges[s].end_word;
    // 40% cpu below
    for (int i = start; i < end; i += 4) {
        __m256i combined = _mm256_load_si256((__m256i*) & P.bits[i]);
        if (_mm256_testz_si256(combined, combined)) continue;

        const int sub_end = (end - i >= 4) ? 4 : (end - i);
        for (int sub = 0; sub < sub_end; ++sub) {
            uint64_t word = ((uint64_t*)&combined)[sub];
            while (word) {
                int v = (i + sub) * 64 + (int)_tzcnt_u64(word);

                // No inline edge check needed! Baked into adjacency and pool filtering.
                if (depth + 1 < K16_SEARCH) {
                    State& next_P = ctx.pool[depth + 1];
                    const uint64_t* rel_adj = &buf.adj[buf.offsets[v]];
                    
                    bool domain_wipeout = false;

                    for (int n = depth + 1; n < K16_SEARCH; n++) {
                        int next_s = ctx.slots[n];
                        int slot_count = 0;
                        const int wStart = buf.ranges[next_s].start_word;
                        const int wEnd = buf.ranges[next_s].end_word;
                        for (int w = wStart; w < wEnd; w += 4) {
                            __m256i p_vec = _mm256_load_si256((__m256i*) & P.bits[w]);
                            __m256i a_vec = _mm256_load_si256((__m256i*) & rel_adj[w]);
                            __m256i res_vec = _mm256_and_si256(p_vec, a_vec);
                            _mm256_store_si256((__m256i*) & next_P.bits[w], res_vec);
                            
                            uint64_t m0 = _mm256_extract_epi64(res_vec, 0);
                            uint64_t m1 = _mm256_extract_epi64(res_vec, 1);
                            uint64_t m2 = _mm256_extract_epi64(res_vec, 2);
                            uint64_t m3 = _mm256_extract_epi64(res_vec, 3);
                            slot_count += (int)(__popcnt64(m0) + __popcnt64(m1) + __popcnt64(m2) + __popcnt64(m3));
                        }
                        ctx.mrv_counts[depth + 1][next_s] = (uint16_t)slot_count;
                        if (slot_count == 0) { domain_wipeout = true; break; }
                    }
                    if (domain_wipeout) { word &= (word - 1); continue; }
                    ctx.counts_valid[depth + 1] = true;
                }

                // --- NEW: Odd Component Pruning ---
                const Mask256_C& fm = buf.edge_masks[v];

                // Safety check: ensure no edge conflict with the prefix
                if ((ctx.used_edges.m[0] & fm.m[0]) | (ctx.used_edges.m[1] & fm.m[1])) {
                    word &= (word - 1); continue;
                }

                Mask256_C old = ctx.used_edges;
                ctx.used_edges.m[0] |= fm.m[0]; ctx.used_edges.m[1] |= fm.m[1];
                clique.push_back(buf.s_to_f[v]);
                internal_solve(depth + 1, clique, ctx, &buf);
                clique.pop_back();
                ctx.used_edges = old;
                word &= (word - 1);
            }
        }
    }
}


K16P1F::PackedAdj K16P1F::pack_factor_adj(const uint8_t* adj) {
    PackedAdj res; memset(&res, 0, sizeof(res));
    memcpy(res.adj, adj, K16_N);
    for (int i = 0; i < K16_N; i++) {
        int v = adj[i];
        if (i < v) {
            int eid = edge_id_table[i][v];
            if (eid != -1) {
                if (eid < 64) res.edge_mask.m[0] |= (1ULL << eid);
                else if (eid < 128) res.edge_mask.m[1] |= (1ULL << (eid - 64));
            }
        }
    }
    return res;
}

bool K16P1F::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    // Fast scalar Hamiltonian cycle check for K16
    uint32_t visited = (1 << 0);
    int count = 1, p = 0, c = adj1[0];
    while (c != 0) {
        visited |= (1 << c);
        int next = (adj1[c] == p) ? adj2[c] : adj1[c];
        p = c; c = next; count++;
        if (count > 16) return false;
    }
    return count == 16;
}

K16P1F::CycleUnion K16P1F::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    uint32_t visited = 0;
    for (int start = 0; start < K16_N; start++) {
        if (visited & (1 << start)) continue;
        int curr = start, prev = -1, len = 0;
        do {
            visited |= (1 << curr);
            cu.cycles[cu.count][len++] = curr;
            int next = (adj1[curr] == prev) ? adj2[curr] : adj1[curr];
            prev = curr; curr = next;
        } while (curr != start);
        cu.lens[cu.count++] = len;
        if (cu.count >= 16) break;
    }
    return cu;
}

void K16P1F::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
    info.count = 0;
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
            p_total.pv = _mm256_loadu_si256((const __m256i*)p_total.p);
            p_total.pinvv = _mm256_loadu_si256((const __m256i*)p_total.p_inv);
        }
    }
    return;
}

void K16P1F::apply_perm_16(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
    // dst[perm.p[i]] = perm.p[src[i]]
    // For N=16, we can use _mm_shuffle_epi8 directly if we are careful.
    __m128i v_src = _mm_loadu_si128((const __m128i*)src_adj);
    __m128i v_p = _mm256_castsi256_si128(perm.pv);
    __m128i v_pinv = _mm256_castsi256_si128(perm.pinvv);
    
    // dst[j] = p[src[pinv[j]]]
    __m128i v_reordered = _mm_shuffle_epi8(v_src, v_pinv);
    __m128i v_final = _mm_shuffle_epi8(v_p, v_reordered);
    _mm_storeu_si128((__m128i*)dst_adj, v_final);
}

K16P1F::FastSortedFactor K16P1F::get_fast_sorted(const uint8_t* adj) {
    alignas(32) uint8_t indices[16] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
    };
    __m128i v_indices = _mm_loadu_si128((const __m128i*)indices);
    __m128i v_adj = _mm_loadu_si128((const __m128i*)adj);
    uint32_t mask = _mm_movemask_epi8(_mm_cmpgt_epi8(v_adj, v_indices)) & 0xFFFF;

    FastSortedFactor sf; memset(&sf, 0, sizeof(sf));
    uint8_t* out = sf.pairs;

    auto process_chunk = [&](uint32_t m, const uint8_t* adj_ptr, int offset) {
        __m128i v_adj_chunk = _mm_loadl_epi64((const __m128i*)adj_ptr);
        __m128i v_table_adj = _mm_load_si128(&gfs_table_adj[m]);
        __m128i v_table_idx = _mm_load_si128(&gfs_table_idx[m]);
        __m128i p_adj = _mm_shuffle_epi8(v_adj_chunk, v_table_adj);
        __m128i p_idx = _mm_add_epi8(v_table_idx, _mm_set1_epi8((char)offset));
        __m128i pairs = _mm_unpacklo_epi8(p_idx, p_adj);
        _mm_storeu_si128((__m128i*)out, pairs);
        out += __popcnt(m) * 2;
    };

    process_chunk(mask & 0xFF, adj, 0);
    process_chunk((mask >> 8) & 0xFF, adj + 8, 8);

    return sf;
}
bool K16P1F::is_canonical(int r4_fid, int r5_fid, const K16P1F::TransInfo* r4_dependent_trans) {
    const uint8_t* r_adj[5] = {
        fixedRows[0].adj, fixedRows[1].adj, fixedRows[2].adj,
        global_pool[r4_fid].adj, global_pool[r5_fid].adj
    };

    K16P1F::FastRowTriplet target_triplet;
    target_triplet.r[0] = fixedRows[2].fs;
    target_triplet.r[1] = global_pool[r4_fid].fs;
    target_triplet.r[2] = global_pool[r5_fid].fs;
    std::sort(target_triplet.r, target_triplet.r + 3, [](const K16P1F::FastSortedFactor& a, const K16P1F::FastSortedFactor& b) {
        return K16P1F::compare_fast_sorted(a, b);
    });

    auto check_trans = [&](const K16P1F::TransInfo& info, int i, int j) {
        for (int k = 0; k < info.count; ++k) {
            const auto& perm = info.perms[k];
            uint8_t tr_fi[16], tr_fj[16];
            apply_perm_16(r_adj[i], perm, tr_fi); 
            apply_perm_16(r_adj[j], perm, tr_fj);

            FastSortedFactor fs_i = get_fast_sorted(tr_fi);
            FastSortedFactor fs_j = get_fast_sorted(tr_fj);

            bool fi_r1 = K16P1F::equal_fast_sorted(fs_i, r1_can);
            bool fi_r2 = K16P1F::equal_fast_sorted(fs_i, r2_can);
            bool fj_r1 = K16P1F::equal_fast_sorted(fs_j, r1_can);
            bool fj_r2 = K16P1F::equal_fast_sorted(fs_j, r2_can);

            if (!((fi_r1 && fj_r2) || (fi_r2 && fj_r1))) continue;

            K16P1F::FastSortedFactor nt_r[3]; int nt_count = 0;
            for (int m = 0; m < 5; ++m) {
                if (m == i || m == j) continue;
                uint8_t tr_m[16]; apply_perm_16(r_adj[m], perm, tr_m);
                nt_r[nt_count++] = get_fast_sorted(tr_m);
            }
            
            if (K16P1F::compare_fast_sorted(nt_r[1], nt_r[0])) std::swap(nt_r[0], nt_r[1]);
            if (K16P1F::compare_fast_sorted(nt_r[2], nt_r[1])) {
                std::swap(nt_r[1], nt_r[2]);
                if (K16P1F::compare_fast_sorted(nt_r[1], nt_r[0])) std::swap(nt_r[0], nt_r[1]);
            }

            for (int m = 0; m < 3; m++) {
                if (K16P1F::compare_fast_sorted(nt_r[m], target_triplet.r[m])) return false;
                if (K16P1F::compare_fast_sorted(target_triplet.r[m], nt_r[m])) break;
            }
        }
        return true;
    };

    static const int pairs[3][2] = { {0,1}, {0,2}, {1,2} };
    for (int idx = 0; idx < 3; ++idx) {
        if (!check_trans(fixed_trans[idx], pairs[idx][0], pairs[idx][1])) return false;
    }
    for (int idx = 0; idx < 3; ++idx) {
        if (!check_trans(r4_dependent_trans[idx], idx, 3)) return false;
    }
    TransInfo infoBuffer;
    for (int i = 0; i < 4; ++i) {
        get_transformations((i < 3 ? fixedRows[i] : global_pool[r4_fid]), global_pool[r5_fid], infoBuffer);
        if (!check_trans(infoBuffer, i, 4)) return false;
    }
    return true;
}

bool K16P1F::is_canonical_stab(int r5_fid, const Permutation* stab, int stab_count) {
    const Factor& f5 = global_pool[r5_fid];

    for (int k = 0; k < stab_count; k++) {
        const auto& perm = stab[k];
        
        // Apply automorphism of {r1, r2, r3, r4} to r5
        uint8_t tr_r5[16];
        apply_perm_16(f5.adj, perm, tr_r5);
        FastSortedFactor fs_r5 = get_fast_sorted(tr_r5);
        
        // If the image of r5 is smaller than original r5, it's not canonical.
        if (compare_fast_sorted(fs_r5, f5.fs)) return false;
    }
    return true;
}

void K16P1F::diagnostic_printout(double current_compr) {
    if (!bTimeSet) {
        std::lock_guard<std::mutex> lock(result_mutex);
        solve_start_time = std::chrono::steady_clock::now();
        last_print_time = solve_start_time;
        bTimeSet = true;
        return;
    }
    if (call_counter++ < 20) return; call_counter = 0; // reduce calls to timer
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::seconds>(now - last_print_time).count() >= 30) {
        std::lock_guard<std::mutex> lock(result_mutex);
    last_print_time = now;
        auto elap = std::chrono::duration_cast<std::chrono::seconds>(now - solve_start_time).count();

        static int print_count = 0;
        print_count++;

        double pct = 0.0, pct_global = 0.0;
        int expected_run_time = 0;
        int total_done = 0;
        for (int t = 0; t < kThreads; t++) total_done += roots_done[t];

        if (total_roots > restart_index) {
            pct = (100.0 * total_done) / (double)(total_roots - restart_index);
            pct_global = (100.0 * (total_done + restart_index)) / (double)total_roots;
            if (pct > 0.0001) expected_run_time = (int)(elap * 100.0 / pct + 0.5);
        }

        int min_idx = 0x7fffffff;
        for (int t = 0; t < kThreads; t++) {
            if (thread_root_idx[t] < restart_index) { printf("internal error 6, exit(6)\n"); exit(6); }
            if (thread_root_idx[t] < min_idx) min_idx = thread_root_idx[t];
        }

        uint64_t sum_total = 0, sum_prop = 0, sum_mask = 0, sum_calls = 0, sum_pm = 0, sum_pp = 0;
    for (int t = 0; t < kThreads; t++) {
        sum_total += thread_buffers[t]->meter.total_cycles;
        sum_prop += thread_buffers[t]->meter.prop_cycles;
        sum_mask += thread_buffers[t]->meter.mask_cycles;
            sum_calls += thread_buffers[t]->meter.calls;
        sum_pm += thread_buffers[t]->meter.masks_checked;
        sum_pp += thread_buffers[t]->meter.props_done;
    }
    
        uint64_t avg_t = sum_calls ? sum_total / sum_calls : 0;
        uint64_t avg_p = sum_pp ? sum_prop / sum_pp : 0;
        uint64_t avg_m = sum_pm ? sum_mask / sum_pm : 0;

        printf("[RUN] M/Restart: %04d/%5d | ", fixed3RowsIndex, min_idx);
        printf("%%:%6.3f | Time/Exp:%4d/%4d | Res(C/N):%d/%d ",
            pct_global, (int)elap, expected_run_time, (int)num_results.load(), (int)num_notCanon.load());
        printf("| Compr:%.1f | Threads:%d\n", current_compr, kThreads);
#if K16_Use_rdtsc
        const char* status = (avg_t < 5000) ? "ok" : (avg_t < 10000 ? "warn" : "slow");
        printf("[RDTSC] AvgCyc Total:%llu (%s) | Prop:%llu | Mask:%llu\n", avg_t, status, avg_p, avg_m);
#endif

#if K16_BENCHMARK_EXIT
        if (elap >= 100) {
            printf("[RUN] Reached 100 seconds. Terminating for benchmarking...\n");
        exit(0);
    }
#endif
    }
}

bool K16P1F::compare_triplets(const K16P1F::FastRowTriplet& a, const K16P1F::FastRowTriplet& b) {
    for (int i = 0; i < 3; i++) {
        if (!K16P1F::equal_fast_sorted(a.r[i], b.r[i])) return K16P1F::compare_fast_sorted(a.r[i], b.r[i]);
    }
    return false;
}

