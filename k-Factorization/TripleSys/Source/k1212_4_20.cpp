#include "k1212_4_20.h"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <format>
#include <iterator>
#include <utility>

K1212_4_20::K1212_4_20(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr, bool bPrint) : KBipartBase(factParam, bPrint) {
    int id_cnt = 0;
    diag_cnt = 0;
    memset(edge_id_table, -1, sizeof(edge_id_table));
    for (int i = 0; i < m_nPlayers; i += 2) { // 0, 2, ... (Even nodes)
        for (int j = 1; j < m_nPlayers; j += 2) { // 1, 3, ... (Odd nodes)
            edge_id_table[i][j] = id_cnt++;
            edge_id_table[j][i] = edge_id_table[i][j];
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

void K1212_4_20::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K(12,12) 4:24, MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
    }
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    start_time = std::chrono::steady_clock::now();
    last_print_time = start_time;

    fixedEdgesMask = { 0, 0, 0 };
    for (int s = 0; s < K1212_FIXED; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * m_nPlayers, m_nPlayers);
        for (int i = 0; i < m_nPlayers; i += 2) {
            const auto a = fixedRows[s].src[i]; const auto b = fixedRows[s].src[i + 1];
            fixedRows[s].adj[a] = b; fixedRows[s].adj[b] = a;
            int eid = edge_id_table[a][b];
            if (eid < 64) fixedEdgesMask.m[0] |= (1ULL << eid);
            else if (eid < 128) fixedEdgesMask.m[1] |= (1ULL << (eid - 64));
            else fixedEdgesMask.m[2] |= (1ULL << (eid - 128));
        }
    }

    target_cu = find_cycles(fixedRows[0].adj, fixedRows[1].adj);
    for (int s = 0; s < K1212_FIXED; s++) fixedRows[s].fs = get_fast_sorted(fixedRows[s].adj);
    r1_can = fixedRows[0].fs;
    r2_can = fixedRows[1].fs;

#if K1212_USE_SIMD_CYCLE_CHECK
    for (int i = 0; i < K1212_FIXED; i++) {
        alignas(16) uint8_t e2o[16] = { 0 };
        alignas(16) uint8_t o2e[16] = { 0 };
        for (int j = 0; j < 12; j++) {
            e2o[j] = (fixedRows[i].adj[j * 2] - 1) / 2;
            o2e[j] = fixedRows[i].adj[j * 2 + 1] / 2;
        }
        fixedRows[i].v_e2o = _mm_loadu_si128((const __m128i*)e2o);
        fixedRows[i].v_o2e = _mm_loadu_si128((const __m128i*)o2e);
        fixed_packed[i] = pack_factor_adj(fixedRows[i].adj);
    }
#else
    for (int i = 0; i < K1212_FIXED; i++) {
        fixed_packed[i] = pack_factor_adj(fixedRows[i].adj);
    }
#endif

    get_transformations(fixedRows[0], fixedRows[1], fixed_trans[0]);
    get_transformations(fixedRows[0], fixedRows[2], fixed_trans[1]);
    get_transformations(fixedRows[1], fixedRows[2], fixed_trans[2]);

    // Reset transient state for reuse
    global_pool.clear();
    packed_pool.clear();
    f_map.clear();
    memset(roots_done, 0, sizeof(roots_done));
    KBase::init();
}

bool K1212_4_20::addRow(int iRow, const unsigned char* source) {
    iRow /= 2;
    if (iRow < m_nFixedRows || iRow >= m_nMatched) return false;
    int slot = iRow - m_nFixedRows;

    Factor f; memset(&f, 0, sizeof(f));
    memcpy(f.src, source, m_nPlayers);
    for (int i = 0; i < m_nPlayers; i += 2) { 
        const auto a = f.src[i]; const auto b = f.src[i + 1]; 
        if (a >= m_nPlayers || b >= m_nPlayers) return false; // Safety check
        f.adj[a] = b; f.adj[b] = a; 
    }

    f.edge_mask = { 0, 0, 0 };
    for (int i = 0; i < m_nPlayers; i += 2) {
        int eid = edge_id_table[f.src[i]][f.src[i+1]];
        if (eid < 64) f.edge_mask.m[0] |= (1ULL << eid);
        else if (eid < 128) f.edge_mask.m[1] |= (1ULL << (eid - 64));
        else f.edge_mask.m[2] |= (1ULL << (eid - 128));
    }

    if ((f.edge_mask.m[0] & fixedEdgesMask.m[0]) || (f.edge_mask.m[1] & fixedEdgesMask.m[1]) || (f.edge_mask.m[2] & fixedEdgesMask.m[2])) {
        printf("\n*** Error: Row has common edges with fixed rows!\n");
        exit(1);
    }

    // Check cycle structure with fixed rows (r1, r2, r3)
    PackedAdj pf = pack_factor_adj(f.adj);
    for (int i = 0; i < K1212_FIXED; i++) {
        if (!is_perfect_packed(pf, fixed_packed[i])) {
            printf("\n*** Error: Row has incorrect cycle structure with fixed row %d!\n", i + 1);
            exit(1);
        }
    }

    f.fs = get_fast_sorted(f.adj);

    std::vector<uint8_t> k(f.adj, f.adj + m_nPlayers);
    int factor_id;
    auto it = f_map.find(k);
    if (it == f_map.end()) {
        factor_id = (int)global_pool.size();
        if (factor_id >= K1212_M_MAX) return false;
        f_map[k] = factor_id;
        global_pool.push_back(f);
        packed_pool.push_back(pack_factor_adj(f.adj));
    } else {
        factor_id = it->second;
    }

    addFactorToSlot(slot, factor_id);
    return true;
}

void K1212_4_20::parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size) {
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

void K1212_4_20::solve(int mode) {
    int Global_M_total = (int)global_pool.size();
    if (m_bPrint) {
        printf("Row slots size: ");
        for (int s = 0; s < K1212_SEARCH; s++) printf("%d:%zd ", s * 2 + 7, temp_slot_ids[s].size());
        printf("\n");
        printf("Iterative Global Filtering (Arc Consistency):\n");
        printf("Total candidates: %d ", Global_M_total);
    }
    std::vector<uint8_t> active(Global_M_total, 1);
    std::vector<uint32_t> f_rows_mask(Global_M_total, 0);
    for (int s = 0; s < K1212_SEARCH; s++) for (int id : temp_slot_ids[s]) f_rows_mask[id] |= (1 << s);

    std::atomic<bool> changed{ true };
    int iter = 0;
    while (changed) {
        changed = false;
        parallel_for(0, Global_M_total, [&](int i, int tid) {
            if (!active[i]) return;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < K1212_SEARCH; r++) {
                if (f_rows_mask[i] & (1 << r)) continue; 
                bool row_found = false;
                const __m128i v_e2o_i = _mm_load_si128(&pi.v_e2o);
                const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
                for (int id : temp_slot_ids[r]) {
                    if (active[id] && check_cycle_simd(v_e2o_i, _mm_load_si128(&packed_pool[id].v_o2e), v_id)) { row_found = true; break; }
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
    }
    if (m_bPrint) { printf("\n"); }

    std::vector<int> old_to_new(Global_M_total, -1);
    std::vector<Factor> f_pool_new; std::vector<PackedAdj> p_pool_new;
    for (int i = 0; i < Global_M_total; i++) if (active[i]) { old_to_new[i] = (int)f_pool_new.size(); f_pool_new.push_back(global_pool[i]); p_pool_new.push_back(packed_pool[i]); }
    global_pool = std::move(f_pool_new); packed_pool = std::move(p_pool_new);
    
    std::vector<std::vector<int>> new_slot_ids(K1212_SEARCH);
    if (m_bPrint) printf("Row: Candidates ");
    for (int s = 0; s < K1212_SEARCH; s++) {
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
        const __m128i v_e2o_i = _mm_load_si128(&pi.v_e2o);
        const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        for (int r = 0; r < K1212_SEARCH; r++) for (int id : new_slot_ids[r]) if (check_cycle_simd(v_e2o_i, _mm_load_si128(&packed_pool[id].v_o2e), v_id)) deg[i]++;
        
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
    for (int s = 0; s < K1212_SEARCH; ++s) {
        std::vector<SortItem> row_items;
        for (int id : new_slot_ids[s]) row_items.push_back({ id, s, deg[id] });
        /*
        std::sort(row_items.begin(), row_items.end(), [](const SortItem& a, const SortItem& b) {
        if (K1212_SORT_INPUT == 1) return a.degree < b.degree;
        if (K1212_SORT_INPUT == 2) return a.degree > b.degree;
            return false; 
        }); */
        while (items_sorted.size() % 256 != 0) items_sorted.push_back({ -1, s, 0 }); 
        items_sorted.insert(items_sorted.end(), row_items.begin(), row_items.end());
    }
    const std::vector<SortItem>& items = items_sorted;

    int MS = (int)items.size(); 
    search_to_factor.assign(MS, -1); 
    for (int i = 0; i < MS; ++i) if (items[i].global_id != -1) search_to_factor[i] = items[i].global_id;
    
    factor_edge_masks.assign(MS, { 0, 0, 0 });
    for (int i = 0; i < MS; i++) {
        if (items[i].global_id == -1) continue;
        const Factor& f = global_pool[items[i].global_id];
        for (int u = 0; u < m_nPlayers; u += 2) { 
            int v = f.adj[u];
            int id = edge_id_table[u][v];
            if (id != -1) {
            if (id < 64) factor_edge_masks[i].m[0] |= (1ULL << id);
            else if (id < 128) factor_edge_masks[i].m[1] |= (1ULL << (id - 64));
            else factor_edge_masks[i].m[2] |= (1ULL << (id - 128));
        }
    }
    }

    row_ranges.assign(K1212_SEARCH, {0, 0});
    int cp = 0;
    for (int s = 0; s < K1212_SEARCH; ++s) {
        row_ranges[s].start_word = cp / 64;
        int count = 0; for (const auto& it : items) if (it.row_id == s) count++;
        cp += count; row_ranges[s].end_word = (cp + 63) / 64;
    }

    factor_offsets.assign(MS, 0); size_t coff = 0;
    for (int i = 0; i < MS; i++) {
        factor_offsets[i] = coff; int d = items[i].row_id;
        if (d + 1 < K1212_SEARCH) { int wn = (row_ranges[K1212_SEARCH - 1].end_word - row_ranges[d + 1].start_word + 3) / 4 * 4; coff += wn; }
    }
    adj_matrix.assign(coff, 0);
    if (m_bPrint) { printf("Populating adj_matrix (%%): "); }
    auto last_matrix_report = std::chrono::steady_clock::now();

    // Create sorted pools for contiguous AVX access
    std::vector<PackedAdj> sorted_packed(MS);
    for (int i = 0; i < MS; i++) {
        if (items[i].global_id != -1) sorted_packed[i] = packed_pool[items[i].global_id];
        else memset(&sorted_packed[i], 0, sizeof(PackedAdj));
    }

    std::atomic<int> matrix_done{ 0 };
    int report_step = MS / 20; if (report_step < 1) report_step = 1;

    parallel_for(0, MS, [&](int i, int tid) {
        if (items[i].global_id == -1) return;
        int d = items[i].row_id; if (d + 1 >= K1212_SEARCH) return;
        const PackedAdj& pi = sorted_packed[i];
        const uint8_t* adj_i = global_pool[items[i].global_id].adj;
        __m256i mi = _mm256_load_si256((const __m256i*)&pi);
        
        uint64_t* row_i = &adj_matrix[factor_offsets[i]];
        int sw = row_ranges[d + 1].start_word;

        int j_start = row_ranges[d + 1].start_word * 64;
        if (j_start < i + 1) j_start = i + 1;

        const __m128i v_e2o_i = _mm_load_si128(&pi.v_e2o);
        const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        for (int j = j_start; j < MS; ++j) {
            if (items[j].global_id == -1) continue;
            const PackedAdj& pj = sorted_packed[j];
            __m256i mj = _mm256_load_si256((const __m256i*)&pj);
            if (_mm256_testz_si256(mi, mj)) {
                if (check_cycle_simd(v_e2o_i, _mm_load_si128(&pj.v_o2e), v_id)) {
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
        const Mask256& m4 = factor_edge_masks[r4_v];

        std::vector<RootPair> local_roots;
        int local_checked = 0, local_kept = 0;
        int local_rej_edge = 0, local_rej_cycle = 0, local_rej_canon = 0;

#if !K1212_DISABLE_IS_CANONICAL_R45_CHECK
#if K1212_USE_STABILIZER_CANON
        // 1. New: Pre-calculate 4-row stabilizer (Automorphisms of prefix r1..r4)
        auto stabilizer_buf = std::make_unique<Permutation24[]>(512); 
        int stab_count = 0;
        
        // Find automorphisms (non-trivial only as requested)
        static const int pairs4[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
        const Factor* r_factors4[4] = { &fixedRows[0], &fixedRows[1], &fixedRows[2], &global_pool[r4_fid] };

        bool r4_canonical = true;
        TransInfo info;
        for (int p_idx = 0; p_idx < 6; p_idx++) {
            get_transformations(*r_factors4[pairs4[p_idx][0]], *r_factors4[pairs4[p_idx][1]], info);
            for (int k = 0; k < info.count; k++) {
                const auto& perm = info.perms[k];
                
                // Image of the 4-row set
                FastSortedFactor img[4];
                for (int m = 0; m < 4; m++) {
                    uint8_t tr_m[24];
                    apply_perm_24(r_factors4[m]->adj, perm, tr_m);
                    img[m] = get_fast_sorted(tr_m);
                }
                std::sort(img, img + 4, compare_fast_sorted);
                
                bool smaller = false, identical = true;
                for (int m = 0; m < 4; m++) {
                    if (compare_fast_sorted(img[m], r_factors4[m]->fs)) { smaller = true; identical = false; break; }
                    if (compare_fast_sorted(r_factors4[m]->fs, img[m])) { identical = false; break; }
                }
                
                if (smaller) { r4_canonical = false; break; }
                if (identical) {
                    // Check if non-trivial
                    bool trivial = true;
                    for (int n = 0; n < 24; n++) if (perm.p[n] != n) { trivial = false; break; }
                    if (!trivial && stab_count < 512) stabilizer_buf[stab_count++] = perm;
                }
            }
            if (!r4_canonical) break;
        }

        if (!r4_canonical) {
            int branch_size = 0;
            for (int r5_v = row_ranges[1].start_word * 64; r5_v < row_ranges[1].end_word * 64; r5_v++) {
                if (r5_v >= MS || items[r5_v].row_id != 1 || items[r5_v].global_id == -1) continue;
                if (r5_v <= r4_v) continue;
                branch_size++;
            }
            g_total_pairs_checked += branch_size;
            g_rejected_canon += branch_size;
            g_total_rejected += branch_size;
            canon_done++;
            return;
        }

        // Update stabilizer stats
        int cur_min = min_stab_size.load();
        while (stab_count < cur_min && !min_stab_size.compare_exchange_weak(cur_min, stab_count));
        int cur_max = max_stab_size.load();
        while (stab_count > cur_max && !max_stab_size.compare_exchange_weak(cur_max, stab_count));
#endif
#endif

        const __m128i v_e2o_4 = _mm_load_si128(&p4.v_e2o);
        const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
        for (int r5_v = row_ranges[1].start_word * 64; r5_v < row_ranges[1].end_word * 64; r5_v++) {
            if (r5_v >= MS || items[r5_v].row_id != 1 || items[r5_v].global_id == -1) continue;
            if (r5_v <= r4_v) continue; // Keep triangle search

            int r5_fid = search_to_factor[r5_v];
            local_checked++;
            
            // 1. Edge compatibility check (Fastest)
            const Mask256& m5 = factor_edge_masks[r5_v];
            if ((m4.m[0] & m5.m[0]) || (m4.m[1] & m5.m[1]) || (m4.m[2] & m5.m[2])) {
                local_rej_edge++;
                continue;
            }

            // 2. Cycle structure check (SIMD)
            if (!is_perfect_scalar(global_pool[r4_fid].adj, global_pool[r5_fid].adj)) {
                local_rej_cycle++;
                continue;
            }

            // 3. Canonicity check
#if !K1212_DISABLE_IS_CANONICAL_R45_CHECK
#if K1212_USE_STABILIZER_CANON
            if (!is_canonical_stab(r5_fid, stabilizer_buf.get(), stab_count)) {
                local_rej_canon++;
                continue;
            }
#else
            if (!is_canonical(r4_fid, r5_fid, r4_trans.get())) {
                local_rej_canon++;
                continue;
            }
#endif
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
#if K1212_USE_ROOT_SORT
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

    parallel_for(0, (int)r4_groups.size(), [&](int group_idx, int tid) {

        const int start_i = r4_groups[group_idx].first;
        const int end_i = r4_groups[group_idx].second;
        const int r4_v = potential_roots[start_i].r4_v;
        const int r4_fid = search_to_factor[r4_v];
        const uint64_t* r4_bits = &adj_matrix[factor_offsets[r4_v]];
        const int r4_col_sw = row_ranges[1].start_word;

        ThreadLocalBuffers& buf = *thread_buffers[tid];
        if (buf.global_to_s4.empty()) buf.global_to_s4.assign(K1212_M_MAX, -1);
        if (buf.g_to_l.empty()) buf.g_to_l.assign(K1212_M_MAX, -1);
        
        buf.s4_to_global.clear();
        buf.active_words.clear(); 
        
        // 1. S4 = {f | f compatible with r4, row(f) >= 2}

        for (int r = 2; r < K1212_SEARCH; r++) {
            while (buf.s4_to_global.size() % 64 != 0) buf.s4_to_global.push_back(-1); // Align row starts for bitwise logic
            const int sw = row_ranges[r].start_word;
            const int ew = row_ranges[r].end_word;
            for (int w = sw; w < ew; w++) {
                if (w < r4_col_sw) continue;
                uint64_t word = r4_bits[w - r4_col_sw];
                if (word) {
                    buf.active_words.push_back({ w, word, (int)buf.s4_to_global.size() });
                    while (word) {
                        int bit = (int)_tzcnt_u64(word);
                        int g_idx = w * 64 + bit;
                        if (g_idx < MS && items[g_idx].row_id == r && items[g_idx].global_id != -1) {
                            buf.global_to_s4[g_idx] = (int)buf.s4_to_global.size();
                            buf.s4_to_global.push_back(g_idx);
                        }
                        word &= (word - 1);
                    }
                }
            }
        }

        const int S4K = (int)buf.s4_to_global.size();
        if (S4K == 0) {
            for (int i = start_i; i < end_i; i++) roots_done[tid]++;
            return;
        }

        // 2. Build Semi-Local Matrix (S4 x S4) with fixed-width rows

        int s4_words = (S4K + 63) / 64 + 1; // PAD with 1 extra word for spill safety
        try {
            buf.s4_adj.assign((size_t)S4K * s4_words, 0);
            buf.s4_offsets.assign(S4K, 0);
        } catch (...) {
            printf("Error: Memory allocation failed for S4 matrix (S4K=%d, words=%d)\n", S4K, s4_words);
            exit(1);
        }
        for (int j = 0; j < S4K; j++) buf.s4_offsets[j] = (size_t)j * s4_words;

        int s4_l_cur = 0;
        for (int r = 2; r < K1212_SEARCH; r++) {
            s4_l_cur = (s4_l_cur + 63) / 64 * 64;
            while (s4_l_cur < S4K && buf.s4_to_global[s4_l_cur] != -1 && items[buf.s4_to_global[s4_l_cur]].row_id == r) 
                s4_l_cur++;
        }


        for (int j = 0; j < S4K; j++) {
            int gj = buf.s4_to_global[j];
            if (gj == -1) continue; // Skip alignment padding
            int r_j = items[gj].row_id;
            uint64_t* s4_row = &buf.s4_adj[buf.s4_offsets[j]];
            if (r_j + 1 >= K1212_SEARCH) continue;
            const uint64_t* g_row = &adj_matrix[factor_offsets[gj]];
            const int g_sw = row_ranges[r_j + 1].start_word;
            
            for (const auto& aw : buf.active_words) {
                if (aw.global_word_idx < g_sw) continue;
                uint64_t g_word = g_row[aw.global_word_idx - g_sw];
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
            int r5_v = potential_roots[i].r5_v;
            int r5_fid = search_to_factor[r5_v];
            const uint64_t* r5_bits = &adj_matrix[factor_offsets[r5_v]];
            const int r5_col_sw = row_ranges[2].start_word; // Row 5 (Slot 1) has columns starting from Slot 2
            
            thread_row4[tid] = (r4_v - row_ranges[0].start_word * 64) + 1;
            thread_row5[tid] = (r5_v - row_ranges[1].start_word * 64) + 1;
            thread_root_idx[tid] = i;

            buf.local_ranges.assign(K1212_SEARCH + 1, { 0, 0 });
            buf.local_to_global.clear();
            buf.local_s_to_f.clear();
            // NEW: Use heap-allocated buffers from ThreadLocalBuffers to avoid stack overflow
            SearchContext& local_ctx = buf.local_ctx;
            for (int r = 2; r < K1212_SEARCH; r++) {
                for (int j = 0; j < buf.dirty_s4_count[r]; j++) buf.r5_row_mask_in_s4[r][buf.dirty_s4_words[r][j]] = 0;
                buf.dirty_s4_count[r] = 0;
                buf.sub_sampling_plan_sizes[r] = 0;
            }
            int l_idx = 0;
            buf.local_to_global.clear();
            buf.local_s_to_f.clear();

            for (int r = 2; r < K1212_SEARCH; r++) {
                l_idx = (l_idx + 255) / 256 * 256;
                buf.local_ranges[r].start_word = l_idx / 64;
                const int g_sw = row_ranges[r].start_word;
                const int g_ew = row_ranges[r].end_word;

                for (int w = g_sw; w < g_ew; w++) {
                    uint64_t g_word = r5_bits[w - r5_col_sw];
                    while (g_word) {
                        int bit = (int)_tzcnt_u64(g_word);
                        int g_idx = w * 64 + bit;
                        if (g_idx < MS && items[g_idx].row_id == r && items[g_idx].global_id != -1) {
                            int s4_idx = buf.global_to_s4[g_idx];
                            if (s4_idx != -1) {
                                int w_idx = s4_idx / 64;
                                if (buf.r5_row_mask_in_s4[r][w_idx] == 0) {
                                    if (buf.dirty_s4_count[r] < K1212_WORDS) buf.dirty_s4_words[r][buf.dirty_s4_count[r]++] = w_idx;
                                }
                                buf.r5_row_mask_in_s4[r][w_idx] |= (1ULL << (s4_idx % 64));
                                buf.g_to_l[g_idx] = l_idx++;
                                buf.local_to_global.push_back(g_idx);
                            }
                        }
                        g_word &= (g_word - 1);
                    }
                }
                buf.local_ranges[r].end_word = (l_idx + 63) / 64;
            }

            if (l_idx == 0) { roots_done[tid]++; continue; }

            buf.local_offsets.assign(l_idx, 0);
            buf.local_edge_masks.assign(l_idx, { 0, 0, 0 });
            buf.local_s_to_f.assign(l_idx, -1);
            size_t l_coff = 0;
            int K = (int)buf.local_to_global.size();
            local_ctx.local_compression = (double)Global_M_total / K;

            for (int j = 0; j < K; j++) {
                int gi = buf.local_to_global[j];
                int li = buf.g_to_l[gi];
                buf.local_edge_masks[li] = factor_edge_masks[gi];
                buf.local_s_to_f[li] = search_to_factor[gi];
                buf.local_offsets[li] = l_coff;
                int d = items[gi].row_id;
                if (d + 1 < K1212_SEARCH) {
                    int wn = (buf.local_ranges[K1212_SEARCH - 1].end_word - buf.local_ranges[d + 1].start_word + 3) / 4 * 4;
                    l_coff += (size_t)wn;
                }
            }
            buf.local_adj.assign(l_coff + 16, 0);
            uint64_t* aligned_base = (uint64_t*)(((uintptr_t)buf.local_adj.data() + 63) & ~63);

            // Sub-sample S4_adj into local_adj
            // Generate sub-sampling plan once for this r5 mate
            for (int r = 2; r < K1212_SEARCH; r++) {
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
                const int l_row_sw = buf.local_ranges[d + 1].start_word;
                const int l_row_ew = buf.local_ranges[K1212_SEARCH - 1].end_word;
                const int l_row_wn = l_row_ew - l_row_sw;

                for (int r = d + 1; r < K1212_SEARCH; r++) {
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

            local_ctx.ranges = buf.local_ranges.data();
            local_ctx.adj = aligned_base;
            local_ctx.adj_size = l_coff;
            local_ctx.offsets = buf.local_offsets.data();
            local_ctx.edge_masks = buf.local_edge_masks.data();
            local_ctx.s_to_f = buf.local_s_to_f.data();
            local_ctx.r4_idx = r4_v; local_ctx.r5_idx = r5_v;

            memset(local_ctx.pool, 0, sizeof(local_ctx.pool));
            buf.cl.clear(); buf.cl.push_back(r4_fid); buf.cl.push_back(r5_fid);
            const Mask256& m4 = factor_edge_masks[r4_v];
            const Mask256& m5 = factor_edge_masks[r5_v];
            local_ctx.used_edges.m[0] = fixedEdgesMask.m[0] | m4.m[0] | m5.m[0];
            local_ctx.used_edges.m[1] = fixedEdgesMask.m[1] | m4.m[1] | m5.m[1];
            local_ctx.used_edges.m[2] = fixedEdgesMask.m[2] | m4.m[2] | m5.m[2];

            for (int idx = 0; idx < K; idx++) local_ctx.pool[2].bits[buf.g_to_l[buf.local_to_global[idx]] / 64] |= (1ULL << (buf.g_to_l[buf.local_to_global[idx]] % 64));

            internal_solve(2, buf.cl, local_ctx);
            roots_done[tid]++;
            thread_root_idx[tid] = 0x7fffffff;
        }
        for (int gi : buf.s4_to_global) if (gi != -1) buf.global_to_s4[gi] = -1;
        for (int gi : buf.local_to_global) if (gi != -1) buf.g_to_l[gi] = -1;
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


void K1212_4_20::internal_solve(int depth, std::vector<int>& clique, SearchContext& ctx) {
    if (depth == 2 && m_bPrint) diagnostic_printout(ctx.local_compression);
    if (depth == K1212_SEARCH) {
        for (int i = 0; i < K1212_FIXED; i++) memcpy(m_pResults + i * m_nPlayers, fixedRows[i].src, m_nPlayers);
        std::vector<Factor> sol_factors;
        for (int fid : clique) sol_factors.push_back(global_pool[fid]);
        std::sort(sol_factors.begin(), sol_factors.end(), [](const Factor& a, const Factor& b) { return a.src[1] < b.src[1]; });
        for (int i = 0; i < K1212_SEARCH; i++) memcpy(m_pResults + (i + K1212_FIXED) * m_nPlayers, sol_factors[i].src, m_nPlayers);
        
        {
            std::lock_guard<std::mutex> lock(result_mutex);
            // Call with mode=1 to check canonicity
            if (resultCallback(cbClass, m_pResults, ctx.r4_idx, ctx.r5_idx, 1)) {
                results_to_sort.push_back(std::vector<unsigned char>(m_pResults, m_pResults + K1212_MATCH * m_nPlayers));
                num_results++;
            }
            else {
                num_notCanon++;
            }
        }
        return;
    }

    State& P = ctx.pool[depth];
    int start = ctx.ranges[depth].start_word;
    int end = ctx.ranges[depth].end_word;

    for (int i = start; i < end; i += 4) {
        __m256i combined = _mm256_load_si256((__m256i*) & P.bits[i]);
        if (_mm256_testz_si256(combined, combined)) continue;

        const int sub_end = (end - i >= 4) ? 4 : (end - i);
        for (int sub = 0; sub < sub_end; ++sub) {
            uint64_t word = ((uint64_t*)&combined)[sub];
            while (word) {
                int v = (i + sub) * 64 + (int)_tzcnt_u64(word);

                const Mask256& fm = ctx.edge_masks[v];
#if K1212_USE_BRANCH_FREE_EDGE
                if ((ctx.used_edges.m[0] & fm.m[0]) | (ctx.used_edges.m[1] & fm.m[1]) | (ctx.used_edges.m[2] & fm.m[2])) {
#else
                if ((ctx.used_edges.m[0] & fm.m[0]) || (ctx.used_edges.m[1] & fm.m[1]) || (ctx.used_edges.m[2] & fm.m[2])) {
#endif
                    word &= (word - 1); continue;
                }
                if (depth + 1 < K1212_SEARCH) {
                    State& next_P = ctx.pool[depth + 1];
                    const uint64_t* rel_adj = &ctx.adj[ctx.offsets[v]];
                    int n_start = ctx.ranges[depth + 1].start_word;
                    int n_end = ctx.ranges[K1212_SEARCH - 1].end_word;
                    for (int j = 0; j < (n_end - n_start); j += 4) {
                        __m256i p_vec = _mm256_load_si256((__m256i*) & P.bits[j + n_start]);
                        __m256i a_vec = _mm256_load_si256((__m256i*) & rel_adj[j]);
                        _mm256_store_si256((__m256i*) & next_P.bits[j + n_start], _mm256_and_si256(p_vec, a_vec));
                    }
                }

                bool domain_wipeout = false;
                for (int d = depth + 1; d < K1212_SEARCH; d++) {
                    bool row_possible = false;
                    int w = ctx.ranges[d].start_word;
                    int wMax = ctx.ranges[d].end_word;
                    for (; w < wMax; w += 4) {
                        const __m256i val = _mm256_load_si256((__m256i*) & ctx.pool[depth + 1].bits[w]);
                        if (!_mm256_testz_si256(val, val)) { row_possible = true; break; }
                    }
                    if (!row_possible) { domain_wipeout = true; break; }
                }
                if (domain_wipeout) { word &= (word - 1); continue; }

                Mask256 old = ctx.used_edges;
                ctx.used_edges.m[0] |= fm.m[0]; ctx.used_edges.m[1] |= fm.m[1]; ctx.used_edges.m[2] |= fm.m[2];
                clique.push_back(ctx.s_to_f[v]);
                internal_solve(depth + 1, clique, ctx);
                clique.pop_back();
                ctx.used_edges = old;
                word &= (word - 1);
            }
        }
    }
}



K1212_4_20::PackedAdj K1212_4_20::pack_factor_adj(const uint8_t* adj) {
    PackedAdj res; memset(&res, 0, sizeof(res));
    for (int i = 0; i < m_nPlayers; i++) {
        int v = adj[i];
        if (i < v) {
            int eid = edge_id_table[i][v];
            if (eid != -1) {
                if (eid < 64) res.m[0] |= (1ULL << eid);
                else if (eid < 128) res.m[1] |= (1ULL << (eid - 64));
                else if (eid < 192) res.m[2] |= (1ULL << (eid - 128));
                else res.m[3] |= (1ULL << (eid - 192));
            }
        }
    }
#if K1212_USE_SIMD_CYCLE_CHECK
    alignas(16) uint8_t e2o[16] = { 0 };
    alignas(16) uint8_t o2e[16] = { 0 };
    for (int j = 0; j < 12; j++) {
        e2o[j] = (adj[j * 2] - 1) / 2;
        o2e[j] = adj[j * 2 + 1] / 2;
    }
    res.v_e2o = _mm_loadu_si128((const __m128i*)e2o);
    res.v_o2e = _mm_loadu_si128((const __m128i*)o2e);
#endif
    return res;
}

bool K1212_4_20::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    // Fast scalar cycle check
    uint32_t visited = (1 << 0);
    int count = 1, p = 0, c = adj1[0];
    while (c != 0) {
        visited |= (1 << c);
        int next = (adj1[c] == p) ? adj2[c] : adj1[c];
        p = c; c = next; count++;
        if (count > 24) return false;
    }
    if (count == 4) {
        uint32_t mask2 = ~visited & 0xFFFFFF;
        if (!mask2) return false;
        int s2 = _tzcnt_u32(mask2);
        int c2 = 1, pp = s2, curr = adj1[s2];
        while (curr != s2) {
            int next = (adj1[curr] == pp) ? adj2[curr] : adj1[curr];
            pp = curr; curr = next; c2++;
            if (c2 > 20) return false;
        }
        return c2 == 20;
    } else if (count == 20) {
        uint32_t mask2 = ~visited & 0xFFFFFF;
        if (!mask2) return false;
        int s2 = _tzcnt_u32(mask2);
        int c2 = 1, pp = s2, curr = adj1[s2];
        while (curr != s2) {
            int next = (adj1[curr] == pp) ? adj2[curr] : adj1[curr];
            pp = curr; curr = next; c2++;
            if (c2 > 4) return false;
        }
        return c2 == 4;
    }
    return false;
}

K1212_4_20::CycleUnion K1212_4_20::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    uint32_t visited = 0;
    for (int start = 0; start < m_nPlayers; start++) {
        if (visited & (1 << start)) continue;
        int curr = start, prev = -1, len = 0;
        do {
            visited |= (1 << curr);
            cu.cycles[cu.count][len++] = curr;
            int next = (adj1[curr] == prev) ? adj2[curr] : adj1[curr];
            prev = curr; curr = next;
        } while (curr != start);
        cu.lens[cu.count++] = len;
        if (cu.count >= m_nPlayers) break;
    }
    return cu;
}

void K1212_4_20::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
    info.count = 0;
    CycleUnion source_cu = find_cycles(fi.adj, fj.adj);
    if (source_cu.count != 2) return;
    if (!((source_cu.lens[0] == 4 && source_cu.lens[1] == 20) || (source_cu.lens[0] == 20 && source_cu.lens[1] == 4))) return;
    int src4 = (source_cu.lens[0] == 4) ? 0 : 1, src20 = (source_cu.lens[0] == 20) ? 0 : 1;
    int tgt4 = (target_cu.lens[0] == 4) ? 0 : 1, tgt20 = (target_cu.lens[0] == 20) ? 0 : 1;
    struct alignas(32) LocalPerm { uint8_t p[32]; uint8_t p_inv[32]; bool swaps; };
    std::vector<LocalPerm> perms4, perms20;
    auto gen_perms = [&](int src_idx, int tgt_idx, int k, std::vector<LocalPerm>& out) {
        for (int dir = 0; dir < 2; dir++) {
            for (int start_offset = 0; start_offset < k; start_offset++) {
                LocalPerm lp; 
                memset(lp.p, 0xFF, sizeof(lp.p));
                memset(lp.p_inv, 0xFF, sizeof(lp.p_inv));
                bool swaps = ((source_cu.cycles[src_idx][start_offset] % 2) != (target_cu.cycles[tgt_idx][0] % 2));
                for (int i = 0; i < k; i++) {
                    int src_pos = (dir == 0) ? (start_offset + i) % k : (start_offset - i + k) % k;
                    uint8_t src_v = (uint8_t)source_cu.cycles[src_idx][src_pos];
                    uint8_t tgt_v = (uint8_t)target_cu.cycles[tgt_idx][i];
                    lp.p[src_v] = tgt_v;
                    lp.p_inv[tgt_v] = src_v;
                }
                lp.swaps = swaps; out.push_back(lp);
            }
        }
    };
    gen_perms(src4, tgt4, 4, perms4); 
    gen_perms(src20, tgt20, 20, perms20);

    __m256i v_ff = _mm256_set1_epi8((char)0xFF);
    for (const auto& p4 : perms4) {
        __m256i v4_p = _mm256_loadu_si256((const __m256i*)p4.p);
        __m256i v4_inv = _mm256_loadu_si256((const __m256i*)p4.p_inv);
        __m256i mask_p = _mm256_cmpeq_epi8(v4_p, v_ff);
        __m256i mask_inv = _mm256_cmpeq_epi8(v4_inv, v_ff);

        for (const auto& p20 : perms20) {
            if (p4.swaps != p20.swaps) continue;
            if (info.count >= 512) break;

            auto& p_total = info.perms[info.count++];
            __m256i v20_p = _mm256_loadu_si256((const __m256i*)p20.p);
            __m256i v20_inv = _mm256_loadu_si256((const __m256i*)p20.p_inv);

            __m256i res_p = _mm256_blendv_epi8(v4_p, v20_p, mask_p);
            __m256i res_inv = _mm256_blendv_epi8(v4_inv, v20_inv, mask_inv);

            _mm256_store_si256((__m256i*)p_total.p, res_p);
            _mm256_store_si256((__m256i*)p_total.p_inv, res_inv);
            p_total.pv = res_p;
            p_total.pinvv = res_inv;
        }
    }

#if K1212_CHECK_BIPARTITE_CONSISTENCY
    for (int k = 0; k < info.count; ++k) {
        const auto& p = info.perms[k].p;
        int par_swap = p[0] % 2;
        for (int i = 0; i < m_nPlayers; i++) {
            if ((p[i] % 2) != (i % 2 ^ par_swap)) {
                printf("FATAL: Bipartite consistency check FAILED! perm[%d]=%d, p[0]=%d\n", i, p[i], p[0]);
                exit(1);
            }
        }
    }
#endif

    return;
}

void K1212_4_20::apply_perm_24(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
    // dst[perm.p[i]] = perm.p[src[i]]
    // Equivalent to lookup: dst[j] = perm.p[src[perm.p_inv[j]]]
    
    __m256i v_src = _mm256_setzero_si256();
    v_src = _mm256_inserti128_si256(v_src, _mm_loadu_si128((const __m128i*)src_adj), 0);
    v_src = _mm256_inserti128_si256(v_src, _mm_loadl_epi64((const __m128i*)(src_adj + 16)), 1);

    __m256i v_reordered = lookup_v(v_src, perm.pinvv);
    __m256i v_final = lookup_v(perm.pv, v_reordered);

    _mm_storeu_si128((__m128i*)dst_adj, _mm256_castsi256_si128(v_final));
    _mm_storel_epi64((__m128i*)(dst_adj + 16), _mm256_extracti128_si256(v_final, 1));
}
#if 0
KBipartBase::FastSortedFactor K1212_4_20::get_fast_sorted(const uint8_t* adj) {
    FastSortedFactor sf; memset(&sf, 0, sizeof(sf));
    int count = 0;
    for (int i = 0; i < m_nPlayers; ++i) {
        if (i < adj[i]) {
            sf.pairs[count * 2] = (uint8_t)i;
            sf.pairs[count * 2 + 1] = adj[i];
            count++;
        }
    }
    return sf;
}
#endif
KBipartBase::FastSortedFactor K1212_4_20::get_fast_sorted(const uint8_t* adj) {
    alignas(32) uint8_t indices[32] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    };
    __m256i v_indices = _mm256_load_si256((const __m256i*)indices);
    __m256i v_adj = _mm256_setzero_si256();
    v_adj = _mm256_inserti128_si256(v_adj, _mm_loadu_si128((const __m128i*)adj), 0);
    v_adj = _mm256_inserti128_si256(v_adj, _mm_loadl_epi64((const __m128i*)(adj + 16)), 1);
    uint32_t mask = _mm256_movemask_epi8(_mm256_cmpgt_epi8(v_adj, v_indices)) & 0x00FFFFFF;

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
    process_chunk((mask >> 16) & 0xFF, adj + 16, 16);

#if 0 // Debug verification
    FastSortedFactor ref; int count = 0;
    for (int i = 0; i < 24; ++i) if (i < adj[i]) {
        ref.pairs[count * 2] = (uint8_t)i; ref.pairs[count * 2 + 1] = adj[i]; count++;
    }
    if (memcmp(sf.pairs, ref.pairs, 24) != 0) {
        printf("GFS Error!\n");
    }
#endif
    return sf;
}
bool K1212_4_20::is_canonical(int r4_fid, int r5_fid, const K1212_4_20::TransInfo* r4_dependent_trans) {
    const uint8_t* r_adj[5] = {
        fixedRows[0].adj, fixedRows[1].adj, fixedRows[2].adj,
        global_pool[r4_fid].adj, global_pool[r5_fid].adj
    };

    FastRowTriplet target_triplet;
    target_triplet.r[0] = fixedRows[2].fs;
    target_triplet.r[1] = global_pool[r4_fid].fs;
    target_triplet.r[2] = global_pool[r5_fid].fs;
    std::sort(target_triplet.r, target_triplet.r + 3, [](const FastSortedFactor& a, const FastSortedFactor& b) {
        return compare_fast_sorted(a, b);
    });

    auto check_trans = [&](const K1212_4_20::TransInfo& info, int i, int j) {
        for (int k = 0; k < info.count; ++k) {
            const auto& perm = info.perms[k];
            uint8_t tr_fi[24], tr_fj[24];
            apply_perm_24(r_adj[i], perm, tr_fi); 
            apply_perm_24(r_adj[j], perm, tr_fj);

            FastSortedFactor fs_i = get_fast_sorted(tr_fi);
            FastSortedFactor fs_j = get_fast_sorted(tr_fj);

            bool fi_r1 = equal_fast_sorted(fs_i, r1_can);
            bool fi_r2 = equal_fast_sorted(fs_i, r2_can);
            bool fj_r1 = equal_fast_sorted(fs_j, r1_can);
            bool fj_r2 = equal_fast_sorted(fs_j, r2_can);

            if (!((fi_r1 && fj_r2) || (fi_r2 && fj_r1))) continue;

            FastSortedFactor nt_r[3]; int nt_count = 0;
            for (int m = 0; m < 5; ++m) {
                if (m == i || m == j) continue;
                uint8_t tr_m[24]; apply_perm_24(r_adj[m], perm, tr_m);
                nt_r[nt_count++] = get_fast_sorted(tr_m);
            }
            
            if (compare_fast_sorted(nt_r[1], nt_r[0])) std::swap(nt_r[0], nt_r[1]);
            if (compare_fast_sorted(nt_r[2], nt_r[1])) {
                std::swap(nt_r[1], nt_r[2]);
                if (compare_fast_sorted(nt_r[1], nt_r[0])) std::swap(nt_r[0], nt_r[1]);
            }

            for (int m = 0; m < 3; m++) {
                if (compare_fast_sorted(nt_r[m], target_triplet.r[m])) return false;
                if (compare_fast_sorted(target_triplet.r[m], nt_r[m])) break;
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

bool K1212_4_20::is_canonical_stab(int r5_fid, const Permutation* stab, int stab_count) {
    const Factor& f5 = global_pool[r5_fid];

    for (int k = 0; k < stab_count; k++) {
        const auto& perm = stab[k];
        
        // Apply automorphism of {r1, r2, r3, r4} to r5
        uint8_t tr_r5[24];
        apply_perm_24(f5.adj, perm, tr_r5);
        FastSortedFactor fs_r5 = get_fast_sorted(tr_r5);
        
        // If the image of r5 is smaller than original r5, it's not canonical.
        if (compare_fast_sorted(fs_r5, f5.fs)) return false;
    }
    return true;
}

void K1212_4_20::diagnostic_printout(double current_compr) {
    if (diag_cnt++ % 16) return;
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::seconds>(now - last_print_time).count() >= 30) {
        last_print_time = now;
        std::lock_guard<std::mutex> lock(result_mutex);
        {
            auto elap = std::chrono::duration_cast<std::chrono::seconds>(now - solve_start_time).count();
            int n_roots_done = 0;
            for (int t = 0; t < kThreads; ++t) {
                n_roots_done += roots_done[t];
            }
            double pct = (total_roots > 0) ? (100.0 * n_roots_done) / total_roots : 0;
            int expected = (pct > 0.001) ? (int)(elap * 100.0 / pct + 0.5) : 0;
            printf("[RUN] M:%05d | %%:%6.3f | T/E:%4d/%4d | Res(C/N):%d/%d | Compr:%7.4f | kThreads:%d\n",
                fixed3RowsIndex, pct, (int)elap, expected, (int)num_results.load(), (int)num_notCanon.load(), current_compr, kThreads);
        }
    }
}

bool K1212_4_20::compare_triplets(const K1212_4_20::FastRowTriplet& a, const K1212_4_20::FastRowTriplet& b) {
    for (int i = 0; i < 3; i++) {
        if (!K1212_4_20::equal_fast_sorted(a.r[i], b.r[i])) return K1212_4_20::compare_fast_sorted(a.r[i], b.r[i]);
    }
    return false;
}
