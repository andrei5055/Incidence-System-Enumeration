#include "k16a2Old.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <intrin.h>

K16A2Old::K16A2Old(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr, bool bPrint) : KBase<Mask256_C>(factParam, bPrint) {
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

void K16A2Old::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K16A2Old solver (Simplified), MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
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
    for (int s = 0; s < m_nFixedRows; s++) {
        std::vector<uint8_t> k(fixedRows[s].adj, fixedRows[s].adj + K16_N);
        f_map[k] = s;
        global_pool.push_back(fixedRows[s]);
        packed_pool.push_back(pack_factor_adj(fixedRows[s].adj));
    }
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }
    KBase::init();
}

bool K16A2Old::addRow(int iRow, const unsigned char* source) {
    iRow--;
    if (iRow < m_nFixedRows || iRow >= K16_MATCH) return false;
    int slot = iRow - m_nFixedRows;

    Factor f; memset(&f, 0, sizeof(f));
    memcpy(f.src, source, K16_N);
    for (int i = 0; i < K16_N; i += 2) { 
        const auto a = f.src[i]; const auto b = f.src[i + 1]; 
        if (a >= K16_N || b >= K16_N) return false;
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
        if (factor_id >= K16_M_MAX) {
            printf("\nError: Number of input candidates > limit(%d), exit(1)\n", factor_id);
            exit(1);
            //return false;
        }
        f_map[k] = factor_id;
        global_pool.push_back(f);
        packed_pool.push_back(pack_factor_adj(f.adj));
    } else {
        factor_id = it->second;
    }

    addFactorToSlot(slot, factor_id);
    return true;
}

bool K16A2Old::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    uint8_t curr = 0;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    if (curr == 0) return false;
    curr = adj1[curr]; curr = adj2[curr];
    return curr == 0;
}

K16A2Old::CycleUnion K16A2Old::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu;
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

void K16A2Old::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
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
}

void K16A2Old::get_transformations_general(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl, TransInfo& info) {
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
            p_total.pv = _mm256_loadu_si256((const __m256i*)p_total.p);
            p_total.pinvv = _mm256_loadu_si256((const __m256i*)p_total.p_inv);
        }
    }
}

void K16A2Old::apply_perm_16(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
    __m128i v_src = _mm_loadu_si128((const __m128i*)src_adj);
    __m128i v_p = _mm256_castsi256_si128(perm.pv);
    __m128i v_pinv = _mm256_castsi256_si128(perm.pinvv);
    
    __m128i v_reordered = _mm_shuffle_epi8(v_src, v_pinv);
    __m128i v_final = _mm_shuffle_epi8(v_p, v_reordered);
    _mm_storeu_si128((__m128i*)dst_adj, v_final);
}

K16A2Old::FastSortedFactor K16A2Old::get_fast_sorted(const uint8_t* adj) {
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

K16A2Old::PackedAdj K16A2Old::pack_factor_adj(const uint8_t* adj) {
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
    
    uint8_t e2o[16];
    uint8_t o2e[16];
    int even_count = 0;
    int odd_count = 0;
    for (int u = 0; u < 16; u++) {
        int v = adj[u];
        if (u % 2 == 0) {
            e2o[even_count++] = v;
        } else {
            o2e[odd_count++] = v;
        }
    }
    pa.v_e2o = _mm_loadu_si128((const __m128i*)e2o);
    pa.v_o2e = _mm_loadu_si128((const __m128i*)o2e);
    return pa;
}

bool K16A2Old::compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b) {
    for (int i = 0; i < 3; i++) {
        if (compare_fast_sorted(a.r[i], b.r[i])) return true;
        if (compare_fast_sorted(b.r[i], a.r[i])) return false;
    }
    return false;
}

// Helper to determine if a permutation of 16 elements has prime order target_p
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

void K16A2Old::solve(int mode) {
    int debug_print_limit = 1000;
    int Global_M_total = (int)global_pool.size();
    if (m_bPrint) {
        printf("K16A2Old solver starting (Simplified on-the-fly pair transition search)...\n");
        printf("Row slots size: ");
        for (int s = 0; s < K16_SEARCH; s++) printf("%d:%zd ", s * 2 + 7, temp_slot_ids[s].size());
        printf("\n");
        printf("Total candidates: %d\n", Global_M_total);
    }
    std::vector<uint8_t> active(Global_M_total, 1);
    std::vector<uint32_t> f_rows_mask(Global_M_total, 0);
    for (int s = 0; s < K16_SEARCH; s++) {
        for (int id : temp_slot_ids[s]) {
            f_rows_mask[id] |= (1 << s);
        }
    }

    bool changed = true;
    int iter = 0;
    while (changed) {
        changed = false;
        for (int i = 0; i < Global_M_total; i++) {
            if (!active[i]) continue;
            if (i < 3) continue;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < K16_SEARCH; r++) {
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
    std::vector<std::vector<int>> new_slot_ids(K16_SEARCH);
    for (int s = 0; s < K16_SEARCH; s++) {
        for (int id : temp_slot_ids[s]) {
            if (active[id]) {
                new_slot_ids[s].push_back(id);
                active_slots[id] = s;
            }
        }
    }

    std::mutex results_mutex;
    std::vector<int> target_primes = {7, 5, 3, 2};
    int pairs_processed = 0;
    uint64_t total_pairs_to_process = 6;
    auto start_time = std::chrono::steady_clock::now();

    auto run_search_for_pair = [&](int type1, int val1, int type2, int val2) {
        Factor T1 = global_pool[val1];
        Factor T2 = global_pool[val2];

        uint64_t total_clique_calls = 0;
        uint64_t total_lookahead_slots = 0;
        uint64_t total_lookahead_orbits = 0;
        uint64_t total_lookahead_failures = 0;

        auto process_permutation = [&](const Permutation& alpha) {
            for (int p : target_primes) {
                if (!has_prime_order_k16(alpha.p, p)) continue;

                // Build orbits under alpha for all factors
                struct Orbit {
                    std::vector<int> fids;
                    std::vector<int> slots;
                    uint16_t slots_mask = 0;
                    Mask256_C edge_mask;
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

                        uint8_t tr_f[16];
                        apply_perm_16(global_pool[curr_fid].adj, alpha, tr_f);

                        std::vector<uint8_t> key_f(tr_f, tr_f + 16);
                        auto it_f = f_map.find(key_f);
                        if (it_f == f_map.end()) {
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
                        orb.edge_mask = Mask256_C();
                        bool conflict = false;
                        for (int f_id : orb.fids) {
                            if ((orb.edge_mask.m[0] & packed_pool[f_id].edge_mask.m[0]) ||
                                (orb.edge_mask.m[1] & packed_pool[f_id].edge_mask.m[1])) {
                                conflict = true;
                                break;
                            }
                            orb.edge_mask.m[0] |= packed_pool[f_id].edge_mask.m[0];
                            orb.edge_mask.m[1] |= packed_pool[f_id].edge_mask.m[1];
                        }
                        if (conflict) {
                            alpha_valid = false;
                            break;
                        }
                        fixed_orbits.push_back(orb);
                    } else {
                        if (orb_ok && (orb.fids.size() == 1 || orb.fids.size() == p) && cycle_ok && !slot_conflict) {
                            // Check internal edge conflict
                            orb.edge_mask = Mask256_C();
                            bool conflict = false;
                            for (int f_id : orb.fids) {
                                  if ((orb.edge_mask.m[0] & packed_pool[f_id].edge_mask.m[0]) ||
                                    (orb.edge_mask.m[1] & packed_pool[f_id].edge_mask.m[1])) {
                                    conflict = true;
                                    break;
                                }
                                orb.edge_mask.m[0] |= packed_pool[f_id].edge_mask.m[0];
                                orb.edge_mask.m[1] |= packed_pool[f_id].edge_mask.m[1];
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
                Mask256_C initial_edges;
                uint16_t initial_slots_mask = 0;
                bool fixed_compatible = true;

                std::vector<int> fixed_fids;
                for (const auto& orb : fixed_orbits) {
                    // Check edge conflicts
                    if ((initial_edges.m[0] & orb.edge_mask.m[0]) || (initial_edges.m[1] & orb.edge_mask.m[1])) {
                        fixed_compatible = false;
                        break;
                    }
                    initial_edges.m[0] |= orb.edge_mask.m[0];
                    initial_edges.m[1] |= orb.edge_mask.m[1];

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
                    if ((initial_edges.m[0] & orb.edge_mask.m[0]) || (initial_edges.m[1] & orb.edge_mask.m[1])) {
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

                // Build a tiny local compatibility mask for filtered_valid orbits
                int num_words = (num_orbits + 63) / 64;
                std::vector<uint64_t> local_compat((size_t)num_orbits * num_words, 0);

                std::vector<uint16_t> o_slots_masks(num_orbits);
                std::vector<Mask256_C> o_edge_masks(num_orbits);
                for (int i = 0; i < num_orbits; i++) {
                    o_slots_masks[i] = filtered_valid[i].slots_mask;
                    o_edge_masks[i] = filtered_valid[i].edge_mask;
                }

                for (int i = 0; i < num_orbits; i++) {
                    uint16_t mask_i = o_slots_masks[i];
                    const auto& emask_i = o_edge_masks[i];
                    for (int j = i + 1; j < num_orbits; j++) {
                        if (mask_i & o_slots_masks[j]) continue;
                        const auto& emask_j = o_edge_masks[j];
                        if ((emask_i.m[0] & emask_j.m[0]) || (emask_i.m[1] & emask_j.m[1])) continue;
                        
                        // Check pairwise compatibility of factors using automorphism symmetry
                        bool compatible = true;
                        int f1 = filtered_valid[i].fids[0];
                        for (int f2 : filtered_valid[j].fids) {
                            if (!is_perfect_scalar(packed_pool[f1].adj, packed_pool[f2].adj)) {
                                compatible = false;
                                break;
                            }
                        }
                        if (compatible) {
                            local_compat[(size_t)i * num_words + (j / 64)] |= (1ULL << (j % 64));
                            local_compat[(size_t)j * num_words + (i / 64)] |= (1ULL << (i % 64));
                        }
                    }
                }

                // Pre-index the valid orbits by slot to enable exact cover backtracking search
                std::vector<int> orbits_by_slot[K16_SEARCH];
                for (int o_idx = 0; o_idx < num_orbits; o_idx++) {
                    for (int s : filtered_valid[o_idx].slots) {
                        if (s >= 0) {
                            orbits_by_slot[s].push_back(o_idx);
                        }
                    }
                }

                // Build slot orbits bitset and non-zero words list
                std::vector<std::vector<uint64_t>> slot_orbits_bitset(K16_SEARCH, std::vector<uint64_t>(num_words, 0));
                std::vector<std::vector<int>> non_zero_words_by_slot(K16_SEARCH);
                for (int s = 0; s < K16_SEARCH; s++) {
                    for (int o_idx : orbits_by_slot[s]) {
                        slot_orbits_bitset[s][o_idx / 64] |= (1ULL << (o_idx % 64));
                    }
                    for (int w = 0; w < num_words; w++) {
                        if (slot_orbits_bitset[s][w] != 0) {
                            non_zero_words_by_slot[s].push_back(w);
                        }
                    }
                }

                // Sort slot indices by number of candidate orbits (static MRV)
                std::vector<int> slot_order(K16_SEARCH);
                for (int i = 0; i < K16_SEARCH; i++) slot_order[i] = i;
                std::sort(slot_order.begin(), slot_order.end(), [&](int a, int b) {
                    return orbits_by_slot[a].size() < orbits_by_slot[b].size();
                });

                if (m_bPrint && debug_print_limit > 0) {
                    debug_print_limit -= 3;
                    printf("  DEBUG [pair %d/%d, p=%d]: fixed=%zd, valid=%zd, filtered=%zd\n", val1, val2, p, fixed_orbits.size(), valid_orbits.size(), filtered_valid.size());
                    printf("  DEBUG [pair %d/%d, p=%d]: slot sizes: ", val1, val2, p);
                    for (int i = 0; i < K16_SEARCH; i++) printf("%d:%zd ", i, orbits_by_slot[i].size());
                    printf("\n  DEBUG [pair %d/%d, p=%d]: slot_order: ", val1, val2, p);
                    for (int i = 0; i < K16_SEARCH; i++) printf("%d ", slot_order[i]);
                    printf("\n");
                }

                // Search state
                std::vector<std::vector<uint64_t>> allowed_at_depth(K16_SEARCH + 1, std::vector<uint64_t>(num_words, ~0ULL));

                // Clique of orbits search
                std::vector<int> chosen_orbits;
                auto find_clique = [&](auto& self, uint16_t slots_mask, Mask256_C edges_mask) -> void {
                    uint64_t calls = ++total_clique_calls;
                    int d = (int)chosen_orbits.size();

                    if (m_bPrint && debug_print_limit > 0) {
                        debug_print_limit--;
                        printf("  CLIQUE_DEBUG [pair %d/%d, p=%d]: depth=%d, slots_mask=0x%03X, calls=%llu\n",
                            val1, val2, p, d, slots_mask, calls);
                    }
                    if (m_bPrint && calls % 1000000 == 0) {
                        auto current_time = std::chrono::steady_clock::now();
                        auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(current_time - last_print_time).count();
                        if (elapsed_sec >= 30) {
                            last_print_time = current_time;
                            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
                            double elapsed_min = (double)elapsed_ms / 60000.0;
                            double pct = (double)pairs_processed * 100.0 / total_pairs_to_process;
                            double est_rem_min = (pairs_processed > 0) ? (elapsed_min * (total_pairs_to_process - pairs_processed) / pairs_processed) : 0.0;
                            printf("Progress: %.1f%% executed. p=%d, calls=%llu. Elapsed: %.1f min. Est. remaining: %.1f min\n",
                                pct, p, calls, elapsed_min, est_rem_min);
                        }
                    }
                    if (slots_mask == 0x0FFF) {
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

                        unsigned char results[K16_MATCH * K16_N];
                        std::vector<Factor> sol_factors;
                        for (int fid : solution_fids) {
                            sol_factors.push_back(global_pool[fid]);
                        }
                        std::sort(sol_factors.begin(), sol_factors.end(), [](const Factor& a, const Factor& b) { return a.src[1] < b.src[1]; });
                        for (int i = 0; i < K16_MATCH; i++) memcpy(results + i * K16_N, sol_factors[i].src, K16_N);

                        std::lock_guard<std::mutex> lock(results_mutex);
                        resultCallback(cbClass, results, val1, val2, 2);
                        return;
                    }

                    // Look-ahead: check if any uncovered slot has 0 compatible candidate orbits (break early)
                    for (int i = 0; i < K16_SEARCH; i++) {
                        int slot = slot_order[i];
                        if (slots_mask & (1 << slot)) continue;
                        
                        total_lookahead_slots++;
                        bool has_compatible = false;
                        for (int w : non_zero_words_by_slot[slot]) {
                            total_lookahead_orbits++;
                            if (allowed_at_depth[d][w] & slot_orbits_bitset[slot][w]) {
                                has_compatible = true;
                                break; 
                            }
                        }
                        if (!has_compatible) {
                            total_lookahead_failures++;
                            return; 
                        }
                    }

                    // Find the first uncovered slot in slot_order (static MRV)
                    int s = -1;
                    for (int i = 0; i < K16_SEARCH; i++) {
                        int slot = slot_order[i];
                        if (!(slots_mask & (1 << slot))) {
                            s = slot;
                            break;
                        }
                    }
                    if (s == -1) return;

                    if (m_bPrint && debug_print_limit > 0) {
                        int compatible_count = 0;
                        for (int o_idx : orbits_by_slot[s]) {
                            if (allowed_at_depth[d][o_idx / 64] & (1ULL << (o_idx % 64))) {
                                compatible_count++;
                            }
                        }
                        debug_print_limit -= 2;
                        printf("  CLIQUE_DEBUG [pair %d/%d, p=%d]: Selected slot s=%d (static candidates=%zd). %d compatible candidate orbits remaining\n",
                            val1, val2, p, s, orbits_by_slot[s].size(), compatible_count);
                        
                        // Print info about the first 3 candidate orbits for slot s
                        int print_cnt = 0;
                        for (int o_idx : orbits_by_slot[s]) {
                            if (print_cnt >= 3) break;
                            const Orbit& orb = filtered_valid[o_idx];
                            bool comp = (allowed_at_depth[d][o_idx / 64] & (1ULL << (o_idx % 64)));
                            debug_print_limit--;
                            printf("    Cand %d: o_idx=%d, size=%zd, slots_mask=0x%03X, compatible=%s\n",
                                print_cnt, o_idx, orb.fids.size(), orb.slots_mask, comp ? "YES" : "NO");
                            print_cnt++;
                        }
                    }

                    // Branch on all candidate orbits covering slot s
                    for (int w : non_zero_words_by_slot[s]) {
                        uint64_t mask = allowed_at_depth[d][w] & slot_orbits_bitset[s][w];
                        while (mask > 0) {
                            unsigned long bit_pos;
                            _BitScanForward64(&bit_pos, mask);
                            int o_idx = w * 64 + bit_pos;
                            const Orbit& orb = filtered_valid[o_idx];

                            // Setup allowed orbits for the next depth
                            for (int wd = 0; wd < num_words; wd++) {
                                allowed_at_depth[d + 1][wd] = allowed_at_depth[d][wd] & local_compat[(size_t)o_idx * num_words + wd];
                            }

                            chosen_orbits.push_back(o_idx);
                            Mask256_C new_edges;
                            new_edges.m[0] = edges_mask.m[0] | orb.edge_mask.m[0];
                            new_edges.m[1] = edges_mask.m[1] | orb.edge_mask.m[1];
                            self(self, slots_mask | orb.slots_mask, new_edges);
                            chosen_orbits.pop_back();

                            mask &= mask - 1; 
                        }
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
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
            double elapsed_min = (double)elapsed_ms / 60000.0;
            double pct = (double)processed * 100.0 / total_pairs_to_process;
            double est_rem_min = (processed > 0) ? (elapsed_min * (total_pairs_to_process - processed) / processed) : 0.0;
            
            char t1_str[64] = "";
            char t2_str[64] = "";
            int p1 = 0, p2 = 0;
            for (int i = 0; i < 16; i += 2) {
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
}
