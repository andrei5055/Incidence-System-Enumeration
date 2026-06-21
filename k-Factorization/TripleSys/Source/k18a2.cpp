#include "k18a2.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <intrin.h>
#include <set>

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

K18A2::K18A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows,
    ResultCallback callback, void* cbClassPtr, bool bPrint) : KBase<Mask18_C>(factParam, bPrint) {
    int id_cnt = 0;
    memset(edge_id_table, -1, sizeof(edge_id_table));
    for (int i = 0; i < K18_N; i++) {
        for (int j = i + 1; j < K18_N; j++) {
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

void K18A2::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    // Initialize static edge mask lookup table
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 18; j++) {
            edge_mask_table[i][j] = Mask18_C();
            if (i != j) {
                int eid = edge_id_table[i][j];
                if (eid != -1) {
                    if (eid < 64) edge_mask_table[i][j].m[0] |= (1ULL << eid);
                    else if (eid < 128) edge_mask_table[i][j].m[1] |= (1ULL << (eid - 64));
                    else if (eid < 192) edge_mask_table[i][j].m[2] |= (1ULL << (eid - 128));
                }
            }
        }
    }

    // Reset lazy candidate cache states
    for (int i = 0; i < K18_SEARCH; i++) {
        slot_generated[i] = false;
    }

    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K18A2 solver (Path Progress version), Compiled: " << __DATE__ << " " << __TIME__ << ", MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
    }
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    start_time = std::chrono::steady_clock::now();
    last_print_time = start_time;
    bTimeSet = false;
    call_counter = 0;

    fixedEdgesMask = Mask18_C();
    for (int s = 0; s < m_nFixedRows; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * K18_N, K18_N);
        for (int i = 0; i < K18_N; i += 2) {
            const auto a = fixedRows[s].src[i]; const auto b = fixedRows[s].src[i + 1];
            fixedRows[s].adj[a] = b; fixedRows[s].adj[b] = a;
            int eid = edge_id_table[a][b];
            if (eid < 64) fixedEdgesMask.m[0] |= (1ULL << eid);
            else if (eid < 128) fixedEdgesMask.m[1] |= (1ULL << (eid - 64));
            else if (eid < 192) fixedEdgesMask.m[2] |= (1ULL << (eid - 128));
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
        std::vector<uint8_t> k(fixedRows[s].adj, fixedRows[s].adj + K18_N);
        f_map[k] = s;
        std::array<uint8_t, 18> k_arr;
        std::copy(fixedRows[s].adj, fixedRows[s].adj + 18, k_arr.begin());
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

        if (!has_prime_order_k18(alpha.p, 17)) {
            continue;
        }

        bool used_v0[18] = { false };
        used_v0[fixedRows[0].adj[0]] = true;
        used_v0[fixedRows[1].adj[0]] = true;
        used_v0[fixedRows[2].adj[0]] = true;

        bool chain_ok = true;
        uint8_t current_adj[18];
        memcpy(current_adj, R_c.adj, 18);

        ValidChain chain;
        chain.alpha = alpha;
        memcpy(chain.factors_adj[1], R_a.adj, 18);
        memcpy(chain.factors_adj[2], R_b.adj, 18);
        memcpy(chain.factors_adj[3], R_c.adj, 18);
        chain.factors_packed[1] = pack_factor_adj(R_a.adj);
        chain.factors_packed[2] = pack_factor_adj(R_b.adj);
        chain.factors_packed[3] = pack_factor_adj(R_c.adj);

        // Check initial compatibility
        if (!is_perfect_packed(chain.factors_packed[1], chain.factors_packed[2]) ||
            !is_perfect_packed(chain.factors_packed[1], chain.factors_packed[3]) ||
            !is_perfect_packed(chain.factors_packed[2], chain.factors_packed[3])) {
            continue;
        }

        for (int k = 4; k <= 17; k++) {
            uint8_t next_adj[18];
            apply_perm_18(current_adj, alpha, next_adj);

            uint8_t v0 = next_adj[0];
            if (used_v0[v0]) {
                chain_ok = false;
                break;
            }
            used_v0[v0] = true;

            memcpy(chain.factors_adj[k], next_adj, 18);
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

            memcpy(current_adj, next_adj, 18);
        }

        if (chain_ok) {
            valid_starter_chains++;
            std::lock_guard<std::mutex> lock(result_mutex);
            if (std::find_if(saved_transitions.begin(), saved_transitions.end(), [&](const Permutation& p) {
                return memcmp(p.p, alpha.p, 18) == 0;
            }) == saved_transitions.end()) {
                saved_transitions.push_back(alpha);
            }
        }
    }

    if (m_bPrint) {
        printf("Transitions: %d candidate transitions, %d valid starter chains recorded.\n", total_candidates, valid_starter_chains);
    }
}

bool K18A2::addRow(int rowNum, const unsigned char* source) {
    Factor f;
    memcpy(f.src, source, 18);
    for (int i = 0; i < 18; i += 2) {
        uint8_t u = f.src[i]; uint8_t v = f.src[i + 1];
        f.adj[u] = v; f.adj[v] = u;
    }

    PackedAdj pf = pack_factor_adj(f.adj);
    
    // Check overlap with fixed rows
    if ((pf.edge_mask.m[0] & fixedEdgesMask.m[0]) ||
        (pf.edge_mask.m[1] & fixedEdgesMask.m[1]) ||
        (pf.edge_mask.m[2] & fixedEdgesMask.m[2])) return false;

    // Check compatibility with the 3 fixed rows
    for (int i = 0; i < K18_FIXED; i++) {
        if (!is_perfect_packed(pf, fixed_packed[i])) return false;
    }

    f.edge_mask = pf.edge_mask;
    f.fs = get_fast_sorted(f.adj);

    std::vector<uint8_t> key(f.adj, f.adj + 18);
    std::array<uint8_t, 18> key_arr;
    std::copy(f.adj, f.adj + 18, key_arr.begin());

    int factor_id;
    {
        std::lock_guard<std::mutex> lock(pool_mutex);
        auto it = f_map.find(key);
        if (it == f_map.end()) {
            factor_id = (int)global_pool.size();
            f_map[key] = factor_id;
            f_map_unordered[key_arr] = factor_id;
            global_pool.push_back(f);
            packed_pool.push_back(pf);
        } else {
            factor_id = it->second;
        }
    }

    int slot = rowNum - 4;
    addFactorToSlot(slot, factor_id);
    return true;
}



K18A2::CycleUnion K18A2::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    bool visited[18] = { false };
    for (int i = 0; i < 18; ++i) {
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

void K18A2::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
    CycleUnion source_cu = find_cycles(fi.adj, fj.adj);
    if (source_cu.count != 1 || source_cu.lens[0] != 18) return;
    if (target_cu.count != 1 || target_cu.lens[0] != 18) return;

    for (int dir = 0; dir < 2; dir++) {
        for (int start_offset = 0; start_offset < 18; start_offset++) {
            if (info.count >= 512) break;
            auto& p_total = info.perms[info.count++];
            memset(p_total.p, 0xFF, sizeof(p_total.p));
            memset(p_total.p_inv, 0xFF, sizeof(p_total.p_inv));
            
            for (int i = 0; i < 18; i++) {
                int src_pos = (dir == 0) ? (start_offset + i) % 18 : (start_offset - i + 18) % 18;
                uint8_t src_v = (uint8_t)source_cu.cycles[0][src_pos];
                uint8_t tgt_v = (uint8_t)target_cu.cycles[0][i];
                p_total.p[src_v] = tgt_v;
                p_total.p_inv[tgt_v] = src_v;
            }
        }
    }
}

void K18A2::get_transformations_general(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl, TransInfo& info) {
    info.count = 0;
    CycleUnion cu_src = find_cycles(fi.adj, fj.adj);
    if (cu_src.count != 1 || cu_src.lens[0] != 18) return;
    CycleUnion cu_tgt = find_cycles(fk.adj, fl.adj);
    if (cu_tgt.count != 1 || cu_tgt.lens[0] != 18) return;

    for (int dir = 0; dir < 2; dir++) {
        for (int start_offset = (dir == 0 ? 0 : 1); start_offset < 18; start_offset += 2) {
            if (info.count >= 512) break;
            auto& p_total = info.perms[info.count++];
            memset(p_total.p, 0xFF, sizeof(p_total.p));
            memset(p_total.p_inv, 0xFF, sizeof(p_total.p_inv));
            
            for (int i = 0; i < 18; i++) {
                int src_pos = (dir == 0) ? (start_offset + i) % 18 : (start_offset - i + 18) % 18;
                uint8_t src_v = (uint8_t)cu_src.cycles[0][src_pos];
                uint8_t tgt_v = (uint8_t)cu_tgt.cycles[0][i];
                p_total.p[src_v] = tgt_v;
                p_total.p_inv[tgt_v] = src_v;
            }
        }
    }
}



K18A2::FastSortedFactor K18A2::get_fast_sorted(const uint8_t* adj) {
    FastSortedFactor sf; memset(&sf, 0, sizeof(sf));
    uint8_t* out = sf.pairs;
    int count = 0;
    for (int u = 0; u < 18; u++) {
        int v = adj[u];
        if (u < v) {
            out[count * 2] = (uint8_t)u;
            out[count * 2 + 1] = (uint8_t)v;
            count++;
        }
    }
    return sf;
}

K18A2::PackedAdj K18A2::pack_factor_adj(const uint8_t* adj) {
    PackedAdj pa;
    memcpy(pa.adj, adj, 18);
    pa.edge_mask = Mask18_C();
    for (int u = 0; u < 18; u++) {
        int v = adj[u];
        if (u < v) {
            int eid = edge_id_table[u][v];
            if (eid < 64) pa.edge_mask.m[0] |= (1ULL << eid);
            else if (eid < 128) pa.edge_mask.m[1] |= (1ULL << (eid - 64));
            else if (eid < 192) pa.edge_mask.m[2] |= (1ULL << (eid - 128));
        }
    }
    return pa;
}

bool K18A2::compare_triplets(const FastRowTriplet& a, const FastRowTriplet& b) {
    for (int i = 0; i < 3; i++) {
        if (compare_fast_sorted(a.r[i], b.r[i])) return true;
        if (compare_fast_sorted(b.r[i], a.r[i])) return false;
    }
    return false;
}

void K18A2::solve(int mode) {
    runTransitionSearch();
}
