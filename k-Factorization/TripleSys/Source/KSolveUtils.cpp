#include "KSolveGen.h"
#include <algorithm>
#include <cstring>
#include <immintrin.h>
#include <nmmintrin.h>

bool KSolveGen::check_cycle_bipartite_simd(const __m128i& v_e2o_i, const __m128i& v_o2e_j) {
    __m128i v_sigma = _mm_shuffle_epi8(v_o2e_j, v_e2o_i);
    __m128i v_curr = v_sigma;
    const __m128i v_id = _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    uint32_t mask_all = (1 << (m_nPlayers / 2)) - 1;
    for (int k = 1; k <= m_nPlayers / 2; k++) {
        __m128i cmp = _mm_cmpeq_epi8(v_curr, v_id);
        uint32_t mask = _mm_movemask_epi8(cmp) & mask_all;
        if (__popcnt(mask) != (uint32_t)m_bipartiteFixedPoints[k]) return false;
        if (k < m_nPlayers / 2) v_curr = _mm_shuffle_epi8(v_sigma, v_curr);
    }
    return true;
}

bool KSolveGen::checkCycleStructure(const uint8_t* adj1, const uint8_t* adj2) {
    if (m_nPlayers > 32) {
        uint8_t found_cycles[32];
        int count = 0;
        uint64_t visited = 0;
        for (int i = 0; i < m_nPlayers; i++) {
            if (!(visited & (1ULL << i))) {
                int len = 0; int curr = i;
                do {
                    visited |= (1ULL << curr);
                    int next_u = adj1[curr];
                    visited |= (1ULL << next_u);
                    curr = adj2[next_u];
                    len += 2;
                } while (curr != i);
                found_cycles[count++] = (uint8_t)len;
            }
        }
        if (count != m_cycleCount) return false;
        if (count > 1) std::sort(found_cycles, found_cycles + count);
        return memcmp(found_cycles, m_targetCycles, count) == 0;
    }

    // Aligned 32-byte path for N <= 32
    alignas(32) uint8_t a1[32] = { 0 }, a2[32] = { 0 };
    memcpy(a1, adj1, 32); memcpy(a2, adj2, 32); // Using 32 consistently
    __m256i v_a1 = _mm256_load_si256((const __m256i*)a1);
    __m256i v_a2 = _mm256_load_si256((const __m256i*)a2);
    __m256i v_a1_high = _mm256_permute2x128_si256(v_a1, v_a1, 0x01);
    __m256i v_a2_high = _mm256_permute2x128_si256(v_a2, v_a2, 0x01);

    uint8_t found_cycles[32];
    int count = 0;
    uint32_t visited = 0;

    for (int i = 0; i < m_nPlayers; i++) {
        if (!(visited & (1U << i))) {
            int len = 0;
            int curr = i;
            do {
                visited |= (1U << curr);
                int next_u;
                if (curr < 16) next_u = ((uint8_t*)&v_a1)[curr];
                else next_u = ((uint8_t*)&v_a1_high)[curr - 16];
                
                visited |= (1U << next_u);
                
                if (next_u < 16) curr = ((uint8_t*)&v_a2)[next_u];
                else curr = ((uint8_t*)&v_a2_high)[next_u - 16];
                
                len += 2;
            } while (curr != i);
            found_cycles[count++] = (uint8_t)len;
        }
    }
    if (count != m_cycleCount) return false;
    if (count > 1) std::sort(found_cycles, found_cycles + count);
    return memcmp(found_cycles, m_targetCycles, count) == 0;
}

KSolveGen::CycleUnion KSolveGen::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    uint64_t visited = 0;
    for (int start = 0; start < m_nPlayers; start++) {
        if (visited & (1ULL << start)) continue;
        int curr = start, prev = -1, len = 0;
        do {
            visited |= (1ULL << curr);
            if (len < 32) cu.cycles[cu.count][len++] = (uint8_t)curr;
            int next = (adj1[curr] == (uint8_t)prev) ? adj2[curr] : adj1[curr];
            prev = curr; curr = next;
        } while (curr != start && len < 40);
        cu.lens[cu.count++] = len;
    }
    return cu;
}

void KSolveGen::get_transformations(const Factor& fi, const Factor& fj, TransInfo& info) {
    info.count = 0;
    CycleUnion src_cu = find_cycles(fi.adj, fj.adj);
    CycleUnion tgt_cu = find_cycles(fixedRows[0].adj, fixedRows[1].adj);

    if (src_cu.count != tgt_cu.count) return;
    if (src_cu.count == 0) return;

    struct Group { std::vector<int> src_indices; std::vector<int> tgt_indices; int len; };
    std::vector<Group> group_list;
    for (int i = 0; i < src_cu.count; i++) {
        bool found = false;
        for (auto& g : group_list) {
            if (g.len == src_cu.lens[i]) { g.src_indices.push_back(i); found = true; break; }
        }
        if (!found) group_list.push_back({{i}, {}, src_cu.lens[i]});
    }
    for (int i = 0; i < tgt_cu.count; i++) {
        bool found = false;
        for (auto& g : group_list) {
            if (g.len == tgt_cu.lens[i]) { g.tgt_indices.push_back(i); found = true; break; }
        }
        if (!found) return; 
    }
    for (auto& g : group_list) if (g.src_indices.size() != g.tgt_indices.size()) return;

    auto get_cycle_syms = [&](int src_idx, int tgt_idx, int L, std::vector<Permutation>& syms) {
        for (int dir = 0; dir < 2; dir++) {
            for (int off = 0; off < L; off++) {
                Permutation p; memset(p.p, 0xFF, 32); memset(p.p_inv, 0xFF, 32);
                for (int i = 0; i < L; i++) {
                    int sv = src_cu.cycles[src_idx][(dir == 0) ? (off + i) % L : (off - i + L) % L];
                    int tv = tgt_cu.cycles[tgt_idx][i];
                    p.p[sv] = (uint8_t)tv;
                }
                syms.push_back(p);
            }
        }
    };

    struct PartialPerm { uint8_t p[32]; uint32_t tgt_mask; };
    std::function<void(int, int, PartialPerm)> backtrack = [&](int g_idx, int c_idx, PartialPerm current_p) {
        if (info.count >= 512) return;
        if (g_idx == (int)group_list.size()) {
            Permutation final_p; memset(final_p.p, 0xFF, 32); memset(final_p.p_inv, 0xFF, 32);
            memcpy(final_p.p, current_p.p, 32);
            for (int i = 0; i < m_nPlayers; i++) {
                if (final_p.p[i] == 0xFF) return;
                final_p.p_inv[final_p.p[i]] = (uint8_t)i;
            }
            final_p.pv = _mm256_loadu_si256((const __m256i*)final_p.p);
            final_p.pinvv = _mm256_loadu_si256((const __m256i*)final_p.p_inv);
            info.perms[info.count++] = final_p;
            return;
        }

        const auto& group = group_list[g_idx];
        int src_c = group.src_indices[c_idx];
        for (int tgt_c : group.tgt_indices) {
            if (current_p.tgt_mask & (1 << tgt_c)) continue; 
            std::vector<Permutation> syms;
            get_cycle_syms(src_c, tgt_c, group.len, syms);
            for (const auto& sym : syms) {
                PartialPerm next_p = current_p;
                next_p.tgt_mask |= (1 << tgt_c);
                bool conflict = false;
                for (int v = 0; v < m_nPlayers; v++) {
                    if (sym.p[v] != 0xFF) {
                        if (next_p.p[v] != 0xFF && next_p.p[v] != sym.p[v]) { conflict = true; break; }
                        next_p.p[v] = sym.p[v];
                    }
                }
                if (!conflict) {
                    if (c_idx + 1 < (int)group.src_indices.size()) backtrack(g_idx, c_idx + 1, next_p);
                    else backtrack(g_idx + 1, 0, next_p);
                }
            }
            if (info.count >= 512) break;
        }
    };

    PartialPerm base_p; memset(base_p.p, 0xFF, 32); base_p.tgt_mask = 0;
    backtrack(0, 0, base_p);
}

FastSortedFactorBase<32> KSolveGen::get_fast_sorted(const uint8_t* adj) {
    FastSortedFactorBase<32> fs;
    memset(&fs, 0, sizeof(fs));
    int count = 0;
    for (int i = 0; i < m_nPlayers; i++) {
        if (i < adj[i]) {
            fs.pairs[count++] = (uint8_t)i;
            fs.pairs[count++] = adj[i];
        }
    }
    for (int gap = count / 4; gap > 0; gap /= 2) {
        for (int i = gap; i < count / 2; i++) {
            uint8_t temp_u = fs.pairs[i * 2];
            uint8_t temp_v = fs.pairs[i * 2 + 1];
            int j;
            for (j = i; j >= gap && fs.pairs[(j - gap) * 2] > temp_u; j -= gap) {
                fs.pairs[j * 2] = fs.pairs[(j - gap) * 2];
                fs.pairs[j * 2 + 1] = fs.pairs[(j - gap) * 2 + 1];
            }
            fs.pairs[j * 2] = temp_u;
            fs.pairs[j * 2 + 1] = temp_v;
        }
    }
    return fs;
}

void KSolveGen::apply_perm_generic(const uint8_t* src_adj, const Permutation& perm, uint8_t* dst_adj) {
    if (m_nPlayers <= 32) {
        __m256i v_src = _mm256_loadu_si256((const __m256i*)src_adj);
        __m256i v_p = perm.pv;
        __m256i v_pinv = perm.pinvv;
        
        auto shuffle32 = []( __m256i v_data, __m256i v_indices) {
            __m256i low_mask = _mm256_set1_epi8(0x0F);
            __m256i high_bit = _mm256_set1_epi8(0x10);
            __m256i v_indices_low = _mm256_and_si256(v_indices, low_mask);
            __m256i is_high = _mm256_and_si256(v_indices, high_bit);
            __m256i data_lo = _mm256_permute2x128_si256(v_data, v_data, 0x00);
            __m256i data_hi = _mm256_permute2x128_si256(v_data, v_data, 0x11);
            __m256i res_lo = _mm256_shuffle_epi8(data_lo, v_indices_low);
            __m256i res_hi = _mm256_shuffle_epi8(data_hi, v_indices_low);
            return _mm256_blendv_epi8(res_lo, res_hi, is_high);
        };
        
        __m256i v_reordered = shuffle32(v_src, v_pinv);
        __m256i v_final = shuffle32(v_p, v_reordered);
        _mm256_storeu_si256((__m256i*)dst_adj, v_final);
    } else {
        for (int i = 0; i < m_nPlayers; i++) dst_adj[perm.p[i]] = perm.p[src_adj[i]];
    }
}

bool KSolveGen::is_canonical(int r4_fid, int r5_fid, const TransInfo* r4_dependent_trans, ThreadLocalBuffers* buf) {
    const uint8_t* rows_adj[5] = {
        fixedRows[0].adj, fixedRows[1].adj, fixedRows[2].adj,
        global_pool[r4_fid].adj, global_pool[r5_fid].adj
    };

    auto check_subset = [&](const TransInfo& info, int idx1, int idx2) {
        for (int k = 0; k < info.count; k++) {
            const auto& perm = info.perms[k];
            alignas(32) uint8_t tr_i[32], tr_j[32];
            apply_perm_generic(rows_adj[idx1], perm, tr_i);
            apply_perm_generic(rows_adj[idx2], perm, tr_j);
            
            FastSortedFactorBase<32> fs_i = get_fast_sorted(tr_i);
            FastSortedFactorBase<32> fs_j = get_fast_sorted(tr_j);

            bool i_r1 = memcmp(fs_i.pairs, r1_can.pairs, 32) == 0; // Using 32 consistently
            bool i_r2 = memcmp(fs_i.pairs, r2_can.pairs, 32) == 0;
            bool j_r1 = memcmp(fs_j.pairs, r1_can.pairs, 32) == 0;
            bool j_r2 = memcmp(fs_j.pairs, r2_can.pairs, 32) == 0;

            if (!((i_r1 && j_r2) || (i_r2 && j_r1))) continue;

            FastSortedFactorBase<32> nt_r[3]; int nt_cnt = 0;
            for (int m = 0; m < 5; m++) {
                if (m == idx1 || m == idx2) continue;
                alignas(32) uint8_t tr_m[32];
                apply_perm_generic(rows_adj[m], perm, tr_m);
                nt_r[nt_cnt++] = get_fast_sorted(tr_m);
            }
            std::sort(nt_r, nt_r + 3, [&](const FastSortedFactorBase<32>& a, const FastSortedFactorBase<32>& b) {
                return memcmp(a.pairs, b.pairs, 32) < 0; 
            });

            FastSortedFactorBase<32> target[3] = { fixedRows[2].fs, global_pool[r4_fid].fs, global_pool[r5_fid].fs };
            std::sort(target, target + 3, [&](const FastSortedFactorBase<32>& a, const FastSortedFactorBase<32>& b) {
                return memcmp(a.pairs, b.pairs, 32) < 0;
            });

            for (int m = 0; m < 3; m++) {
                int cmp = memcmp(nt_r[m].pairs, target[m].pairs, 32);
                if (cmp < 0) return false;
                if (cmp > 0) break;
            }
        }
        return true;
    };

    if (!check_subset(fixed_trans[0], 0, 1)) return false;
    if (!check_subset(fixed_trans[1], 0, 2)) return false;
    if (!check_subset(fixed_trans[2], 1, 2)) return false;

    if (r4_dependent_trans) {
        for (int k = 0; k < 3; k++) if (!check_subset(r4_dependent_trans[k], k, 3)) return false;
    }

    TransInfo r5_trans;
    for (int k = 0; k < 4; k++) {
        get_transformations((k < 3) ? fixedRows[k] : global_pool[r4_fid], global_pool[r5_fid], r5_trans);
        if (!check_subset(r5_trans, k, 4)) return false;
    }
    
    return true;
}

bool KSolveGen::is_canonical_stab(int r5_fid, const Permutation* stab, int stab_count) {
    const uint8_t* r5_adj = global_pool[r5_fid].adj;
    const uint8_t* r5_can_pairs = global_pool[r5_fid].fs.pairs;
    for (int k = 0; k < stab_count; k++) {
        alignas(32) uint8_t tr_r5[32];
        apply_perm_generic(r5_adj, stab[k], tr_r5);
        FastSortedFactorBase<32> fs = get_fast_sorted(tr_r5);
        if (memcmp(fs.pairs, r5_can_pairs, 32) < 0) return false;
    }
    return true;
}
