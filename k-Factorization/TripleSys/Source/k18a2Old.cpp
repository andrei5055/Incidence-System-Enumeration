#include "k18a2Old.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <functional>
#include <mutex>
#include <thread>
#include <atomic>

// =============================================================================
// K18A2Old — see k18a2Old.h and MD/k18a2old_spec.md §2.
// Scalar NP=18 adaptation of K16A2Old's pool + arc-consistency + orbit-clique engine.
// =============================================================================

K18A2Old::K18A2Old(const FactorParams& factParam, int fixed3RowsIndex, int kThreads,
                   const unsigned char* first3Rows, ResultCallback callback,
                   void* cbClassPtr, bool bPrint)
    : KBase<Mask18O_C>(factParam, bPrint) {
    int id_cnt = 0;
    memset(edge_id_table, -1, sizeof(edge_id_table));
    for (int i = 0; i < K18O_N; i++) {
        for (int j = i + 1; j < K18O_N; j++) {
            int eid = id_cnt++;                 // 0..152
            edge_id_table[i][j] = eid;
            edge_id_table[j][i] = eid;
        }
    }
    init(fixed3RowsIndex, kThreads, first3Rows, callback, cbClassPtr);
}

void K18A2Old::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows,
                    ResultCallback callback, void* cbClassPtr) {
    if (kThreads < 1) kThreads = 1;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    start_time = std::chrono::steady_clock::now();
    last_print_time = start_time;
    bTimeSet = false;

    if (m_bPrint)
        std::cout << "Init: K18A2Old solver (starter pool + arc-consistency + orbit-clique), "
                  << "Compiled: " << __DATE__ << " " << __TIME__ << ", kThreads: " << kThreads << std::endl;

    fixedEdgesMask = Mask18O_C();
    for (int s = 0; s < m_nFixedRows; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * K18O_N, K18O_N);
        fixedRows[s].edge_mask = Mask18O_C();
        for (int i = 0; i < K18O_N; i += 2) {
            const auto a = fixedRows[s].src[i];
            const auto b = fixedRows[s].src[i + 1];
            fixedRows[s].adj[a] = b;
            fixedRows[s].adj[b] = a;
            int eid = edge_id_table[a][b];
            set_edge_mask(fixedRows[s].edge_mask, eid);
            set_edge_mask(fixedEdgesMask, eid);
        }
        fixed_packed[s] = pack_factor_adj(fixedRows[s].adj);
    }

    // pool index 0..2 are the fixed rows themselves
    global_pool.clear();
    packed_pool.clear();
    f_map.clear();
    temp_slot_ids.assign(K18O_SEARCH, {});
    for (int s = 0; s < m_nFixedRows; s++) {
        std::vector<uint8_t> k(fixedRows[s].adj, fixedRows[s].adj + K18O_N);
        f_map[k] = s;
        global_pool.push_back(fixedRows[s]);
        packed_pool.push_back(fixed_packed[s]);
    }
    KBase::init();
}

K18A2Old::CycleUnion K18A2Old::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    bool visited[K18O_N] = { false };
    for (int i = 0; i < K18O_N; ++i) {
        if (visited[i]) continue;
        int curr = i, cycle_pos = 0;
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

K18A2Old::PackedAdj K18A2Old::pack_factor_adj(const uint8_t* adj) {
    PackedAdj pa;
    memcpy(pa.adj, adj, K18O_N);
    pa.edge_mask = Mask18O_C();
    for (int u = 0; u < K18O_N; u++) {
        int v = adj[u];
        if (u < v) set_edge_mask(pa.edge_mask, edge_id_table[u][v]);
    }
    return pa;
}

bool K18A2Old::has_order(const uint8_t* p, int target_order) {
    // order(p) == lcm of cycle lengths; require it to equal target_order exactly
    bool visited[K18O_N] = { false };
    long ord = 1;
    auto lcm = [](long a, long b) { long g = a; long bb = b; while (bb) { long t = a % bb; a = bb; bb = t; } return g / a * b; };
    for (int i = 0; i < K18O_N; i++) {
        if (visited[i]) continue;
        int len = 0, curr = i;
        while (!visited[curr]) { visited[curr] = true; curr = p[curr]; len++; }
        if (len > 1) ord = lcm(ord, len);
    }
    return ord == target_order;
}

bool K18A2Old::addRow(int iRow, const unsigned char* source) {
    auto _t0 = std::chrono::steady_clock::now();
    struct _Acc { K18A2Old* s; std::chrono::steady_clock::time_point t0;
        ~_Acc() { s->m_addRowSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); s->m_addRowCalls++; } } _acc{ this, _t0 };
    // append one candidate factor row (compatible with the 3 fixed rows) to the pool.
    iRow--;
    if (iRow < m_nFixedRows || iRow >= K18O_MATCH) return false;
    int slot = iRow - m_nFixedRows;

    Factor f; memset(&f, 0, sizeof(f));
    memcpy(f.src, source, K18O_N);
    f.edge_mask = Mask18O_C();
    for (int i = 0; i < K18O_N; i += 2) {
        const auto a = f.src[i], b = f.src[i + 1];
        if (a >= K18O_N || b >= K18O_N) return false;
        f.adj[a] = b; f.adj[b] = a;
        set_edge_mask(f.edge_mask, edge_id_table[a][b]);
    }
    if ((f.edge_mask.m[0] & fixedEdgesMask.m[0]) ||
        (f.edge_mask.m[1] & fixedEdgesMask.m[1]) ||
        (f.edge_mask.m[2] & fixedEdgesMask.m[2])) {
        printf("\n*** K18A2Old: row shares an edge with a fixed row!\n");
        return false;
    }
    PackedAdj pf = pack_factor_adj(f.adj);
    for (int i = 0; i < m_nFixedRows; i++)
        if (!is_perfect_packed(pf, fixed_packed[i])) {
            printf("\n*** K18A2Old: row not Hamiltonian with fixed row %d!\n", i + 1);
            return false;
        }

    std::vector<uint8_t> k(f.adj, f.adj + K18O_N);
    int fid;
    auto it = f_map.find(k);
    if (it == f_map.end()) {
        fid = (int)global_pool.size();
        if (fid >= K18O_M_MAX) { printf("\n*** K18A2Old: pool > %d, exit\n", K18O_M_MAX); exit(1); }
        f_map[k] = fid;
        global_pool.push_back(f);
        packed_pool.push_back(pf);
    } else {
        fid = it->second;
    }
    temp_slot_ids[slot].push_back(fid);
    return true;
}

// ---- Arc-consistency: drop any candidate that has no compatible partner in some slot ----
void K18A2Old::arcConsistencyPrune(std::vector<uint8_t>& active) {
    int M = (int)global_pool.size();
    std::vector<uint32_t> f_rows_mask(M, 0);
    for (int s = 0; s < K18O_SEARCH; s++)
        for (int id : temp_slot_ids[s]) f_rows_mask[id] |= (1u << s);

    bool changed = true;
    int iter = 0;
    while (changed) {
        changed = false;
        for (int i = 0; i < M; i++) {
            if (!active[i] || i < m_nFixedRows) continue;
            const PackedAdj& pi = packed_pool[i];
            for (int r = 0; r < K18O_SEARCH; r++) {
                if (f_rows_mask[i] & (1u << r)) continue;      // i already serves slot r
                bool found = false;
                for (int id : temp_slot_ids[r])
                    if (active[id] && is_perfect_packed(pi, packed_pool[id])) { found = true; break; }
                if (!found) { active[i] = 0; changed = true; break; }
            }
        }
        if (m_bPrint) {
            int rem = 0; for (int i = 0; i < M; i++) if (active[i]) rem++;
            printf("  [K18A2Old] arc-consistency iter %d: remaining %d\n", ++iter, rem);
            fflush(stdout);
        }
    }
}

// ---- Per-alpha: group the pruned pool into alpha-orbits and clique exact-cover ----
// Faithful adaptation of K16A2Old::process_permutation + find_clique, with:
//   * order-4 (composite): valid orbit size must DIVIDE K18O_AUT_ORDER (i.e. 1,2,4),
//   * NP=18 -> full slots mask = (1<<K18O_SEARCH)-1, 17 result rows,
//   * 256-bit edge masks (m[0..2]), scalar is_perfect.
void K18A2Old::completeForAlpha(const uint8_t* alpha, const std::vector<uint8_t>& active,
                                const std::vector<int>& active_slots) {
    const int M = (int)global_pool.size();
    if (!has_order(alpha, K18O_AUT_ORDER)) return;

    struct Orbit {
        std::vector<int> fids;
        std::vector<int> slots;
        uint32_t slots_mask = 0;
        Mask18O_C edge_mask;
    };
    std::vector<Orbit> fixed_orbits, valid_orbits;
    std::vector<uint8_t> visited(M, 0);
    bool alpha_valid = true;

    for (int fid = 0; fid < M && alpha_valid; fid++) {
        if (fid >= m_nFixedRows && !active[fid]) continue;
        if (visited[fid]) continue;

        Orbit orb;
        int curr = fid;
        bool orb_ok = true;
        while (true) {
            visited[curr] = 1;
            orb.fids.push_back(curr);
            orb.slots.push_back(active_slots[curr]);
            uint8_t tr[K18O_N];
            apply_perm(global_pool[curr].adj, alpha, tr);
            std::vector<uint8_t> key(tr, tr + K18O_N);
            auto it = f_map.find(key);
            if (it == f_map.end()) { orb_ok = false; break; }     // image not a candidate
            int nxt = it->second;
            if (nxt >= m_nFixedRows && !active[nxt]) { orb_ok = false; break; }
            if (nxt == fid) break;                                  // orbit closed
            if (visited[nxt]) { orb_ok = false; break; }
            curr = nxt;
        }

        // mutual compatibility (edge-disjoint Hamiltonian) within the orbit
        bool cycle_ok = true;
        for (size_t j = 1; j < orb.fids.size() && cycle_ok; j++)
            if (!is_perfect_packed(packed_pool[orb.fids[0]], packed_pool[orb.fids[j]])) cycle_ok = false;

        // slot conflicts + slots_mask
        bool slot_conflict = false;
        orb.slots_mask = 0;
        for (int sl : orb.slots) if (sl >= 0) {
            if (orb.slots_mask & (1u << sl)) { slot_conflict = true; break; }
            orb.slots_mask |= (1u << sl);
        }
        const int d = (int)orb.fids.size();
        const bool size_ok = (K18O_AUT_ORDER % d == 0);            // size divides aut order (1,2,4)

        bool contains_fixed = false;
        for (int f : orb.fids) if (f < m_nFixedRows) contains_fixed = true;

        // internal edge-disjointness + accumulate orbit edge mask
        auto build_edges = [&](Orbit& o) -> bool {
            o.edge_mask = Mask18O_C();
            for (int f : o.fids) {
                const auto& em = packed_pool[f].edge_mask;
                if ((o.edge_mask.m[0] & em.m[0]) || (o.edge_mask.m[1] & em.m[1]) ||
                    (o.edge_mask.m[2] & em.m[2])) return false;
                o.edge_mask.m[0] |= em.m[0]; o.edge_mask.m[1] |= em.m[1]; o.edge_mask.m[2] |= em.m[2];
            }
            return true;
        };

        if (contains_fixed) {
            if (!orb_ok || !size_ok || !cycle_ok || slot_conflict || !build_edges(orb)) { alpha_valid = false; break; }
            fixed_orbits.push_back(orb);
        } else if (orb_ok && size_ok && cycle_ok && !slot_conflict && build_edges(orb)) {
            valid_orbits.push_back(orb);
        }
    }
    if (!alpha_valid) return;

    // assemble fixed orbits: edge/slot compatibility + collect fixed fids
    Mask18O_C init_edges; uint32_t init_slots = 0;
    std::vector<int> fixed_fids;
    for (const auto& o : fixed_orbits) {
        if ((init_edges.m[0] & o.edge_mask.m[0]) || (init_edges.m[1] & o.edge_mask.m[1]) ||
            (init_edges.m[2] & o.edge_mask.m[2]) || (init_slots & o.slots_mask)) return;
        init_edges.m[0] |= o.edge_mask.m[0]; init_edges.m[1] |= o.edge_mask.m[1]; init_edges.m[2] |= o.edge_mask.m[2];
        init_slots |= o.slots_mask;
        for (int f : o.fids) fixed_fids.push_back(f);
    }
    for (size_t i = 0; i < fixed_fids.size(); i++)
        for (size_t j = i + 1; j < fixed_fids.size(); j++)
            if (!is_perfect_packed(packed_pool[fixed_fids[i]], packed_pool[fixed_fids[j]])) return;

    // filter free orbits to those compatible with the fixed selection
    std::vector<Orbit> cand;
    for (const auto& o : valid_orbits) {
        if ((init_edges.m[0] & o.edge_mask.m[0]) || (init_edges.m[1] & o.edge_mask.m[1]) ||
            (init_edges.m[2] & o.edge_mask.m[2]) || (init_slots & o.slots_mask)) continue;
        bool ok = true;
        for (int ff : fixed_fids) { for (int nf : o.fids) if (!is_perfect_packed(packed_pool[ff], packed_pool[nf])) { ok = false; break; } if (!ok) break; }
        if (ok) cand.push_back(o);
    }
    const int N = (int)cand.size();
    // NOTE: do NOT early-return on N==0. When the fixed orbits already tile every search slot
    // (e.g. an order-NM automorphism whose single orbit IS the whole P1F), there are no free
    // orbits, but cover() must still emit. With N==0 the structures below are empty and
    // cover(init_slots, .) emits immediately because init_slots == FULL.

    // pairwise orbit-compatibility (edge-disjoint + Hamiltonian); slot-disjoint implicit
    std::vector<uint8_t> compat((size_t)N * N, 0);
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            if (cand[i].slots_mask & cand[j].slots_mask) continue;
            if ((cand[i].edge_mask.m[0] & cand[j].edge_mask.m[0]) ||
                (cand[i].edge_mask.m[1] & cand[j].edge_mask.m[1]) ||
                (cand[i].edge_mask.m[2] & cand[j].edge_mask.m[2])) continue;
            bool ok = true;
            for (int a : cand[i].fids) { for (int b : cand[j].fids) if (!is_perfect_scalar(packed_pool[a].adj, packed_pool[b].adj)) { ok = false; break; } if (!ok) break; }
            if (ok) { compat[(size_t)i * N + j] = 1; compat[(size_t)j * N + i] = 1; }
        }

    // orbits-by-slot index for exact-cover branching
    std::vector<std::vector<int>> by_slot(K18O_SEARCH);
    for (int o = 0; o < N; o++) for (int s : cand[o].slots) if (s >= 0) by_slot[s].push_back(o);
    const uint32_t FULL = (K18O_SEARCH >= 32) ? 0xFFFFFFFFu : ((1u << K18O_SEARCH) - 1u);

    std::vector<int> chosen;
    std::vector<uint8_t> allowed(N, 1);                 // depth-0 allowed set
    // recursive exact cover over slots (MRV: fewest candidate orbits first)
    std::function<void(uint32_t, std::vector<uint8_t>&)> cover =
        [&](uint32_t slots_mask, std::vector<uint8_t>& allow) {
        if (slots_mask == FULL) {
            // emit: fixed_fids + chosen orbits' fids => 17 rows
            std::vector<int> sol(fixed_fids);
            for (int o : chosen) for (int f : cand[o].fids) sol.push_back(f);
            if ((int)sol.size() != K18O_MATCH) return;
            for (size_t i = 0; i < sol.size(); i++)
                for (size_t j = i + 1; j < sol.size(); j++)
                    if (!is_perfect_scalar(global_pool[sol[i]].adj, global_pool[sol[j]].adj)) return;
            unsigned char res[K18O_MATCH * K18O_N];
            std::vector<Factor> sf; for (int f : sol) sf.push_back(global_pool[f]);
            std::sort(sf.begin(), sf.end(), [](const Factor& a, const Factor& b) { return a.src[1] < b.src[1]; });
            for (int i = 0; i < K18O_MATCH; i++) memcpy(res + i * K18O_N, sf[i].src, K18O_N);
            std::lock_guard<std::mutex> lk(result_mutex);
            solutions++;
            if (resultCallback) resultCallback(cbClass, res, 0, 0, 2);
            return;
        }
        // pick uncovered slot with fewest live candidate orbits (MRV)
        int best = -1, bestcnt = 1 << 30;
        for (int s = 0; s < K18O_SEARCH; s++) {
            if (slots_mask & (1u << s)) continue;
            int c = 0; for (int o : by_slot[s]) if (allow[o]) c++;
            if (c == 0) return;                          // dead slot -> prune
            if (c < bestcnt) { bestcnt = c; best = s; }
        }
        if (best < 0) return;
        for (int o : by_slot[best]) {
            if (!allow[o]) continue;
            std::vector<uint8_t> na(allow);
            na[o] = 0;
            for (int k = 0; k < N; k++) if (na[k] && k != o && !compat[(size_t)o * N + k]) na[k] = 0;
            chosen.push_back(o);
            cover(slots_mask | cand[o].slots_mask, na);
            chosen.pop_back();
        }
    };
    cover(init_slots, allowed);
}

// =============================================================================
// Candidate order-4 automorphism generator — port of K18A2::CycleBacktrackState (L=4).
// 4 does NOT divide 18, so order-4 alphas cannot be Hamiltonian-cycle rotations
// (K16A2Old's get_transformations). We enumerate them by backtracking exactly as the
// proven cyclic K18A2 does, but WITHOUT the decomposeMissingEdges existence check:
// every structurally valid alpha is saved, then completed by completeForAlpha (orbit
// clique) instead. The two search types mirror K18A2 (fixed vs transposed out-of-cycle).
// =============================================================================

// Prune alphas that would force a fixed-row edge onto another fixed-row edge within L/2
// steps of the orbit (i.e. no invariant completion is possible). F0[0] is the first row.
static bool k18o_validateDefinedOrbits(const uint8_t* alpha_p, int L, const bool* defined,
                                       const uint8_t F0[][K18O_N]) {
    int limit = L / 2;
    for (int u = 0; u < K18O_N; u++) {
        if (!defined[u]) continue;
        int v = F0[0][u];
        if (!defined[v] || u > v) continue;
        uint8_t curr_u = (uint8_t)u, curr_v = (uint8_t)v;
        for (int j = 1; j <= limit; j++) {
            curr_u = alpha_p[curr_u];
            curr_v = alpha_p[curr_v];
            if (F0[0][curr_u] == curr_v) return false;
        }
    }
    return true;
}

void K18A2Old::CycleBacktrackState::apply_perm(const uint8_t* src_adj, const uint8_t* perm, uint8_t* dst_adj) {
    for (int i = 0; i < K18O_N; i++) dst_adj[perm[i]] = perm[src_adj[i]];
}

bool K18A2Old::CycleBacktrackState::is_perfect_scalar(const uint8_t* adj1, const uint8_t* adj2) {
    uint8_t curr = 0;
    for (int i = 0; i < K18O_HALF; i++) {
        curr = adj1[curr];
        curr = adj2[curr];
        if (curr == 0 && i < K18O_HALF - 1) return false;
    }
    return curr == 0;
}

void K18A2Old::CycleBacktrackState::backtrack(int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    if (depth == L - 1) { processCandidate(c, used, L); return; }
    for (int v = 0; v < K18O_N; v++) {
        if (used[v]) continue;
        tryVertexForCycle(v, depth, pairs_visited, c, used, L);
    }
}

void K18A2Old::CycleBacktrackState::tryVertexForCycle(int v, int depth, int pairs_visited, uint8_t* c, bool* used, int L) {
    int p_idx = vertex_to_pair[v];
    if (p_idx == -2) { recurseWithVertex(v, depth, pairs_visited, c, used, L); return; }
    bool pair_visited = used[pair_elements[p_idx][0]] || used[pair_elements[p_idx][1]];
    if (pair_visited) {
        recurseWithVertex(v, depth, pairs_visited, c, used, L);
    } else if (v == pair_elements[pairs_visited][0]) {
        recurseWithVertex(v, depth, pairs_visited + 1, c, used, L);
    }
}

void K18A2Old::CycleBacktrackState::recurseWithVertex(int v, int depth, int next_pairs, uint8_t* c, bool* used, int L) {
    used[v] = true;
    c[depth] = v;
    backtrack(depth + 1, next_pairs, c, used, L);
    used[v] = false;
}

void K18A2Old::CycleBacktrackState::buildPermutation(uint8_t* c, uint8_t* alpha_p, int L) {
    for (int i = 0; i < K18O_N; i++) alpha_p[i] = (uint8_t)i;
    if (L > 1) {
        alpha_p[v0] = c[0];
        for (int i = 0; i < L - 2; i++) alpha_p[c[i]] = c[i + 1];
        alpha_p[c[L - 2]] = v0;
    }
}

void K18A2Old::CycleBacktrackState::processCandidate(uint8_t* c, bool* used, int L) {
    total_generated++;

    // For even powers-of-2 L, a valid completion needs an invariant factor, which requires
    // the main-cycle antipodal edges to be missing in F[0]. Prune otherwise.
    if (L == 2 || L == 4 || L == 8 || L == 16) {
        int half = L / 2;
        for (int i = 0; i < half; i++) {
            uint8_t u = (i == 0) ? v0 : c[i - 1];
            uint8_t v = c[i + half - 1];
            if (F[0][u] == v) return;
        }
    }

    uint8_t alpha_p[K18O_N];
    buildPermutation(c, alpha_p, L);
    if (search_type == 2) { alpha_p[0] = 1; alpha_p[1] = 0; }

    bool defined[K18O_N] = { false };
    defined[v0] = true;
    for (int i = 0; i < L - 1; i++) defined[c[i]] = true;
    if (!k18o_validateDefinedOrbits(alpha_p, L, defined, F)) return;

    bool in_main_cycle[K18O_N] = { false };
    in_main_cycle[v0] = true;
    for (int i = 0; i < L - 1; i++) in_main_cycle[c[i]] = true;

    uint8_t rem[K18O_N];
    int rem_size = 0;
    int start_u = (search_type == 2) ? 2 : 1;
    for (int u = start_u; u < K18O_N; u++)
        if (!in_main_cycle[u]) rem[rem_size++] = u;

    bool rem_used[K18O_N] = { false };
    generate_remaining_cycles(0, rem, rem_size, rem_used, alpha_p, L, 0);
}

void K18A2Old::CycleBacktrackState::generate_remaining_cycles(int start_idx, const uint8_t* rem, int rem_size,
                                                              bool* rem_used, uint8_t* alpha_p, int L, int depth) {
    if (enum_aborted) return;
    int first_unused = -1;
    for (int i = start_idx; i < rem_size; i++)
        if (!rem_used[i]) { first_unused = i; break; }

    if (first_unused == -1) {
        if (checkPermutationPassed(alpha_p, rem_used, L)) saveAlpha(alpha_p);
        return;
    }

    uint8_t v0_rem = rem[first_unused];
    rem_used[first_unused] = true;

    for (int d = 1; d <= rem_size; d++) {
        if (L % d != 0) continue;        // remaining cycle length must divide L (order stays L)

        if (d == 1) {
            alpha_p[v0_rem] = v0_rem;     // fixed point
            bool defined[K18O_N];
            for (int i = 0; i < K18O_N; i++) defined[i] = true;
            for (int i = 0; i < rem_size; i++) if (!rem_used[i]) defined[rem[i]] = false;
            if (k18o_validateDefinedOrbits(alpha_p, L, defined, F))
                generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L, depth + 1);
        } else {
            // Build the d-cycle node-by-node with SOUND incremental edge-disjointness pruning.
            // This engine assumes rem[i]==i (index == vertex), which generateOrder4Alphas guarantees.
            // When we extend the cycle by a directed edge u->w (alpha(u)=w), the F[0]-edge
            // {u, F[0][u]} maps to {w, alpha(F[0][u])}. If F[0][u] is already placed and that image
            // equals the F[0]-edge {w, F[0][w]}, then alpha(F[0]) shares an edge with F[0] (step-1
            // j=1 fails) -> prune immediately. Each F[0]-edge is checked when its 2nd endpoint is
            // assigned, so coverage is complete; it never prunes a valid alpha (sound). Replaces the
            // old build-whole-cycle-then-validate arrange_stack, which exploded (~1e8/first L-cycle).
            int unused_idx[K18O_N];
            int unused_count = 0;
            for (int j = first_unused + 1; j < rem_size; j++)
                if (!rem_used[j]) unused_idx[unused_count++] = j;     // rem[i]==i -> these are vertices

            if (unused_count >= d - 1) {
                uint8_t cyc[K18O_N];
                cyc[0] = v0_rem;
                bool used_local[K18O_N] = { false };

                auto build = [&](auto& self_fn, int pos) -> void {
                    if (enum_aborted) return;
                    if (pos == d) {                                   // close cycle: edge cyc[d-1] -> cyc[0]
                        const uint8_t u = cyc[d - 1], w = cyc[0], fu = F[0][u];
                        if (rem_used[fu] && alpha_p[fu] == F[0][w]) return;   // edge-disjoint prune
                        alpha_p[u] = w;
                        bool defined[K18O_N];
                        for (int i = 0; i < K18O_N; i++) defined[i] = rem_used[i];
                        if (k18o_validateDefinedOrbits(alpha_p, L, defined, F))
                            generate_remaining_cycles(first_unused + 1, rem, rem_size, rem_used, alpha_p, L, depth + 1);
                        alpha_p[u] = u;                               // undo closing edge
                        return;
                    }
                    const uint8_t u = cyc[pos - 1], fu = F[0][u];
                    for (int t = 0; t < unused_count; t++) {
                        if (used_local[t]) continue;
                        const uint8_t w = (uint8_t)unused_idx[t];
                        if (rem_used[fu] && alpha_p[fu] == F[0][w]) continue;  // edge-disjoint prune
                        if (++enum_steps > 200000000LL) { enum_aborted = true; return; }   // safety cap
                        used_local[t] = true; rem_used[w] = true; cyc[pos] = w; alpha_p[u] = w;
                        self_fn(self_fn, pos + 1);
                        alpha_p[u] = u; rem_used[w] = false; used_local[t] = false;
                        if (enum_aborted) return;
                    }
                };
                build(build, 1);
            }
        }
    }

    rem_used[first_unused] = false;
}

bool K18A2Old::CycleBacktrackState::checkPermutationPassed(const uint8_t* alpha_p, bool* /*used*/, int L) {
    n_reach_leaf++;
    // (1) Fast pre-filter: each power alpha^j(F[0]) (j=1..L/2) must be edge-disjoint
    //     Hamiltonian with F[0].
    uint8_t G[K18O_HALF][K18O_N];
    memcpy(G[0], F[0], K18O_N);
    int limit = L / 2;
    for (int j = 1; j <= limit; j++) {
        apply_perm(G[j - 1], alpha_p, G[j]);
        if (!is_perfect_scalar(G[0], G[j])) return false;
    }
    n_pass_step1++;

    // (2) alpha must respect ALL 3 fixed rows, not just F[0]. For alpha to be an automorphism
    //     of a P1F containing the fixed rows, the image of every fixed row must itself be a
    //     factor of that P1F, hence pairwise-Hamiltonian with all fixed rows, hence present in
    //     the candidate pool (which holds exactly those rows). Necessary AND complete (the pool
    //     is the full set of fixed-row-compatible rows). Replaces the cyclic decompose check and
    //     collapses the candidate-alpha flood that F[0]-only generation produces.
    for (int s = 0; s < K18O_FIXED; s++) {
        uint8_t img[K18O_N];
        apply_perm(self->fixedRows[s].adj, alpha_p, img);
        std::vector<uint8_t> key(img, img + K18O_N);
        if (self->f_map.find(key) == self->f_map.end()) return false;
    }
    n_pass_step2++;
    return true;
}

void K18A2Old::CycleBacktrackState::saveAlpha(const uint8_t* alpha_p) {
    std::array<uint8_t, K18O_N> a;
    std::copy(alpha_p, alpha_p + K18O_N, a.begin());
    valid_alphas.push_back(a);
    total_passed_p1f++;
}

// ---- backtrack-state setup (port of K18A2 setupBacktrackState/setupPairsTable/...) ----
void K18A2Old::setupBacktrackState(CycleBacktrackState& state, int search_type) {
    memcpy(state.F[0], fixedRows[0].adj, K18O_N);
    state.search_type = search_type;
    state.v0 = (search_type == 2) ? 3 : fixedRows[0].adj[0];
    setupPairsTable(state, search_type);
}

void K18A2Old::setupPairsTable(CycleBacktrackState& state, int search_type) {
    for (int i = 0; i < K18O_N; i++) { state.vertex_to_pair[i] = -1; state.vertex_to_pos[i] = -1; }
    fillPairsTable(state, search_type);
    uint8_t partner = fixedRows[0].adj[state.v0];
    state.vertex_to_pair[partner] = -2;
}

void K18A2Old::fillPairsTable(CycleBacktrackState& state, int search_type) {
    int pair_count = 0;
    for (int u = 1; u < K18O_N; u++) {
        if (search_type == 2 && u == 1) continue;
        if (u == state.v0) continue;
        uint8_t v = fixedRows[0].adj[u];
        if (v == state.v0) continue;
        if (state.vertex_to_pair[u] != -1) continue;
        storePair(state, pair_count++, u, v);
    }
}

void K18A2Old::storePair(CycleBacktrackState& state, int pair_idx, uint8_t u, uint8_t v) {
    state.pair_elements[pair_idx][0] = u;
    state.pair_elements[pair_idx][1] = v;
    state.vertex_to_pair[u] = pair_idx; state.vertex_to_pos[u] = 0;
    state.vertex_to_pair[v] = pair_idx; state.vertex_to_pos[v] = 1;
}

// Scalar NP=18 port of K16A2Old::get_transformations_general (no __m128i): align the
// Hamiltonian cycle of (fi,fj) onto that of (fk,fl). Each (dir, offset) alignment gives a
// permutation p with p[fi-vertex] = corresponding fk-vertex (maps fi->fk, fj->fl as factors).
void K18A2Old::alignmentAlphas(const Factor& fi, const Factor& fj, const Factor& fk, const Factor& fl,
                               std::vector<std::array<uint8_t, K18O_N>>& out) {
    CycleUnion src = find_cycles(fi.adj, fj.adj);
    if (src.count != 1 || src.lens[0] != K18O_N) return;   // source union must be one 18-cycle
    CycleUnion tgt = find_cycles(fk.adj, fl.adj);
    if (tgt.count != 1 || tgt.lens[0] != K18O_N) return;   // target union must be one 18-cycle

    // offsets keep cycle-position parity (even with dir=0, odd with dir=1) so fi-edges map to
    // fk-edges and fj-edges to fl-edges (the union cycles alternate the two factors' edges).
    for (int dir = 0; dir < 2; dir++) {
        for (int off = (dir == 0 ? 0 : 1); off < K18O_N; off += 2) {
            std::array<uint8_t, K18O_N> p;
            for (int i = 0; i < K18O_N; i++) {
                int src_pos = (dir == 0) ? (off + i) % K18O_N
                                         : (((off - i) % K18O_N) + K18O_N) % K18O_N;
                p[(uint8_t)src.cycles[0][src_pos]] = (uint8_t)tgt.cycles[0][i];
            }
            out.push_back(p);
        }
    }
}

// Enumerate candidate order-L (=K18O_AUT_ORDER) alphas via factor alignment and complete each
// into P1Fs. OPTION 1 (K16A2Old-faithful): align the source pair (F[0],F[1]) onto each of the 6
// ordered target pairs drawn from the 3 fixed rows, i.e. alphas that permute {F[0],F[1],F[2]}
// among themselves. Fast, but INCOMPLETE: misses automorphisms that map a fixed row to a
// non-fixed factor (same limitation as K16A2Old's order-5/7 gap). Completeness is a follow-up.
void K18A2Old::generateOrder4Alphas(const std::vector<uint8_t>& active,
                                    const std::vector<int>& active_slots) {
    const int L = K18O_AUT_ORDER;

    static const int tpairs[6][2] = { {0,1},{0,2},{1,0},{1,2},{2,0},{2,1} };
    std::vector<std::array<uint8_t, K18O_N>> raw;
    for (auto& pr : tpairs)
        alignmentAlphas(fixedRows[0], fixedRows[1], fixedRows[pr[0]], fixedRows[pr[1]], raw);

    std::vector<std::array<uint8_t, K18O_N>> all_alphas;
    for (auto& a : raw)
        if (has_order(a.data(), L)) all_alphas.push_back(a);

    std::sort(all_alphas.begin(), all_alphas.end());
    all_alphas.erase(std::unique(all_alphas.begin(), all_alphas.end()), all_alphas.end());
    if (m_bPrint) {
        printf("  [K18A2Old] L=%d alignment-gen: %zu raw -> %zu distinct order-%d alphas; completing...\n",
               L, raw.size(), all_alphas.size(), L);
        fflush(stdout);
    }

    std::atomic<size_t> next_alpha{ 0 };
    int nthreads = (kThreads > 1) ? kThreads : 1;
    std::vector<std::thread> workers;
    for (int t = 0; t < nthreads; t++) {
        workers.emplace_back([&]() {
            while (true) {
                size_t i = next_alpha.fetch_add(1);
                if (i >= all_alphas.size()) break;
                completeForAlpha(all_alphas[i].data(), active, active_slots);
            }
        });
    }
    for (auto& w : workers) w.join();
}

void K18A2Old::solve(int mode) {
    const int M = (int)global_pool.size();
    if (m_bPrint) {
        double total_feed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
        printf("K18A2Old: pool=%d candidates. [PERF] feed-phase total=%.1fs | addRow body=%.1fs (%lld calls) | gen(main SW)=%.1fs\n",
               M, total_feed, m_addRowSeconds, m_addRowCalls, total_feed - m_addRowSeconds);
        printf("K18A2Old: arc-consistency...\n");
        fflush(stdout);
    }
    std::vector<uint8_t> active(M, 1);
    arcConsistencyPrune(active);

    std::vector<int> active_slots(M, -1);
    for (int s = 0; s < K18O_SEARCH; s++)
        for (int id : temp_slot_ids[s]) if (active[id]) active_slots[id] = s;

    generateOrder4Alphas(active, active_slots);

    if (m_bPrint) {
        printf("K18A2Old: done. solutions emitted: %llu\n", (unsigned long long)solutions.load());
        fflush(stdout);
    }
}
