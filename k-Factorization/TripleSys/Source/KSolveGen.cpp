#include "KSolveGen.h" 
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>

int g_printCount = 0;

KSolveGen::KSolveGen(const FactorParams& factParam,  bool bipartite, const unsigned char* targetCycles, int cycleCount,
                     int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, 
                     ResultCallback callback, void* cbClassPtr, bool bPrint,
                     const uint8_t* v4Row, const uint8_t* v4Val) 
    : KBase<Mask512>(factParam, bPrint),
      m_bipartite(bipartite), m_cycleCount(cycleCount),
      m_useOptimized(!bipartite && cycleCount == 1),
      m_generateOnTheFly(false)
{
    if (v4Row && v4Val) {
        memcpy(m_v4Row, v4Row, 4);
        memcpy(m_v4Val, v4Val, 4);
    } else {
        memset(m_v4Row, 0, 4);
        memset(m_v4Val, 0, 4);
    }
    // Configure Cycle Structure
    memcpy(m_targetCycles, targetCycles, cycleCount);
    std::sort(m_targetCycles, m_targetCycles + cycleCount);
    
    memset(m_bipartiteFixedPoints, 0, sizeof(m_bipartiteFixedPoints));
    if (m_bipartite) {
        for (int k = 1; k <= m_nPlayers / 2; k++) {
            int fix = 0;
            for (int i = 0; i < m_cycleCount; i++) {
                int L_perm = m_targetCycles[i] / 2;
                if (k % L_perm == 0) fix += L_perm;
            }
            m_bipartiteFixedPoints[k] = fix;
        }
    }

    // Setup Edge ID Table
    memset(edge_id_table, -1, sizeof(edge_id_table));
    int id_cnt = 0;
    if (bipartite) {
        for (int i = 0; i < m_nPlayers; i += 2) { // Even nodes (set A)
            for (int j = 1; j < m_nPlayers; j += 2) { // Odd nodes (set B)
                edge_id_table[i][j] = id_cnt++;
                edge_id_table[j][i] = edge_id_table[i][j];
            }
        }
    } else {
        for (int i = 0; i < m_nPlayers; i++) {
            for (int j = i + 1; j < m_nPlayers; j++) {
                edge_id_table[i][j] = id_cnt++;
                edge_id_table[j][i] = edge_id_table[i][j];
            }
        }
    }

    total_edges_mask.clear();
    for (int i = 0; i < id_cnt; i++) total_edges_mask.set_bit(i);

    init(fixed3RowsIndex, kThreads, first3Rows, callback, cbClassPtr);
}

void KSolveGen::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    this->start_time = std::chrono::steady_clock::now();
    this->last_report_time = start_time;

    fixedEdgesMask.clear();
    for (int s = 0; s < 3; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * m_nPlayers, m_nPlayers);
        for (int i = 0; i < m_nPlayers; i += 2) {
            uint8_t u = fixedRows[s].src[i];
            uint8_t v = fixedRows[s].src[i + 1];
            fixedRows[s].adj[u] = v;
            fixedRows[s].adj[v] = u;
            int eid = edge_id_table[u][v];
            if (eid != -1) fixedEdgesMask.set_bit(eid);
        }
        fixedRows[s].fs = get_fast_sorted(fixedRows[s].adj);
        fixed_packed[s] = pack_factor_adj(fixedRows[s].adj);
        if (0 && m_bPrint) {
            printf("KSolveGen: FixedRow %d: adj[0]=%d, adj[1]=%d, adj[2]=%d\n", s, fixedRows[s].adj[0], fixedRows[s].adj[1], fixedRows[s].adj[2]);
        }
    }
    r1_can = fixedRows[0].fs;
    r2_can = fixedRows[1].fs;

    m_presolveDone = false;
    
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    
    // Reset transient state
    m_resultsFound = 0;
    m_rootsDone = 0;
    global_pool.clear();
    packed_pool.clear();
    f_map.clear();
    results_to_sort.clear();
    m_referenceSets.clear();
    m_parityReportDone = false;
    KBase::init();

    // Precalculate transformations for fixed row pairs
    get_transformations(fixedRows[0], fixedRows[1], fixed_trans[0]);
    get_transformations(fixedRows[0], fixedRows[2], fixed_trans[1]);
    get_transformations(fixedRows[1], fixedRows[2], fixed_trans[2]);
}

KSolveGen::PackedAdj KSolveGen::pack_factor_adj(const uint8_t* adj) {
    PackedAdj p;
    p.edge_mask.clear();
    memcpy(p.adj, adj, 32);
    for (int i = 0; i < m_nPlayers; i++) {
        if (i < adj[i]) {
            int eid = edge_id_table[i][adj[i]];
            if (eid != -1) p.edge_mask.set_bit(eid);
        }
    }
    if (m_bipartite) {
        alignas(16) uint8_t e2o[16] = { 0 }, o2e[16] = { 0 };
        for (int j = 0; j < m_nPlayers / 2; j++) {
            e2o[j] = (adj[j * 2] - 1) / 2;
            o2e[j] = adj[j * 2 + 1] / 2;
        }
        p.v_e2o = _mm_loadu_si128((const __m128i*)e2o);
        p.v_o2e = _mm_loadu_si128((const __m128i*)o2e);
    }
    return p;
}

bool KSolveGen::addRow(int iRow, const unsigned char* source) {
    // Sanity Check: Ensure Node 0's partner in 'source' matches iRow
    uint8_t partnerOf0 = 0xFF;
    for (int i = 0; i < m_nPlayers; i += 2) {
        if (source[i] == 0) partnerOf0 = source[i+1];
        else if (source[i+1] == 0) partnerOf0 = source[i];
        if (partnerOf0 != 0xFF) break;
    }

    if (partnerOf0 != (uint8_t)iRow) {
        // Only error if it's the legacy capture phase (we want to know the legacy format)
        // If it's on-the-fly generation, we trust rowK
        if (!m_generateOnTheFly) {
            printf("CRITICAL ERROR: addRow(iRow=%d) called, but Node 0's partner in source is %d!\n", 
                   iRow, partnerOf0);
            printf("Source Data (first 6 bytes): %d %d %d %d %d %d\n",
                   source[0], source[1], source[2], source[3], source[4], source[5]);
            exit(1);
        }
    }

    int slot = m_bipartite ? (iRow / 2 - 3) : (iRow - 4);
    if (slot < 0 || slot >= 32) return true; // Skip rows outside search range (fixed rows)

    Factor f;
    memcpy(f.src, source, m_nPlayers);
    for (int i = 0; i < m_nPlayers; i += 2) {
        uint8_t u = f.src[i]; uint8_t v = f.src[i + 1];
        f.adj[u] = v; f.adj[v] = u;
    }

    PackedAdj pf = pack_factor_adj(f.adj);
    if (pf.edge_mask.has_overlap(fixedEdgesMask)) return false;

    for (int i = 0; i < 3; i++) {
        if (!is_perfect_packed(pf, fixed_packed[i])) return false;
    }

    // Apply v4-Row constraints (Player 1 partner)
    for (int i = 0; i < 4; i++) {
        if (m_v4Row[i] == iRow && m_v4Val[i] != 0) {
            if (f.adj[1] != m_v4Val[i]) return false;
        }
    }

    f.edge_mask = pf.edge_mask;
    f.fs = get_fast_sorted(f.adj);

    std::vector<uint8_t> key(f.adj, f.adj + m_nPlayers);
    int factor_id;
    {
        std::lock_guard<std::mutex> lock(pool_mutex);
        auto it = f_map.find(key);
        if (it == f_map.end()) {
            factor_id = (int)global_pool.size();
            if (factor_id >= KGEN_M_MAX) return false;
            f_map[key] = factor_id;
            global_pool.push_back(f);
            packed_pool.push_back(pf);
        } else {
            factor_id = it->second;
        }

        addFactorToSlot(slot, factor_id);
    }
    return true;
}

void KSolveGen::parallel_for(int start, int end, std::function<void(int, int)> task, int grain_size) {
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
                for (int j = i; j < limit; ++j) task(j, t);
            }
        });
    }
    for (auto& w : workers) w.join();
}

void KSolveGen::generateSlotCandidates(int slotIndex, int startIdx, int endIdx) {
    int P0;
    if (m_bipartite) {
        P0 = (slotIndex + 3) * 2 + 1;
    } else {
        P0 = slotIndex + 4;
    }
    
    //if (m_bPrint) printf("KSolveGen: Slot %d, partner P0=%d starting...\n", slotIndex, P0);

    uint8_t factor[32];
    memset(factor, 0, 32);
    factor[0] = 0;
    factor[1] = (uint8_t)P0;
    
    uint64_t usedMask = (1ULL << 0) | (1ULL << P0);
    int start_depth = 0;

    // Check for v4-Row constraints (Partner of Node 1)
    for (int i = 0; i < 4; i++) {
        if (m_v4Row[i] == P0 && m_v4Val[i] != 0) {
            uint8_t v1 = m_v4Val[i];
            if (!(usedMask & (1ULL << 1)) && !(usedMask & (1ULL << v1))) {
                factor[2] = 1; factor[3] = v1;
                usedMask |= (1ULL << 1) | (1ULL << v1);
                start_depth = 1;
            }
            break;
        }
    }
    
    // The generator will try all valid matchings and call addRow
    generateRecursive(start_depth, usedMask, factor, slotIndex, startIdx, endIdx, P0);
}

void KSolveGen::generateRecursive(int depth, uint64_t usedMask, uint8_t* factor, int slotIndex, int startPairIdx, int endPairIdx, int rowK) {
    if (depth == m_nPlayers / 2 - 1) {
        addRow(rowK, factor);
        return;
    }

    unsigned long u;
    if (!_BitScanForward64(&u, ~usedMask)) return;

    for (int v = 0; v < m_nPlayers; ++v) {
        if (usedMask & (1ULL << v)) continue;
        
        if (m_bipartite) {
            if ((u % 2) == (v % 2)) continue;
        } else {
            if (v <= (int)u) continue;
        }

        int eid = edge_id_table[u][v];
        if (eid != -1 && fixedEdgesMask.m[eid >> 6] & (1ULL << (eid & 63))) continue;

        factor[(depth + 1) * 2] = (uint8_t)u;
        factor[(depth + 1) * 2 + 1] = (uint8_t)v;
        generateRecursive(depth + 1, usedMask | (1ULL << u) | (1ULL << v), factor, slotIndex, startPairIdx, endPairIdx, rowK);
    }
}
