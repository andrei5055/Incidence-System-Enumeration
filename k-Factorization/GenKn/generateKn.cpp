#include "generateKn.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <set>
#include <mutex>
#include <thread>
#include <atomic>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <filesystem>
#include <immintrin.h>

// Toggle to run only a fast subset (excluding slow p=2 cases where k < 4)
const bool FAST_SUBSET_ONLY = false;

// Search limit configurations for benchmarking baselines
const bool LIMIT_SEARCH = false;
const int SEARCH_LIMIT_CONFIGS = 1000;
static std::atomic<bool> stop_search(false);

// Matrix status definitions
const int UNASSIGNED = -1;

struct FactorOrbit {
    int representative;
    int size;
    std::vector<int> factors;
};

struct EdgePlacement {
    int u, w, f;
};

struct SymmetryComb {
    int p, k, m;
};

struct ThreadProgress {
    int p = 0, k = 0, m = 0;
    int current_t = 0;
    std::chrono::steady_clock::time_point start_time;
    bool active = false;
};
static std::vector<ThreadProgress> global_thread_progress;
static thread_local int current_thread_idx = -1;
static int completed_combs_count = 0;
static double total_completed_combs_time = 0.0;

// Global mutexes for synchronization
std::mutex unique_matrices_mutex;

// Helper to get elapsed time in [HH:MM:SS] format
static std::string get_timestamp() {
    static const auto start_time = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    auto elapsed_secs = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
    long long hours = elapsed_secs / 3600;
    long long minutes = (elapsed_secs % 3600) / 60;
    long long seconds = elapsed_secs % 60;
    std::ostringstream oss;
    oss << "[" << std::setfill('0')
        << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":"
        << std::setw(2) << seconds << "] ";
    return oss.str();
}

// ============================================================================
// PART 1: FASTER ISOMORPHISM CHECKER
// ============================================================================
bool backtrack_isomorphism(
    int pair_idx,
    int num_pairs,
    int N,
    int num_factors,
    const int M1[32][32],
    const int M2[32][32],
    const std::pair<int, int> pairs2[32],
    int pi[32],
    bool pair2_used[32],
    int sigma[32],
    int sigma_inv[32]
) {
    if (pair_idx == num_pairs) {
        return true;
    }

    int u1 = 2 * pair_idx;
    int v1 = 2 * pair_idx + 1;

    for (int j = 0; j < num_pairs; ++j) {
        if (pair2_used[j]) continue;

        int u2 = pairs2[j].first;
        int v2 = pairs2[j].second;

        // Try mapping u1 -> u2, v1 -> v2
        for (int swap = 0; swap < 2; ++swap) {
            int map_u2 = (swap == 0) ? u2 : v2;
            int map_v2 = (swap == 0) ? v2 : u2;

            pi[u1] = map_u2;
            pi[v1] = map_v2;

            // Check compatibility with all vertices mapped in pairs < pair_idx
            bool ok = true;
            int sigma_changes_g[32];
            int sigma_changes_h[32];
            int num_changes = 0;

            auto check_edge = [&](int src, int dest) {
                int g = M1[src][dest];
                int h = M2[pi[src]][pi[dest]];
                if (sigma[g] != -1) {
                    if (sigma[g] != h) return false;
                } else {
                    if (sigma_inv[h] != -1) return false;
                    sigma[g] = h;
                    sigma_inv[h] = g;
                    sigma_changes_g[num_changes] = g;
                    sigma_changes_h[num_changes] = h;
                    num_changes++;
                }
                return true;
            };

            for (int prev_u = 0; prev_u < u1; ++prev_u) {
                if (!check_edge(prev_u, u1) || !check_edge(prev_u, v1)) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                pair2_used[j] = true;
                if (backtrack_isomorphism(pair_idx + 1, num_pairs, N, num_factors, M1, M2, pairs2, pi, pair2_used, sigma, sigma_inv)) {
                    return true;
                }
                pair2_used[j] = false;
            }

            // Undo sigma changes
            for (int c = 0; c < num_changes; ++c) {
                sigma[sigma_changes_g[c]] = -1;
                sigma_inv[sigma_changes_h[c]] = -1;
            }
        }
    }
    return false;
}

bool are_isomorphic(int N, const std::vector<std::vector<int>>& F1, const std::vector<std::vector<int>>& F2) {
    int num_factors = N - 1;
    int num_pairs = N / 2;
    int M1[32][32];
    int M2[32][32];
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            M1[i][j] = -1;
            M2[i][j] = -1;
        }
    }

    for (int f = 0; f < num_factors; ++f) {
        for (int i = 0; i < N; i += 2) {
            int u1 = F1[f][i], v1 = F1[f][i+1];
            M1[u1][v1] = f; M1[v1][u1] = f;
            int u2 = F2[f][i], v2 = F2[f][i+1];
            M2[u2][v2] = f; M2[v2][u2] = f;
        }
    }

    // We try to map factor 0 of F1 to factor f2 of F2
    for (int f2 = 0; f2 < num_factors; ++f2) {
        std::pair<int, int> pairs2[32];
        for (int i = 0; i < N; i += 2) {
            pairs2[i/2] = {F2[f2][i], F2[f2][i+1]};
        }

        int pi[32];
        bool pair2_used[32];
        int sigma[32];
        int sigma_inv[32];
        for (int i = 0; i < 32; ++i) {
            pi[i] = -1;
            sigma[i] = -1;
            sigma_inv[i] = -1;
            pair2_used[i] = false;
        }

        sigma[0] = f2;
        sigma_inv[f2] = 0;

        if (backtrack_isomorphism(0, num_pairs, N, num_factors, M1, M2, pairs2, pi, pair2_used, sigma, sigma_inv)) {
            return true;
        }
    }
    return false;
}

bool is_isomorphic_to_any(int N, const std::vector<std::vector<int>>& candidate, const std::set<std::vector<std::vector<int>>>& unique_matrices) {
    for (const auto& existing : unique_matrices) {
        if (are_isomorphic(N, candidate, existing)) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// PART 2: CANONICALIZATION HELPER
// ============================================================================
std::vector<std::vector<int>> canonicalize_matrix(int N, const int matrix[32][32]) {
    int num_factors = N - 1;
    std::vector<std::vector<int>> sorted_rows;
    for (int j = 0; j < num_factors; j++) {
        std::vector<int> row_pairs;
        std::vector<bool> pair_used(N, false);
        for (int vertex = 0; vertex < N; ++vertex) {
            if (pair_used[vertex]) continue;
            int neighbor = matrix[j][vertex];

            if (vertex < neighbor) {
                row_pairs.push_back(vertex);
                row_pairs.push_back(neighbor);
            }
            else {
                row_pairs.push_back(neighbor);
                row_pairs.push_back(vertex);
            }
            pair_used[vertex] = true;
            pair_used[neighbor] = true;
        }
        sorted_rows.push_back(row_pairs);
    }

    std::sort(sorted_rows.begin(), sorted_rows.end());

    std::vector<int> pi(N);
    for (int i = 0; i < N; ++i) {
        pi[sorted_rows[0][i]] = i;
    }

    std::vector<std::vector<int>> permuted_rows;
    for (const auto& row : sorted_rows) {
        std::vector<std::pair<int, int>> pairs;
        for (size_t i = 0; i < row.size(); i += 2) {
            int u = pi[row[i]];
            int v = pi[row[i + 1]];
            if (u > v) std::swap(u, v);
            pairs.push_back({u, v});
        }
        std::sort(pairs.begin(), pairs.end());

        std::vector<int> flat_row;
        for (const auto& p : pairs) {
            flat_row.push_back(p.first);
            flat_row.push_back(p.second);
        }
        permuted_rows.push_back(flat_row);
    }

    std::sort(permuted_rows.begin(), permuted_rows.end());
    return permuted_rows;
}

// ============================================================================
// PART 3: BACKTRACKING ORBIT SOLVER AND LOOK-AHEAD SUBCYCLE PRUNING
// ============================================================================
bool forms_subcycle_general(int N, int u, int w, int f, const int matrix[32][32], const int factor_assigned_mask[32]) {
    int num_factors = N - 1;
    const int* row_f = &matrix[f][0];
    int mask_uw = (1 << u) | (1 << w);
    for (int g = 0; g < num_factors; ++g) {
        if (g == f) continue;
        if ((factor_assigned_mask[g] & mask_uw) != mask_uw) continue;
        
        const int* row_g = &matrix[g][0];
        int curr = row_g[w];
        int len = 2;
        while (curr >= 0) {
            curr = row_f[curr];
            if (curr < 0) break;
            curr = row_g[curr];
            len += 2;
            if (curr == u) {
                if (len < N) return true;
                break;
            }
        }
    }
    return false;
}

void backtrack_matching(
    int orbit_idx,
    int v,
    int N,
    int p,
    const std::vector<FactorOrbit>& orbits,
    const int alpha_pow_col[32][8],
    const int beta_pow_col_shl[32][8],
    int matrix[32][32],
    int edge_used[32][32],
    int factor_assigned_mask[32],
    std::atomic<long long>& configs_checked,
    std::set<std::vector<std::vector<int>>>& unique_matrices
) {
    if (stop_search) return;

    if (orbit_idx == (int)orbits.size()) {
        long long current_val = ++configs_checked;
        if (LIMIT_SEARCH && current_val >= SEARCH_LIMIT_CONFIGS) {
            stop_search = true;
        }
        std::vector<std::vector<int>> candidate = canonicalize_matrix(N, matrix);
        
        std::lock_guard<std::mutex> lock(unique_matrices_mutex);
        if (!is_isomorphic_to_any(N, candidate, unique_matrices)) {
            unique_matrices.insert(candidate);
            int current_idx = (int)unique_matrices.size();
            
            // Print status line with timestamp
            std::cout << get_timestamp() << "[PROGRESS] Found unique P1F #" << current_idx << std::endl;
            
            // Append to file results/N/P0000000001.txt (exact current format, no timestamps)
            std::string file_path = "results/" + std::to_string(N) + "/P0000000001.txt";
            std::FILE* result_file = std::fopen(file_path.c_str(), "a");
            if (result_file) {
                std::fprintf(result_file, "        %5d: | Aut(M) | = 1\n", current_idx);
                for (const auto& row : candidate) {
                    std::fprintf(result_file, " \"");
                    for (size_t i = 0; i < row.size(); i += 2) {
                        std::fprintf(result_file, " %3d%3d", row[i], row[i + 1]);
                    }
                    std::fprintf(result_file, " \"\n");
                }
                std::fprintf(result_file, "\n\n\n\n");
                std::fflush(result_file);
                std::fclose(result_file);
            }
            
            // Duplicate on screen (WITHOUT timestamps on result lines)
            printf("        %5d: | Aut(M) | = 1\n", current_idx);
            for (const auto& row : candidate) {
                printf(" \"");
                for (size_t i = 0; i < row.size(); i += 2) {
                    printf(" %3d%3d", row[i], row[i + 1]);
                }
                printf(" \"\n");
            }
            printf("\n\n\n\n");
            fflush(stdout);
        }
        return;
    }

    int f = orbits[orbit_idx].representative;

    if (v == N) {
        backtrack_matching(orbit_idx + 1, 0, N, p, orbits, alpha_pow_col, beta_pow_col_shl, matrix, edge_used, factor_assigned_mask, configs_checked, unique_matrices);
        return;
    }

    if (matrix[f][v] != UNASSIGNED) {
        backtrack_matching(orbit_idx, v + 1, N, p, orbits, alpha_pow_col, beta_pow_col_shl, matrix, edge_used, factor_assigned_mask, configs_checked, unique_matrices);
        return;
    }

    const __m256i v_unassigned = _mm256_set1_epi32(-1);
    const __m256i v_zero       = _mm256_setzero_si256();
    const int target_mask = (1 << (4 * p)) - 1;

    // Load symmetry mapping matrices for vertex mapping u and factor f (constant in the t loop)
    __m256i u_map = _mm256_load_si256((const __m256i*)&alpha_pow_col[v][0]);
    __m256i f_map_shl = _mm256_load_si256((const __m256i*)&beta_pow_col_shl[f][0]);

    for (int t = v + 1; t < N; ++t) {
        if (orbit_idx == 0 && v == 0) {
            std::lock_guard<std::mutex> lock(unique_matrices_mutex);
            if (current_thread_idx >= 0 && current_thread_idx < (int)global_thread_progress.size()) {
                global_thread_progress[current_thread_idx].current_t = t;
            }
        }
        // Load symmetry mapping matrices for vertex mapping w
        __m256i w_map = _mm256_load_si256((const __m256i*)&alpha_pow_col[t][0]);

        // Standardize edge endpoints u_min <= w_max (equivalent to standard std::swap in scalar)
        __m256i u_min = _mm256_min_epi32(u_map, w_map);
        __m256i w_max = _mm256_max_epi32(u_map, w_map);

        // Compute flat 1D offsets for matrix[32][32] and edge_used[32][32]
        // Offset = Row * 32 + Col
        __m256i idx_u_f  = _mm256_add_epi32(f_map_shl, u_min);
        __m256i idx_w_f  = _mm256_add_epi32(f_map_shl, w_max);
        __m256i idx_edge = _mm256_add_epi32(_mm256_slli_epi32(u_min, 5), w_max);

        // Vectorized gather reads
        __m256i val_u_f  = _mm256_i32gather_epi32((const int*)matrix, idx_u_f, 4);
        __m256i val_w_f  = _mm256_i32gather_epi32((const int*)matrix, idx_w_f, 4);
        __m256i val_edge = _mm256_i32gather_epi32((const int*)edge_used, idx_edge, 4);

        // Check if endpoints are unassigned in factor f_idx (val_u_f == UNASSIGNED && val_w_f == UNASSIGNED)
        __m256i cmp_u = _mm256_cmpeq_epi32(val_u_f, v_unassigned);
        __m256i cmp_w = _mm256_cmpeq_epi32(val_w_f, v_unassigned);
        // Check if edge (u_min, w_max) is unused (val_edge == 0)
        __m256i cmp_edge = _mm256_cmpeq_epi32(val_edge, v_zero);

        // Combine logic (we want all of them to be true)
        __m256i combined = _mm256_and_si256(_mm256_and_si256(cmp_u, cmp_w), cmp_edge);
        int mask = _mm256_movemask_epi8(combined);

        // Prune search tree early if any lane for the active prime order is invalid
        if ((mask & target_mask) != target_mask) {
            continue;
        }

        // Store the SIMD results to aligned stack arrays for final deduplication
        alignas(32) int u_arr[8];
        alignas(32) int w_arr[8];
        alignas(32) int f_arr[8];
        _mm256_store_si256((__m256i*)u_arr, u_min);
        _mm256_store_si256((__m256i*)w_arr, w_max);
        _mm256_store_si256((__m256i*)f_arr, f_map_shl);

        EdgePlacement placements[8];
        int num_placements = 0;
        for (int i = 0; i < p; ++i) {
            // Recover factor index by dividing f_arr[i] by 32
            EdgePlacement ep = { u_arr[i], w_arr[i], f_arr[i] >> 5 };
            
            bool already_exists = false;
            for (int k = 0; k < num_placements; ++k) {
                if (placements[k].u == ep.u && placements[k].w == ep.w && placements[k].f == ep.f) {
                    already_exists = true;
                    break;
                }
            }
            if (!already_exists) {
                placements[num_placements++] = ep;
            }
        }

        bool valid = true;
        // Check conflicts within the placements (self-conflicts)
        for (int i = 0; i < num_placements; ++i) {
            for (int j = i + 1; j < num_placements; ++j) {
                if (placements[i].u == placements[j].u && placements[i].w == placements[j].w) {
                    valid = false;
                    break;
                }
                if (placements[i].f == placements[j].f) {
                    if (placements[i].u == placements[j].u || placements[i].u == placements[j].w ||
                        placements[i].w == placements[j].u || placements[i].w == placements[j].w) {
                        valid = false;
                        break;
                    }
                }
            }
            if (!valid) break;
        }

        if (valid) {
            // Temporarily place ALL edges of the orbit
            for (int i = 0; i < num_placements; ++i) {
                const auto& ep = placements[i];
                matrix[ep.f][ep.u] = ep.w;
                matrix[ep.f][ep.w] = ep.u;
                factor_assigned_mask[ep.f] |= (1 << ep.u) | (1 << ep.w);
            }

            // Check look-ahead subcycles on the fly using forms_subcycle_general
            bool subcycle = false;
            for (int i = 0; i < num_placements; ++i) {
                const auto& ep = placements[i];
                if (forms_subcycle_general(N, ep.u, ep.w, ep.f, matrix, factor_assigned_mask)) {
                    subcycle = true;
                    break;
                }
            }

            // Restore ALL edges of the orbit
            for (int i = 0; i < num_placements; ++i) {
                const auto& ep = placements[i];
                matrix[ep.f][ep.u] = UNASSIGNED;
                matrix[ep.f][ep.w] = UNASSIGNED;
                factor_assigned_mask[ep.f] &= ~((1 << ep.u) | (1 << ep.w));
            }

            if (subcycle) continue;

            // Place edges
            for (int i = 0; i < num_placements; ++i) {
                const auto& ep = placements[i];
                matrix[ep.f][ep.u] = ep.w;
                matrix[ep.f][ep.w] = ep.u;
                edge_used[ep.u][ep.w] = 1;
                factor_assigned_mask[ep.f] |= (1 << ep.u) | (1 << ep.w);
            }

            backtrack_matching(orbit_idx, v + 1, N, p, orbits, alpha_pow_col, beta_pow_col_shl, matrix, edge_used, factor_assigned_mask, configs_checked, unique_matrices);

            for (int i = 0; i < num_placements; ++i) {
                const auto& ep = placements[i];
                matrix[ep.f][ep.u] = UNASSIGNED;
                matrix[ep.f][ep.w] = UNASSIGNED;
                edge_used[ep.u][ep.w] = 0;
                factor_assigned_mask[ep.f] &= ~((1 << ep.u) | (1 << ep.w));
            }
        }
    }
}

// ============================================================================
// PART 4: SYMMETRY COMBINATION GENERATION VALIDATION
// ============================================================================
bool is_valid_combination(int N, int p, int k, int m) {
    int F = N - k * p;
    if (F < 0) return false;
    if (F > 1 && F % 2 != 0) return false;
    int s = (N - 1) - m * p;
    if (s < 0) return false;
    if (F == 0) {
        if (s != 0) return false;
    } else if (F == 1) {
        if (s > 1) return false;
    } else {
        if (s > F - 1) return false;
    }
    return true;
}

// ============================================================================
// PART 5: MULTITHREADED WORKER FUNCTION
// ============================================================================
void worker(
    int thread_idx,
    int N,
    const std::vector<SymmetryComb>& combinations,
    std::atomic<size_t>& next_comb_idx,
    std::atomic<long long>& total_configurations_checked,
    std::set<std::vector<std::vector<int>>>& unique_matrices
) {
    current_thread_idx = thread_idx;
    while (true) {
        size_t idx = next_comb_idx++;
        if (idx >= combinations.size()) break;
        
        SymmetryComb comb = combinations[idx];
        int p = comb.p;
        int k = comb.k;
        int m = comb.m;
        
        {
            std::lock_guard<std::mutex> lock(unique_matrices_mutex);
            global_thread_progress[thread_idx].p = p;
            global_thread_progress[thread_idx].k = k;
            global_thread_progress[thread_idx].m = m;
            global_thread_progress[thread_idx].current_t = 0;
            global_thread_progress[thread_idx].start_time = std::chrono::steady_clock::now();
            global_thread_progress[thread_idx].active = true;

            std::cout << get_timestamp() << "[INFO] Running search for p=" << p << ", k=" << k << ", m=" << m 
                      << " on thread " << thread_idx << std::endl;
        }

        // Setup vertex permutation alpha
        std::vector<int> alpha(N);
        for (int j = 0; j < k; ++j) {
            for (int i = 0; i < p; ++i) {
                alpha[j * p + i] = j * p + (i + 1) % p;
            }
        }
        for (int v = k * p; v < N; ++v) {
            alpha[v] = v;
        }

        // Setup factor permutation beta
        int num_factors = N - 1;
        std::vector<int> beta(num_factors);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < p; ++i) {
                beta[j * p + i] = j * p + (i + 1) % p;
            }
        }
        for (int f = m * p; f < num_factors; ++f) {
            beta[f] = f;
        }

        // Precompute powers of alpha and beta in column-major with 32-byte alignment
        alignas(32) int alpha_pow_col[32][8] = {0};
        alignas(32) int beta_pow_col_shl[32][8] = {0};

        for (int v = 0; v < N; ++v) {
            int curr = v;
            for (int i = 0; i < p; ++i) {
                alpha_pow_col[v][i] = curr;
                curr = alpha[curr];
            }
            // Pad remaining lanes with safe zero-values
            for (int i = p; i < 8; ++i) {
                alpha_pow_col[v][i] = 0;
            }
        }
        for (int f = 0; f < num_factors; ++f) {
            int curr = f;
            for (int i = 0; i < p; ++i) {
                beta_pow_col_shl[f][i] = curr * 32;
                curr = beta[curr];
            }
            for (int i = p; i < 8; ++i) {
                beta_pow_col_shl[f][i] = 0;
            }
        }

        // Construct factor orbits and sort them (size 1 first, then size p)
        std::vector<FactorOrbit> orbits;
        std::vector<bool> factor_visited(num_factors, false);
        for (int f = 0; f < num_factors; ++f) {
            if (factor_visited[f]) continue;
            std::vector<int> orb;
            int curr = f;
            do {
                orb.push_back(curr);
                factor_visited[curr] = true;
                curr = beta[curr];
            } while (curr != f);
            orbits.push_back({f, (int)orb.size(), orb});
        }

        std::sort(orbits.begin(), orbits.end(), [](const FactorOrbit& a, const FactorOrbit& b) {
            return a.size < b.size;
        });

        alignas(32) int matrix[32][32];
        alignas(32) int edge_used[32][32];
        alignas(32) int factor_assigned_mask[32];
        for (int i = 0; i < 32; ++i) {
            factor_assigned_mask[i] = 0;
            for (int j = 0; j < 32; ++j) {
                matrix[i][j] = UNASSIGNED;
                edge_used[i][j] = 0;
            }
        }

        backtrack_matching(0, 0, N, p, orbits, alpha_pow_col, beta_pow_col_shl, matrix, edge_used, factor_assigned_mask, total_configurations_checked, unique_matrices);
        
        auto comb_end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration<double>(comb_end - global_thread_progress[thread_idx].start_time).count();
        {
            std::lock_guard<std::mutex> lock(unique_matrices_mutex);
            completed_combs_count++;
            total_completed_combs_time += duration;
        }
    }
    {
        std::lock_guard<std::mutex> lock(unique_matrices_mutex);
        global_thread_progress[thread_idx].active = false;
    }
}

// ============================================================================
// PIPELINE WRAPPER (MAIN EXECUTION ENTRY)
// ============================================================================
int generateKn(int N, int target_p, int target_k, int target_m, int num_threads, int verbose) {
    // Initialize start time for elapsed tracking
    get_timestamp();

    {
        std::lock_guard<std::mutex> lock(unique_matrices_mutex);
        global_thread_progress.clear();
        global_thread_progress.resize(num_threads);
        completed_combs_count = 0;
        total_completed_combs_time = 0.0;
    }

    // Open on start, truncate to create a new file, fault if can't open
    std::string dir_path = "results/" + std::to_string(N);
    try {
        std::filesystem::create_directories(dir_path);
    } catch (...) {
        std::cerr << get_timestamp() << "[ERROR] Could not create directory " << dir_path << "! Exiting." << std::endl;
        exit(1);
    }
    std::string file_path = dir_path + "/P0000000001.txt";

    {
        std::FILE* test_file = std::fopen(file_path.c_str(), "w");
        if (!test_file) {
            std::cerr << get_timestamp() << "[ERROR] Could not open " << file_path << " for writing! Exiting." << std::endl;
            exit(1);
        }
        std::fclose(test_file);
    }

    std::cout << get_timestamp() << "--- STARTING PROGRAMMATIC FACTORIZATION SCHEME FOR K_" << N << " ---" << std::endl;

    stop_search = false;
    std::atomic<long long> total_configurations_checked(0);
    std::set<std::vector<std::vector<int>>> unique_matrices;

    // Collect all valid combinations of (p, k, m)
    std::vector<SymmetryComb> combinations;
    std::vector<int> primes = { 7, 5, 3, 2 };
    for (int p : primes) {
        for (int k = 1; k * p <= N; ++k) {
            for (int m = 0; m * p <= N - 1; ++m) {
                if (is_valid_combination(N, p, k, m)) {
                    if (target_p != 0) {
                        if (p == target_p && k == target_k && m == target_m) {
                            combinations.push_back({ p, k, m });
                        }
                    } else {
                        if (FAST_SUBSET_ONLY && p == 2 && k < 4) {
                            continue;
                        }
                        combinations.push_back({ p, k, m });
                    }
                }
            }
        }
    }

    if (target_p != 0 && combinations.empty()) {
        std::cerr << get_timestamp() << "[ERROR] Specified target combination (p=" << target_p 
                  << ", k=" << target_k << ", m=" << target_m << ") is invalid for N=" << N << "." << std::endl;
        return 1;
    }

    std::cout << get_timestamp() << "[INFO] Found " << combinations.size() << " valid prime-order symmetry targets." << std::endl;

    std::atomic<size_t> next_comb_idx(0);
    std::cout << get_timestamp() << "[INFO] Spawning " << num_threads << " worker threads." << std::endl;

    std::atomic<bool> monitor_stop(false);
    std::thread monitor_thread([&]() {
        while (!monitor_stop) {
            std::this_thread::sleep_for(std::chrono::seconds(10));
            if (monitor_stop) break;
            
            std::lock_guard<std::mutex> lock(unique_matrices_mutex);
            
            auto now = std::chrono::steady_clock::now();
            double min_left_pct = 100.0;
            bool any_active = false;
            bool has_pct = false;
            
            for (size_t i = 0; i < global_thread_progress.size(); ++i) {
                const auto& tp = global_thread_progress[i];
                if (tp.active) {
                    any_active = true;
                    int curr_t = tp.current_t;
                    if (curr_t > 0) {
                        double pct_left = (double)(N - curr_t) * 100.0 / (N - 1);
                        if (!has_pct || pct_left < min_left_pct) {
                            min_left_pct = pct_left;
                            has_pct = true;
                        }
                    }
                }
            }
            
            std::cout << get_timestamp() << "[MONITOR] Configs checked: " << total_configurations_checked 
                      << " | Unique P1Fs found: " << unique_matrices.size();
            
            if (any_active) {
                if (has_pct) {
                    std::cout << " | Left: " << std::fixed << std::setprecision(1) << min_left_pct << "%";
                } else {
                    std::cout << " | Left: 100.0%";
                }
            }
            std::cout << std::endl;
            
            if (verbose) {
                for (size_t i = 0; i < global_thread_progress.size(); ++i) {
                    const auto& tp = global_thread_progress[i];
                    if (tp.active) {
                        int curr_t = tp.current_t;
                        if (curr_t > 0) {
                            double progress_pct = 100.0 * (curr_t - 1) / (N - 1);
                            auto elapsed_secs = std::chrono::duration_cast<std::chrono::seconds>(now - tp.start_time).count();
                            
                            std::cout << "  Thread " << i << " (p=" << tp.p << ",k=" << tp.k << ",m=" << tp.m 
                                      << "): branch " << curr_t << "/" << (N - 1) 
                                      << " (" << std::fixed << std::setprecision(1) << progress_pct << "% done)";
                            
                            if (curr_t > 1) {
                                double avg_sec_per_branch = (double)elapsed_secs / (curr_t - 1);
                                double remaining_secs = avg_sec_per_branch * (N - curr_t);
                                
                                long long eta_h = (long long)remaining_secs / 3600;
                                long long eta_m = ((long long)remaining_secs % 3600) / 60;
                                long long eta_s = (long long)remaining_secs % 60;
                                
                                std::cout << " | Elapsed: " << elapsed_secs << "s | ETA: "
                                          << std::setfill('0') << std::setw(2) << eta_h << ":"
                                          << std::setfill('0') << std::setw(2) << eta_m << ":"
                                          << std::setfill('0') << std::setw(2) << eta_s;
                            } else {
                                std::cout << " | Elapsed: " << elapsed_secs << "s | ETA: Estimating...";
                            }
                            std::cout << std::endl;
                        } else {
                            std::cout << "  Thread " << i << " (p=" << tp.p << ",k=" << tp.k << ",m=" << tp.m 
                                      << "): Starting..." << std::endl;
                        }
                    }
                }
            }
        }
    });

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(worker, i, N, std::ref(combinations), std::ref(next_comb_idx), 
                                      std::ref(total_configurations_checked), std::ref(unique_matrices)));
    }

    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }

    monitor_stop = true;
    if (monitor_thread.joinable()) {
        monitor_thread.join();
    }

    // Print all unique matrices found
    int matrix_idx = 1;
    for (const auto& matrix : unique_matrices) {
        printf("        %5d: | Aut(M) | = 1\n", matrix_idx++);
        //std::cout << "[SUCCESS] Valid P1F Matrix #" << matrix_idx++ << " discovered!" << std::endl;
        for (const auto& row : matrix) {
            printf(" \"");
            for (size_t i = 0; i < row.size(); i += 2) {
                printf(" %3d%3d", row[i], row[i + 1]);
            }
            printf(" \"\n");
        }
        printf("\n\n\n\n");
    }

    std::cout << get_timestamp() << "[INFO] Search complete. Total configurations processed: " << total_configurations_checked << std::endl;
    std::cout << get_timestamp() << "[INFO] Total unique P1F matrices found: " << unique_matrices.size() << std::endl;
    return 0;
}
