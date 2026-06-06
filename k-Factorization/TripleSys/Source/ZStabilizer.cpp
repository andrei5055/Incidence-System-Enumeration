#include "triplesys.h"
#include <new>
#include <map>

#if USE_INTRINSIC
extern void transform32_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut, int count);
extern void sortGroupsI(ctchar* pIn, tchar* pOut, int nPairs);
#endif

/**
 * Name: init
 * Purpose: Initializes the ZStabilizer component with graph settings, memory allocations, and search constraints.
 * Parameters:
 *   - size: Number of players/vertices in the graph (must be <= 32).
 *   - min_aut: The minimum target automorphism group size (|Aut(M)| >= min_aut).
 *   - bipartite: True if searching on a bipartite graph, false for complete graphs.
 *   - print: Flag to enable verbose debug console logging.
 */
void ZStabilizer::init(int size, int min_aut, bool bipartite, bool print) {
	if (size > 32) {
		printfRed("*** ZStabilizer: number of vertices(%d) > limit(32), exit 1.\n", size);
		exit(1);
	}
	num_rows = bipartite ? size / 2 : size - 1;
	player_index = 0;
	row_size = size;
	min_aut_val = min_aut;
	bBipartite = bipartite;
	bPrint = print;
	current_row = 0;
	transition_stack.clear();
	transition_stack.assign(num_rows, std::vector<std::vector<unsigned char>>());
	candidates.clear();
	candidates.assign(num_rows, std::vector<std::vector<unsigned char>>());
	memset(current_index, 0, sizeof(current_index));
	survivors.clear();
}

/**
 * Name: clearCurrent
 * Purpose: Clears candidates, transition stack, and indexes for a given row during backtracking to reset search states.
 * Parameters:
 *   - iRow: Index of the row to clear.
 */
void ZStabilizer::clearCurrent(int iRow) {
	if (iRow < num_rows - 1 && candidates[iRow + 1].size() > 0)
		candidates[iRow + 1].clear();
	if (iRow < num_rows) {
		transition_stack[iRow].clear();
		candidates[iRow].clear();
		current_index[iRow] = 0;
	}
}

/**
 * Name: sortCandidates
 * Purpose: Sorts and removes duplicate candidate rows at the current row level.
 */
void ZStabilizer::sortCandidates() {
	int iRow = current_row;
	if (iRow < 3) {
		printfRed("*** Request to sort row %d. Exit 1\n", iRow);
		exit(1);
	}
	std::sort(candidates[iRow].begin(), candidates[iRow].end());
	auto it = std::unique(candidates[iRow].begin(), candidates[iRow].end());
	candidates[iRow].erase(it, candidates[iRow].end());
#if 0
	if (bPrint && (int)candidates[iRow].size() > 0)
		printf("Number of candidates for row %d from automorphism of rows 1-%d: %d\n", 
			iRow + 1, bUseKSolve ? 3 : iRow - zebra_size + 1, (int)candidates[iRow].size());
	//if (iRow > 5)exit(1);
#endif
}

/**
 * Name: generateCandidates
 * Purpose: Orchestrates candidate generation for the next row by applying prefix automorphisms, calling the solver, and performing early backtracking checks.
 * Parameters:
 *   - rows: Pointer to the prefix matrix row array.
 *   - iRow: The index of the row being generated.
 * Returns: True if candidates were successfully generated, false if the branch should be pruned.
 */
bool ZStabilizer::generateCandidates(ctchar* rows, int iRow) {
	if (iRow < 3) {
		printfRed("*** Request to generate row %d. Exit 1\n", iRow);
		exit(1);
	}
	// Backtrack the solver's row pointer up to match the requested iRow
	while (iRow < current_row)
		goUp();
	// Move the solver's row pointer down to match the requested iRow
	while (iRow > current_row)
		goDown();
	int iTrToUse = iRow - zebra_size;

	size_t initial_size = candidates[iRow].size();

	// Iterate over all active prefix automorphisms to generate candidate rows
	for (const auto& tr : transition_stack[iTrToUse]) {
		// 1. Trace transition stack for prefix automorphism compatibility.
		bool is_case_a = false;
		if (iRow >= 3) {
			auto* prev_row = rows + ((iRow - 1) * row_size);
			alignas(32) tchar transformedPre[32 * 2];
			apply_permutation(transformedPre, prev_row, tr.data(), row_size);
			is_case_a = (memcmp(transformedPre, prev_row, row_size) == 0);
		}

		// 2. Generate new row candidates for each compatible automorphism.
		size_t size_before = candidates[iRow].size();
		generateRow(rows, tr.data(), iRow);
		size_t size_after = candidates[iRow].size();

		// 3. Early backtrack check: Prune the branch if a Case A automorphism fails and there are no initial candidates.
		if (is_case_a && size_after == size_before) {
			if (initial_size == 0) {
				return false;
			}
		}
	}

	sortCandidates();
	return candidates[iRow].size() > 0;
}
// Add automorphisms from starter
/**
 * Name: addTr
 * Purpose: Adds a permutation to the transition stack if it is a non-identity automorphism, and handles initial row 3 starter candidate extensions.
 * Parameters:
 *   - tr: Pointer to the permutation array.
 *   - rows: Pointer to the prefix matrix rows array.
 *   - iRow: The row index target for the starter translation.
 */
void ZStabilizer::addTr(ctchar* tr, ctchar* rows, int iRow) {
	ctchar* row3 = rows + row_size * 2;
	if (iRow < 3 - zebra_size || iRow > 2) {
		printfRed("*** ZStabilizer internal error: addTr iRow(%d) is incorrect for zebra size %d. Exit(1)\n", iRow, zebra_size);
		exit(1);
	}
	// Align current_row up to target iRow
	while (iRow < current_row)
		goUp();
	// Align current_row down to target iRow
	while (iRow > current_row)
		goDown();
	bool bOk = false;
	// Check if the permutation is different from the identity permutation
	for (int i = 0; i < row_size; i++) {
		if (tr[i] != i) {
			bOk = true;
			break;
		}
	}
	if (!bOk)
		return;
	alignas(32) tchar transformed[32 * 2]; // Buffer for (row_size*2 max)
	alignas(32) tchar transformedPre[32 * 2]; // Buffer for (row_size*2 max)
	if (iRow == 1) {
		apply_permutation(transformed, row3, tr, row_size);
		auto v3 = transformed[1];
		if (v3 == row3[1]) {
			if (memcmp(transformed, row3, row_size) != 0)
				return;
			//return; // temporary allow only automorphism with rows n,n+1 swap
		}
		else {
			apply_permutation(transformedPre, transformed, tr, row_size);
			if (memcmp(transformedPre, row3, row_size) != 0)
				return;
			tchar tRow = (bBipartite ? transformed[1] / 2 : transformed[1]) - 1;

			memset(active_cand_adj, 0, 32);
			ctchar* row = transformed;
			// Reconstruct the candidate adjacency list representation for the transformed row
			for (int i = 0; i < row_size; i += 2) {
				active_cand_adj[row[i]] = row[i + 1];
				active_cand_adj[row[i+1]] = row[i]; 
			}
			if (!isP1F(rows, 3))
				return;
			// without this candidate we get 22 results for K(11,11) aut > 1 (instead of 23)
			std::vector<unsigned char> candidate(transformed, transformed + row_size);
			candidates[tRow].push_back(candidate);
#if 0
			if (bPrint)
				printTable("c", transformed, 1, row_size, 2);
#endif
		}
	}
	std::vector<unsigned char> perm(tr, tr + row_size);
	transition_stack[current_row].push_back(perm);
}

/**
 * Name: getCandidate
 * Purpose: Returns the next candidate matching from the candidate list for a given row index.
 * Parameters:
 *   - iRow: The row index of the candidates to retrieve.
 * Returns: Pointer to the candidate matching array, or NULL if list is exhausted.
 */
tchar* ZStabilizer::getCandidate(int iRow) {
	if (iRow >= candidates.size() || current_index[iRow] >= candidates[iRow].size())
		return NULL;
	return candidates[iRow][current_index[iRow]++].data();
}

/**
 * Name: getTr
 * Purpose: Returns a specific automorphism permutation from the transition stack of a given row level.
 * Parameters:
 *   - iRow: The row index of the transition stack.
 *   - iTr: The index of the automorphism in the stack.
 * Returns: Pointer to the permutation array, or NULL if indices are out of bounds.
 */
tchar* ZStabilizer::getTr(int iRow, int iTr) {
	if (iRow >= num_rows || iTr >= transition_stack[iRow].size())
		return NULL;
	return transition_stack[iRow][iTr].data();
}

/**
 * Name: goUp
 * Purpose: Climbs up the backtracking search tree (decrements current_row and resets current level data).
 */
void ZStabilizer::goUp() {
	// Backtrack
	clearCurrent(current_row);
	if (current_row > 0)
		current_row--;
}

/**
 * Name: goDown
 * Purpose: Descends down the backtracking search tree (increments current_row).
 */
void ZStabilizer::goDown() {
	if (current_row < num_rows - 1)
		current_row++;
}

/**
 * Name: sortGroups2CPU
 * Purpose: Standardizes and sorts edge matchings on the CPU so that the vertex pairs in each match group are ordered lexicographically.
 * Parameters:
 *   - mi: Input pointer to vertex matches.
 *   - mo: Output pointer to sorted matches.
 *   - nGroups: Number of matching groups (n / 2).
 */
void sortGroups2CPU(ctchar* mi, tchar* mo, int nGroups)
{
	ctchar* pmi[32];
	memset(pmi, 0, sizeof(pmi));
	auto* mic = mi;

	// Map each group edge matching to an array index for sorting
	for (tchar ir = 0; ir < nGroups; ir++)
	{
		pmi[*mic] = mic;
		mic += 2;
	}
	const auto iMax = nGroups * 2 - 1;
	auto* mos = (short int*)mo;
	// Write out sorted matching pairs sequentially to the output buffer
	for (tchar i = 0; i < iMax; i++)
	{
		if (auto* mis = (short int*)pmi[i])
			*mos++ = *mis;
	}
}

/**
 * Name: apply_permutation
 * Purpose: Applies an automorphism permutation to a row matching and standardizes the result.
 * Parameters:
 *   - trRow: Output buffer for the permuted row.
 *   - row: Input row matching array.
 *   - tr: Automorphism permutation mapping table.
 *   - n: Vertex count (row_size).
 */
void ZStabilizer::apply_permutation(tchar* trRow, ctchar* row, ctchar* tr, int n) {
	alignas(32) tchar pTmp[32];
#if USE_INTRINSIC
	transform32_and_sort_each_pair(row, tr, pTmp, n);
	sortGroupsI(pTmp, trRow, n / 2);
#else
	// 1. Apply permutation and sort individual edge vertices
	for (int i = 0; i < n; i += 2) {
		unsigned char u = tr[row[i]];
		unsigned char v = tr[row[i + 1]];
		ASSERT_IF(u > 31 || v > 31);
		if (u < v) { pTmp[i] = u; pTmp[i + 1] = v; }
		else { pTmp[i] = v; pTmp[i + 1] = u; }
	}

	// 2. Internal Sort of Pairs (Standardize)
	sortGroups2CPU(pTmp, trRow, n / 2)
#endif
}

/**
 * Name: setPlayerIndexByPos
 * Purpose: Finds the maximum player position index that has changed after mapping an automorphism.
 * Parameters:
 *   - tr: Automorphism permutation table.
 *   - co: Candidate row matching.
 *   - ciFrom: Reference row matching.
 *   - ip: Prefix match index.
 * Returns: Calculated player index to backtrack to.
 */
int ZStabilizer::setPlayerIndexByPos(ctchar* tr, ctchar* co, ctchar* ciFrom, int ip)
{
	tchar ttr[32];
	tchar tpr[32];
	tchar i, j = 0, k;
	// Initialize permutation inverse mappings for quick index lookups
	for (i = 0; i < row_size; i++)
	{
		tpr[ciFrom[i]] = ttr[tr[i]] = i;
	}
	// Find the maximum index of the player that was permuted
	for (i = 0; i <= ip; i++)
	{
		k = ttr[co[i]];
		if (j < tpr[k])
			j = tpr[k];
	}
	if (player_index > j) {
		player_index = j;
	}
	return player_index;
}
/**
 * Name: setPlayerIndex
 * Purpose: Finds the first index where the candidate row differs from the prefix canonical row.
 * Parameters:
 *   - tr: Automorphism permutation table.
 *   - co: Candidate row matching.
 *   - ci: Prefix canonical row.
 *   - ciFrom: Reference row matching.
 * Returns: Index of the player to backtrack to.
 */
int ZStabilizer::setPlayerIndex(ctchar* tr, ctchar* co, ctchar* ci, ctchar* ciFrom) {
	const auto iMax = row_size - 3;
	int i = 0;
	// Scan matching elements to locate the first mismatch between the candidate and canonical row
	for (; i < iMax; i++) {
		if (co[i] < ci[i])
			break;
	}
	return setPlayerIndexByPos(tr, co, ciFrom, i);
}

/**
 * Name: check
 * Purpose: Performs canonicalization check on the current generated candidate row under the active stabilizer group action, filtering survivors.
 * Parameters:
 *   - row: The candidate row being checked.
 *   - preRow: The previous row of the prefix.
 *   - playerIndex: Output parameter to receive the backtracking player index in case of non-canonicity.
 *   - iRow: The row index of the candidate.
 * Returns: 1 if canonical and valid, 0 if not canonical, -1 if fails automorphism order requirement, -2 if prefix stabilizer group size is too small.
 */
int ZStabilizer::check(unsigned char* row, unsigned char* preRow, int* playerIndex, int iRow) {
	player_index = *playerIndex;
	if (iRow < 3) {
		printfRed("*** ZStabilizer internal error: call to check() with iRow(%d) < 3. Exit(1)\n", iRow);
		exit(1);
	}
	// Align current_row up to match iRow
	while (iRow < current_row)
		goUp();
	// Align current_row down to match iRow
	while (iRow > current_row)
		goDown();

	if ((int)transition_stack[current_row - zebra_size].size() + 1 < min_aut_val) {
		*playerIndex = 0;
		return -2;
	}

	survivors.clear();
	survivors_candidates.clear();
	int cmp;
	int iTr = 0;
	alignas(32) tchar transformed[32 * 2]; // Buffer for (row_size*2 max)
	alignas(32) tchar transformedPre[32 * 2]; // Buffer for (row_size*2 max)
	// Iterate over prefix automorphism group elements to perform symmetry checks on the candidate
	for (const auto& tr : transition_stack[current_row - zebra_size]) {
		apply_permutation(transformed, row, tr.data(), row_size);
		if (zebra_size == 1) {
			if (transformed[1] != row[1])
				continue;
			cmp = std::memcmp(transformed, row, row_size);
			if (cmp < 0) {
				*playerIndex = setPlayerIndex(tr.data(), transformed, row, row);
				return 0;      // Not canonical
			}
			if (cmp == 0) survivors.push_back(tr); // Symmetry survives
		}
		else if (transformed[1] == preRow[1]) { // swap with prev row
			cmp = std::memcmp(transformed, preRow, row_size);
			if (cmp < 0) {
				*playerIndex = setPlayerIndex(tr.data(), transformed, preRow, row);
				return 0;      // Not canonical
			}
			if (cmp == 0) {
				apply_permutation(transformedPre, preRow, tr.data(), row_size);
				cmp = std::memcmp(transformedPre, row, row_size);
				if (cmp < 0) {
					// can't modify *playerIndex
					return 0;      // Not canonical
				}
				if (cmp == 0) survivors.push_back(tr); // Symmetry survives
			}
		}
		else { // same row
			apply_permutation(transformedPre, preRow, tr.data(), row_size);
			if (transformedPre[1] == preRow[1]) {
				cmp = std::memcmp(transformedPre, preRow, row_size);
				if (cmp < 0) {
					return 0;
					printfYellow("*** ZStabilizer input error: matrix (%d-rows) not canonical, exit(1)\n", iRow);
					printTable("Row ", preRow, 1, row_size, 2);
					printTable("Tr  ", tr.data(), 1, row_size, 2);
					printTable("Rslt", transformedPre, 1, row_size, 2);
					exit(1);
				}
				if (cmp == 0) {
					cmp = std::memcmp(transformed, row, row_size);
					if (cmp < 0) {
						*playerIndex = setPlayerIndex(tr.data(), transformed, row, row);
						return 0;      // Not canonical
					}
					if (cmp == 0) {
						survivors.push_back(tr); // Symmetry survives
						//printf("iTr = %d size=%d ", iTr, (int)survivors.size());
						//printTable("", tr.data(), 1, row_size, 2);
					}
				}
			}
		}/**
		else if (current_row + zebra_size < num_rows && transformed[1] == row[1] + bBipartite ? 2 : 1) {
			apply_permutation(transformedPre, transformed, tr.data(), row_size);
			cmp = std::memcmp(transformedPre, row, row_size);
			if (cmp == 0) survivors.push_back(tr); // Symmetry survives
		}**/
	}
	// 4. Threshold check: Survivors + Identity
	if ((int)survivors.size() + 1 >= min_aut_val) {
		if (current_row < num_rows - 1) {
			transition_stack[current_row] = survivors;
			if (zebra_size == 2) {
				// For zebra_size == 2, find swapping symmetries with the next row level
				for (const auto& tr : transition_stack[current_row - 1]) {
					apply_permutation(transformed, row, tr.data(), row_size);
					if (transformed[1] == row[1] + (bBipartite ? 2 : 1)) { // swap with next row?
						apply_permutation(transformedPre, transformed, tr.data(), row_size);
						cmp = std::memcmp(transformedPre, row, row_size);
						if (cmp == 0) {
							memset(active_cand_adj, 0, 32);
							// Convert the candidate row matching into an adjacency representation
							for (int i = 0; i < row_size; i += 2) {
								active_cand_adj[row[i]] = row[i + 1];
								active_cand_adj[row[i + 1]] = row[i]; 
							}
							tchar* r = preRow - row_size * (iRow - 1);
							if (isP1F(r, iRow)) {
								memset(active_cand_adj, 0, 32);
								// Set up the permuted matching adjacency for Hamiltonian cycle verification
								for (int i = 0; i < row_size; i += 2) {
									active_cand_adj[transformed[i]] = transformed[i + 1]; 
									active_cand_adj[transformed[i + 1]] = transformed[i];
								}
								if (isP1F(r, iRow) && isP1F(row, 1)) {
									std::vector<unsigned char> candidate(transformed, transformed + row_size);
									candidates[current_row + 1].push_back(candidate);
#if 0
									if (bPrint)
										printTable("c", transformed, 1, row_size, 2);
#endif
								}
							}
						}
					}
				}
			}
		}
		return 1;
	}
	return -1;
}
// --- GLOBAL ABSOLUTE P1F VERIFICATION (NO LEAKS) ---
/**
 * Name: isP1F
 * Purpose: Verifies whether the active candidate matching forms a Hamiltonian cycle with all prefix rows.
 * Parameters:
 *   - rows: Pointer to the prefix matrix rows array.
 *   - iRow: The row index of the candidate.
 * Returns: True if all prefix union graphs are valid Hamiltonian cycles, false otherwise.
 */
bool ZStabilizer::isP1F(ctchar* rows, int iRow) {
	// --- GLOBAL ABSOLUTE P1F VERIFICATION (NO LEAKS) ---
	// Iterate over all prefix rows to check for Hamiltonian cycle validity
	for (int r = 0; r < iRow; ++r) {
		auto* base_row = rows + (r * row_size);
		uint8_t base_adj[32] = { 0 };
		// Set up the adjacency representation of the current prefix row
		for (int i = 0; i < row_size; i += 2) {
			base_adj[base_row[i]] = base_row[i + 1];
			base_adj[base_row[i + 1]] = base_row[i];
		}

		uint8_t curr = 0;
		int visited_count = 0;
		// Chase the alternating matching pointers starting from 0 to verify cycle length
		do {
			curr = base_adj[curr];
			curr = active_cand_adj[curr];
			visited_count += 2;
		} while (curr != 0 && visited_count <= row_size);

		if (curr != 0 || visited_count != row_size) 
			return false;
	}
	return true;
}

#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>


/**
 * Name: checkCandidate
 * Purpose: Verifies a generated matching candidate, standardizes it, and adds it to the candidate list for the row.
 * Parameters:
 *   - rows: Pointer to the prefix rows.
 *   - candidate: The generated matching candidate.
 *   - iRow: The row index.
 */
void ZStabilizer::checkCandidate(ctchar* rows, ctchar* candidate, int iRow)
{
	alignas(32) tchar pTmp[32 * 2];
	//printTable("Tr", candidate, 1, row_size, 2);
#if USE_INTRINSIC
	sortGroupsI(candidate, pTmp, row_size / 2);
#elif
	sortGroups2CPU(pTmp, candidate, row_size / 2);
#endif
	std::vector<unsigned char> vcandidate(pTmp, pTmp + row_size);
#if 0
	static tchar r24[24] =
		//{ 0,4,1,8,2,6,3,10,5,12,7,14,9,16,11,18,13,20,15,22,17,23,19,21 };
	{ 0,5,1,9,2,11,3,7,4,13,6,15,8,17,10,19,12,21,14,23,16,22,18,20 };
	int icmp;
	if (current_row == 4) {
		icmp = memcmp(r24, candidate, 15);
		if (icmp == 0) {
			printTable("t", candidate, 1, row_size, 2);
			icmp = icmp;
		}
	}
#endif
	candidates[iRow].push_back(vcandidate);
	int iExp = bBipartite ? iRow * 2 + 1: iRow + 1;
	if (pTmp[1] != iExp) {
		printfRed("*** Candidate not for Row %d (expected 0 %d ...), exit 1\n", iRow + 1, iExp);
		printTable("Row candidate", pTmp, 1, row_size, 2);
		exit(1);
	}
}

// Full self-contained generation controller - Functionally Safe & Stable for All Starters
/**
 * Name: generateRow
 * Purpose: Generates candidate row matchings under a specific prefix automorphism.
 * Parameters:
 *   - rows: Pointer to the prefix matrix row array.
 *   - automorphism: Automorphism permutation array.
 *   - iRow: Index of the row being generated.
 * Returns: The row index.
 */
int ZStabilizer::generateRow(ctchar* rows, ctchar* automorphism, int iRow) {
	size_t initial_size = candidates[iRow].size();

	// Check if this automorphism preserves prev_row (Case A)
	bool is_case_a = false;
	if (iRow >= 3) {
		auto* prev_row = rows + ((iRow - 1) * row_size);
		alignas(32) tchar transformedPre[32 * 2];
		apply_permutation(transformedPre, prev_row, automorphism, row_size);
		is_case_a = (memcmp(transformedPre, prev_row, row_size) == 0);
	}

	// Pre-compute the adjacency lists for all prefix rows 0 to iRow - 1
	alignas(32) uint8_t precomputed_adj[32 * 32];
	// Precompute adjacency representations for prefix rows for fast O(1) cycle checks
	for (int r = 0; r < iRow; ++r) {
		auto* base_row = rows + (r * row_size);
		auto* adj = precomputed_adj + r * 32;
		memset(adj, 0, 32);
		// Build adjacency entries for each edge in the prefix row
		for (int i = 0; i < row_size; i += 2) {
			adj[base_row[i]] = base_row[i + 1];
			adj[base_row[i + 1]] = base_row[i];
		}
	}

	// 1. Mark edges already consumed by all previous rows
	bool edge_taken[32][32] = { false };
	// Mark edges consumed by prefix rows as unavailable
	for (int r = 0; r < iRow; ++r) {
		auto* row = rows + (r * row_size);
		for (int i = 0; i < row_size; i += 2) {
			uint8_t u = row[i];
			uint8_t v = row[i + 1];
			// Mark both directions of the consumed edge
			edge_taken[u][v] = true;
			edge_taken[v][u] = true;
		}
	}

	// 2. Gather all remaining available edges
	std::vector<Edge> remaining_edges;
	// Loop through vertices to discover and build remaining available edges
	for (uint8_t u = 0; u < row_size; ++u) {
		// Identify edge endpoints that are not yet taken and satisfy bipartite restrictions if active
		for (uint8_t v = u + 1; v < row_size; ++v) {
			if (!edge_taken[u][v]) {
				if (bBipartite && ((u % 2) == (v % 2))) continue;
				remaining_edges.push_back(Edge(u, v));
			}
		}
	}

	// 3. Group remaining edges into orbits under the passed automorphism
	std::vector<std::vector<Edge>> all_orbits;
	std::vector<bool> visited_edges(remaining_edges.size(), false);
	// Construct edge orbits under the automorphism action by tracing cycles of edges
	for (size_t i = 0; i < remaining_edges.size(); ++i) {
		if (visited_edges[i]) continue;
		std::vector<Edge> orbit;
		Edge current = remaining_edges[i];

		// Cycle through edge mappings until the orbit cycle closes
		while (true) {
			auto it = std::find(orbit.begin(), orbit.end(), current);
			if (it != orbit.end()) break;
			orbit.push_back(current);

			auto rem_it = std::find(remaining_edges.begin(), remaining_edges.end(), current);
			if (rem_it != remaining_edges.end()) {
				visited_edges[std::distance(remaining_edges.begin(), rem_it)] = true;
			}

			uint8_t next_u = automorphism[current.u];
			uint8_t next_v = automorphism[current.v];
			current = Edge(next_u, next_v);
		}
		all_orbits.push_back(orbit);
	}

	// 4. Filter orbits down to only those mathematically capable of fitting inside Row
	std::vector<std::vector<Edge>> valid_orbits;
	std::vector<uint32_t> orbit_masks;

	// Filter out edge orbits that have vertex conflicts (e.g. self-intersecting orbits)
	for (const auto& orbit : all_orbits) {
		if (orbit.size() > (size_t)(row_size / 2)) continue;
		uint32_t mask = 0;
		bool vertex_conflict = false;
		// Calculate vertex coverage bitmask for the orbit and check for overlap
		for (const auto& edge : orbit) {
			uint32_t edge_mask = (1UL << edge.u) | (1UL << edge.v);
			if ((mask & edge_mask) != 0) {
				vertex_conflict = true;
				break;
			}
			mask |= edge_mask;
		}
		if (!vertex_conflict) {
			valid_orbits.push_back(orbit);
			orbit_masks.push_back(mask);
		}
	}

	// FIXED: Return iRow unconditionally if empty to let your top-level handle initialization naturally
	if (valid_orbits.empty()) {
		return iRow;
	}

	// Orbit Coverage Check: Verify if all vertices can be covered by the union of all valid orbits.
	uint32_t all_covered = 0;
	// Accumulate union of all vertex masks covered by the orbits
	for (uint32_t mask : orbit_masks) {
		all_covered |= mask;
	}
	if (all_covered != ((1ULL << row_size) - 1)) {
		return iRow; // Reject early: cannot cover all vertices
	}

	// --- MANDATORY EDGE ORBIT FILTERING & PRE-SELECTION ---
	uint8_t neighbor0 = bBipartite ? (iRow * 2 + 1) : (iRow + 1);
	int target_orbit_idx = -1;

	// Scan valid orbits to find the target orbit containing the baseline edge (0, neighbor0)
	for (size_t i = 0; i < valid_orbits.size(); ++i) {
		// Search edges within the current orbit
		for (const auto& edge : valid_orbits[i]) {
			if (edge.u == 0 && edge.v == neighbor0) {
				target_orbit_idx = (int)i;
				break;
			}
		}
		if (target_orbit_idx != -1) break;
	}

	// If no orbit contains (0, neighbor0), then no candidate row is possible.
	if (target_orbit_idx == -1) {
		return iRow;
	}

	// Filter out any orbit that covers vertex 0 or neighbor0 but is NOT the target orbit.
	std::vector<std::vector<Edge>> filtered_orbits;
	std::vector<uint32_t> filtered_masks;
	int new_target_orbit_idx = -1;

	// Discard orbits that conflict with the pre-selected target orbit (i.e. cover 0 or neighbor0)
	for (size_t i = 0; i < valid_orbits.size(); ++i) {
		const auto& orbit = valid_orbits[i];
		bool covers_0 = false;
		bool covers_neighbor0 = false;
		// Check vertex coverage within each orbit
		for (const auto& edge : orbit) {
			if (edge.u == 0 || edge.v == 0) covers_0 = true;
			if (edge.u == neighbor0 || edge.v == neighbor0) covers_neighbor0 = true;
		}

		if ((covers_0 || covers_neighbor0) && (int)i != target_orbit_idx) {
			continue; // Discard conflicting orbit
		}

		if ((int)i == target_orbit_idx) {
			new_target_orbit_idx = (int)filtered_orbits.size();
		}
		filtered_orbits.push_back(orbit);
		filtered_masks.push_back(orbit_masks[i]);
	}

	valid_orbits = std::move(filtered_orbits);
	orbit_masks = std::move(filtered_masks);
	target_orbit_idx = new_target_orbit_idx;
	
	// Fast bitmask solver variables
	memset(active_cand_adj, 0, 32);

	auto get_match = [&](uint8_t v) -> int {
		if (v == 0) {
			return active_cand_adj[0] == 0 ? -1 : active_cand_adj[0];
		}
		if (active_cand_adj[v] != 0) return active_cand_adj[v];
		if (active_cand_adj[0] == v) return 0;
		return -1;
	};

	auto check_new_subcycles = [&](const std::vector<Edge>& orbit, int num_r) -> bool {
		if (bBipartite) return false; // Bypass look-ahead pruning for bipartite graphs to avoid over-filtering valid matchings
		auto* base_adj = precomputed_adj;
		// Check cycle closure against each prefix row independently
		for (int r = 0; r < num_r; ++r, base_adj += 32) {
			// Evaluate cycles generated by placing the proposed edge orbit
			for (const auto& edge : orbit) {
				uint8_t u = edge.u;
				uint8_t v = edge.v;
				uint8_t curr = v;
				uint32_t path_mask = (1UL << u) | (1UL << v);
				bool is_cycle = false;
				// Chase alternating pointers until we close a cycle, hit a path end, or detect an infinite loop
				while (true) {
					curr = base_adj[curr];
					if (curr == u) {
						is_cycle = true;
						break;
					}
					path_mask |= (1UL << curr);
					int next = get_match(curr);
					if (next == -1) {
						break; // Path ends
					}
					if (path_mask & (1UL << next)) {
						break; // Already visited
					}
					curr = next;
				}
				if (is_cycle) {
					int path_len = __popcnt(path_mask);
					if (path_len < row_size) {
						return true; // Found a subcycle!
					}
				}
			}
		}
		return false;
	};

	// Pre-select the target orbit
	std::vector<size_t> chosen_indices;
	uint32_t initial_vertices = 0;
	if (target_orbit_idx != -1) {
		chosen_indices.push_back(target_orbit_idx);
		initial_vertices |= orbit_masks[target_orbit_idx];
		// Copy target orbit edges to the active adjacency
		for (const auto& edge : valid_orbits[target_orbit_idx]) {
			active_cand_adj[edge.u] = edge.v;
			active_cand_adj[edge.v] = edge.u;
		}
		if (check_new_subcycles(valid_orbits[target_orbit_idx], iRow)) {
			memset(active_cand_adj, 0, 32);
			return iRow;
		}
	}

	// 5. Recursive Bitmask Solver (Algorithm X / Exact Cover style)
	auto solve = [&](auto& self, uint32_t current_vertices, std::vector<size_t>& chosen_indices) -> void {
		if (current_vertices == ((1ULL << row_size) - 1)) {
			// Fast O(1) early check: Candidate must match the expected neighbor of 0.
			// This is guaranteed since the target orbit containing (0, neighbor0) was pre-selected,
			// but we keep it as a lightweight defensive guard.
			uint8_t neighbor0_check = bBipartite ? (iRow * 2 + 1) : (iRow + 1);
			if (active_cand_adj[0] != neighbor0_check)
				return;

			// Fast P1F Cycle Check:
			// Verifies if the candidate forms a Perfect 1-Factorization (Hamiltonian cycle) with each prefix row.
			// Traversing the alternate matching paths is O(L) time and O(1) space.
			// If the cycle starting at 0 returns to 0 after exactly row_size steps, it is mathematically 
			// guaranteed to be a single Hamiltonian cycle covering all vertices, meaning no other disjoint
			// cycles exist. Thus, we can completely omit visited-element tracking and redundant loops.
			auto isP1F_fast = [&](int num_r) -> bool {
				auto* base_adj = precomputed_adj;
				// Verify cycle connectivity for each prefix row
				for (int r = 0; r < num_r; ++r, base_adj += 32) {
					uint8_t curr = 0;
					int visited_count = 0;
					// Fast trace of alternate paths starting from vertex 0
					do {
						curr = base_adj[curr];
						curr = active_cand_adj[curr];
						visited_count += 2;
					} while (curr != 0 && visited_count <= row_size);

					if (curr != 0 || visited_count != row_size) 
						return false;
				}
				return true;
			};

			if (!isP1F_fast(iRow))
				return;

			// --- CASE B: SWAPPING SYMMETRY VALIDATOR ---
			uint8_t swap_adj[32] = { 0 };
			bool is_case_b = false;

			if (zebra_size == 2 && iRow >= 3) {
				auto* prev_row = rows + ((iRow - 1) * row_size);
				bool candidate_matches_swap = true;

				// Verify if candidate matching maps exactly to the swapped previous matching
				for (int i = 0; i < row_size; i += 2) {
					uint8_t mapped_u = automorphism[prev_row[i]];
					uint8_t mapped_v = automorphism[prev_row[i + 1]];
					if (active_cand_adj[mapped_u] != mapped_v) {
						candidate_matches_swap = false;
						break;
					}
				}

				if (!candidate_matches_swap) {
					is_case_b = true;
					// Set up swapping symmetry adjacency relations
					for (int i = 0; i < row_size; i += 2) {
						uint8_t mu = automorphism[prev_row[i]];
						uint8_t mv = automorphism[prev_row[i + 1]];
						swap_adj[mu] = mv;
						swap_adj[mv] = mu;
					}
				}
			}

			// --- CRITICAL CRISS-CROSS CHECK FOR CASE B ---
			if (is_case_b) {
				bool twin_visited[32] = { false };
				uint8_t curr = 0;
				int twin_visited_count = 0;
				// Verify Hamiltonian cycle condition for the union of the swapped matchings
				do {
					twin_visited[curr] = true;
					curr = swap_adj[curr];
					twin_visited[curr] = true;
					curr = active_cand_adj[curr];
					twin_visited_count += 2;
				} while (curr != 0);

				if (twin_visited_count < row_size) return;
				// Ensure every single vertex was visited in the criss-cross traversal
				for (int v = 0; v < row_size; ++v) {
					if (!twin_visited[v]) return;
				}
			}

			// PERFORMANCE HOTSPOT OPTIMIZATION: Defer Edge List Reconstruction & Sorting
			// 99.9% of generated candidates fail the P1F cycle or Case B symmetry checks.
			// By checking P1F and symmetries first, we avoid constructing the Edge vector,
			// memory allocations, copying, and calling std::sort for the discarded matchings.
			std::vector<Edge> candidate_edges;
			// Collect all matched edges from the chosen orbits
			for (size_t idx : chosen_indices) {
				// Verify active edge assignments against the current solver configuration
				for (const auto& edge : valid_orbits[idx]) {
					if (active_cand_adj[edge.u] == edge.v) {
						candidate_edges.push_back(edge);
					}
				}
			}

			std::sort(candidate_edges.begin(), candidate_edges.end(), [](const Edge& a, const Edge& b) {
				return a.u < b.u;
				});

			unsigned char candidate[32];
			int byte_idx = 0;
			// Populate the candidate array using the sorted edge vertices
			for (const auto& edge : candidate_edges) {
				candidate[byte_idx++] = edge.u;
				candidate[byte_idx++] = edge.v;
			}
			checkCandidate(rows, candidate, iRow);
			return;
		}

		// 1. Choose the first uncovered vertex
		int c = -1;
		// Search bitmask to choose the first uncovered vertex (Exact Cover column selection)
		for (int i = 0; i < row_size; ++i) {
			if (!(current_vertices & (1UL << i))) {
				c = i;
				break;
			}
		}
		if (c == -1) return;

		// 2. Branch only on compatible orbits covering vertex 'c'
		// Branch on all candidate orbits that cover the selected vertex 'c'
		for (size_t idx = 0; idx < valid_orbits.size(); ++idx) {
			if ((orbit_masks[idx] & (1UL << c)) && (current_vertices & orbit_masks[idx]) == 0) {
				// Temporarily place the candidate orbit edge matches into the active adjacency
				for (const auto& edge : valid_orbits[idx]) {
					active_cand_adj[edge.u] = edge.v;
					active_cand_adj[edge.v] = edge.u;
				}

				if (!check_new_subcycles(valid_orbits[idx], iRow)) {
					chosen_indices.push_back(idx);
					self(self, current_vertices | orbit_masks[idx], chosen_indices);
					chosen_indices.pop_back();
				}

				// Backtrack: Remove candidate orbit edge matches from the active adjacency
				for (const auto& edge : valid_orbits[idx]) {
					active_cand_adj[edge.u] = 0;
					active_cand_adj[edge.v] = 0;
				}
			}
		}
	};

	solve(solve, initial_vertices, chosen_indices);

	// Clean up active_cand_adj after solver finishes
	memset(active_cand_adj, 0, 32);

	return iRow;
}
