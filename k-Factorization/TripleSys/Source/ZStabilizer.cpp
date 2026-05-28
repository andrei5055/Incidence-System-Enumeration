#include "triplesys.h"
#include <new>

#if USE_INTRINSIC
extern void transform32_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut, int count);
extern void sortGroupsI(ctchar* pIn, tchar* pOut, int nPairs);
#endif

void ZStabilizer::init(int size, int min_aut, bool useKSolve, bool bipartite, bool print) {
	if (size > 32) {
		printfRed("*** ZStabilizer: number of vertices(%d) > limit(32), exit 1.\n", size);
		exit(1);
	}
	bUseKSolve = useKSolve;
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

void ZStabilizer::clearCurrent(int iRow) {
	if (iRow < num_rows) {
		transition_stack[iRow].clear();
		candidates[iRow].clear();
		current_index[iRow] = 0;
	}
}
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

void ZStabilizer::generateCandidates(ctchar* rows, int iRow) {
	if (iRow < 3) {
		printfRed("*** Request to generate row %d. Exit 1\n", iRow);
		exit(1);
	}
	while (iRow < current_row)
		goUp();
	while (iRow > current_row)
		goDown();
	//int iTrToUse = bUseKSolve ? 2 - (iRow % zebra_size) : iRow - zebra_size;
	int iTrToUse = bUseKSolve ? iRow - zebra_size : iRow - zebra_size;
	for (const auto& tr : transition_stack[iTrToUse]) {
		generateRow(rows, tr.data(), iRow);
	}/**
	for (const auto& tr : transition_stack[iTrToUse + 1]) {
		generateRow(rows, tr.data(), iRow);
	}**/
	sortCandidates();
}
// Add automorphisms from starter
void ZStabilizer::addTr(ctchar* tr, ctchar* row3, int iRow) {
	if (iRow < 3 - zebra_size || iRow > 2) {
		printfRed("*** ZStabilizer internal error: addTr iRow(%d) is incorrect for zebra size %d. Exit(1)\n", iRow, zebra_size);
		exit(1);
	}
	while (iRow < current_row)
		goUp();
	while (iRow > current_row)
		goDown();
	bool bOk = false;
	for (int i = 0; i < row_size; i++) {
		if (tr[i] != i) {
			bOk = true;
			break;
		}
	}
	if (!bOk)
		return;
	alignas(32) tchar transformed[MAX_PLAYER_NUMBER * 2]; // Buffer for (row_size*2 max)
	alignas(32) tchar transformedPre[MAX_PLAYER_NUMBER * 2]; // Buffer for (row_size*2 max)
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
			//if ((bUseKSolve && tRow < 3) || (!bUseKSolve && tRow != 3))
				//return;
			std::vector<unsigned char> candidate(transformed, transformed + row_size);
			candidates[tRow].push_back(candidate);
		}
	}
	std::vector<unsigned char> perm(tr, tr + row_size);
	transition_stack[current_row].push_back(perm);
}

tchar* ZStabilizer::getCandidate(int iRow) {
	if (iRow >= candidates.size() || current_index[iRow] >= candidates[iRow].size())
		return NULL;
	return candidates[iRow][current_index[iRow]++].data();
}

tchar* ZStabilizer::getTr(int iRow, int iTr) {
	if (iRow >= num_rows || iTr >= transition_stack[iRow].size())
		return NULL;
	return transition_stack[iRow][iTr].data();
}

void ZStabilizer::goUp() {
	// Backtrack
	clearCurrent(current_row);
	if (current_row > 0)
		current_row--;
}

void ZStabilizer::goDown() {
	if (current_row < num_rows - 1)
		current_row++;
}

void sortGroups2CPU(ctchar* mi, tchar* mo, int nGroups)
{
	ctchar* pmi[32];
	memset(pmi, 0, sizeof(pmi));
	auto* mic = mi;

	for (tchar ir = 0; ir < nGroups; ir++)
	{
		pmi[*mic] = mic;
		mic += 2;
	}
	const auto iMax = nGroups * 2 - 1;
	auto* mos = (short int*)mo;
	for (tchar i = 0; i < iMax; i++)
	{
		if (auto* mis = (short int*)pmi[i])
			*mos++ = *mis;
	}
}

void ZStabilizer::apply_permutation(tchar* trRow, ctchar* row, ctchar* tr, int n) {
	alignas(32) tchar pTmp[32];
#if USE_INTRINSIC
	transform32_and_sort_each_pair(row, tr, pTmp, n);
	sortGroupsI(pTmp, trRow, n / 2);
#else
	// 1. Apply permutation
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
int ZStabilizer::setPlayerIndexByPos(ctchar* tr, ctchar* co, ctchar* ciFrom, int ip)
{
	tchar ttr[MAX_PLAYER_NUMBER];
	tchar tpr[MAX_PLAYER_NUMBER];
	tchar i, j = 0, k;
	for (i = 0; i < row_size; i++)
	{
		tpr[ciFrom[i]] = ttr[tr[i]] = i;
	}
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
int ZStabilizer::setPlayerIndex(ctchar* tr, ctchar* co, ctchar* ci, ctchar* ciFrom) {
	const auto iMax = row_size - 3;
	int i = 0;
	for (; i < iMax; i++) {
		if (co[i] < ci[i])
			break;
	}
	return setPlayerIndexByPos(tr, co, ciFrom, i);
}
int ZStabilizer::check(unsigned char* row, unsigned char* preRow, int* playerIndex, int iRow) {
	player_index = *playerIndex;
	if (iRow < 3) {
		printfRed("*** ZStabilizer internal error: call to check() with iRow(%d) < 3. Exit(1)\n", iRow);
		exit(1);
	}
	while (iRow < current_row)
		goUp();
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
	alignas(32) tchar transformed[MAX_PLAYER_NUMBER * 2]; // Buffer for (row_size*2 max)
	alignas(32) tchar transformedPre[MAX_PLAYER_NUMBER * 2]; // Buffer for (row_size*2 max)
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
				for (const auto& tr : transition_stack[current_row - 1]) {
					apply_permutation(transformed, row, tr.data(), row_size);
					if (transformed[1] == row[1] + (bBipartite ? 2 : 1)) { // swap with next row?
						apply_permutation(transformedPre, transformed, tr.data(), row_size);
						cmp = std::memcmp(transformedPre, row, row_size);
						if (cmp == 0) {
							std::vector<unsigned char> candidate(transformed, transformed + row_size);
							candidates[current_row + 1].push_back(candidate);
						}
					}
				}
			}
		}
		return 1;
	}
	return -1;
}
#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>


void ZStabilizer::checkCandidate(ctchar* rows, ctchar* candidate, int iRow)
{
	alignas(32) tchar pTmp[32 * 2];
	//printTable("Tr", candidate, 1, row_size, 2);
#if USE_INTRINSIC
	sortGroupsI(candidate, pTmp, row_size / 2);
#elif
	sortGroups2CPU(pTmp, candidate, row_size / 2);
#endif
#if 0
	alignas(32) tchar n0[32], nj[32];
	for (auto m = 0; m < row_size; m += 2) {
		n0[pTmp[m]] = pTmp[m + 1];
		n0[pTmp[m + 1]] = pTmp[m];
	}
	for (auto j = 0; j < iRow; j++) {
		for (auto m = 0; m < row_size; m += 2) {
			nj[rows[j * row_size + m]] = rows[j * row_size + m + 1];
			nj[rows[j * row_size + m + 1]] = rows[j * row_size + m];
		}
		tchar k = 0;
		for (auto i = 2; i < row_size; i += 2)
		{
			if ((k = nj[n0[k]]) == 0)
				return;
		}
	}
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

// Main generation controller
void ZStabilizer::generateRow(ctchar* rows, ctchar* automorphism, int iRow) {
	// 1. Mark edges already consumed by the first iRow rows
	// Every row is represented as row_size bytes (row_size / 2 pairs of adjacent bytes)
	bool edge_taken[MAX_PLAYER_NUMBER][MAX_PLAYER_NUMBER] = { false };
	int nRows = bUseKSolve ? 3 : iRow;
	for (int r = 0; r < nRows; ++r) {
		auto* row = rows + (r * row_size);
		for (int i = 0; i < row_size; i += 2) {
			uint8_t u = row[i];
			uint8_t v = row[i + 1];
			edge_taken[u][v] = true;
			edge_taken[v][u] = true;
		}
	}
	// 2. Gather all remaining available edges
	std::vector<Edge> remaining_edges;
	for (uint8_t u = 0; u < row_size; ++u) {
		for (uint8_t v = u + 1; v < row_size; ++v) {
			if (!edge_taken[u][v]) {
				// If bipartite, enforce that one vertex is even and the other is odd
				if (bBipartite) {
					if ((u % 2) == (v % 2)) {
						continue; // Skip illegal same-parity (intra-partition) edge
					}
				}
				remaining_edges.push_back(Edge(u, v));
			}
		}
	}

	// 3. Group remaining edges into orbits under the automorphism
	std::vector<std::vector<Edge>> all_orbits;
	std::vector<bool> visited_edges(remaining_edges.size(), false);

	for (size_t i = 0; i < remaining_edges.size(); ++i) {
		if (visited_edges[i]) continue;

		std::vector<Edge> orbit;
		Edge current = remaining_edges[i];

		// Trace the full orbit cyclic pathway
		while (true) {
			// Check if this edge is already inside our current orbit tracking loop
			auto it = std::find(orbit.begin(), orbit.end(), current);
			if (it != orbit.end()) break;

			orbit.push_back(current);

			// Mark this edge globally as visited so we don't restart a new orbit from it
			auto rem_it = std::find(remaining_edges.begin(), remaining_edges.end(), current);
			if (rem_it != remaining_edges.end()) {
				visited_edges[std::distance(remaining_edges.begin(), rem_it)] = true;
			}

			// Apply automorphism to both vertices to step to the next edge in the orbit
			uint8_t next_u = automorphism[current.u];
			uint8_t next_v = automorphism[current.v];
			current = Edge(next_u, next_v);
		}
		all_orbits.push_back(orbit);
	}

	// 4. Filter orbits down to only those mathematically capable of fitting inside Row
	std::vector<std::vector<Edge>> valid_orbits;
	std::vector<uint32_t> orbit_masks;

	for (const auto& orbit : all_orbits) {
		if (orbit.size() > row_size / 2) continue; // Row can only hold row_size / 2 edges total

		uint32_t mask = 0;
		bool vertex_conflict = false;

		for (const auto& edge : orbit) {
			uint32_t edge_mask = (1 << edge.u) | (1 << edge.v);
			// If a vertex is already hit within the same orbit, it's an illegal overlapping matching
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
	// 5. Recursive Bitmask Solver with Mandatory (0,v) Sorting Constraint
	auto solve = [&](auto& self, size_t orbit_idx, uint32_t current_vertices, std::vector<size_t>& chosen_indices) -> void {
		if (current_vertices == ((1UL << row_size) - 1)) {
			std::vector<Edge> candidate_edges;
			bool vertex_written[MAX_PLAYER_NUMBER] = { false };
			uint8_t neighbor0 = bBipartite ? iRow * 2 + 1 : iRow + 1;
	
			for (size_t idx : chosen_indices) {
				for (const auto& edge : valid_orbits[idx]) {
					if (!vertex_written[edge.u] && !vertex_written[edge.v]) {
						candidate_edges.push_back(edge);
						vertex_written[edge.u] = true;
						vertex_written[edge.v] = true;
					}
				}
			}

			// FORCE CANONICAL SORTING: Enforces that the edge containing 0 is placed at the front
			std::sort(candidate_edges.begin(), candidate_edges.end(), [](const Edge& a, const Edge& b) {
				return a.u < b.u;
				});

			// Double check that it meets lexicographical starter criteria
			if (candidate_edges[0].u == 0 && candidate_edges[0].v == neighbor0) {
				unsigned char candidate[MAX_PLAYER_NUMBER];
				int byte_idx = 0;
				for (const auto& edge : candidate_edges) {
					candidate[byte_idx++] = edge.u;
					candidate[byte_idx++] = edge.v;
				}
				checkCandidate(rows, candidate, iRow);
			}
			return;
		}

		if (orbit_idx >= orbit_masks.size()) return;

		// Choice Branch A: Skip the current orbit entirely
		self(self, orbit_idx + 1, current_vertices, chosen_indices);

		// Choice Branch B: Take the current orbit (only if it has zero vertex overlaps)
		if ((current_vertices & orbit_masks[orbit_idx]) == 0) {
			chosen_indices.push_back(orbit_idx);
			self(self, orbit_idx + 1, current_vertices | orbit_masks[orbit_idx], chosen_indices);
			chosen_indices.pop_back(); // Backtrack step
		}
		};

	std::vector<size_t> chosen_indices;
	solve(solve, 0, 0, chosen_indices);
}
#if 0
// Structure to track orbit metrics across different row symmetries
struct OrbitInfo {
	std::vector<Edge> edges;
	uint32_t vertex_mask;
	uint64_t edge_mask[4]; // 4 * 64 = 256 bits, enough to map all 240 remaining edges
};
// Call this ONCE per starter - DO NOT LOOP!
std::string genR2 = findGenerator(automorphismsR2);
std::string genR3 = findGenerator(automorphismsR3);

generateRow4And5(rows123, genR2.c_str(), genR3.c_str());
#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>

// Global definitions from your project architecture
#define MAX_PLAYER_NUMBER 24
typedef const char ctchar;

struct Edge {
	uint8_t u, v;
	Edge(uint8_t a, uint8_t b) {
		if (a < b) { u = a; v = b; }
		else { u = b; v = a; }
	}
	bool operator==(const Edge& other) const {
		return u == other.u && v == other.v;
	}
};

struct OrbitInfo {
	std::vector<Edge> edges;
	uint32_t vertex_mask;
	uint64_t edge_mask[4]; // 256 bits to track up to 240 remaining edges
};

// External callback provided by your validation pipeline
void checkCandidate(unsigned char* candidate);

// Helper function to build unified orbit pools across ALL group permutations
std::vector<OrbitInfo> getValidOrbits(const std::vector<Edge>& remaining_edges,
	int edge_id[MAX_PLAYER_NUMBER][MAX_PLAYER_NUMBER],
	const std::vector<std::string>& group_permutations,
	int row_size) {
	std::vector<std::vector<Edge>> all_orbits;
	std::vector<bool> visited_edges(remaining_edges.size(), false);

	// If the group has no non-trivial symmetries, every edge is its own size-1 orbit
	if (group_permutations.empty()) {
		for (const auto& edge : remaining_edges) {
			all_orbits.push_back({ edge });
		}
	}
	else {
		for (size_t i = 0; i < remaining_edges.size(); ++i) {
			if (visited_edges[i]) continue;

			std::vector<Edge> orbit;
			std::vector<Edge> queue;

			queue.push_back(remaining_edges[i]);
			orbit.push_back(remaining_edges[i]);

			// Breadth-First Search (BFS) to map the orbit under ALL active group elements
			size_t head = 0;
			while (head < queue.size()) {
				Edge current = queue[head++];

				for (const auto& perm_str : group_permutations) {
					const unsigned char* perm = reinterpret_cast<const unsigned char*>(perm_str.data());

					// Direct binary index mapping (0 to row_size - 1)
					uint8_t next_u = perm[current.u];
					uint8_t next_v = perm[current.v];
					Edge mapped_edge(next_u, next_v);

					if (std::find(orbit.begin(), orbit.end(), mapped_edge) == orbit.end()) {
						orbit.push_back(mapped_edge);
						queue.push_back(mapped_edge);

						auto rem_it = std::find(remaining_edges.begin(), remaining_edges.end(), mapped_edge);
						if (rem_it != remaining_edges.end()) {
							visited_edges[std::distance(remaining_edges.begin(), rem_it)] = true;
						}
					}
				}
			}
			all_orbits.push_back(orbit);
		}
	}

	// Filter orbits for structural validity (disjoint matchings)
	std::vector<OrbitInfo> valid_orbits;
	for (const auto& orbit : all_orbits) {
		if (orbit.size() > (size_t)(row_size / 2)) continue;

		uint32_t v_mask = 0;
		bool vertex_conflict = false;
		uint64_t e_mask[4] = { 0, 0, 0, 0 };

		for (const auto& edge : orbit) {
			uint32_t edge_v_mask = (1UL << edge.u) | (1UL << edge.v);
			if ((v_mask & edge_v_mask) != 0) {
				vertex_conflict = true;
				break;
			}
			v_mask |= edge_v_mask;

			int id = edge_id[edge.u][edge.v];
			if (id != -1) {
				e_mask[id / 64] |= (1ULL << (id % 64));
			}
		}

		if (!vertex_conflict) {
			OrbitInfo info;
			info.edges = orbit;
			info.vertex_mask = v_mask;
			std::copy(std::begin(e_mask), std::end(e_mask), std::begin(info.edge_mask));
			valid_orbits.push_back(info);
		}
	}
	return valid_orbits;
}

// Main Generation Controller (Call EXACTLY ONCE per 3-row starter - do not loop)
void ZStabilizer::generateRow4And5(ctchar* rows123,
	const std::vector<std::string>& automorphismsR2,
	const std::vector<std::string>& automorphismsR3) {
	// 1. Mark edges already consumed by the first 3 rows
	bool edge_taken[MAX_PLAYER_NUMBER][MAX_PLAYER_NUMBER] = { false };
	for (int r = 0; r < 3; ++r) {
		auto* row = rows123 + (r * row_size);
		for (int i = 0; i < row_size; i += 2) {
			uint8_t u = row[i];
			uint8_t v = row[i + 1];
			edge_taken[u][v] = true;
			edge_taken[v][u] = true;
		}
	}

	// 2. Gather all remaining available edges and assign unique global IDs
	std::vector<Edge> remaining_edges;
	int edge_id[MAX_PLAYER_NUMBER][MAX_PLAYER_NUMBER];
	std::fill(&edge_id[0][0], &edge_id[0][0] + sizeof(edge_id) / sizeof(int), -1);

	for (uint8_t u = 0; u < row_size; ++u) {
		for (uint8_t v = u + 1; v < row_size; ++v) {
			if (!edge_taken[u][v]) {
				if (bBipartite && ((u % 2) == (v % 2))) continue;

				edge_id[u][v] = edge_id[v][u] = (int)remaining_edges.size();
				remaining_edges.push_back(Edge(u, v));
			}
		}
	}

	// 3. Generate independent orbit streams using full group vectors
	std::vector<OrbitInfo> valid_orbits_R4 = getValidOrbits(remaining_edges, edge_id, automorphismsR2, row_size);
	std::vector<OrbitInfo> valid_orbits_R5 = getValidOrbits(remaining_edges, edge_id, automorphismsR3, row_size);

	// 4. Set up target masks and dynamic starter anchors
	uint32_t target_vertices_mask = (1UL << row_size) - 1;
	uint8_t neighbor0_r4 = bBipartite ? 7 : 4;

	// STAGE 2 BITMASK SOLVER: Verifies if Row 5 can exist in the leftover edge space
	auto solveRow5 = [&](auto& self, size_t orbit_idx, uint32_t current_vertices, const uint64_t* r4_edges) -> bool {
		if (current_vertices == target_vertices_mask) return true;
		if (orbit_idx >= valid_orbits_R5.size()) return false;

		// Option A: Skip this Row 5 orbit
		if (self(self, orbit_idx + 1, current_vertices, r4_edges)) return true;

		const auto& orbit = valid_orbits_R5[orbit_idx];

		// Option B: Take orbit if vertices don't intersect AND edges weren't taken by Row 4
		if ((current_vertices & orbit.vertex_mask) == 0) {
			if ((r4_edges[0] & orbit.edge_mask[0]) == 0 && (r4_edges[1] & orbit.edge_mask[1]) == 0 &&
				(r4_edges[2] & orbit.edge_mask[2]) == 0 && (r4_edges[3] & orbit.edge_mask[3]) == 0) {

				if (self(self, orbit_idx + 1, current_vertices | orbit.vertex_mask, r4_edges)) return true;
			}
		}
		return false;
		};

	// STAGE 1 BITMASK SOLVER: Generates valid matchings for Row 4
	auto solveRow4 = [&](auto& self, size_t orbit_idx, uint32_t current_vertices,
		uint64_t* r4_edges, std::vector<size_t>& chosen_indices) -> void {
			if (current_vertices == target_vertices_mask) {
				// Prune instantly if Row 5 can't find a complementary structure in the remaining edge space
				if (solveRow5(solveRow5, 0, 0, r4_edges)) {
					std::vector<Edge> candidate_edges;
					bool vertex_written[MAX_PLAYER_NUMBER] = { false };

					for (size_t idx : chosen_indices) {
						for (const auto& edge : valid_orbits_R4[idx].edges) {
							if (!vertex_written[edge.u] && !vertex_written[edge.v]) {
								candidate_edges.push_back(edge);
								vertex_written[edge.u] = true;
								vertex_written[edge.v] = true;
							}
						}
					}

					// Sort pairs sequentially to enforce your canonical format
					std::sort(candidate_edges.begin(), candidate_edges.end(), [](const Edge& a, const Edge& b) {
						return a.u < b.u;
						});

					if (candidate_edges[0].u == 0 && candidate_edges[0].v == neighbor0_r4) {
						unsigned char candidate[MAX_PLAYER_NUMBER];
						int byte_idx = 0;
						for (const auto& edge : candidate_edges) {
							candidate[byte_idx++] = edge.u;
							candidate[byte_idx++] = edge.v;
						}
						// Streams the unique Row 4 candidate straight to your P1F and canonization engines
						checkCandidate(candidate);
					}
				}
				return;
			}

			if (orbit_idx >= valid_orbits_R4.size()) return;

			// Choice Branch A: Skip the current orbit entirely
			self(self, orbit_idx + 1, current_vertices, r4_edges, chosen_indices);

			// Choice Branch B: Take the current orbit (only if vertices are free)
			const auto& orbit = valid_orbits_R4[orbit_idx];
			if ((current_vertices & orbit.vertex_mask) == 0) {
				r4_edges[0] |= orbit.edge_mask[0]; r4_edges[1] |= orbit.edge_mask[1];
				r4_edges[2] |= orbit.edge_mask[2]; r4_edges[3] |= orbit.edge_mask[3];
				chosen_indices.push_back(orbit_idx);

				self(self, orbit_idx + 1, current_vertices | orbit.vertex_mask, r4_edges, chosen_indices);

				chosen_indices.pop_back(); // Backtrack
				r4_edges[0] &= ~orbit.edge_mask[0]; r4_edges[1] &= ~orbit.edge_mask[1];
				r4_edges[2] &= ~orbit.edge_mask[2]; r4_edges[3] &= ~orbit.edge_mask[3];
			}
		};

	std::vector<size_t> chosen_indices;
	uint64_t r4_edges[4] = { 0, 0, 0, 0 };
	solveRow4(solveRow4, 0, 0, r4_edges, chosen_indices);
}
#endif