#pragma once
#include "TripleSys.h"
#include <vector>
class ZStabilizer {
private:
	int row_size;
	int min_aut_val;
	// transition_stack[i] = list of automorphisms of the submatrix (rows 1 to i)
	std::vector<std::vector<std::vector<unsigned char>>> transition_stack;
	std::vector<std::vector<std::vector<unsigned char>>> candidates;
	// Fast bitmask solver variables
	uint8_t active_cand_adj[32] = { 0 };
	int current_index[32];
	int num_rows = 0;
	int current_row;
	int player_index;
	std::vector<std::vector<unsigned char>> survivors;
	std::vector<std::vector<unsigned char>> survivors_candidates;
	tchar* transformed = NULL;
	tchar* transformedPre = NULL;
	bool bBipartite = false;
	bool bPrint = true;
	int zebra_size = 2;

public:
	ZStabilizer() {}

	void init(int size, int min_aut, bool bipartite, bool print);
	void clearCurrent(int iRow);
	void addTr(ctchar* tr, ctchar* rows, int iRow);
	void goUp();
	void goDown();
	tchar* getTr(int iRow, int iTr);
	tchar* getCandidate(int iRow);
	void sortCandidates();
	void apply_permutation(tchar* trRow, ctchar* row, ctchar* tr, int n);
	int check(unsigned char* row, unsigned char* preRow, int* playerIndex, int iRow);
	int setPlayerIndex(ctchar* tr, ctchar* co, ctchar* ci, ctchar* ciFrom);
	int setPlayerIndexByPos(ctchar* tr, ctchar* co, ctchar* ciFrom, int ip);
	void checkCandidate(ctchar* rows, ctchar* candidate, int iRow);
	int generateRow(ctchar* rows, ctchar* automorphism, int iRow);
	bool isP1F(ctchar* rows, int iRow);
	bool generateCandidates(ctchar* rows, int iRow);
	// Structure to hold a single edge
	struct Edge {
		uint8_t u, v;

		// Normalize edge representation so u < v for easy comparison
		Edge(uint8_t a, uint8_t b) {
			if (a < b) { u = a; v = b; }
			else { u = b; v = a; }
		}

		bool operator==(const Edge& other) const {
			return u == other.u && v == other.v;
		}
	};
};
#if 1
// Each array represents the mapping: target_permutation[input_element] = output_element

// 1. Shape A: Order-7 Pair (7^2 1^2)
// Cycles: (0 1 2 3 4 5 6)(7 8 9 10 11 12 13)(14)(15)
// Expected Result: Exactly 1 primary master family
const char perm_shape_A[16] = { 1, 2, 3, 4, 5, 6, 0, 8, 9, 10, 11, 12, 13, 7, 14, 15 };

// 2. Shape B: Order-5 Triple (5^3 1^1)
// Cycles: (0 1 2 3 4)(5 6 7 8 9)(10 11 12 13 14)(15)
// Expected Result: Exactly 1 primary master family
const char perm_shape_B[16] = { 1, 2, 3, 4, 0, 6, 7, 8, 9, 5, 11, 12, 13, 14, 10, 15 };

// 3. Shape C: Order-3 Quintuple (3^5 1^1)
// Cycles: (0 1 2)(3 4 5)(6 7 8)(9 10 11)(12 13 14)(15)
// Expected Result: Exactly 4 primary master families
const char perm_shape_C[16] = { 1, 2, 0, 4, 5, 3, 7, 8, 6, 10, 11, 9, 13, 14, 12, 15 };

// 4. Shape D: Order-2 Septuple (2^7 1^2)
// Cycles: (0 1)(2 3)(4 5)(6 7)(8 9)(10 11)(12 13)(14)(15)
// Expected Result: Exactly 23 primary master families
const char perm_shape_D[16] = { 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 14, 15 };
#endif