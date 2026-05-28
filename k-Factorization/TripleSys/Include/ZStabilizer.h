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
	bool bUseKSolve = false;
	int zebra_size = 2;

public:
	ZStabilizer() {}

	void init(int size, int min_aut, bool useKSolve, bool bipartite, bool print);
	void clearCurrent(int iRow);
	void addTr(ctchar* tr, ctchar* row3, int iRow);
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
	void generateRow(ctchar* rows, ctchar* automorphism, int iRow);
	void generateCandidates(ctchar* rows, int iRow);
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