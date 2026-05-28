#include "TripleSys.h"

extern const char* arrayParamNames[];
/**
 * Returns false if the given cycle structure is theoretically impossible
 * for a Uniform 1-Factorization (U1F).
 *
 * @param cycles     Array of cycle lengths (e.g., {4, 6, 16})
 * @param nCycles    Number of cycles in the array
 * @param bBipartite True if the graph is K_{n,n}, false if it is K_{2n}
 */
bool isTheoreticalPossible(unsigned char* cycles, int nCycles, bool bBipartite) {
	int totalVertices = 0;

	for (int i = 0; i < nCycles; ++i) {
		// 1. Every cycle in the union of two 1-factors must be EVEN.
		// The union of two matchings is bipartite, so it cannot have odd cycles.
		if (cycles[i] % 2 != 0) return false;

		totalVertices += cycles[i];
	}

	// 2. A 1-factorization requires an even number of total vertices.
	if (totalVertices % 2 != 0) return false;

	// 3. Kotzig's Theorem for Complete Bipartite Graphs (K_{n,n})
	// A P1F (one cycle of length 2n) exists in K_{n,n} ONLY if n is odd.
	if (bBipartite) {
		int n = totalVertices / 2;
		// If it's a Perfect 1-Factorization (only one cycle)
		if (nCycles == 1) {
			// If n is even, K_{n,n} cannot have a P1F.
			if (n % 2 == 0) return false;
		}
		/*
		 * Checks if a cycle structure is theoretically possible for K_{n,n}.
		 * In bipartite graphs, the union of two 1-factors must correspond
		 * to an EVEN permutation if the 1-factorization is to be consistent.
		 */
		int evenPermutationCycleCount = 0;
		for (int k = 0; k < nCycles; k++) {
			// A graph cycle of length G actually represents 
			// a permutation cycle of length L = G / 2.
			int L = cycles[k] / 2;
			if (L % 2 == 0) {
				evenPermutationCycleCount++;
			}
		}
		if (evenPermutationCycleCount % 2 != 0) {
			return false;
		}
	}
	return true;
}
bool checkInputParam(const kSysParam &param, const char** paramNames) {
	const auto& val = param.val;
#ifndef USE_CUDA
    if (val[t_useGPU]) {
		printfRed("*** Program is compiled without GPU support, but %s=%d is used\n", paramNames[t_useGPU], val[t_useGPU]);
		return false;
	}
#endif

	const auto& strVal = param.strVal;
	const auto numPlayers = val[t_numPlayers];
	const auto vCBMPgraph = val[t_CBMP_Graph];
	const auto bCBMPgraph = vCBMPgraph > 1;
	const auto nRowStart = val[t_nRowsInStartMatrix];
	const auto nPrecalcRows = val[t_useRowsPrecalculation];
	const auto pCycles = param.u1fCycles[0];
	bool bP1f = (!pCycles || pCycles[1] == numPlayers) && !val[t_allowUndefinedCycles];
	int nCycles = pCycles ? pCycles[0] : 1;
	if (numPlayers > MAX_PLAYER_NUMBER) {
		printfRed("*** Program is compiled for %d players maximum, but %s=%d is used\n", MAX_PLAYER_NUMBER, paramNames[t_numPlayers], numPlayers);
		return false;
	}

	if (numPlayers < 4) {
		printfRed("*** Program is not designed for less then 4 players, but %s=%d is used\n", paramNames[t_numPlayers], numPlayers);
		return false;
	}

	const auto groupSize = val[t_groupSize];
	if (groupSize < 2) {
		printfRed("*** %s(%d) must be > 1, Exit\n", paramNames[t_groupSize], groupSize);
		return false;
	}

	const auto numGroups = numPlayers / groupSize;
	if (numGroups < groupSize) {
		printfRed("*** %s/%s (%d/%d) must be >= %s, Exit\n",
			paramNames[t_numPlayers], paramNames[t_groupSize], numPlayers, groupSize, paramNames[t_groupSize]);
		return false;
	}
	if (numGroups > MAX_GROUP_NUMBER) {
		printfRed("*** Program is compiled for %d groups of players maximum, but %s/%s = %d/%d = %d is used\n",
			MAX_GROUP_NUMBER, paramNames[t_numPlayers], paramNames[t_groupSize], numPlayers, groupSize, numGroups);
		return false;
	}

	int numDays;
	
	if (val[t_generateMatrixExample]) {
		if (val[t_MultiThreading] > 0) {
			printfRed("*** Program can't generate matrix(GenerateMatrixExample=%d) in multithread mode, Exit\n",
				val[t_generateMatrixExample]);
			return false;
		}

		if (!bCBMPgraph) {
			printfRed("*** Program can generate (GenerateMatrixExample=%d) only n-partite graph, but CBMP_Graph=%d, Exit\n",
				val[t_generateMatrixExample], vCBMPgraph);
			return false;
		}
		if (val[t_generateMatrixExample] > 1 && vCBMPgraph != 2) {
			printfRed("*** With GenerateMatrixExample > 1 this program supports only 2-partite graphs (CBMP_Graph=2). Exit\n");
			return false;
		}
	}

	if (strVal[t_InputDataFileName]->length() && (vCBMPgraph != 2  || nRowStart != numPlayers / 2)) {
		printfRed("*** Incorrect parameters: InputDataFileName('%s') can be used only with %s=2; NRowsInStartMatrix=%d. Exit\n",
			strVal[t_InputDataFileName]->c_str(), paramNames[t_CBMP_Graph], numPlayers / 2);
		return false;
	}

	if (bCBMPgraph) {
		numDays = numPlayers / groupSize;
		if (numDays < 1 || numDays * groupSize != numPlayers || groupSize < 0)
		{
			printfRed("*** Incorrect parameters: %s=%d %s=%d %s=%d, Exit\n", 
				paramNames[t_numPlayers], numPlayers, paramNames[t_groupSize], groupSize, paramNames[t_CBMP_Graph], vCBMPgraph);
			return false;
		}
	}
	else {
		numDays = (numPlayers - 1) / (groupSize - 1);
		if (numDays < 1 || numDays * (groupSize - 1) != numPlayers - 1 ||
			groupSize < 0 || numPlayers / groupSize * groupSize != numPlayers)
		{
			printfRed("*** Incorrect parameters: %s=%d %s=%d, Exit\n", paramNames[t_numPlayers], numPlayers, paramNames[t_groupSize], groupSize);
			return false;
		}
	}

	if (val[t_nFirstIndexOfStartMatrices] < 1) {
		printfRed("*** %s(%d) must be > 0\n", paramNames[t_nFirstIndexOfStartMatrices], val[t_nFirstIndexOfStartMatrices]);
		return false;
	}

	if (bP1f && val[t_useCompatibilityCheck] > 1) {
		printfRed("*** %s(%d) > 1 can be used only with P1F\n", paramNames[t_useCompatibilityCheck], val[t_useCompatibilityCheck]);
		return false;
	}
	if (pCycles) {
		auto pU1F = pCycles + 1;
		for (tchar j = 0; j < *pCycles; j++) {
			// Iterating through all prescribed combinations of cycle versions. 
			int nElem = 0;
			for (int i = 0; i < MAX_CYCLES_PER_SET; i++) {
				if (pU1F[i])
				{
					if (i > 0 && pU1F[i - 1] > pU1F[i]) {
						printfRed("*** Incorrect parameter '%s': cycle length(%d) is less than next cycle length(%d), Exit\n",
							arrayParamNames[0], pU1F[i - 1], pU1F[i]);
						return false;
					}
					nElem += pU1F[i];
				}
				else
					break;
			}

			if (nElem != numPlayers) {
				printfRed("*** Incorrect parameters: Sum of all cycle lengths %d is not equal to %s=%d, Exit\n",
					nElem, paramNames[t_numPlayers], numPlayers);
				return false;
			}
			pU1F += MAX_CYCLES_PER_SET;
		}
	}

	if (val[t_use2RowsCanonization]) {
		if (groupSize > 3 && vCBMPgraph <= 1)
		{
			printfRed("*** %s with CBMP_Graph < 2 cannot be used with %s > 3. Exit\n", paramNames[t_use2RowsCanonization], paramNames[t_groupSize]);
			return false;
		}
		if (!val[t_u1f])
		{
			printfRed("*** %s cannot be used with %s=%d. Exit\n", paramNames[t_use2RowsCanonization], paramNames[t_u1f], val[t_u1f]);
			return false;
		}
	}

	const auto multiThreading = val[t_MultiThreading];
	if (nRowStart) {
		const auto nRowRes = val[t_nRowsInResultMatrix];	
		if (nRowRes && nRowRes < nRowStart) {
			printfRed("*** When it's not 0, %s=%d should be > %s=%d. Exit\n", paramNames[t_nRowsInResultMatrix], nRowRes, paramNames[t_nRowsInStartMatrix], nRowStart);
			return false;
		}
		
		if (!val[t_orderMatrices]) {
			if (!multiThreading) {
				printfRed("*** %s=%d should be 0 if %s=0. Exit\n", paramNames[t_nRowsInStartMatrix], nRowStart, paramNames[t_MultiThreading]);
				return false;
			}
		}
		else if (val[t_nMaxNumberOfStartMatrices] < MaxNumberOfStartMatrices) {
			// When we order the set of matrices, we need to read all of them.
			printfYellow("*** Because matrix ordering requires considering all previously found configurations, the input value\n"
				"%s=%d is ignored and replaced with the default maximum value (%d)\n", 
				paramNames[t_nMaxNumberOfStartMatrices], val[t_nMaxNumberOfStartMatrices], MaxNumberOfStartMatrices);
			*(int *)(val + t_nMaxNumberOfStartMatrices) = MaxNumberOfStartMatrices - val[t_nFirstIndexOfStartMatrices] + 1;
		}
	}
	// current version does not work with nRowStart != nPrecalcRows
	if (multiThreading && nPrecalcRows && nRowStart != nPrecalcRows) {
		printfRed("*** With %s=%d, %s(%d) should be equal %s(%d). Exit\n", paramNames[t_MultiThreading], multiThreading,
			paramNames[t_nRowsInStartMatrix], nRowStart, paramNames[t_useRowsPrecalculation], nPrecalcRows);
		return false;
	}

	if (nPrecalcRows) {
		if (groupSize == 2 && nPrecalcRows != 0 && !(nPrecalcRows >= 3 && nPrecalcRows <= MAX_PRECALC_ROWS)) {
			printfRed("*** With GroupSize = 2 the value of %s(%d), can be 0, or from 3 to %d. Exit\n",
				paramNames[t_useRowsPrecalculation], nPrecalcRows, MAX_PRECALC_ROWS);
			return false;
		}
		if (groupSize == 3 && nPrecalcRows != 3) {
			printfRed("*** With GroupSize = 3 the value of %s(%d), can be 0 or 3. Exit\n",
				paramNames[t_useRowsPrecalculation], nPrecalcRows);
			return false;
		}

		if (USE_GROUP_4_2_ROWS && multiThreading == 2) {
			printfRed("*** With %s=%d the use of the Aut(M) on 2 rows (USE_GROUP_4_2_ROWS) is not implemented. Exit\n", 
				paramNames[t_MultiThreading], multiThreading);
			return false;
		}
	}

	if (multiThreading && (nRowStart < 2 || nRowStart > numDays))
	{
		printfRed("*** %s(%d) with %s must be in range 2:%d\n",
			paramNames[t_nRowsInStartMatrix], nRowStart, paramNames[t_MultiThreading], numDays);
		return false;
	}
	if (val[t_groupSize] == 2 && val[t_any2RowsConvertToFirst2]) {
		for (int i = 0; i < nCycles; i++) {
			tchar cycle[2] = { 0 }, * cycles = cycle;
			cycle[0] = numPlayers;
			int iCycleLength = 1;
			if (pCycles) {
				cycles = pCycles + 1 + i * MAX_CYCLES_PER_SET;
				for (iCycleLength = 1; iCycleLength < MAX_CYCLES_PER_SET && cycles[iCycleLength] != 0; iCycleLength++)
					if (cycles[iCycleLength] == 0)
						break;
			}
			if (!isTheoreticalPossible(cycles, iCycleLength, vCBMPgraph == 2)) {
				printfRed("*** Cycle {%d", cycles[0]);
				for (int j = 1; j < iCycleLength; j++) printfRed(",%d", cycles[j]);
				printfRed("} is theoretically impossible\n");
				return false;
			}
		}
	}

	return true;
}
