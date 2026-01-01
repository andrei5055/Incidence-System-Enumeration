#include "TripleSys.h"

extern const char* arrayParamNames[];

bool checkInputParam(const kSysParam &param, const char** paramNames) {
	const auto& val = param.val;
	const auto& strVal = param.strVal;
	const auto numPlayers = val[t_numPlayers];
	const auto vCBMPgraph = val[t_CBMP_Graph];
	const auto bCBMPgraph = vCBMPgraph> 1;
	const auto nRowStart = val[t_nRowsInStartMatrix];
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

	const auto nRowsOut = val[t_nRowsInResultMatrix];
	if (nRowsOut == 1 || nRowsOut > numDays) {
		printfRed("*** %s(%d) must be either 0 or within the range 2:%d\n", paramNames[t_nRowsInResultMatrix], nRowsOut, numDays);
		return false;
	}

	const auto pCycles = param.u1fCycles[0];
	if (/*val[t_u1f] */ pCycles) {
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
		else {
			// When we order the set of matrices, we need to read all of them.
			printfRed("*** Because matrix ordering requires considering all previously found configurations, the input value\n"
				"\"%s\"=%d is ignored and replaced with the maximum integer value (%d)\n", 
				paramNames[t_nMaxNumberOfStartMatrices], val[t_nMaxNumberOfStartMatrices], MaxNumberOfStartMatrices);
			*(int *)(val + t_nMaxNumberOfStartMatrices) = MaxNumberOfStartMatrices - val[t_nFirstIndexOfStartMatrices];
		}
	}
	// current version does not work with nRowStart != val[t_useRowsPrecalculation]
	if (multiThreading && val[t_useRowsPrecalculation] && nRowStart != val[t_useRowsPrecalculation]) {
		printfRed("*** With %s=%d, %s(%d) should be equal %s(%d). Exit\n", paramNames[t_MultiThreading], multiThreading,
			paramNames[t_nRowsInStartMatrix], nRowStart, paramNames[t_useRowsPrecalculation], val[t_useRowsPrecalculation]);
		return false;
	}

	if (val[t_useRowsPrecalculation]) {
		if (groupSize == 2 && val[t_useRowsPrecalculation] != 3 && val[t_useRowsPrecalculation] != 4) {
			printfRed("*** With GroupSize = 2 the value of %s(%d), can be 0,3 or 4. Exit\n",
				paramNames[t_useRowsPrecalculation], val[t_useRowsPrecalculation]);
			return false;
		}
		if (groupSize == 3 && val[t_useRowsPrecalculation] != 3) {
			printfRed("*** With GroupSize = 3 the value of %s(%d), can be 0 or 3. Exit\n",
				paramNames[t_useRowsPrecalculation], val[t_useRowsPrecalculation]);
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

	return true;
}
