#include "TripleSys.h"

extern const char* arrayParamNames[];

bool checkInputParam(const kSysParam &param, const char** paramNames) {
	const auto& val = param.val;
	const auto numPlayers = val[t_numPlayers];
	if (numPlayers > MAX_PLAYER_NUMBER) {
		printfRed("*** Program is compiled for %d players maximum, but %s=%d is used\n", MAX_PLAYER_NUMBER, paramNames[t_numPlayers], numPlayers);
		return false;
	}

	if (numPlayers < 4) {
		printfRed("*** Program is not designed for less then 4 players, but %s=%d is used\n", paramNames[t_numPlayers], numPlayers);
		return false;
	}

	const auto groupSize = val[t_groupSize];
	const auto numGroups = numPlayers / groupSize;
	if (numGroups > MAX_GROUP_NUMBER) {
		printfRed("*** Program is compiled for %d groups of players maximum, but %s/%s = %d/%d = %d is used\n",
			MAX_GROUP_NUMBER, paramNames[t_numPlayers], paramNames[t_groupSize], numPlayers, groupSize, numGroups);
		return false;
	}

	const auto numDays = (numPlayers - 1) / (groupSize - 1);
	if (numDays < 1 || numDays * (groupSize - 1) != numPlayers - 1 ||
		groupSize < 0 || numPlayers / groupSize * groupSize != numPlayers)
	{
		printfRed("*** Incorrect parameters: %s=%d %s=%d, Exit\n", paramNames[t_numPlayers], numPlayers, paramNames[t_groupSize], groupSize);
		return false;
	}

	const auto pCycles = param.u1fCycles[0];
	if (/*val[t_u1f] */ pCycles) {
		auto pU1F = pCycles + 1;
		if (0) {//??? leo if (*pCycles != 1) {
			printfRed("*** Incorrect parameter '%s': this version supports only one cycles set definition, Exit\n", arrayParamNames[0]);
			return false;
		}
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
		if (groupSize > 3)
		{
			printfRed("*** %s cannot be used with %s=%d. Exit\n", paramNames[t_use2RowsCanonization], paramNames[t_groupSize], groupSize);
			return false;
		}
		if (!val[t_u1f])
		{
			printfRed("*** %s cannot be used with %s=%d. Exit\n", paramNames[t_use2RowsCanonization], paramNames[t_u1f], val[t_u1f]);
			return false;
		}
	}

	const auto nRowStart = val[t_nRowsInStartMatrix];
	if (nRowStart && !val[t_MultiThreading]) {
		printfRed("*** %s=%d should be 0 if %s=0. Exit\n", paramNames[t_nRowsInStartMatrix], nRowStart, paramNames[t_MultiThreading]);
		return false;
	}

	if (nRowStart) {
		const auto nRowRes = val[t_nRowsInResultMatrix];
		if (nRowRes && nRowRes < nRowStart) {
			printfRed("*** When it's not 0, %s=%d should be > %s=%d. Exit\n", paramNames[t_nRowsInResultMatrix], nRowRes, paramNames[t_nRowsInStartMatrix], nRowStart);
			return false;
		}
	}
	if (val[t_MultiThreading] == 2 && val[t_nRowsInStartMatrix] != val[t_useRowsPrecalculation]) {
		printfRed("*** With %s=%d, %s(%d) should be equal %s(%d). Exit\n", paramNames[t_MultiThreading], val[t_MultiThreading],
			paramNames[t_nRowsInStartMatrix], val[t_nRowsInStartMatrix], paramNames[t_useRowsPrecalculation], val[t_useRowsPrecalculation]);
		return false;
	}

	if (val[t_useRowsPrecalculation]) {
		if (groupSize == 3 && val[t_useRowsPrecalculation] != 3) {
			printfRed("*** With GroupSize = 3 the value of %s(%d), can be 0 or 3. Exit\n",
				paramNames[t_useRowsPrecalculation], val[t_useRowsPrecalculation]);
			return false;
		}

		if (groupSize == 3 && val[t_useAutForPrecRows] != 3) {
			printfRed("*** With GroupSize = 3 the value of %s(%d), must be 3. Exit\n", 
				paramNames[t_useAutForPrecRows], val[t_useAutForPrecRows]);
			return false;
		}
		if (groupSize == 2 && val[t_useAutForPrecRows] != 2 && val[t_useAutForPrecRows] != 3) {
			printfRed("*** With GroupSize = 2 the value of %s(%d), can be 2 or 3. Exit\n", 
				paramNames[t_useAutForPrecRows], val[t_useAutForPrecRows]);
			return false;
		}
		if (USE_GROUP_4_2_ROWS && val[t_MultiThreading] == 2) {
			printfRed("*** With %s=%d the use of the Aut(M) on 2 rows (USE_GROUP_4_2_ROWS) is not implemented. Exit\n", 
				paramNames[t_MultiThreading], val[t_MultiThreading]);
			return false;
		}
	}

	return true;
}
