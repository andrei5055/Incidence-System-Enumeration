#include "TripleSys.h"

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

	if (/*val[t_u1f] */ param.u1f[0]) {
		auto pU1F = param.u1f[0] + 1;
		for (tchar i = 0; i < *param.u1f[0]; i++) {
			// Iterating through all prescribed combinations of cycle versions. 
			int nElem = 0;
			for (int i = 0; i < MAX_UNIFOM_CONF_LENGTH; i++) {
				if (pU1F[i])
					nElem += pU1F[i];
				else
					break;
			}

			if (nElem != numPlayers) {
				printfRed("*** Incorrect parameters: Sum of all cycle lengths %d is not equal to %s=%d, Exit\n",
					nElem, paramNames[t_numPlayers], numPlayers);
				return false;
			}
			pU1F += MAX_UNIFOM_CONF_LENGTH;
		}
	}

	if (val[t_p1f]) {
		if (groupSize > 3)
		{
			printfRed("*** %s cannot be used with %s=%d. Exit\n", paramNames[t_p1f], paramNames[t_groupSize], groupSize);
			return false;
		}
	}

	const auto nRowStart = val[t_nRowsInStartMatrix];
	if (nRowStart && !val[t_MultiThreading]) {
		printfRed("*** %s=%d should be 0 if UseMultiThreading=0. Exit\n", paramNames[t_nRowsInStartMatrix], nRowStart);
		return false;
	}

	if (nRowStart) {
		const auto nRowRes = val[t_nRowsInResultMatrix];
		if (nRowRes && nRowRes < nRowStart) {
			printfRed("*** When it's not 0, %s=%d should be > %s=%d. Exit\n", paramNames[t_nRowsInResultMatrix], nRowRes, paramNames[t_nRowsInStartMatrix], nRowStart);
			return false;
		}
	}

	return true;
}
