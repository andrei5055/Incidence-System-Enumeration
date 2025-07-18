#include "TripleSys.h"
#include "CheckCanon.h"

bool alldata::improveMatrix(int improveResult, tchar* bResults, const int lenResult, tchar **pbRes1)
{
	const auto nDays = iDay;
	const auto inputMatrix = (unsigned char*)result();
	addCanonCall(0);

	m_pCheckCanon->setAllData(this);
	if (m_pCheckCanon->CheckCanonicity(inputMatrix, nDays, &m_groupIndex, bResults))
		return false;

	const int improveResultMax = IMPROVE_RESULT_MAX;
	if (improveResult > 1 || improveResult && PrintImprovedResults) {
		int cntr = 0;
		auto* bRes1 = bResults;
		auto* bRes2 = bResults + lenResult;
		unsigned char* bRes3 = NULL;
		bool flag;
		do {
			flag = m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainMatr);
			if (PrintImprovedResults) {
				if (!cntr) {
					// Output of initial results
					addCanonCall(1);
					outputResults(nDays, inputMatrix);
				}

				outputResults(nDays, bRes1, ++cntr);
			}

			if (improveResult == 1 || !flag || ++improveResult > improveResultMax)
				break;   // No need OR further improvement is impossible

#if CHECK_PERMUTATIONS
			int errLine, errGroup, dubLine;
			char lnks[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
			if (!_CheckMatrix((char *)bRes1, nDays, numPlayers(), m_groupSize, lnks, true, &errLine, &errGroup, &dubLine, sysParam()->numFactors()))
				outputError();
#endif
			// Swap the the best results buffers
			auto* bRes = bRes1;
			bRes1 = bRes2;
			bRes3 = bRes2 = bRes;
			addCanonCall(0);
		} while (!m_pCheckCanon->CheckCanonicity(bRes2, nDays, &m_groupIndex, bRes1));

		const auto bestResult = improveResult == 1 ? (flag ? bRes1 : NULL)
								: improveResult > improveResultMax? (flag? bRes1 : bRes3)
								: bRes3; // Best improved result
		if (pbRes1)
			*pbRes1 = bestResult;

#if CHECK_PERMUTATIONS && 0
		if (!m_pCheckCanon->CheckPermutations(inputMatrix, bestResult, nDays, sysParam()->numFactors())) {
			outputError();
			abort();
		}
#endif
	}

	return true;
}

