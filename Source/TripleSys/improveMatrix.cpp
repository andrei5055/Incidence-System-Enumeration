
#include <iostream>
#include <assert.h>
#include "TripleSys.h"

#ifdef CD_TOOLS
#include "../CanonicityChecker.h"
#else
#include "CheckCanon.h"
#endif

bool alldata::improveMatrix(int improveResult, unsigned char* bResults, const int lenResult, unsigned char **pbRes1)
{
	const auto nDays = iDay;
	const auto inputMatrix = (unsigned char*)result();
	addCanonCall(0);
	if (m_pCheckCanon->CheckCanonicity(inputMatrix, nDays, bResults))
		return false;

	//improveResult = 1;
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
			char lnks[21 * 21];
			if (!_CheckMatrix((char *)bRes1, nDays, numPlayers(), lnks, true, &errLine, &errGroup, &dubLine))
				outputError();
#endif
			// Swap the the best results buffers
			auto* bRes = bRes1;
			bRes1 = bRes2;
			bRes3 = bRes2 = bRes;
			addCanonCall(0);
		} while (!m_pCheckCanon->CheckCanonicity(bRes2, nDays, bRes1));

		const auto bestResult = improveResult == 1 ? (flag ? bRes1 : NULL)
								: improveResult > improveResultMax? (flag? bRes1 : bRes3)
								: bRes3; // Best improved result
		if (pbRes1)
			*pbRes1 = bestResult;

#if CHECK_PERMUTATIONS && 0
		if (!m_pCheckCanon->CheckPermutations(inputMatrix, bestResult, nDays)) {
			outputError();
			abort();
		}
#endif
	}

	return true;
}

