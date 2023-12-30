
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
	addCanonCall(0);
	if (m_pCheckCanon->CheckCanonicity((unsigned char*)result(), iDay + 1, bResults))
		return false;

	improveResult = 2;
	const int improveResultMax = 3;
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
					outputResults(iDay, (unsigned char*)result());
				}

				outputResults(iDay, bRes1, ++cntr);
			}

			if (improveResult == 1 || !flag || ++improveResult > improveResultMax)
				break;   // No need OR further improvement is impossible

			// Swap the the best results buffers
			auto* bRes = bRes1;
			bRes1 = bRes2;
			bRes3 = bRes2 = bRes;
			addCanonCall(0);
		} while (!m_pCheckCanon->CheckCanonicity(bRes2, iDay + 1, bRes1));

		if (pbRes1) {
			*pbRes1 = improveResult == 1 ? (flag ? bRes1 : NULL) 
				    : flag ? bRes1 : bRes3; // Best improved result
		}
	}

	return true;
}

