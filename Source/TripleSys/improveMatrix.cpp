
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


	if (improveResult > 1 || improveResult && PrintImprovedResults) {
		int cntr = 0;
		auto* bRes1 = bResults;
		auto* bRes2 = bResults + lenResult;
		bool flag;
		while (true) {
			flag = m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainMatr);
			if (PrintImprovedResults) {
				if (!cntr) {
					// Output of initial results
					addCanonCall(1);
					outputResults(iDay, (unsigned char*)result());
				}

				outputResults(iDay, bRes1, ++cntr);
			}

			if (improveResult == 1)
				break;   // No need OR impossible to go deeper.

			if (!flag)
				break;  // Further improvement is impossible

			// Swap the the best results buffers
			auto* bRes = bRes1;
			bRes1 = bRes2;
			bRes2 = bRes;
			addCanonCall(0);
			if (m_pCheckCanon->CheckCanonicity(bRes2, iDay + 1, bRes1))
				break; // Current bRes2 result is canonical

		}

		if (pbRes1)
		    *pbRes1 = improveResult == 1? (flag? bRes1 : NULL) : bRes2; // Best improved result
	}

	return true;
}
