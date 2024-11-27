#include "TripleSys.h"

CC void CompSolStorage::addCompatibleSolutions(uint jBase, tmask mask, int kMax, const CRowStorage* pRowStorage, long long* pSolMask)
{
	tmask solMask;
	unsigned long iBit = 0;
	do {
		const auto retVal = _BitScanForward64(&iBit, mask);
		ASSERT(!retVal);

		const auto solID = (uint)(((tmask)jBase << SHIFT) + iBit);
		solMask = (tmask)1 << iBit;
		int k = 0;
		CompSol* pCompSol = kMax? NULL : compatibleSolutions(solID, 0);
		for (; k < kMax; k++) {         // for all previously constructed groups of solutions
			auto& solCur = m_solDB[k];
			for (auto& sol : solCur) {	// for all solutions of current group
				const auto solIdx = sol->solIdx();
				if (solMask & pRowStorage->getSolutionMask(solIdx)[jBase]) {
					if (!pCompSol)
						pCompSol = compatibleSolutions(solID, kMax);

					pCompSol->addNeighbor(solIdx, k);
				}
			}
			if (!pCompSol || pCompSol->noNeighborsOnLevel(k)) {
				// The current solution is NOT compatible with any solution of the current group.
				// This solution cannot be used in the matrix being built.
				if (pCompSol) 
					releaseCompatibleSolutions(kMax);

				pSolMask[jBase] ^= (tmask)1 << iBit;
				break;
			}
		}

		if (k >= kMax)
			m_solDB[kMax].push_back(pCompSol);

	} while (mask ^= solMask);
}

void CompSolStorage::removeUnreachableVertices(int rowIdx) {
	for (int i = rowIdx; i; i--) {
		auto& compSol = m_solDB[i];
		auto& compSolPrev = m_solDB[i-1];
		// Mark used solutions of previous set 
		for (const auto sol : compSol) {
			const auto& compSolList = sol->compatibleSolutions(i - 1);
			uint j = 0;
			for (auto id : compSolList[0]) {
				while (id > compSolPrev[j]->solIdx())
					j++;

				if (id == compSolPrev[j]->solIdx())
					compSolPrev[j]->setUsed();
			}
		}

		// Remove unused elements of previous list
		compSolPrev.erase(std::remove_if(compSolPrev.begin(), compSolPrev.end(),
			[](CompSol *pSol) { return !pSol->isUsed(); }), compSolPrev.end());

		// Reset used flag for remaining elements
		for (auto it = compSolPrev.begin(), end = compSolPrev.end(); it != end; ++it) {
			(*it)->setUnused();
		}

		/*for (const auto sol : compSolPrev) {
			if (sol->isUsed())
		}*/
/*		for (int j = i; j--;) {
			;
		}
		*/
	}
}