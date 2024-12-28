#include "TripleSys.h"

CC CompSolStorage::CompSolStorage(const CRowStorage* const pRowStorage, int lenGroup) :
	m_pRowStorage(pRowStorage) {
	m_nRowMax = pRowStorage->numPlayers() - 2;
	const int nRowsRes = pRowStorage->sysParam()->val[t_nRowsInResultMatrix];
	if (nRowsRes && m_nRowMax > nRowsRes - 1)
		m_nRowMax = nRowsRes - 1;

	m_nGroups = m_nRowMax + 2 - pRowStorage->useCliquesAfterRow();
	if (m_nGroups > 3) {
		m_nGroups -= 3;
		m_ppCompSol = new CompSolSet * [m_nGroups];
		m_ppCompSol[0] = new CompSolSet(lenGroup);
		for (int i = 1; i < m_nGroups; i++) {
			m_ppCompSol[i] = new CompSolSet(lenGroup);
			for (int j = 0; j < lenGroup; j++)
				m_ppCompSol[i]->getObject(j)->Init(i);
		}
		m_solDB = new std::vector<CompSol*>[m_nGroups];
		m_idxUsedSol = new int[m_nGroups];
	}
	else
		m_nGroups = 0;
}

CC CompSolStorage::~CompSolStorage() {
	for (int i = 0; i < m_nGroups; i++)
		delete m_ppCompSol[i];

	delete[] m_ppCompSol;
	delete[] m_solDB;
	delete[] m_idxUsedSol;
}


void out64bits(FILE* f, const char* prefix, const void *pntr, const char* postFix) {
	if (prefix)
		fprintf(f, prefix);
#if 0
	fprintf(f, "%016llx", *(tmask*)pntr);
#else
	const auto pChar = (unsigned char*)pntr;
	for (int i = 0; i < 8; i++)
		fprintf(f, "%02x", pChar[i]);
#endif
	if (postFix)
		fprintf(f, postFix);
}


CC void CompSolStorage::addCompatibleSolutions(uint jBase, tmask& mask, int kMax)
{
	long long solMask;
	unsigned long iBit = 0;
	do {
		_BitScanForward64(&iBit, mask);
		const auto solID = (uint)(((long long)jBase << 6) + iBit);
		solMask = (long long)1 << iBit;
		CompSol* pCompSol = kMax? NULL : compatibleSolutions(solID, 0);
		int k = 0;
		for (; k < kMax; k++) {         // for all previously constructed groups of solutions
			for (const auto& sol : m_solDB[k]) {	// for all solutions of current group
				if (solMask & m_pRowStorage->getSolutionMask(sol->solIdx())[jBase]) {
					if (!pCompSol)
						pCompSol = compatibleSolutions(solID, kMax);

					pCompSol->addNeighbor(sol, k);
				}
			}

			if (!pCompSol)
				break;

			if (pCompSol->noNeighborsOnLevel(k)) {
				// The current solution is NOT compatible with any solution of the current group.
				// This solution cannot be used in the matrix being built.
				releaseCompatibleSolutions(kMax);
				break;
			}
		}

		if (k >= kMax)
			m_solDB[kMax].push_back(pCompSol);

	} while (mask ^= solMask);
}

bool CompSolStorage::removeUnreachableVertices(int rowIdx) {
	for (int i = rowIdx; i;) {
		const auto& compSol = m_solDB[i--];
		auto& compSolPrev = m_solDB[i];
		// Mark used solutions of previous set 
		for (const auto sol : compSol) {
			const auto& compSolList = sol->compatibleSolutions(i);
			uint j = 0;
			for (const auto solPrev : compSolList) {
				const auto id = solPrev->solIdx();
				while (id > compSolPrev[j]->solIdx())
					j++;

				if (id == compSolPrev[j]->solIdx())
					compSolPrev[j++]->setUsed();
			}
		}

		// Remove unused elements of previous list
		if (i) {
			for (int j = i + 2; j <= rowIdx; j++) {
				auto& compSolNext = m_solDB[j];
				for (const auto sol : compSolNext) {
					auto compSolList = sol->compatibleSolutions(i);
					compSolList.erase(std::remove_if(compSolList.begin(), compSolList.end(),
						[](CompSol* pSol) { return !pSol->isUsed(); }), compSolList.end());
				}
			}
		}

		compSolPrev.erase(std::remove_if(compSolPrev.begin(), compSolPrev.end(),
			[](CompSol *pSol) { return !pSol->isUsed(); }), compSolPrev.end());

		if (!compSolPrev.size())  // There are no solutions on that level
			return false;

		// Reset used flag for remaining elements
		for (auto it = compSolPrev.begin(), end = compSolPrev.end(); it != end; ++it) {
			(*it)->setUnused();
		}
	}

	return true;
}

bool CompSolStorage::ConstructCompatibleSolutionGraph(tmask* pToA, int iRow)
{
	auto pRowSolutionMasksIdx = m_pRowStorage->rowSolutionMasksIdx();
	auto pRowSolutionMasks = m_pRowStorage->rowSolutionMasks();
	releaseSolDB();
	int kMax = 0;
	auto i = iRow + 1;

	auto jMax = pRowSolutionMasksIdx[iRow];
	for (; i <= m_nRowMax; i++, kMax++) {
		auto j = jMax;
		jMax = pRowSolutionMasksIdx[i];

		// Check left, middle and right parts of the solution interval for i-th row
		auto mask = pRowSolutionMasks[i - 1];
		if (mask && (mask &= pToA[j++])) {
			// at least one solution masked by left part of the interval is still valid
			addCompatibleSolutions(j - 1, mask, kMax);
		}

		// middle part
		while (true) {
			while (j < jMax && !pToA[j])
				j++;

			if (j >= jMax)
				break;

			addCompatibleSolutions(j, pToA[j], kMax);
		}

		mask = pRowSolutionMasks[i];
		if (mask && (mask = (~mask) & pToA[jMax]))
			addCompatibleSolutions(j, mask, kMax);

		if (kMax > 1 && !removeUnreachableVertices(kMax))
			return false;
	}

	return true;
}

CC bool CompSolStorage::completeMatrix(tchar* row, tchar* neighbors, int nRows, int iRow) {
	const auto numPlayers = m_pRowStorage->numPlayers();
	const auto iMax = nRows - iRow - 1;
	int i, j = m_idxUsedSol[i = iMax];
	if (j < 0)
		i = 0;

	do {
		const CompSol* sol;
		do {
			// Looking for next solution to be used on i-th level
			while (m_solDB[i].size() <= ++j) {
				if (!i--)
					return false;

				j = m_idxUsedSol[i];
			}

			// Verify that the current solution is compatible 
			// with all previously chosen solutions.
			sol = m_solDB[i][j];
			int l, k = i;
			while (k--) {
				l = k;
				do {
					const auto& compSol = sol->compatibleSolutions(l);
					const auto it = find(compSol.begin(), compSol.end(), currentSolution(l));
					if (it == compSol.end())
						break;
				} while (l--);

				if (l >= 0)
					break; // At least one of previous solution is not compatible with a current one.
			}

			if (k < 0)
				break;     // All previous solution are compatible with a current one.

			m_idxUsedSol[i] = -1;
		} while (true);

		m_idxUsedSol[i] = j;
		j = -1; // When moving to the next matrix row, always start with the first solution.
		const auto* pObj = m_pRowStorage->getObject(sol->solIdx() + m_pRowStorage->numRecAdj());
		const auto shift = (iRow + i) * numPlayers;
		memcpy(row + shift, pObj, numPlayers);
		memcpy(neighbors + shift, pObj + numPlayers, numPlayers);
	} while (++i <= iMax);

	return true;
}