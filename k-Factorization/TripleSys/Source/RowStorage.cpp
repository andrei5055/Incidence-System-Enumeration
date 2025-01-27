#include "TripleSys.h"

CC CRowStorage::CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects, const alldata* pAllData) :
	m_pSysParam(pSysParam), m_numPlayers(numPlayers),
	m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]),
	m_numDaysResult(pAllData ? pAllData->numDaysResult() : numPlayers - 1), m_pAllData(pAllData),
	m_bUseCombinedSolutions(pSysParam->val[t_useCombinedSolutions]),
	m_step(pSysParam->val[t_MultiThreading] == 2 ? pSysParam->val[t_numThreads] : 1),
	CStorage<tchar>(numObjects, 2 * numPlayers) {
	m_numObjectsMax = numObjects;
	m_pPlayerSolutionCntr = new uint[2 * numPlayers];
	m_pNumLongs2Skip = m_pPlayerSolutionCntr + numPlayers;
	initMaskStorage(numObjects);
	m_lenMask = m_pMaskStorage->lenObject();
	const auto useCliquesAfterRow = pSysParam->val[t_useSolutionCliquesAfterRow];
	m_useCliquesAfterRow = useCliquesAfterRow ? useCliquesAfterRow : m_numDaysResult;
	memset(m_pRowsCompatMasks, 0, sizeof(m_pRowsCompatMasks));
	m_fRowToBitmask = m_pAllData ? &CRowStorage::rowToBitmask3 : &CRowStorage::rowToBitmask2;
	m_fSolutionInterval = m_pAllData ? &CRowStorage::solutionInterval3 : &CRowStorage::solutionInterval2;
	m_lenDayResults = m_numDaysResult + 1;
}

CC CRowStorage::~CRowStorage() {
	delete[] m_pPlayerSolutionCntr;
	for (int i = countof(m_pRowsCompatMasks); i--;)
		delete[] m_pRowsCompatMasks[i];

	delete m_pMaskStorage;
	delete[] m_pRowSolutionMasks;
	delete[] m_pRowSolutionMasksIdx;
}

CC bool CRowStorage::p1fCheck2(ctchar* neighborsi, ctchar* neighborsj) const {
	uint checked = 0;
	for (tchar m = 0; m < m_numPlayers; m++)
	{
		if (!(checked & (1 << m))) {
			tchar k = m;
			for (int i = 2; i <= m_numPlayers; i += 2)
			{
				if ((k = neighborsj[neighborsi[k]]) == m)
					return i == m_numPlayers;

				checked |= 1 << k;
			}
		}
	}
	return false;
}

CC void CRowStorage::addRow(ctchar* pRow, ctchar* pNeighbors) {
	if (m_numObjects == m_numObjectsMax) {
		reallocStorageMemory(m_numObjectsMax <<= 1);
		m_pMaskStorage->reallocStorageMemory(m_numObjectsMax);
	}
#if 0
	FOPEN_F(f, "aaa.txt", m_numObjects ? "a" : "w");
	char buf[256], * pBuf = buf;
	for (int i = 0; i < m_numPlayers; i++)
		SPRINTFD(pBuf, buf, "%2d ", pRow[i]);

	fprintf(f, "%3d: %s\n", m_numObjects, buf);
	FCLOSE_F(f);
#endif
	(this->*m_fRowToBitmask)(pRow, (tmask*)(m_pMaskStorage->getObject(m_numObjects)));
	auto* pntr = getObject(m_numObjects++);
#if !USE_CUDA
	const auto ptr = pntr - 2 * m_numPlayers + 1;
#endif
#if 1
	ASSERT_(m_numObjects > 1 && (pRow[1] != *ptr && pRow[1] != *ptr + 1) && (m_playersMask && (ll) 1 << *ptr),
		printfRed("\nError in code: and pRow[1](%d) != %d, and pRow[1] != %d and m_playersMask = %llx & (1 << %d) != 0\n",
			pRow[1], *ptr, *ptr + 1, m_playersMask, *ptr);
		exit(1)
	);
#else
	// NOTE: This conditions are valid for groupSize == 2 or groupSize == 3 and m_numPreconstructedRows == 3.
	//       For the other cases they are more complicated.
	if (m_numObjects > 1 && (pRow[1] != *ptr && pRow[1] != *ptr + 1) && (m_playersMask && (ll)1 << *ptr) {
		printfRed("\nError in code: and pRow[1](%d) != %d, and pRow[1] != %d and m_playersMask = %x & (1 << %d) != 0\n",
			pRow[1], *ptr, *ptr + 1, m_playersMask, *ptr);
		exit(1);
	}
#endif
#if 0
	static int cntr = 0;
	FOPEN_F(f, "aaa.txt", cntr? "a" : "w");
	char buf[256], *pBuf = buf;
	for (int i = 0; i < m_numPlayers; i++)
		SPRINTFD(pBuf, buf, "%2d ", pRow[i]);

	fprintf(f, " %2d:  %s\n", ++cntr, buf);
	FCLOSE_F(f);
#endif
	memcpy(pntr, pRow, m_numPlayers);
	memcpy(pntr + m_numPlayers, pNeighbors, m_numPlayers);
	m_pPlayerSolutionCntr[pRow[1] - 1]++;  // Increasing the number of solutions for player pRow[1] 
	                                       // NOTE: for groupSize = 2, this number is equal to the number of solutions for (pRow[1]-1)-th row
}

CC bool CRowStorage::checkCompatibility(ctchar* neighborsi, const ll* rm, uint idx) const {
	// Let's check if the masks are mutually compatible
	auto* pMask = (const ll*)(m_pMaskStorage->getObject(idx));
	int j = m_lenMask >> 3;
	while (j-- && !(rm[j] & pMask[j]));

	if (j >= 0)
		return false;

	const auto pObj = getObject(idx);
	if (m_pAllData)
		return m_pAllData->p1fCheck3(neighborsi - m_numPlayers, pObj, neighborsi, pObj + m_numPlayers);

	return p1fCheck2(neighborsi, pObj + m_numPlayers);
}

CC bool CRowStorage::maskForCombinedSolutions(tmask* pMaskOut, uint & solIdx) const {
	// Constructing a mask for the "combined" solution, which represents 
	// the combination of solutions for the first two "non-predefined" rows.
	const auto n = m_numRec[1];
	do {
		auto* rm1 = getSolutionMask(solIdx / n - m_numRecAdj);
		const auto idx2 = solIdx % n;
		if (CHECK_MASK_BIT(rm1, idx2)) {
			const auto numLongs2Skip = m_pNumLongs2Skip[numPreconstructedRows()];
			const auto numBites2Skip = numLongs2Skip << 3;
			const auto len = m_numSolutionTotalB - numBites2Skip;
			auto pMaskOutStart = (tmask *)((tchar *)pMaskOut + numBites2Skip);
			auto pCompSol = (tmask*)((tchar *)pMaskOut - m_numSolutionTotalB);
			memcpy(pMaskOutStart, ((const ll*)getSolutionMask(solIdx % n)) + numLongs2Skip, len);
			memcpy((long long*)pCompSol + numLongs2Skip, pMaskOutStart, len);

			const auto last = m_numSolutionTotal - m_numRecAdj;
			auto first = m_numRecAdj2 - m_numRecAdj;
			const auto lastB = IDX(last);
			while (true) {
				// Skip all bytes/longs equal to 0
				auto firstB = first >> SHIFT;
				while (firstB < lastB && !pCompSol[firstB])
					firstB++;

				if (firstB >= lastB)
					break;
#if USE_64_BIT_MASK
				unsigned long iBit;
				_BitScanForward64(&iBit, pCompSol[firstB]);
#else
				const auto iBit = this->firstOnePosition(pCompSol[firstB]);
#endif
				if ((first = (firstB << SHIFT) + iBit) >= last)
					break;

				pCompSol[firstB] ^= (tmask)1 << iBit;
				if (!CHECK_MASK_BIT(rm1, first))
					RESET_MASK_BIT(pMaskOut, first);     //  reset bit
			}

			return true;
		}
	} while ((solIdx += m_step) < m_lastInFirstSet);

	return false;
}

CC void CRowStorage::generateCompatibilityMasks(tmask* pMaskOut, uint solIdx, uint idx) const {
	auto* rm = (const ll*)m_pMaskStorage->getObject(solIdx);
	ctchar* pSolution;
	const auto pNeighbors = (pSolution = getObject(solIdx)) + m_numPlayers;
	do {
		if (checkCompatibility(pNeighbors, rm, idx)) {
			const auto newIdx = idx - m_numRecAdj;
			SET_MASK_BIT(pMaskOut, newIdx);     // 1 - means OK
		}
	} while (++idx < m_numSolutionTotal);

	if (m_pAllData) {
		// If groupSize > 2, we also need to create mask which will keep  
		// the information regarding the players used by current solution.
		// We will store it as 0's of corresponding bites.
		auto* pMaskOutLong = (ll*)pMaskOut + m_lenSolutionMask - 1;
		*pMaskOutLong = m_playersMask;
		for (auto i = m_pAllData->groupSize(); --i;)
			*pMaskOutLong ^= (ll)1 << pSolution[i];
	}
}

CC void CRowStorage::initCompatibilityMasks(ctchar* u1fCycles) {
	m_u1fCycles = u1fCycles;
	for (int i = 1; i < m_numPlayers; i++)
		m_pPlayerSolutionCntr[i] += m_pPlayerSolutionCntr[i - 1];

	// Define the number of first long long's we don't need to copy to the next row.
	memset(m_pNumLongs2Skip, 0, m_numPlayers * sizeof(m_pNumLongs2Skip[0]));
	int i = m_numPreconstructedRows;
	m_pNumLongs2Skip[i] = m_pPlayerSolutionCntr[i] >> 6;
	m_lastInFirstSet = m_numRecAdj = m_pPlayerSolutionCntr[i];
	while (++i < m_numPlayers)
		m_pNumLongs2Skip[i] = ((m_pPlayerSolutionCntr[i] - m_numRecAdj) >> 6);

	m_numSolutionTotal = m_pPlayerSolutionCntr[m_numPlayers - 1];

	const auto useCombinedSolutions = sysParam()->val[t_useCombinedSolutions];
	if (useCombinedSolutions) {
		m_lastInFirstSet *= (m_numRec[1] = ((m_numRecAdj2 = m_pPlayerSolutionCntr[m_numPreconstructedRows + 1]) - m_numRecAdj));
		delete [] m_pRowsCompatMasks[0];
	}

	// Adding additional long long when we use groupSize > 2
	m_numSolutionTotalB = ((m_numSolutionTotal - m_numRecAdj + 7) / 8 + 7) / 8 * 8 + (m_pAllData ? 8 : 0);
	m_lenSolutionMask = m_numSolutionTotalB / sizeof(tmask);

	m_solAdj = useCombinedSolutions ? m_numRecAdj : 0;
	const auto len = m_numSolutionTotal - (m_numRecAdj - m_solAdj);
	delete[] m_pRowsCompatMasks[1];
	m_pRowsCompatMasks[1] = new tmask[len * m_lenSolutionMask];
	memset(m_pRowsCompatMasks[1], 0, len * m_numSolutionTotalB);

	delete[] m_pRowSolutionMasksIdx;
	delete[] m_pRowSolutionMasks;

	if (m_pAllData) {
		// Create a mask to manage players utilized in the predefined rows of the matrix.
		const auto groupSize = m_pAllData->groupSize();
		// Excluding players of the first group from ...
		m_playersMask = (ll)(-1) << groupSize;			// ...first row
		m_playersMask ^= (ll)(-1) << numPlayers();
		auto const* pSolution = m_pAllData->result();
		for (int j = numPreconstructedRows(); --j;) {   // ... remaining pre-constructed rows
			pSolution += m_numPlayers;
			for (auto i = groupSize; --i;)
				m_playersMask ^= (ll)1 << pSolution[i];
		}
	}
	else {
		const auto numDays = numDaysResult();
		m_pRowSolutionMasksIdx = new uint[numDays];

		m_pRowSolutionMasks = new tmask[numDays];
		memset(m_pRowSolutionMasks, 0, numDays * sizeof(m_pRowSolutionMasks[0]));
		m_pRowSolutionMasksIdx[0] = 0;
		m_playersMask = -1;
	}

#if !USE_64_BIT_MASK
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

	tmask* pRowsCompatMasks[] = { m_pRowsCompatMasks[0], m_pRowsCompatMasks[1] };
	tmask* pCompatMask = pRowsCompatMasks[1];
	unsigned int first, last = 0;
	i = m_numPreconstructedRows - 1;
#if 1
	auto availablePlayers = m_playersMask;
	unsigned int rem;
	const auto iMax = m_numPlayers - 1;
	while (++i < iMax && availablePlayers) {
		first = last;
		if (m_pAllData) {
			unsigned long iBit;
			_BitScanForward64(&iBit, availablePlayers);
			last = m_pPlayerSolutionCntr[iBit - 1];
			availablePlayers ^= (ll)1 << iBit;
		}
		else {
			last = m_pPlayerSolutionCntr[i];
		}
		if (m_pRowSolutionMasksIdx) {
			const auto lastAdj = last - m_numRecAdj;
			m_pRowSolutionMasksIdx[i] = lastAdj >> SHIFT;
			if (i == m_numPreconstructedRows || m_pRowSolutionMasksIdx[i] > m_pRowSolutionMasksIdx[i - 1]) {
				if (rem = REM(lastAdj))
					m_pRowSolutionMasks[i] = (tmask)(-1) << rem;
			}
			else {
				// We cannot use code for UseSolutionMasks, because now our code is not ready for such cases
				delete[] m_pRowSolutionMasksIdx;
				delete[] m_pRowSolutionMasks;
				m_pRowSolutionMasksIdx = NULL;
				m_pRowSolutionMasks = NULL;
			}
		}

		if (!first && !useCombinedSolutions) {
			// Skip construction of masks for the first set of solutions.
			// The threads will do this latter.
			continue;
		}

		if (m_pAllData)
			pCompatMask = getSolutionMask(first);

		while (first < last) {
			generateCompatibilityMasks(pCompatMask, first++, last);
			pCompatMask += m_lenSolutionMask;
		}
	}

	if (m_numRecAdj) {
		for (int i = m_numPreconstructedRows; i < m_numPlayers; i++)
			m_pPlayerSolutionCntr[i] -= m_numRecAdj;
	}

#else
	// Calculate the number of mutually compatible pairs of solutions
	int cntrs[8];
	memset(cntrs, 0, sizeof(cntrs));
	unsigned ll fff = 0;
	const auto jMax = m_lenMask >> 3;
	int a = 0;
	while (i < iMax) {
		auto first = last;
		last = m_pPlayerSolutionCntr[i];
		i++;
		while (first < last) {
			auto* rm = (const ll*)m_pMaskStorage->getObject(first);
			const auto pRow = getObject(first++);
			ASSERT(pRow[1] != i);
			const auto pNeighbors = pRow + m_numPlayers;
			auto idx = last - 1;
			while (++idx < m_pPlayerSolutionCntr[i]) {
				// Let's check if the masks are mutually compatible
				auto* pMask = (const ll*)(m_pMaskStorage->getObject(idx));
				int j = jMax;
				while (j-- && !(rm[j] & pMask[j]));

				if (j < 0 && p1fCheck2(pNeighbors, getObject(idx) + m_numPlayers)) {
					cntrs[i / 2]++;
					if (first <= m_numRecAdj && idx < m_numRecAdj2)
						a++;
				}
			}
		}
		fff += cntrs[i / 2];
		last = m_pPlayerSolutionCntr[i++];
	}
#endif

	if (useCombinedSolutions) {
		delete m_pMaskStorage;
		m_pMaskStorage = NULL;
	}
}

CC void CRowStorage::getMatrix(tchar* row, tchar* neighbors, int nRows, uint* pRowSolutionIdx) const {
	auto iRow = numPreconstructedRows();
	uint savedIdx;
	if (m_bUseCombinedSolutions) {
		const auto ind = (savedIdx = pRowSolutionIdx[iRow]) - m_step;
		pRowSolutionIdx[iRow] = ind / numRec(1) + m_step;
		pRowSolutionIdx[iRow + 1] = ind % numRec(1) + 1;
	}

	size_t shift = iRow * m_numPlayers;
	const int adj = numRecAdj() - 1;
	auto* pObj = getObject(pRowSolutionIdx[iRow] - m_step);
	while (true) {
		memcpy(row + shift, pObj, m_numPlayers);
		memcpy(neighbors + shift, pObj + m_numPlayers, m_numPlayers);
		if (++iRow == nRows)
			break;

		pObj = getObject(pRowSolutionIdx[iRow] + adj);
		shift += m_numPlayers;
	}

	if (m_bUseCombinedSolutions)
		pRowSolutionIdx[numPreconstructedRows()] = savedIdx;
}

CC int CRowStorage::initRowUsage(tchar** ppCompatibleSolutions, bool *pUsePlayersMask) const {
	const auto len = (numDaysResult() - numPreconstructedRows()) * m_numSolutionTotalB;
	*ppCompatibleSolutions = new tchar[len];
	*pUsePlayersMask = m_pAllData != NULL;
	return m_numSolutionTotalB;
}

CC uint& CRowStorage::solutionInterval2(uint* pRowSolutionIdx, int iRow, uint* pLast, ll availablePlayers) const {
	*pLast = pRowSolutionIdx[iRow + 1] = m_pPlayerSolutionCntr[iRow];
	if (iRow == numPreconstructedRows())
		*pLast = lastInFirstSet();

	return pRowSolutionIdx[iRow];
}

CC uint& CRowStorage::solutionInterval3(uint* pRowSolutionIdx, int iRow, uint* pLast, ll availablePlayers) const {
	pRowSolutionIdx[iRow + 1] = 0;
	if (pRowSolutionIdx[iRow]) {
		*pLast = pRowSolutionIdx[iRow + m_lenDayResults];
		return pRowSolutionIdx[iRow];
	}

#if USE_64_BIT_MASK
	unsigned long iBit;
	_BitScanForward64(&iBit, availablePlayers);
#else
	const auto iBit = this->firstOnePosition(availablePlayers);
#endif

	*pLast = pRowSolutionIdx[iRow + m_lenDayResults] = m_pPlayerSolutionCntr[iBit - 1];
	return pRowSolutionIdx[iRow] = m_pPlayerSolutionCntr[iBit - 2];
}
