#include "TripleSys.h"
#include "Table.h"

CC CRowStorage::CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects, const alldata* pAllData) :
	m_pSysParam(pSysParam), m_numPlayers(numPlayers),
	m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]),
	m_numDaysResult(pAllData->numDaysResult()), m_pAllData(pAllData),
	m_bGroupSize2(pAllData->groupSize() == 2),
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
	m_fRowToBitmask = m_bGroupSize2 ? &CRowStorage::rowToBitmask2 : &CRowStorage::rowToBitmask3;
	m_fSolutionInterval = m_bGroupSize2 ? &CRowStorage::solutionInterval2 : &CRowStorage::solutionInterval3;
	m_lenDayResults = m_numDaysResult + 1;
	m_pSolMemory = new tchar[2 * numPlayers];
}

CC CRowStorage::~CRowStorage() {
	delete[] m_pPlayerSolutionCntr;
	for (int i = countof(m_pRowsCompatMasks); i--;)
		delete[] m_pRowsCompatMasks[i];

	delete m_pMaskStorage;
	delete[] m_pSolMemory;
	releaseSolMaskInfo();
}

CC void CRowStorage::initMaskStorage(uint numObjects) {
	m_pMaskStorage = new CStorage<tchar>(numObjects, (((m_numPlayers * (m_numPlayers - 1) / 2) + 63) / 64) * 8);
	memset(m_pPlayerSolutionCntr, 0, m_numPlayers * sizeof(m_pPlayerSolutionCntr[0]));
	reset();
}

CC void CRowStorage::initPlayerMask(bool groupSize_2) {
	ll playersMask = -1;
	if (!groupSize_2) {
		// Create a mask to manage players utilized in the predefined rows of the matrix.
		const auto groupSize = m_pAllData->groupSize();
		// Excluding players of the first group from ...
		playersMask <<= groupSize;						// ...first row
		auto const* pSolution = m_pAllData->result();
		for (int j = numPreconstructedRows(); --j;) {   // ... remaining pre-constructed rows
			pSolution += m_numPlayers;
			for (auto i = groupSize; --i;)
				playersMask ^= (ll)1 << pSolution[i];
		}
	}
	else {
		// Predefined rows for groupSuze = 2, are assumed to be correct.
		playersMask <<= (numPreconstructedRows() + 1);
	}

	m_playersMask[0] = m_playersMask[1] = playersMask ^ ((ll)(-1) << m_numPlayers);
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

#define USE_EXIT	0
#if USE_EXIT
#define EXIT	exit(1)
#else
#define EXIT	return false
#endif

#define USE_PRINT_RED	0
#if USE_PRINT_RED
#define PRINT_RED(format, ...) printfRed(format, __VA_ARGS__)
#else
#define PRINT_RED(format, ...)
#endif

CC bool CRowStorage::addRow(ctchar* pRow, ctchar* pNeighbors) {
	if (m_numObjects == m_numObjectsMax) {
		reallocStorageMemory(m_numObjectsMax <<= 1);
		m_pMaskStorage->reallocStorageMemory(m_numObjectsMax);
	}
#if 0 && !USE_CUDA
	FOPEN_F(f, "aaa.txt", m_numObjects ? "a" : "w");
	char buf[32];
	sprintf_s(buf, "%3d: ", m_numObjects);
	outMatrix(pRow, 1, m_numPlayers, m_pAllData->groupSize(), 0, f, false, false, buf);
	FCLOSE_F(f);
#endif
	(this->*m_fRowToBitmask)(pRow, (tmask*)(m_pMaskStorage->getObject(m_numObjects)));
	auto* pntr = getObject(m_numObjects++);
#if !USE_CUDA
	const auto ptr = pntr - 2 * m_numPlayers + 1;

	auto i = m_pAllData->groupSize();
	while (--i) {
		if (!(getPlayersMask() & ((ll)1 << pRow[i]))) {
			PRINT_RED("\nSolution rejected : Player %d was already assigned to one of predefined rows, violating constraints.\n", pRow[i]);
			EXIT;
		}
	}

	if (!m_bGroupSize2 && m_playersMask[1] && pRow[1] > (i = minPlayer(m_playersMask[1]))) {
		PRINT_RED("\nSolution rejected : The solution involving player #%d, should precede the solution for player #%d\n", i, pRow[1]);
		EXIT;
	}

	if (m_numObjects > 1) {
		if (pRow[1] != *ptr && pRow[1] != *ptr + 1) {
			char buf[256], *pBuf = buf;
			SPRINTFD(pBuf, buf, "Error in code: pRow[1](%d) != %d, and pRow[1] != %d", pRow[1], *ptr, *ptr + 1);
			if (!m_bGroupSize2) {
				// Calculate the number of non-referenced players
				int counter = 0;
				for (auto i = (*ptr + 1); i < pRow[1]; i++) {
					if (m_playersMask[1] & ((ll)1 << i))
						counter++;
				}

				if (counter) {
					SPRINTFD(pBuf, buf, ".\nPlayer%s#", counter > 1 ? "s" : "");
					for (auto i = (*ptr + 1); i < pRow[1]; i++) {
						if (m_playersMask[1] & ((ll)1 << i))
							SPRINTFD(pBuf, buf, " %d", i);
					}

					SPRINTFD(pBuf, buf, " %s not referenced in any previous solutions", counter > 1 ? "were" : "was");
				}
				else
					pBuf = buf;
			}

			if (pBuf != buf) {
				PRINT_RED("\n%s\n", buf);
				EXIT;
			}
		}
	}

	if (!m_bGroupSize2) {
		// Mark referenced players:
		for (int i = m_pAllData->groupSize(); --i;)
			m_playersMask[1] &= ((ll)-1 ^ ((ll)1 << pRow[i]));
	}
#endif

	memcpy(pntr, pRow, m_numPlayers);
	memcpy(pntr + m_numPlayers, pNeighbors, m_numPlayers);
	m_pPlayerSolutionCntr[pRow[1] - 1]++;  // Increasing the number of solutions for player pRow[1] 
	                                       // NOTE: for groupSize = 2, this number is equal to the number of solutions for (pRow[1]-1)-th row
	return true;
}

CC bool CRowStorage::checkCompatibility(ctchar* neighborsi, const ll* rm, uint idx) const {
	// Let's check if the masks are mutually compatible
	auto* pMask = (const ll*)(m_pMaskStorage->getObject(idx));
	int j = m_lenMask >> 3;
	while (j-- && !(rm[j] & pMask[j]));

	if (j >= 0)
		return false;

	const auto pObj = getObject(idx);
	if (m_bGroupSize2)
		return p1fCheck2(neighborsi, pObj + m_numPlayers);

	return m_pAllData->p1fCheck3(neighborsi - m_numPlayers, pObj, neighborsi, pObj + m_numPlayers);
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
	bool compSolFound = false;
	const auto pNeighbors = getObject(solIdx) + m_numPlayers;
	do {
		if (checkCompatibility(pNeighbors, rm, idx)) {
			const auto newIdx = idx - m_numRecAdj;
			SET_MASK_BIT(pMaskOut, newIdx);     // 1 - means OK
			compSolFound = true;
		}
	} while (++idx < m_numSolutionTotal);

	if (!m_bGroupSize2 && compSolFound) {
		// If groupSize > 2, we also need to create mask which will keep  
		// the information regarding the players used by current solution.
		// We will store it as 0's of corresponding bites.
		ctchar* pSolution = pNeighbors - m_numPlayers;
		auto* pMaskOutLong = (ll*)pMaskOut + m_lenSolutionMask - 1;
		*pMaskOutLong = getPlayersMask();
		for (auto i = m_pAllData->groupSize(); --i;)
			*pMaskOutLong ^= (ll)1 << pSolution[i];
	}
}

CC void CRowStorage::updateMasksByAut(uint idxMax, const CGroupInfo* pGroupInfo) const {
	uint last = 0;
	if (!m_bGroupSize2) {
		auto availablePlayers = getPlayersMask();
		getSolutionRange(last, availablePlayers, m_numPreconstructedRows);
	}
	else
		last = m_lastInFirstSet;

	const auto solIdxLast = last - 1;
	auto pMask = m_pRowsCompatMasks[1];
	for (uint solIdx = 0; solIdx < solIdxLast; solIdx++) {
		const auto pSolution = getObject(solIdx);
		auto pPermSolution = m_pSolMemory + numPlayers();
		// For all non-trivial automorphisms:
		for (auto i = pGroupInfo->numObjects(); --i;) { 
			auto* pntr = pGroupInfo->getObject(i);
			for (auto j = numPlayers(); j--;)
				m_pSolMemory[j] = pntr[pSolution[j]];


			(m_pAllData->sortGroupsFn)(m_pSolMemory);
			m_pAllData->kmSortGroupsByFirstValue(m_pSolMemory, pPermSolution);
			const auto idx = findSolution(pPermSolution, 0, idxMax);
			if (idx < last) {
				if (idx == solIdx)
					continue;

				// The solution # idx is not canonical, because there 
				// is an automorphism converting it to solution # < solIdx.
				if (idx < solIdx)
					break;

				memset(m_pRowsCompatMasks[1] + m_lenSolutionMask * idx, 0, m_numSolutionTotalB);
			} else {
				if (idx == UINT_MAX)
					continue;

				// All solutions: solIdx < solutionIdx <= solIdxLast cannot use the solution #idx
				// otherwise, the matrix will be non-canonical.
				auto pMaskOut = pMask + (idx >> SHIFT);
				const auto exMask = ~MASK_BIT(idx);
				for (auto j = solIdx; ++j < last;)
					*(pMaskOut += m_lenSolutionMask) &= exMask;
			}
		}
		pMask += m_lenSolutionMask;
	}
}

CC uint CRowStorage::findSolution(ctchar* tr, uint low, uint high) const {
	// search for element 
	high--;
	while (low <= high) {
		const auto itr = low + ((high - low) >> 1);
		const auto cmp = MEMCMP(getObject(itr), tr, m_numPlayers);
		if (!cmp)
			return itr;

		if (cmp < 0)
			low = itr + 1;  // ignore left half
		else
			high = itr - 1; // ignore right half
	}

	return UINT_MAX;		// not found 
}

CC void CRowStorage::initCompatibilityMasks(const CGroupInfo *pGroupInfo) {
	for (int i = 1; i < m_numPlayers; i++)
		m_pPlayerSolutionCntr[i] += m_pPlayerSolutionCntr[i - 1];

	// Define the number of first long long's we don't need to copy to the next row.
	memset(m_pNumLongs2Skip, 0, m_numPlayers * sizeof(m_pNumLongs2Skip[0]));
	int i = m_numPreconstructedRows;
	m_pNumLongs2Skip[i] = m_pPlayerSolutionCntr[i] >> 6;
	m_lastInFirstSet = m_pPlayerSolutionCntr[i];
	// Trivial groups will not be used.
	m_bUseAut = pGroupInfo && pGroupInfo->numObjects() > 1;
	m_numRecAdj = !m_bUseAut ? m_lastInFirstSet : 0;
	while (++i < m_numPlayers)
		m_pNumLongs2Skip[i] = ((m_pPlayerSolutionCntr[i] - m_numRecAdj) >> 6);

	m_numSolutionTotal = m_pPlayerSolutionCntr[m_numPlayers - 1];



	const auto useCombinedSolutions = sysParam()->val[t_useCombinedSolutions];
	if (useCombinedSolutions) {
		m_lastInFirstSet *= (m_numRec[1] = ((m_numRecAdj2 = m_pPlayerSolutionCntr[m_numPreconstructedRows + 1]) - m_numRecAdj));
		delete [] m_pRowsCompatMasks[0];
	}

	// Adding additional long long when we use groupSize > 2
	m_numSolutionTotalB = ((m_numSolutionTotal - m_numRecAdj + 7) / 8 + 7) / 8 * 8 + (m_bGroupSize2 ? 0 : 8);
	m_lenSolutionMask = m_numSolutionTotalB / sizeof(tmask);

	m_solAdj = useCombinedSolutions || m_bUseAut ? m_numRecAdj : 0;
	const auto len = m_numSolutionTotal - (m_numRecAdj - m_solAdj);
	delete[] m_pRowsCompatMasks[1];
	m_pRowsCompatMasks[1] = new tmask[len * m_lenSolutionMask];
	memset(m_pRowsCompatMasks[1], 0, len * m_numSolutionTotalB);
	if (m_bGroupSize2) {
		releaseSolMaskInfo();
		const auto numDays = numDaysResult();
		m_pRowSolutionMasksIdx = new uint[numDays];
		memset(m_pRowSolutionMasksIdx, 0, numDays * sizeof(m_pRowSolutionMasksIdx[0]));

		m_pRowSolutionMasks = new tmask[numDays];
		memset(m_pRowSolutionMasks, 0, numDays * sizeof(m_pRowSolutionMasks[0]));
		m_pRowSolutionMasksIdx[0] = 0;

#if SAME_MASK_IDX
		m_pMaskTestingCompleted = new bool[numDays];
		memset(m_pMaskTestingCompleted, 0, numDays * sizeof(m_pMaskTestingCompleted[0]));
#endif
	}

#if !USE_64_BIT_MASK
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

	tmask* pRowsCompatMasks[] = { m_pRowsCompatMasks[0], m_pRowsCompatMasks[1] };
	tmask* pCompatMask = pRowsCompatMasks[1];

	const bool skipFirstAllowed = m_bGroupSize2 && !(useCombinedSolutions || m_bUseAut);
	unsigned int first, last = 0;
	i = m_numPreconstructedRows - 1;
#if 1
	auto availablePlayers = getPlayersMask();
	unsigned int rem;
	const auto iMax = m_numPlayers - 1;
	while (++i < iMax && availablePlayers) {
		first = getSolutionRange(last, availablePlayers, i);
		if (m_pRowSolutionMasksIdx) {
			const auto lastAdj = last - m_numRecAdj;
			m_pRowSolutionMasksIdx[i] = lastAdj >> SHIFT;
			if (i == m_numPreconstructedRows || m_pMaskTestingCompleted || m_pRowSolutionMasksIdx[i] > m_pRowSolutionMasksIdx[i - 1]) {
				if (rem = REM(lastAdj)) {
					m_pRowSolutionMasks[i] = (tmask)(-1) << rem;
					if (m_pRowSolutionMasksIdx[i - 1] == m_pRowSolutionMasksIdx[i]) {
						// At this point, the solutions for three consecutive rows are masked by the same `tmask` element
						// Clear the bits in the previous mask that do not correspond to the solutions of the previous row
						m_pRowSolutionMasks[i - 1] ^= m_pRowSolutionMasksIdx[i];
#if SAME_MASK_IDX
						m_pMaskTestingCompleted[i - 1] = 1;
#endif
					}
				}
			}
			else {
				// We cannot use code for UseSolutionMasks, because now our code is not ready for such cases
				// Actually, in that case m_pRowSolutionMasks for the start and the end of the interval are stored in the same element of the array m_pRowSolutionMasks. 
				releaseSolMaskInfo();
				m_pRowSolutionMasksIdx = NULL;
				m_pRowSolutionMasks = NULL;
				m_pMaskTestingCompleted = NULL;
			}

#if 0
			FOPEN_F(f, "aaa.txt", "w");
			char buf[256], * pBuf = buf;
			for (int j = 0; j <= i; j++)
				SPRINTFD(pBuf, buf, "%18d", m_pRowSolutionMasksIdx[j]);

			fprintf(f, "%s\n", buf);
			pBuf = buf;
			for (int j = 0; j <= i; j++)
				SPRINTFD(pBuf, buf, "  %016llx", m_pRowSolutionMasks[j]);

			fprintf(f, "%s\n", buf);
			FCLOSE_F(f);
#endif
		}

		if (skipFirstAllowed && !first) {
			// Skip construction of masks for the first set of solutions.
			// The threads will do this latter.
			continue;
		}

		if (!m_bGroupSize2)
			pCompatMask = getSolutionMask(first);

		for (; first < last; first++) {
			generateCompatibilityMasks(pCompatMask, first, last);
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

	if (m_bUseAut)
		updateMasksByAut(last, pGroupInfo);

	if (useCombinedSolutions) {
		delete m_pMaskStorage;
		m_pMaskStorage = NULL;
	}
}

CC uint CRowStorage::getSolutionRange(uint& last, ll &availablePlayers, int i) const {
	const uint first = last;
	if (!m_bGroupSize2) {
		unsigned long iBit = 0;
#if !USE_CUDA
		_BitScanForward64(&iBit, availablePlayers);
#else
		#pragma message("A GPU-equivalent function similar to `_BitScanForward64` needs to be implemented.")
#endif
		last = m_pPlayerSolutionCntr[iBit - 1];
			availablePlayers ^= (ll)1 << iBit;
	}
	else
		last = m_pPlayerSolutionCntr[i];

	return first;
}

CC bool CRowStorage::checkSolutionByMask(int iRow, const tmask* pToASol) const {
	auto jMax = *(rowSolutionMasksIdx() + iRow);
	const auto iMax = numPlayers() - 2;
	for (int i = iRow + 1; i <= iMax; i++) {
		auto j = jMax;
		jMax = *(rowSolutionMasksIdx() + i);

		// Check left, middle and right parts of the solution interval for i-th row
		auto mask = *(rowSolutionMasks() + i - 1);
		if (mask && (mask & pToASol[j++]))
			continue;  // at least one solution masked by left part of the interval is still valid

#if NEW
		if (m_pMaskTestingCompleted[i - 1])
			return false;
#endif

		// middle part
		while (j < jMax && !pToASol[j])
			j++;

		if (j < jMax)
			continue;   // at least one solution masked by middle part of the interval is still valid

		// There are no valid solutions with the indices inside 
		// the interval defined by set of long longs
		mask = *(rowSolutionMasks() + i);
		// If mask != 0, we need to check the right side of the intervals.
		if (!mask || !((~mask) & pToASol[jMax]))
			return false;
	}

	return true;
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
	*pUsePlayersMask = !m_bGroupSize2;
	return m_numSolutionTotalB;
}

CC uint& CRowStorage::solutionInterval2(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const {
	const auto iRow = *pLast;
	*pLast = pRowSolutionIdx[1] = m_pPlayerSolutionCntr[iRow];
	if (iRow == numPreconstructedRows())
		*pLast = lastInFirstSet();

	return *pRowSolutionIdx;
}

CC uint& CRowStorage::solutionInterval3(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const {
	pRowSolutionIdx[1] = 0;
	if (*pRowSolutionIdx) {
		*pLast = pRowSolutionIdx[m_lenDayResults];
		return *pRowSolutionIdx;
	}

	if (!availablePlayers) {
		*pLast = UINT_MAX;
		return m_pPlayerSolutionCntr[0];  // Dummy return, it will not be used.
	}

	const auto iBit = minPlayer(availablePlayers);
	*pLast = pRowSolutionIdx[m_lenDayResults] = m_pPlayerSolutionCntr[iBit - 1];
	return *pRowSolutionIdx = m_pPlayerSolutionCntr[iBit - 2];
}

CC void CRowStorage::passCompatibilityMask(tmask* pCompatibleSolutions, uint first, uint last) const {
	if (!m_bUseAut) {
		memset(pCompatibleSolutions, 0, m_numSolutionTotalB);
		generateCompatibilityMasks((tmask*)pCompatibleSolutions, first, last);
	}
	else {
		memcpy(pCompatibleSolutions, m_pRowsCompatMasks[1] + first * m_lenSolutionMask, m_numSolutionTotalB);
	}
}
