#include "Table.h"

#define LATEST_IMPROVEMENT_FOR_TRIPLES   1

CC CRowStorage::CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects, const alldata* pAllData) :
	m_pSysParam(pSysParam), m_numPlayers(numPlayers),
	m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]),
	m_numDaysResult(pAllData->numDaysResult()), m_pAllData(pAllData),
	m_bGroupSize2(pAllData->groupSize() == 2),
	m_bSelectPlayerByMask(pAllData->groupSize() > 2 || pAllData->param(t_CBMP_Graph) > 1),
	m_bUseCombinedSolutions(pSysParam->val[t_useCombinedSolutions]),
	m_step(pSysParam->val[t_MultiThreading] == 2 ? pSysParam->val[t_numThreads] : 1),
	CStorage<tchar>(numObjects, 3 * numPlayers) {
	m_numObjectsMax = numObjects;
	m_pPlayerSolutionCntr = new uint[numPlayers + m_numDaysResult];
	m_pNumLongs2Skip = m_pPlayerSolutionCntr + numPlayers;
	initMaskStorage(numObjects);
	m_lenMask = m_pMaskStorage->lenObject();
	const auto useCliquesAfterRow = pSysParam->val[t_useSolutionCliquesAfterRow];
	m_useCliquesAfterRow = useCliquesAfterRow ? useCliquesAfterRow : m_numDaysResult;
	memset(m_pRowsCompatMasks, 0, sizeof(m_pRowsCompatMasks));
	m_fRowToBitmask = groupSize2() ? &CRowStorage::rowToBitmask2 : &CRowStorage::rowToBitmask3;
	m_fSolutionInterval = !selectPlayerByMask() ? &CRowStorage::solutionInterval2 : &CRowStorage::solutionInterval3;
	m_lenDayResults = m_numDaysResult + 1;
	m_pSolMemory = new tchar[2 * numPlayers];
	setLenCompare(m_numPlayers);
}

CC CRowStorage::~CRowStorage() {
	delete[] m_pPlayerSolutionCntr;
	for (int i = countof(m_pRowsCompatMasks); i--;)
		delete[] m_pRowsCompatMasks[i];

	delete m_pMaskStorage;
	delete[] m_pSolMemory;
	releaseSolMaskInfo();
	delete m_pTRTSN_Storage;
}

CC void CRowStorage::initMaskStorage(uint numObjects) {
	m_pMaskStorage = new CStorage<tchar>(numObjects, (((m_numPlayers * (m_numPlayers - 1) / 2) + 63) / 64) * 8);
	memset(m_pPlayerSolutionCntr, 0, m_numPlayers * sizeof(m_pPlayerSolutionCntr[0]));
	reset();
}

CC void CRowStorage::initPlayerMask(ctchar* pFirstMatr, ctchar lastNeighborOfPlayer0) {
	m_pFirstMatr = pFirstMatr;
	m_stepCombSolution = 0;
	ll playersMask = -1;
	int shift = numPlayers();
	const auto groupSize = m_pAllData->groupSize();
	if (!groupSize2()) {
		// Create a mask to manage players utilized in the predefined rows of the matrix.
		// Excluding players of the first group from ...
		playersMask <<= groupSize;						// ...first row
		auto const* pSolution = m_pAllData->result();
		for (int j = numPreconstructedRows(); --j;) {   // ... remaining pre-constructed rows
			pSolution += numPlayers();
			for (auto i = groupSize; --i;)
				playersMask ^= (ll)1 << pSolution[i];
		}
	}
	else {
		if (selectPlayerByMask()) {
			m_stepCombSolution = 2;
			playersMask <<= 2 * numPreconstructedRows() + 1;
		}
		else {
			// Predefined rows for groupSize = 2, are assumed to be correct.
			m_stepCombSolution = 1;
			playersMask <<= (numPreconstructedRows() + 1);
			shift = numDaysResult() + 1;
		}
	}

	if (m_pAllData->param(t_CBMP_Graph) > 1) {
		// Remove players flagged as ineligible from all generated solutions
		for (int i = groupSize; i < numPlayers(); i += groupSize)
			playersMask &= ((ll)-1 ^ ((ll)1 << i));
	}

	for (int i = lastNeighborOfPlayer0 + 1; i < numPlayers(); i++)
		playersMask &= ((ll)-1 ^ ((ll)1 << i));

	m_playersMask[0] = m_playersMask[1] = playersMask ^ ((ll)(-1) << shift);
}

CC bool CRowStorage::p1fCheck2P1F(ctchar* neighborsi, ctchar* neighborsj) const {
	tchar k = 0;
	tchar i = 2;
	tchar np = (tchar)m_numPlayers;
	if (np >= 16)
	{
		if ((k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0 ||
		    (k = neighborsj[neighborsi[k]]) == 0)
			return false;
		if (np == 16)
			return true;
		i = np - 14;
	}
	for (; i < np; i += 2)
	{
		if ((k = neighborsj[neighborsi[k]]) == 0)
			return false;
	}
	return true;
}

#define USE_EXIT	0
#if USE_EXIT
#define EXIT	exit(1)
#else
#define EXIT	return false
#endif

#define USE_PRINT_RED	1
#if USE_PRINT_RED
#define PRINT_RED(format, ...) printfRed(format, __VA_ARGS__)
#else
#define PRINT_RED(format, ...)
#endif

CC bool CRowStorage::addRow(ctchar* pRow, ctchar* pNeighbors, ctchar* pNeighbors2) {
	if (m_numObjects == m_numObjectsMax) {
		reallocStorageMemory(m_numObjectsMax <<= 1);
		m_pMaskStorage->reallocStorageMemory(m_numObjectsMax);
	}
#if 0 && !USE_CUDA
	FOPEN_F(f, "aaa.txt", m_numObjects ? "a" : "w");
	char buf[32];
	if (!m_numObjects) {
		// Output of the last precomputed row
		sprintf_s(buf, "%3d: ", -1);
		outMatrix(m_pAllData->result() + (m_numPreconstructedRows-1) * m_numPlayers, 1, m_numPlayers, m_pAllData->groupSize(), 0, f, false, false, buf);
	}
	sprintf_s(buf, "%3d: ", m_numObjects);
	outMatrix(pRow, 1, m_numPlayers, m_pAllData->groupSize(), 0, f, false, false, buf);
	FCLOSE_F(f);
#endif
	(this->*m_fRowToBitmask)(pRow, (tmask*)(m_pMaskStorage->getObject(m_numObjects)));
	auto* pntr = getObject(m_numObjects++);
#if !USE_CUDA
	const auto ptr = pntr - lenObject() + 1;

	auto i = m_pAllData->groupSize();
	while (--i) {
		if (!(getPlayersMask() & ((ll)1 << pRow[i]))) {
			PRINT_RED("\nSolution rejected : Player %d was already assigned to one of predefined rows, violating constraints.\n", pRow[i]);
			EXIT;
		}
	}

	if (!groupSize2() && m_playersMask[1] && pRow[1] > (i = minPlayer(m_playersMask[1]))) {
		//PRINT_RED(".");
		//printTable("Solution", pRow, 1, m_numPlayers, m_pAllData->groupSize(), 0, true);
		//PRINT_RED("\nSolution rejected : The solution involving player #%d, should precede the solution for player #%d\n", i, pRow[1]);
		//EXIT;
		return false; // this can happen if a row with requested player was skipped in main code because of some other reason
	}

	if (m_numObjects > 1) {
		if (pRow[1] != *ptr && pRow[1] != *ptr + m_stepCombSolution) {
			char buf[256], *pBuf = buf;
			SPRINTFD(pBuf, buf, "Error in code: pRow[1](%d) != %d, and pRow[1] != %d", pRow[1], *ptr, *ptr + m_stepCombSolution);
			if (selectPlayerByMask()) {
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

	if (selectPlayerByMask()) {
		// Mark referenced players:
		for (int i = m_pAllData->groupSize(); --i;)
			m_playersMask[1] &= ((ll)-1 ^ ((ll)1 << pRow[i]));
	}
#endif

	memcpy(pntr, pRow, m_numPlayers);
	memcpy(pntr + m_numPlayers, pNeighbors, m_numPlayers);
	memcpy(pntr + m_numPlayers * 2, pNeighbors2, m_numPlayers);
	m_pPlayerSolutionCntr[pRow[1] - 1]++;  // Increasing the number of solutions for player pRow[1] 
	                                       // NOTE: for groupSize = 2, this number is equal to the number of solutions for (pRow[1]-1)-th row
	return true;
}

CC bool CRowStorage::checkCompatibility(ctchar* neighborsi, const ll* rm, uint idx) const {
	// Let's check if the masks are mutually compatible
	auto* pMask = (const ll*)(m_pMaskStorage->getObject(idx));
	int j = m_lenMask >> 3;
	while (j-- && !(rm[j] & pMask[j]));
	// (j-- && !(rm[j] & pMask[j]));

	if (j >= 0)
		return false;
	ASSERT_IF(idx >= m_numSolutionTotal);
	const auto pObj = getObject(idx);
	ASSERT_IF(!pObj);
	TrCycles* tc = &m_pAllData->m_TrCyclesFirst2Rows[0];

	if (!m_pAllData->m_allRowPairsSameCycles)
		return true;

	if (groupSize2() && tc->length[0] == m_numPlayers)
		return p1fCheck2P1F(neighborsi + m_numPlayers, pObj + m_numPlayers * 2);

	TrCycles tcs;
	tcs.irow1 = tcs.irow2 = 0;
	const auto iret = m_pAllData->u1fGetCycleLength(&tcs, neighborsi, pObj + m_numPlayers,
		neighborsi - m_numPlayers, pObj, eCheckErrors);
	if (iret <= 0)
		return false;
	for (int itr0 = 0; itr0 < MAX_CYCLE_SETS; itr0++, tc++)
	{
		if (tc->counter == 0)
			break;
		if (!MEMCMP(tc->length, tcs.length, MAX_CYCLES_PER_SET))
			return true;
	}
	return false;
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

CC bool CRowStorage::generateCompatibilityMasks(tmask* pMaskOut, uint solIdx, uint idx, ll * pUsedPlayers) const {
	auto* rm = (const ll*)m_pMaskStorage->getObject(solIdx);
	bool compSolFound = false;
	const auto pSolution = getObject(solIdx);
	const auto pNeighbors = pSolution + m_numPlayers;
	for (;  idx < m_numSolutionTotal; idx++) {
		if (checkCompatibility(pNeighbors, rm, idx)) {
			const auto newIdx = idx - m_numRecAdj;
			SET_MASK_BIT(pMaskOut, newIdx);     // 1 - means OK
			compSolFound = true;
		}
	}

	if (compSolFound && selectPlayerByMask()) {
		// If groupSize > 2, we also need to create mask which will keep  
		// the information regarding the players used by current solution.
		// We will store it as 0's of corresponding bites.
		auto* pMaskOutLong = (ll*)pMaskOut + m_lenSolutionMask - 1;
		*pMaskOutLong = getPlayersMask();
		for (auto i = m_pAllData->groupSize(); --i;)
			*pMaskOutLong ^= (ll)1 << pSolution[i];

		if (pUsedPlayers)
			*pUsedPlayers &= *pMaskOutLong;
	}

	return compSolFound;
}

CC void CRowStorage::updateMasksByAutForSolution(ctchar *pSolution, const CGroupInfo* pGroupInfo, tmask *pMask, uint solIdx, uint last, uint idxMin) const {
	// For all non-trivial automorphisms:
	int i = 0;
	while (++i < pGroupInfo->orderOfGroup()) {
		const auto idx = getTransformerSolIndex(pSolution, pGroupInfo->getObject(i), m_numSolutionTotal, idxMin);
		if (idx == solIdx || idx == UINT_MAX)
			continue;

		if (idx < last) {
			// The solution #idx is not canonical, because there 
			// is an automorphism converting it to solution # < solIdx.
			if (idx < solIdx)
				break;

			// Should never happen if all first non-fixed row solutions are canonical. 
	//		ASSERT_IF(true);
			// Just in case, we'll reset the mask.
			memset(m_pRowsCompatMasks[1] + m_lenSolutionMask * idx, 0, m_numSolutionTotalB);
		}
		else {
			// All solutions: solIdx < solutionIdx <= solIdxLast cannot use 
			// the solution #idx, otherwise, the matrix will be non-canonical.
			// Mark the solution #idx as noncompatible with these solutions 
			auto pMaskOut = pMask + (idx >> SHIFT);
			const auto exMask = ~MASK_BIT(idx);
			for (auto j = solIdx; ++j < last;)
				*(pMaskOut += m_lenSolutionMask) &= exMask;
		}
	}
}

CC void CRowStorage::updateMasksByAut(const CGroupInfo* pGroupInfo) const {
	uint last = 0;
	if (selectPlayerByMask()) {
		auto availablePlayers = getPlayersMask();
		getSolutionRange(last, availablePlayers, numPreconstructedRows());
	}
	else
		last = m_lastInFirstSet;

	const auto solIdxLast = last - 1;
	auto pMask = m_pRowsCompatMasks[1];
	for (uint solIdx = 0; solIdx < solIdxLast; solIdx++) {
		const auto pSolution = getObject(solIdx);
		updateMasksByAutForSolution(pSolution, pGroupInfo, pMask, solIdx, last);
		pMask += m_lenSolutionMask;
	}
}

#define USE_PREVIOUS_SOLUTIONS	0
#if USE_PREVIOUS_SOLUTIONS
int my_counter = 0;
#endif

#define TRACE_MASKS					0
#if !USE_CUDA && TRACE_MASKS
int ggg = 0;
#define REPORT_REGECTION(reason)	fprintf(f, "  Rejected %d\n", reason); \
									FCLOSE_F(f)
#else
#define REPORT_REGECTION(reason)
#endif


CC bool CRowStorage::initCompatibilityMasks(CStorageSet<tchar>** ppSolRecast) {
#if !USE_CUDA && TRACE_MASKS
	FOPEN_F(f, "aaa.txt", ggg++ ? "a" : "w");
	fprintf(f, "Matrix #%4d", ggg);
#endif
	const auto* param = sysParam()->val;
	const auto pGroupInfo = m_pAllData->groupInfo(param[t_useAutForPrecRows]);
	for (int i = 1; i < m_numPlayers; i++)
		m_pPlayerSolutionCntr[i] += m_pPlayerSolutionCntr[i - 1];

	// Define the number of first long long's we don't need to copy to the next row.
	memset(m_pNumLongs2Skip, 0, m_numDaysResult * sizeof(m_pNumLongs2Skip[0]));
	int i = numPreconstructedRows();
	m_pNumLongs2Skip[i] = (m_lastInFirstSet = m_pPlayerSolutionCntr[i]) >> 6;

	// Trivial groups will not be used.
	m_bUseAut = pGroupInfo && pGroupInfo->numObjects() > 1;
	m_numRecAdj = !m_bUseAut ? m_lastInFirstSet : 0;
	while (++i < m_numDaysResult)
		m_pNumLongs2Skip[i] = ((m_pPlayerSolutionCntr[i] - m_numRecAdj) >> 6);

	m_numSolutionTotal = m_pPlayerSolutionCntr[m_numPlayers - 1];
	if (m_numSolutionTotal != m_numObjects)
		return false;

#if LATEST_IMPROVEMENT_FOR_TRIPLES
	int numRowToConstruct = m_numDaysResult - numPreconstructedRows();
	if (numRowToConstruct > (int)m_numSolutionTotal) {
		REPORT_REGECTION(1);
		return false;
	}
#endif

#if 0   // Output of table with total numbers of solutions for different input matrices 
	static int ccc = 0;
	FOPEN_F(f, "aaa.txt", ccc ? "a" : "w");
	if (!ccc++)
		fprintf(f, "Matrix#:     Num.Solutions:\n");
	fprintf(f, " %5d:          %6d\n", ccc, m_numSolutionTotal);
	FCLOSE_F(f);
#endif

	const auto useCombinedSolutions = param[t_useCombinedSolutions];
	if (useCombinedSolutions) {
		m_lastInFirstSet *= (m_numRec[1] = ((m_numRecAdj2 = m_pPlayerSolutionCntr[numPreconstructedRows() + 1]) - m_numRecAdj));
		delete [] m_pRowsCompatMasks[0];
	}

	const auto flg = !selectPlayerByMask();
	// Adding additional long long when we use groupSize > 2
	m_numSolutionTotalB = ((m_numSolutionTotal - m_numRecAdj + 7) / 8 + 7) / 8 * 8 + (flg ? 0 : 8);
	m_lenSolutionMask = m_numSolutionTotalB / sizeof(tmask);

	m_solAdj = useCombinedSolutions || m_bUseAut ? m_numRecAdj : 0;
	const auto len = m_numSolutionTotal - (m_numRecAdj - m_solAdj);
	delete[] m_pRowsCompatMasks[1];
	m_pRowsCompatMasks[1] = new tmask[len * m_lenSolutionMask];
	memset(m_pRowsCompatMasks[1], 0, len * m_numSolutionTotalB);

	releaseSolMaskInfo();
	const auto lenMasks = flg ? numDaysResult() : m_numPlayers;
	m_pRowSolutionMasksIdx = new uint[lenMasks];
	memset(m_pRowSolutionMasksIdx, 0, lenMasks * sizeof(m_pRowSolutionMasksIdx[0]));

	m_pRowSolutionMasks = new tmask[lenMasks];
	memset(m_pRowSolutionMasks, 0, lenMasks * sizeof(m_pRowSolutionMasks[0]));
	m_pRowSolutionMasksIdx[0] = 0;

#if SAME_MASK_IDX
	m_pMaskTestingCompleted = new bool[lenMasks];
	memset(m_pMaskTestingCompleted, 0, lenMasks * sizeof(m_pMaskTestingCompleted[0]));
#endif


#if !USE_64_BIT_MASK
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

	tmask* pRowsCompatMasks[] = { m_pRowsCompatMasks[0], m_pRowsCompatMasks[1] };
	tmask* pCompatMask = pRowsCompatMasks[1];

	bool skipAllowed = flg && !(useCombinedSolutions || m_bUseAut);
	unsigned int first, last = 0;
	i = numPreconstructedRows() - 1;

#if LATEST_IMPROVEMENT_FOR_TRIPLES
	int numRemainingSolution = m_numSolutionTotal;
#endif
	ll playerMask;
	auto availablePlayers = playerMask = getPlayersMask();
	ll* pUsedPlayers = flg ? NULL : &playerMask;
	while (availablePlayers) {
		first = getSolutionRange(last, availablePlayers, ++i);
		if (first >= last) {
			if (flg)
				break;

			continue;
		}

		if (m_pRowSolutionMasksIdx) {
			const auto lastAdj = last - m_numRecAdj;
			m_pRowSolutionMasksIdx[i] = lastAdj >> SHIFT;
			if (i == numPreconstructedRows() || m_pMaskTestingCompleted || m_pRowSolutionMasksIdx[i] > m_pRowSolutionMasksIdx[i - 1]) {
				unsigned int rem;
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

#if LATEST_IMPROVEMENT_FOR_TRIPLES
		numRemainingSolution -= last - first;
		if (--numRowToConstruct > numRemainingSolution) {
			REPORT_REGECTION(2);
			return false;
		}
#endif
		if (skipAllowed) {
			// Skip construction of masks for the first set of solutions.
			// The threads will do this latter.
			skipAllowed = false;
			continue;
		}

		if (last == m_numSolutionTotal) {
			if (pUsedPlayers) {
				// We need to mark players who will be potentially used with the last portion of solutions
				for (; first < last; first++) {
					auto* pSolution = getObject(first);
					for (auto i = m_pAllData->groupSize(); --i;)
						*pUsedPlayers &= (ll)(-1) ^ ((ll)1 << pSolution[i]);
				}
			}
			break;
		}

		if (!flg)
			pCompatMask = getSolutionMask(first);

		bool flag = false;
		for (; first < last; first++) {
			flag |= generateCompatibilityMasks(pCompatMask, first, last, pUsedPlayers);
			pCompatMask += m_lenSolutionMask;
		}

		if (LATEST_IMPROVEMENT_FOR_TRIPLES && !flag && pUsedPlayers) {
			// Any solutions from [first, last] is not compatible with any remaining solutions
			// Let's check if player #1 of any solution from that intarval was used in one of previous solutions
			const auto pSolution = getObject(last - 1);
			if (*pUsedPlayers & ((ll)1 << pSolution[1])) {
				REPORT_REGECTION(3);
				return false;
			}
		}
	}

	if (m_numRecAdj) {
		for (int i = numPreconstructedRows(); i < m_numPlayers; i++)
			m_pPlayerSolutionCntr[i] -= m_numRecAdj;
	}


	if (LATEST_IMPROVEMENT_FOR_TRIPLES && pUsedPlayers && *pUsedPlayers) {
		REPORT_REGECTION(4);
		return false;    // At least one player was not present with player 0 in any matrix row solution. 
	}

	ASSERT_IF(last != m_numSolutionTotal);
	if (USE_GROUP_4_2_ROWS) {
		const auto pGroupInfo = m_pAllData->groupInfo(2);
		const auto groupOrder = pGroupInfo->orderOfGroup();
		if (groupOrder > 1) {
			m_pTRTSN_Storage = new CRepository<uint>(2 * sizeof(uint), 32);
			const auto pRow = m_pAllData->result(numPreconstructedRows() - 1);

			// Populate the database with all solutions that the current 
			// solution of the 3rd row transforms into.
			// For all non-trivial automorphisms:
			uint solInfo[2];
			// This will be used for triples
			//auto availablePlayers = getPlayersMask();
			//first = getSolutionRange(m_lastInFirstSet, availablePlayers, 2);
			for (auto i = 0; ++i < groupOrder;) {
				// Transform the 3rd row solution and search results among all solutions.
				solInfo[0] = getTransformerSolIndex(pRow, pGroupInfo->getObject(i), m_numSolutionTotal, m_lastInFirstSet);

				if (solInfo[0] != UINT_MAX) {
					solInfo[1] = i;
					m_pTRTSN_Storage->updateRepo(solInfo);
				}
			}

#if USE_PREVIOUS_SOLUTIONS
			const auto shift = (numPreconstructedRows() - 1) * numPlayers();
			const auto pCurrentMatr = m_pAllData->result() + shift;

			// The transformed solutions of the 3rd row which are smaller than 
			// the current one cannot be compatible with any solution of the 4th row
			const auto lenRecord = shift + numPlayers();
			auto const* pSolution = m_pFirstMatr - numPlayers();
			my_counter = 0;
			while (memcmp(pSolution += lenRecord, pCurrentMatr, numPlayers()) < 0)
				updateMasksByAutForSolution(pSolution, pGroupInfo, m_pRowsCompatMasks[1], 0, m_lastInFirstSet, m_lastInFirstSet);
#endif
		}
	}

	if (m_bUseAut)
		updateMasksByAut(pGroupInfo);

	if (useCombinedSolutions) {
		delete m_pMaskStorage;
		m_pMaskStorage = NULL;
	}

#if !USE_CUDA && TRACE_MASKS
	fprintf(f, "  EXIT OK\n");
	FCLOSE_F(f);

	FOPEN_F(f1, "bbb.txt", "a");
	fprintf(f1, "Matrix #%4d is OK\n", ggg);
	FCLOSE_F(f1);
#endif

	if (ppSolRecast) {
#define COUNT_MASK_WEIGHT	1
#if COUNT_MASK_WEIGHT
		const auto prevWeight = countMaskFunc();
#endif
		modifyMask(ppSolRecast);
#if COUNT_MASK_WEIGHT
		countMaskFunc(prevWeight);
#endif
	}

	return true;
}

size_t CRowStorage::countMaskFunc(size_t prevWeight) const {
	static size_t matrWeight = 0;  // Max number of ones located above the upper main diagonal
	static uint numFirstSol = 0;   // Number of solutions for the first non-preconstruted row 
	static tchar table[256];
	const auto numMatrRow = m_numSolutionTotal - (m_numRecAdj - m_solAdj);
	if (!matrWeight) {
		// Preparing the table to calculate the weight of the mask matrix.
		table[0] = 0;
		for (int i = 1; i < sizeof(table); i++)
			table[i] = table[i >> 1] + i % 2;

		matrWeight = numMatrRow * (m_numSolutionTotal - m_numRecAdj);
		auto availablePlayers = getPlayersMask();
		unsigned int first, last = 0;
		int i = numPreconstructedRows() - 1;		
		while (availablePlayers) {
			first = getSolutionRange(last, availablePlayers, ++i);
			const auto portion = last - first;
			matrWeight -= portion * portion;
			if (!numFirstSol)
				numFirstSol = portion;
		}
	}

	ctchar* pMask = (ctchar * )getSolutionMask(0);
	size_t changingWeight, totalWeight = 0;
	const auto ind = m_numSolutionTotalB * numFirstSol;
	size_t i = m_numSolutionTotalB * numMatrRow;
	while (i--) {
		totalWeight += table[pMask[i]];
		if (i == ind)
			changingWeight = totalWeight;
	}

	changingWeight = totalWeight - changingWeight;

	printf("\nMask matrix total weight: %zd  %5.2f%%\n", totalWeight, double(totalWeight) / matrWeight * 100);
	printf("\nMask matrix changing weight: %zd  %5.2f%%\n", changingWeight, double(changingWeight) / matrWeight * 100);
	if (prevWeight)
		printf("Weight change is %zd:  %5.2f%%\n", prevWeight - changingWeight, double(prevWeight - changingWeight)/ prevWeight * 100);

	return changingWeight;
}

int CRowStorage::findIndexInRange(int left, int right, ctchar* pSol) const {
	while (left <= right) {
		const int mid = left + (right - left) / 2;
		const auto cmp = memcmp(getObject(mid), pSol, numPlayers());
		if (!cmp)
			return mid;

		if (cmp < 0)
			left = mid + 1;
		else
			right = mid - 1;
	}

	return -1;
}

CC void CRowStorage::modifyMask(CStorageSet<tchar>** ppSolRecast) {
#define OUT_RECASTED_SOLUTIONS  1
	// Update the masks to align with the recast solutions database.
	char buf[64], * pBuf = buf;
	SPRINTFD(pBuf, buf, "\nRow4 and corresponding Row #");
	const auto len = sizeof(buf) - (pBuf - buf);
	const auto groupSize = m_pAllData->groupSize();
	int f1, f2, f3, f4;	f1 = f2 = f3 = f4 = 0;
	int cntr = 0;
	int i = numPreconstructedRows() + 1;
	auto right = m_pPlayerSolutionCntr[numPreconstructedRows()];
	while (++i < numPlayers()) {
		int cmp = -1;
		uint last = 0;
		auto availablePlayers = getPlayersMask();
		auto first = getSolutionRange(last, availablePlayers, numPreconstructedRows());
		const auto left = right;
		right = m_pPlayerSolutionCntr[i-1];
		auto* pSolRecast = ppSolRecast[i];
		printf("\ni = %2d: numObj = %3d", i, pSolRecast->numObjects());
		for (int j = 0; j < pSolRecast->numObjects(); j++) {
			const auto sol = pSolRecast->getObject(j);
#if OUT_RECASTED_SOLUTIONS
			sprintf_s(pBuf, len, "%d", ++cntr);
			printTable(buf, sol, 2, numPlayers(), groupSize);
#endif
			int a = 0;
			while (cmp < 0 && first < last) {
				cmp = memcmp(getObject(first++), sol, numPlayers());
				a++;
			}

			//ASSERT_IF(first >= last/* || cmp > 0 */);
			// 
			// Let's try to find index of the second solution
			const auto idx = findIndexInRange(left, right, sol + numPlayers());
			if (!cmp) {
				f2++;
				// First vector from current pair of vectors was found,
				if (idx < 0)
					f4++;

				auto k = first;
				auto *pSolMask = getSolutionMask(first) + ((idx) >> SHIFT);
				const auto maskBit = (tmask)(-1) ^ MASK_BIT(idx);
				while (++k < last) {
					*(pSolMask += m_lenSolutionMask) &= maskBit;
				}
			}
			else {
				f1++;
				// First vector from current pair of vectors was NOT found
				if (idx < 0)
					f3++;
				first--; // for comparing the previous solution
			}

			cmp = -1;
		}
		pSolRecast->releaseAllObjects();
	}
}

CC uint CRowStorage::getSolutionRange(uint& last, ll &availablePlayers, int i) const {
	const uint first = last;
	if (selectPlayerByMask()) {
		const auto iBit = minPlayer(availablePlayers);
		last = m_pPlayerSolutionCntr[iBit - 1];
		availablePlayers ^= (ll)1 << iBit;
	}
	else {
		last = m_pPlayerSolutionCntr[i];
		availablePlayers &= (availablePlayers - 1);
	}

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

#if SAME_MASK_IDX
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

CC int CRowStorage::initRowUsage(tmask** ppCompatibleSolutions, bool *pUsePlayersMask) const {
	const auto lenMask = m_numSolutionTotalB >> (SHIFT - 3);
	const auto len = (numDaysResult() - numPreconstructedRows()) * lenMask;
	*ppCompatibleSolutions = new tmask[len];
	*pUsePlayersMask = selectPlayerByMask();
	return lenMask;
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
		generateCompatibilityMasks(pCompatibleSolutions, first, last);
	}
	else {
		memcpy(pCompatibleSolutions, m_pRowsCompatMasks[1] + first * m_lenSolutionMask, m_numSolutionTotalB);
	}

	if (m_pTRTSN_Storage && first) {
		// Using the group of automorphisms on 2 rows of matrix.
		const auto pGroupInfo = m_pAllData->groupInfo(2);
		// For all previous solution of 4-th row which were NOT yet used for 2-row Aut elimination 
		for (uint solIdx = first - m_step; solIdx < first; solIdx++) {
			auto* pSol = getObject(solIdx);
			uint minIdxTr;
			for (int i = 0; i < m_pTRTSN_Storage->numObjects(); i++) {
				const auto *p = (const uint*)m_pTRTSN_Storage->getObject(i);
				uint solIdxTr = getTransformerSolIndex(pSol, pGroupInfo->getObject(p[1]), m_numSolutionTotal, m_lastInFirstSet);
				if (solIdxTr == UINT_MAX)
					continue;

				if (solIdxTr < (minIdxTr = p[0])) {
					minIdxTr = solIdxTr;
					solIdxTr = p[0];
				}

				auto pMaskOut = m_pRowsCompatMasks[1] + m_lenSolutionMask * minIdxTr + (solIdxTr >> SHIFT);
				*pMaskOut &= ~MASK_BIT(solIdxTr);
			}
		}
	}
}

CC uint CRowStorage::getTransformerSolIndex(ctchar* pSol, ctchar *pPerm, uint last, uint first) const {
	auto* pPermSolution = m_pSolMemory + numPlayers();
	PERMUTATION_OF_PLAYERS(numPlayers(), pSol, pPerm, m_pSolMemory);
	(m_pAllData->sortGroupsFn)(m_pSolMemory);
	m_pAllData->kmSortGroupsByFirstValue(m_pSolMemory, pPermSolution);
#if USE_PREVIOUS_SOLUTIONS && !USE_CUDA
	FOPEN_F(f, "bbb.txt", my_counter ? "a" : "w");
	char buf[32];
	sprintf_s(buf, "%3d: ", ++my_counter);
	outMatrix(pPermSolution, 1, m_numPlayers, m_pAllData->groupSize(), 0, f, false, false, buf);
	FCLOSE_F(f);
#endif
	return findObject(pPermSolution, first, last);
}

#if !USE_CUDA
CC void CRowStorage::outSelectedSolution(int iRow, uint first, uint last, int threadID) const {
	char buf[32];
	sprintf_s(buf, "aaa_%2d.txt", threadID);
	FOPEN_F(f, buf, "a");
	sprintf_s(buf, "row = %2d: (%2d - %2d) ", iRow, first, last);
	const auto& param = sysParam()->val;
	outMatrix(getObject(first), 1, param[t_numPlayers], param[t_groupSize], 0, f, false, false, buf);
	FCLOSE_F(f);
}
#endif
