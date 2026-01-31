#include "Table.h"
//#include "CompatMasks.h"

CCompatMasks::CCompatMasks(const kSysParam* pSysParam, const alldata* pAllData) :
	m_numPlayers(pSysParam->paramVal(t_numPlayers)),
	m_numPreconstructedRows(pSysParam->paramVal(t_useRowsPrecalculation)),
	m_numDaysResult(pAllData->numDaysResult()),
	m_bSelectPlayerByMask(pAllData->groupSize() > 2 || pAllData->param(t_CBMP_Graph) > 1) {
	m_pPlayerSolutionCntr = new uint[numPlayers() + numDaysResult()];
	m_pNumLongs2Skip = m_pPlayerSolutionCntr + numPlayers();
	memset(numPlayerSolutionsPtr(), 0, numDaysResult() * sizeof(m_pPlayerSolutionCntr[0]));
	m_lenDayResults = numDaysResult() + 1;
	m_fSolutionInterval = !selectPlayerByMask() ? &CCompatMasks::solutionInterval2 : &CCompatMasks::solutionInterval3;
	initSolMaskIndices();
}

CCompatMasks::~CCompatMasks() {
	releaseCompatMaskMemory();
	releaseSolMaskInfo();
	delete[] numPlayerSolutionsPtr();
}

CC void CCompatMasks::releaseSolMaskInfo() {
	delete[] rowSolutionMasks();
	delete[] rowSolutionMasksIdx();
	m_pRowSolutionMasksIdx = NULL;
	m_pRowSolutionMasks = NULL;
#if SAME_MASK_IDX
	delete[] maskTestingCompleted();
	m_pMaskTestingCompleted = NULL;
#endif
}

void CCompatMasks::initSolMaskIndices() {
	const auto lenMaskIdx = selectPlayerByMask() ? numPlayers() : numDaysResult();
	releaseSolMaskInfo();
	m_pRowSolutionMasksIdx = new uint[lenMaskIdx];
	memset(m_pRowSolutionMasksIdx, 0, lenMaskIdx * sizeof(m_pRowSolutionMasksIdx[0]));

	m_pRowSolutionMasks = new tmask[lenMaskIdx];
	memset(m_pRowSolutionMasks, 0, lenMaskIdx * sizeof(m_pRowSolutionMasks[0]));
	m_pRowSolutionMasksIdx[0] = 0;
#if SAME_MASK_IDX
	m_pMaskTestingCompleted = new bool[lenMaskIdx];
	memset(maskTestingCompleted(), 0, lenMaskIdx * sizeof(m_pMaskTestingCompleted[0]));
#endif
}

void CCompatMasks::initMaskMemory(uint numSolutions, int lenUsedMask, int numRecAdj, int numSolAdj)
{
	// Adding additional long long when we use groupSize > 2
	setNumSolutions(numSolutions);
	m_numSolutionTotalB = ((numSolutions - numRecAdj + 7) / 8 + 7) / 8 * 8 + lenUsedMask;
	m_lenSolutionMask = numSolutionTotalB() / sizeof(tmask);

	m_numMasks = m_numSolutionTotal - (numRecAdj - (m_solAdj = numSolAdj));
	releaseCompatMaskMemory();
	m_pRowsCompatMasks = new tmask[numMasks() * lenSolutionMask()];
	memset(rowsCompatMasks(), 0, numMasks() * numSolutionTotalB());
}

CC void CCompatMasks::initRowUsage(tmask** ppCompatibleSolutions, bool* pUsePlayersMask) const {
	const auto lenMask = m_numSolutionTotalB >> (SHIFT - 3);
	if (!*ppCompatibleSolutions) {
		const auto len = (numDaysResult() - numPreconstructedRows()) * lenMask;
		*ppCompatibleSolutions = new tmask[len];
	}

	*pUsePlayersMask = selectPlayerByMask();
}

CC uint& CCompatMasks::solutionInterval2(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const {
	const auto iRow = *pLast;
	*pLast = pRowSolutionIdx[1] = m_pPlayerSolutionCntr[iRow];
	if (iRow == numPreconstructedRows())
		*pLast = lastInFirstSet();

	return *pRowSolutionIdx;
}

CC uint& CCompatMasks::solutionInterval3(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const {
	pRowSolutionIdx[1] = 0;
	if (pRowSolutionIdx[0]) {
		*pLast = pRowSolutionIdx[lenDayResults()];
		return *pRowSolutionIdx;
	}

	if (!availablePlayers) {
		*pLast = UINT_MAX;
		return m_pPlayerSolutionCntr[0];  // Dummy return, it will not be used.
	}

	const auto iBit = minPlayer(availablePlayers);
	*pLast = pRowSolutionIdx[lenDayResults()] = m_pPlayerSolutionCntr[iBit - 1];
	return *pRowSolutionIdx = m_pPlayerSolutionCntr[iBit - 2];
}

int CCompatMasks::setNumLongs2Skip(bool adjustRecCounter)
{
	// Define the number of first long long's we don't need to copy to the next row.
	int i = numPreconstructedRows();
	setLastInFirstSet(m_pPlayerSolutionCntr[i]);
	int recAdj = adjustRecCounter ? lastInFirstSet() : 0;

	m_pNumLongs2Skip[i] = lastInFirstSet() >> 6;
	while (++i < numDaysResult())
		m_pNumLongs2Skip[i] = ((m_pPlayerSolutionCntr[i] - recAdj) >> 6);

	return recAdj;
}

CC uint CCompatMasks::getSolutionRange(uint& last, ll& availablePlayers, int i) const {
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

void CCompatMasks::defineMask4SolutionIntervals(int nRow, unsigned int last)
{
	if (m_pRowSolutionMasksIdx) {
		const auto lastAdj = last - numRecAdj();
		m_pRowSolutionMasksIdx[nRow] = lastAdj >> SHIFT;
		if (nRow == numPreconstructedRows() || maskTestingCompleted() || m_pRowSolutionMasksIdx[nRow] > m_pRowSolutionMasksIdx[nRow - 1]) {
			unsigned int rem;
			if (rem = REM(lastAdj)) {
				m_pRowSolutionMasks[nRow] = (tmask)(-1) << rem;
				if (m_pRowSolutionMasksIdx[nRow - 1] == m_pRowSolutionMasksIdx[nRow]) {
					// At this point, the solutions for three consecutive rows are masked by the same `tmask` element
					// Clear the bits in the previous mask that do not correspond to the solutions of the previous row
					m_pRowSolutionMasks[nRow - 1] ^= m_pRowSolutionMasksIdx[nRow];
#if SAME_MASK_IDX
					maskTestingCompleted()[nRow - 1] = 1;
#endif
				}
			}
		}
		else {
			// We cannot use the code for UseSolutionMasks because our code is not currently ready to handle such cases.
			// Actually, in that case m_pRowSolutionMasks for the start and the end of the interval are stored in the same element of the array m_pRowSolutionMasks. 
			releaseSolMaskInfo();
		}
#if 0
		FOPEN_F(f, "aaa.txt", "w");
		char buf[256], * pBuf = buf;
		for (int j = 0; j <= nRow; j++)
			SPRINTFD(pBuf, buf, "%18d", m_pRowSolutionMasksIdx[j]);

		fprintf(f, "%s\n", buf);
		pBuf = buf;
		for (int j = 0; j <= nRow; j++)
			SPRINTFD(pBuf, buf, "  %016llx", m_pRowSolutionMasks[j]);

		fprintf(f, "%s\n", buf);
		FCLOSE_F(f);
#endif
	}
}
#if COUNT_MASK_WEIGHT || OUT_MASK_WEIGHT
size_t totalWeighChange = 0;
#endif

uint CCompatMasks::countMaskFunc(const tmask* pCompSolutions, size_t numMatrRow, uint prevWeight) const {
	static size_t matrWeight = 0;  // Max number of ones located above the upper main diagonal
	static uint numFirstSol = 0;   // Number of solutions for the first non-preconstruted row 
	static tchar table[256];
	if (!matrWeight) {
		// Preparing the table to calculate the weight of the mask matrix.
		table[0] = 0;
		for (int i = 1; i < sizeof(table); i++)
			table[i] = table[i >> 1] + i % 2;

		auto availablePlayers = getPlayersMask();
		unsigned int last = 0;
		int i = numPreconstructedRows() - 1;
		unsigned int first = getSolutionRange(last, availablePlayers, ++i);
		numFirstSol = last - first;
		matrWeight = numMatrRow * (m_numSolutionTotal - m_numRecAdj);
#if OUT_MASK_WEIGHT
		matrWeight -= numFirstSol * numFirstSol;
		while (availablePlayers) {
			first = getSolutionRange(last, availablePlayers, ++i);
			const auto portion = last - first;
			matrWeight -= portion * portion;
		}
#endif
	}

	ctchar* pMask = (ctchar*)pCompSolutions;
	uint changingWeight = 0, totalWeight = 0;
	const auto ind = numSolutionTotalB() * numFirstSol;
	size_t i = numSolutionTotalB() * numMatrRow;
	while (i--) {
		totalWeight += table[pMask[i]];
		if (i == ind)
			changingWeight = totalWeight;
	}

	changingWeight = totalWeight - changingWeight;

#if OUT_MASK_WEIGHT
	printf("\nMask matrix total weight: %d  %5.2f%%", totalWeight, double(totalWeight) / matrWeight * 100);
	printf("\nMask matrix changing weight: %d  %5.2f%%", changingWeight, double(changingWeight) / matrWeight * 100);
	if (prevWeight) {
		const auto weightChange = prevWeight - changingWeight;
		totalWeighChange += weightChange;
		printf("\nWeight change is %d:  %5.2f%%\n", weightChange, double(weightChange) / prevWeight * 100);

	}
#endif
	return changingWeight;
}

void CCompressedMask::compressCompatMasks(tmask* pCompSol, const CCompatMasks* pCompMask)
{
	// Count the number of valid solutions in the compatibility mask 
	// for the current solution of the first non-predefined row (including itself)
	const auto nValidSol = pCompMask->countMaskFunc(pCompSol) + 1; 

	initMaskMemory(nValidSol, (selectPlayerByMask() ? 8 : 0));
	releaseSolIndices();
	m_pSolIdx = new uint[nValidSol];

	auto const *pSolMasksIniIdx = pCompMask->numPlayerSolutionsPtr() + numPreconstructedRows();
	auto pSolMasksCompIdx = numPlayerSolutionsPtr() + numPreconstructedRows();
	auto first = *pSolMasksIniIdx;
	const auto lastB = IDX(pCompMask->numMasks());
	uint idxSol = 0;
	m_pSolIdx[0] = 0;
	//const auto n = m_numRec[1];
	while (idxSol < nValidSol) {
		// Skip all bytes/longs equal to 0
		auto firstB = first >> SHIFT;
		while (firstB < lastB && !pCompSol[firstB])
			firstB++;

#if USE_64_BIT_MASK
		unsigned long iBit;
		_BitScanForward64(&iBit, pCompSol[firstB]);
#else
		const auto iBit = this->firstOnePosition(pCompSol[firstB]);
#endif
		first = (firstB << SHIFT) + iBit;
		pCompSol[firstB] ^= (tmask)1 << iBit;

		if ((m_pSolIdx[++idxSol] = first) >= *pSolMasksIniIdx) {
			*pSolMasksCompIdx++ = idxSol;
			pSolMasksIniIdx++;
		}

		for (uint j = 1; j < idxSol; j++) {
			const auto idx = m_pSolIdx[j];
			/*
			for (auto k = j + 1; k < idxSol; k++) {
				const auto solIdx = m_pSolIdx[k];
				ASSERT_IF(solIdx >= idxSol);
				if (CHECK_MASK_BIT(pCompMask->getSolutionMask(solIdx), idxIniSol)) {
					SET_MASK_BIT(rowsCompatMasks() + lenSolutionMask() * j, idxSol);
				}
			}*/

			if (CHECK_MASK_BIT(pCompMask->getSolutionMask(idx), first)) {
				SET_MASK_BIT(rowsCompatMasks() + lenSolutionMask() * j, idxSol);
			}
		}
	}

	*pSolMasksCompIdx = idxSol;
	setNumLongs2Skip();
}
