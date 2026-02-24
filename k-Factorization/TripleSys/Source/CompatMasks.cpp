#include "Table.h"
//#include "CompatMasks.h"

void initTripleSysData() {
	CCompatMasks::initBitsInByteTable();
}

tchar CCompatMasks::m_pBitsInByte[256];

void CCompatMasks::initBitsInByteTable() {
	// Preparing the table to calculate the weight of the mask matrix.
	m_pBitsInByte[0] = 0;
	for (int i = 1; i < sizeof(m_pBitsInByte); i++)
		m_pBitsInByte[i] = m_pBitsInByte[i >> 1] + i % 2;
}

CCompatMasks::CCompatMasks(const kSysParam* pSysParam, const alldata* pAllData) :
	m_numPlayers(pSysParam->paramVal(t_numPlayers)),
	m_numPreconstructedRows(pSysParam->paramVal(t_useRowsPrecalculation)),
	m_numDaysResult(pAllData->numDaysResult()),
	m_bSelectPlayerByMask(pAllData->groupSize() > 2 || pAllData->param(t_CBMP_Graph) > 1),
	m_pSysParam(pSysParam), m_pAllData(pAllData) {
	m_pPlayerSolutionCntr = new uint[numPlayers() + numDaysResult()];
	m_pNumLongs2Skip = numPlayerSolutionsPtr() + numPlayers();
	memset(numPlayerSolutionsPtr(), 0, numPlayers() * sizeof(m_pPlayerSolutionCntr[0]));
	m_lenDayResults = numDaysResult() + 1;
	m_fSolutionInterval = selectPlayerByMask() ? &CCompatMasks::solutionInterval3 : &CCompatMasks::solutionInterval2;
	initSolMaskIndices();

#if !USE_64_BIT_MASK
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif
}

CCompatMasks::~CCompatMasks() {
	releaseCompatMaskMemory();
	releaseSolMaskInfo();
	delete[] numPlayerSolutionsPtr();
	delete[] m_pCompatSolutions;
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
	const auto lenMaskIdx = numPlayers();
	releaseSolMaskInfo();
	m_pRowSolutionMasksIdx = new uint[lenMaskIdx];

	m_pRowSolutionMasks = new tmask[lenMaskIdx];
	memset(m_pRowSolutionMasks, 0, lenMaskIdx * sizeof(m_pRowSolutionMasks[0]));
	*rowSolutionMasksIdx() = 0;
#if SAME_MASK_IDX
	m_pMaskTestingCompleted = new bool[lenMaskIdx];
	memset(maskTestingCompleted(), 0, lenMaskIdx * sizeof(m_pMaskTestingCompleted[0]));
#endif
}

void CCompatMasks::resetSolMaskIndices(bool resetSolutionForPlayers) {
	if (!rowSolutionMasksIdx())	 // Previously allocated memory for thee indices could be released 
		initSolMaskIndices();    // Let's allocate it again 

	const auto lenMaskIdx = numPlayers();
	memset(rowSolutionMasksIdx(), 0, lenMaskIdx * sizeof(m_pRowSolutionMasksIdx[0]));
	memset(m_pRowSolutionMasks, 0, lenMaskIdx * sizeof(m_pRowSolutionMasks[0]));
	if (resetSolutionForPlayers)
		memset(numPlayerSolutionsPtr(), 0, lenMaskIdx * sizeof(m_pPlayerSolutionCntr[0]));
}

void CCompatMasks::initMaskMemory(uint nSolutions, int nMasks, int numRecAdj, int numSolAdj)
{
	// Adding additional long long when we use groupSize > 2
	setNumSolutions(nSolutions);
	m_numSolutionTotalB = ((numSolutions() - numRecAdj + 7) / 8 + 7) / 8 * 8 + (selectPlayerByMask() ? AVALABLE_PLAYER_MASK_LENGTH : 0);
	m_lenSolutionMask = numSolutionTotalB() / sizeof(tmask);

	m_solAdj = numSolAdj;
	releaseCompatMaskMemory();

	m_pRowsCompatMasks = new tmask[nMasks * lenSolutionMask()];
	memset(rowsCompatMasks(), 0, nMasks * numSolutionTotalB());
}

CC void CCompatMasks::initRowUsage(tmask** ppCompatibleSolutions, bool fullMatrix, bool* pUsePlayersMask, ll* pAvailablePlayerMask) const {
	if (!*ppCompatibleSolutions) {
		auto lenMask = numSolutionTotalB() >> (SHIFT - 3);
		if (fullMatrix)
			lenMask *= (numDaysResult() - numPreconstructedRows());

		*ppCompatibleSolutions = new tmask[lenMask];
	}

	*pUsePlayersMask = selectPlayerByMask();
	if (pAvailablePlayerMask) {
		auto pCompatibleSolutions = *ppCompatibleSolutions;
		// NOTE; Let's make a trivial mask for now and improve it later
		memset(pCompatibleSolutions, 0xff, numSolutionTotalB());
		pCompatibleSolutions[0]--;
		if (selectPlayerByMask())
			pCompatibleSolutions[lenSolutionMask() - 1] = *pAvailablePlayerMask;
	}
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

	const auto iBit = minBit(availablePlayers);
	*pLast = pRowSolutionIdx[lenDayResults()] = m_pPlayerSolutionCntr[iBit - 1];
	return *pRowSolutionIdx = m_pPlayerSolutionCntr[iBit - 2];
}

int CCompatMasks::setNumLongs2Skip(bool adjustRecCounter)
{
	// Define the number of first long long's we don't need to copy to the next row.
	int i = numPreconstructedRows();
	setLastInFirstSet(m_pPlayerSolutionCntr[i]);
	const int recAdj = adjustRecCounter ? lastInFirstSet() : 0;

	m_pNumLongs2Skip[i] = lastInFirstSet() >> 6;
	while (++i < numDaysResult())
		m_pNumLongs2Skip[i] = ((m_pPlayerSolutionCntr[i] - recAdj) >> 6);

	return recAdj;
}

CC uint CCompatMasks::getSolutionRange(uint& last, ll& availablePlayers, int i) const {
	const uint first = last;
	if (selectPlayerByMask()) {
		const auto iBit = minBit(availablePlayers);
		last = m_pPlayerSolutionCntr[iBit - 1];
		availablePlayers ^= (ll)1 << iBit;
	}
	else {
		last = m_pPlayerSolutionCntr[i];
		availablePlayers &= (availablePlayers - 1);
	}

	return first;
}

void CCompatMasks::defineSolutionIntervals() {
	for (int i = 1; i < numPlayers(); i++)
		m_pPlayerSolutionCntr[i] += m_pPlayerSolutionCntr[i - 1];

	setNumSolution(m_pPlayerSolutionCntr[numPlayers() - 1]);
}

void CCompatMasks::defineMask4SolutionIntervals(int nRow, uint last)
{
#if PRINT_MASK_SOLUTION_INTERVALS
	static int fff = 0;
	static int ggg = 0;
	static int ccc = 0;
	FOPEN_F(f, "ccc.txt", ccc++? "a" : "w");
	if (fff > nRow) {
		fprintf(f, "\n*** %3d: Compressed Matrix info ***\n", ++ggg);
		if (ggg >= 428)
			ggg += 0;
	}

	fff = nRow;
	fprintf(f, "nRow = %2d  last = %4d\n", nRow, last);
	FCLOSE_F(f);
#endif
	auto* pRowSolutionMasksIdx = rowSolutionMasksIdx();
	if (!pRowSolutionMasksIdx)
		return;

	const auto lastAdj = last - numRecAdj();
	m_pRowSolutionMasksIdx[nRow] = lastAdj >> SHIFT;
	if (nRow == numPreconstructedRows() || maskTestingCompleted() || m_pRowSolutionMasksIdx[nRow] > m_pRowSolutionMasksIdx[nRow - 1]) {
		unsigned int rem;
		if (rem = REM(lastAdj)) {
			m_pRowSolutionMasks[nRow] = (tmask)(-1) << rem;
			if (pRowSolutionMasksIdx[nRow - 1] == pRowSolutionMasksIdx[nRow]) {
				// At this point, the solutions for three consecutive rows are masked by the same `tmask` element
				// Clear the bits in the previous mask that do not correspond to the solutions of the previous row
				m_pRowSolutionMasks[nRow - 1] ^= pRowSolutionMasksIdx[nRow]; // ???
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
		SPRINTFD(pBuf, buf, "%18d", *(rowSolutionMasksIdx()+j));

	fprintf(f, "%s\n", buf);
	pBuf = buf;
	for (int j = 0; j <= nRow; j++)
		SPRINTFD(pBuf, buf, "  %016llx", m_pRowSolutionMasks[j]);

	fprintf(f, "%s\n", buf);
	FCLOSE_F(f);
#endif
}

#if COUNT_MASK_WEIGHT || OUT_MASK_WEIGHT
size_t totalWeighChange = 0;
#endif

uint CCompatMasks::countMaskFunc(const tmask* pCompSolutions, size_t numMatrRow, int bytes2skeep, uint prevWeight) const {
	ctchar* pMask = (ctchar*)pCompSolutions + bytes2skeep;
	size_t i = numSolutionTotalB() - bytes2skeep;
	if (selectPlayerByMask())
		i -= AVALABLE_PLAYER_MASK_LENGTH;

	uint weight = 0;
	if (numMatrRow == 1) {
		while (i--)
			weight += m_pBitsInByte[pMask[i]];
	}
#if OUT_MASK_WEIGHT
	else {
		static size_t matrWeight = 0;  // Max number of ones located above the upper main diagonal
		static uint numFirstSol = 0;   // Number of solutions for the first non-preconstruted row 
		if (!matrWeight) {
			auto availablePlayers = getPlayersMask();
			unsigned int last = 0;
			int i = numPreconstructedRows() - 1;
			unsigned int first = getSolutionRange(last, availablePlayers, ++i);
			numFirstSol = last - first;
			matrWeight = numMatrRow * (numSolutions() - m_numRecAdj);

			matrWeight -= numFirstSol * numFirstSol;
			while (availablePlayers) {
				first = getSolutionRange(last, availablePlayers, ++i);
				const auto portion = last - first;
				matrWeight -= portion * portion;
			}
		}

		uint totalWeight = 0;
		for (int j = 0; j < numMatrRow; j++) {
			auto k = i;
			while (k--)
				totalWeight += m_pBitsInByte[pMask[k]];

			pMask += numSolutionTotalB();
			if (j == numFirstSol)
				weight = totalWeight;
		}

		printf("\nMask matrix total weight: %d  %5.2f%%", totalWeight, double(totalWeight) / matrWeight * 100);
		printf("\nMask matrix changing weight: %d  %5.2f%%", weight, double(weight) / matrWeight * 100);
		if (prevWeight) {
			const auto weightChange = prevWeight - weight;
			totalWeighChange += weightChange;
			printf("\nWeight change is %d:  %5.2f%%\n", weightChange, double(weightChange) / prevWeight * 100);

		}
	}
#endif

	return weight;
}

int CCompatMasks::setNumMasks() {
	// Define the number of rows in mask table.
    // It is equal to number of solutions without the last set.
	auto i = numPlayers();
	while (m_pPlayerSolutionCntr[--i] == m_numObjects);
	return m_numMasks = m_pPlayerSolutionCntr[i] - (numRecAdj() - solAdj());
}

void CCompatMasks::outputPrecalcSolution(int num_obj, ctchar* pSol, const char* pFileName) const {
	FOPEN_F(f, pFileName, num_obj ? "a" : "w");
	char buf[32];
	if (!num_obj) {
		// Output of the last precomputed row
		sprintf_s(buf, "%5d: ", -1);
		outMatrix(allData()->result() + (numPreconstructedRows() - 1) * numPlayers(), 1, numPlayers(), allData()->groupSize(), 0, f, false, false, buf);
	}
	sprintf_s(buf, "%5d: ", num_obj);
	outMatrix(pSol, 1, numPlayers(), allData()->groupSize(), 0, f, false, false, buf, -1, NULL, false);
	FCLOSE_F(f);
}

void CCompressedMask::compressCompatMasks(tmask* pCompSol, const CCompatMasks* pCompMask, uint solIdx)
{
	// solIdx - index of the solution used for first non-predefined row

	// Count the number of valid solutions in the compatibility mask 
	// for the current solution of the first non-predefined row, including itself.
	const auto nValidSol = pCompMask->countMaskFunc(pCompSol) + 1;
	
	auto first = pCompMask->numMasks();
	const auto idxFirst = first >> 3;
	auto nn = pCompMask->countMaskFunc(pCompSol, 1, idxFirst);
	tchar bnd;
	if ((first -= (idxFirst << 3)) && (bnd = *((tchar*)pCompSol + idxFirst))) {
		// The bnd byte may contain units corresponding to two consecutive sets of solutions
		// We need to decrease nn by the number of solutions from the previous set. 
		do {
			if (bnd & 1)
				nn--;
		} while (--first && (bnd >>= 1));	
	}

	const auto nMasks = nValidSol - nn;
	initMaskMemory(nValidSol, nMasks);
	resetSolMaskIndices();
	releaseSolIndices();
	m_pSolIdx = new uint[nValidSol];

	auto const* pSolMasksIniIdx = pCompMask->numPlayerSolutionsPtr() - 1;
	auto pSolMasksCompIdx = numPlayerSolutionsPtr() - 1;

	uint idxSol = 0;
	m_pSolIdx[0] = 0;

	auto* pUnusedPlaeyerMask = selectPlayerByMask()? rowsCompatMasks() - 1 : NULL;
	const auto offsetSrc = pCompMask->lenSolutionMask() - 1;
	auto nRow = numPreconstructedRows();
	pSolMasksCompIdx[nRow] = 1;
	OUTPUT_PRECALC_SOLUTION(idxSol, pCompMask->getSolution(solIdx), "bbb.txt");

	bool maskDefined = true;
	auto availablePlayers = pCompMask->getPlayersMask();
	auto playerIdx = minBit(availablePlayers);
	auto last = pSolMasksIniIdx[playerIdx];
	const auto recAdj = pCompMask->numRecAdj();
	uint jMax = idxSol;
	while (availablePlayers ^= (ll)1 << playerIdx) {
		defineMask4SolutionIntervals(nRow++, ++jMax);
		auto first = last;
		last = pSolMasksIniIdx[playerIdx = minBit(availablePlayers)];
		const auto lastB = IDX(last);
		while (true) {
			// Skip all bytes/longs equal to 0
			auto firstB = first >> SHIFT;
			while (firstB < lastB && !pCompSol[firstB])
				firstB++;

			if (firstB >= lastB)
				break;

			const auto iBit = minBit(pCompSol[firstB]);
			if ((first = (firstB << SHIFT) + iBit) >= last)
				break;

			pCompSol[firstB] ^= (tmask)1 << iBit;


			pSolMasksCompIdx[playerIdx]++;
			m_pSolIdx[++idxSol] = first;
			OUTPUT_PRECALC_SOLUTION(idxSol, pCompMask->getSolution(first + recAdj), "bbb.txt");

			//if (idxSol >= nMasks)
			//	continue;

			if (pUnusedPlaeyerMask && idxSol < nMasks) {
				// Copying the mask of available players.
				memcpy(pUnusedPlaeyerMask += lenSolutionMask(), pCompMask->getSolutionMask(first) + offsetSrc,
					AVALABLE_PLAYER_MASK_LENGTH);
			}

			for (uint j = 1; j < jMax; j++) {
				if (CHECK_MASK_BIT(pCompMask->getSolutionMask(m_pSolIdx[j]), first)) {
					ASSERT_IF(j >= nMasks);
					SET_MASK_BIT(rowsCompatMasks() + lenSolutionMask() * j, idxSol);
				}
			}
		}

		jMax = idxSol;
	}

	defineMask4SolutionIntervals(nRow++, idxSol++);

	ASSERT_IF(idxSol != nValidSol);
	if (recAdj) {
		// Adjust the solution indices.
		for (uint j = 1; j < nValidSol; j++)
			m_pSolIdx[j] += recAdj;
	}

	defineSolutionIntervals();
	setNumLongs2Skip();
}
