#include "Table.h"
//#include "CompatMasks.h"

CCompatMasks::CCompatMasks(const kSysParam* pSysParam, const alldata* pAllData) :
	m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]),
	m_numDaysResult(pAllData->numDaysResult()),
	m_bSelectPlayerByMask(pAllData->groupSize() > 2 || pAllData->param(t_CBMP_Graph) > 1) {
	const auto numPlayers = pSysParam->paramVal(t_numPlayers);
	m_pPlayerSolutionCntr = new uint[numPlayers + m_numDaysResult];
	m_pNumLongs2Skip = m_pPlayerSolutionCntr + numPlayers;
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

void CCompressedMask::compressCompatMasks(tmask* pCompSol, uint nValidSol, uint last, const CCompatMasks* pCompMask)
{
	initMaskMemory(nValidSol, (selectPlayerByMask() ? 8 : 0));

	auto first = last;
	last = pCompMask->numMasks();
	const auto lastB = IDX(last);
	while (true) {
		// Skip all bytes/longs equal to 0
		auto firstB = first >> SHIFT;
		while (firstB < lastB && !pCompSol[firstB])
			firstB++;

		if (firstB >= lastB)
			return;

#if USE_64_BIT_MASK
		unsigned long iBit;
		_BitScanForward64(&iBit, pCompSol[firstB]);
#else
		const auto iBit = this->firstOnePosition(pCompSol[firstB]);
#endif
		if ((first = (firstB << SHIFT) + iBit) >= last)
			break;

		pCompSol[firstB] ^= (tmask)1 << iBit;
	}
}
