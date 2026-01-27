#include "Table.h"
//#include "CompatMasks.h"

CCompatMasks::CCompatMasks(const kSysParam* pSysParam, const alldata* pAllData) :
	m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]),
	m_numDaysResult(pAllData->numDaysResult()),
	m_bSelectPlayerByMask(pAllData->groupSize() > 2 || pAllData->param(t_CBMP_Graph) > 1) {
	const auto numPlayers = pSysParam->paramVal(t_numPlayers);
	m_pPlayerSolutionCntr = new uint[numPlayers + m_numDaysResult];
	m_pNumLongs2Skip = m_pPlayerSolutionCntr + numPlayers;
	memset(numPlayerSolutionsPtr(), 0, m_numDaysResult * sizeof(m_pPlayerSolutionCntr[0]));
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

CC void CCompatMasks::initRowUsage(tmask** ppCompatibleSolutions, bool* pUsePlayersMask) const {
	const auto lenMask = m_numSolutionTotalB >> (SHIFT - 3);
	if (!*ppCompatibleSolutions) {
		const auto len = (numDaysResult() - numPreconstructedRows()) * lenMask;
		*ppCompatibleSolutions = new tmask[len];
	}

	*pUsePlayersMask = selectPlayerByMask();
}

void CCompressedMask::compressCompatMasks(tmask* pCompSol, uint nValidSol, const CCompatMasks* pCompMask)
{
	initMaskMemory(nValidSol, (selectPlayerByMask() ? 8 : 0));
	releaseSolIndices();
	m_pSolIdx = new uint[nValidSol];

	auto const *pSolMasksIniIdx = pCompMask->numPlayerSolutionsPtr() + numPreconstructedRows();
	auto pSolMasksCompIdx = numPlayerSolutionsPtr() + numPreconstructedRows();
	auto first = *pSolMasksIniIdx;
	const auto lastB = IDX(pCompMask->numMasks());
	uint idxSol = 0;
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

		if ((m_pSolIdx[idxSol++] = first) >= *pSolMasksIniIdx) {
			*pSolMasksCompIdx++ = idxSol;
			pSolMasksIniIdx++;
		}

		for (uint j = 1; j < idxSol; j++) {
			const auto idx = m_pSolIdx[j];
			if (CHECK_MASK_BIT(pCompMask->getSolutionMask(idx), idx)) {
				SET_MASK_BIT(rowsCompatMasks() + lenSolutionMask() * j, idxSol);
			}
		}
	}

	*pSolMasksCompIdx = idxSol;
}
