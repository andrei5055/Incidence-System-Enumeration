#pragma once
#include <vector>
#include "k-SysSupport.h"
#include "Storage.h"

#define NEW					1

#define NEW_GET_ROW			1
#define USE_64_BIT_MASK		!USE_CUDA
#define UseSolutionMasks	1
#define UseSolutionCliques	!USE_CUDA	// The graph whose vertices are the remaining solutions must have a maximum 
										// clique whose size is equal to the number of unconstructed rows of the matrix.

#if USE_64_BIT_MASK
typedef long long tmask;
#define SHIFT					6 
#define IDX(n)					(n + 63) >> 6
#define SET_MASK_BIT(mask, idx)	mask[(idx) >> 6] |= (tmask)1 << ((idx) & 0x3f)
#else
typedef tchar tmask;
#define SHIFT					3 
#define IDX(n)					(n + 7) >> 3
#define SET_MASK_BIT(mask, idx)	mask[(idx) >> 3] |= (tmask)1 << ((idx) & 0x07)
#endif
#define REM(n)					(n % ((tmask)1<<SHIFT))	// remainder from division

#include "CompSolGraph.h"

class CRowStorage : public CStorage<tchar> {
public:
	CC CRowStorage(const kSysParam *pSysParam, int numPlayers, int numObjects = 1000) : m_pSysParam(pSysParam), m_numPlayers(numPlayers),
		m_numPreconstructedRows(pSysParam->val[t_useRowsPrecalculation]), CStorage<tchar>(numObjects, 2 * numPlayers) {
		m_numObjectsMax = numObjects;
		m_pRowSolutionCntr = new uint[2*numPlayers];
		m_pNumLongs2Skip = m_pRowSolutionCntr + numPlayers;
		initMaskStorage(numObjects);
		m_lenMask = m_pMaskStorage->lenObject();
		const auto useCliquesAfterRow = pSysParam->val[t_useSolutionCliquesAfterRow];
		m_useCliquesAfterRow = useCliquesAfterRow ? useCliquesAfterRow : numPlayers;
	}
	CC ~CRowStorage() {
		delete[] m_pRowSolutionCntr;
		delete[] m_fullExcludeTable;
		delete m_pMaskStorage;
		delete[] m_pRowSolutionMasks;
		delete[] m_pRowSolutionMasksIdx;
	}
	CC void init() {
		initMaskStorage(m_numObjectsMax);
	}
	CC void addRow(ctchar* pRow, ctchar* pNeighbors) {
		if (m_numObjects == m_numObjectsMax) {
			reallocStorageMemory(m_numObjectsMax <<= 1);
			m_pMaskStorage->reallocStorageMemory(m_numObjectsMax);
		}

		row2bitmask(pRow, (tmask*)(m_pMaskStorage->getObject(m_numObjects)), false);
		auto* pntr = getObject(m_numObjects++);
		ASSERT(m_numObjects > 1 && pRow[1] != *(pntr - 2 * m_numPlayers + 1) &&
			pRow[1] != *(pntr - 2 * m_numPlayers + 1) + 1);
		memset(pntr, 0, m_lenObj);
		memcpy(pntr, pRow, m_numPlayers);
		memcpy(pntr + m_numPlayers, pNeighbors, m_numPlayers);
		m_pRowSolutionCntr[pRow[1] - 1]++;  // Increasing the number of solutions for (pRow[1]-1)-th row
	}
	CC inline auto numPlayers() const					{ return m_numPlayers; }
	CC void initCompatibilityMasks(ctchar* u1fCycles = NULL);

	CC inline auto numPreconstructedRows() const		{ return m_numPreconstructedRows; }
	CC inline auto numSolutionTotalB() const			{ return m_numSolutionTotalB; }
	CC inline auto numRowSolutions(int nRow) const		{ return m_pRowSolutionCntr[nRow]; }
	CC inline auto getSolutionMask(uint solNumb) const	{ return (const long long *)(m_fullExcludeTable + solNumb * m_numSolutionTotalB); }
	CC inline auto numLongs2Skip(int iRow) const		{ return m_pNumLongs2Skip[iRow]; }
	CC inline const auto rowSolutionMasksIdx() const	{ return m_pRowSolutionMasksIdx; }
	CC inline const auto rowSolutionMasks() const		{ return m_pRowSolutionMasks; }
	CC inline auto useCliquesAfterRow() const			{ return m_useCliquesAfterRow; }
	CC inline auto useCliques(int iRow) const			{ return iRow > m_useCliquesAfterRow; }
	CC inline const kSysParam* sysParam() const			{ return m_pSysParam; }
	CC inline const auto numRecAdj() const				{ return m_numRecAdj; }
	CC void getMatrix(tchar* row, tchar* neighbors, int nRows, int iStep, const uint* pRowSolutionIdx) const {
		size_t shift = numPreconstructedRows() * m_numPlayers;
		const int adj = numRecAdj();
		for (int iRow = numPreconstructedRows(); iRow < nRows; iRow++) {
			const auto* pObj = getObject(pRowSolutionIdx[iRow] - iStep);
			iStep = 1 - adj;
			memcpy(row + shift, pObj, m_numPlayers);
			memcpy(neighbors + shift, pObj + m_numPlayers, m_numPlayers);
			shift += m_numPlayers;
		}
	}
#if !(USE_64_BIT_MASK && NEW_GET_ROW)
	CC inline auto firstOnePosition(tchar byte) const	{ return m_FirstOnePosition[byte]; }
private:
	tchar m_FirstOnePosition[256]; // Table for fast determination of the first 1's position in byte.
#endif

private:
	CC void row2bitmask(ctchar* pRow, tmask* bm, bool bAdd)
	{
		if (!bAdd)
			memset(bm, 0, m_lenMask);
		for (int i = 0; i < m_numPlayers; i += 2) {
			auto ir = pRow[i];
			auto ic = pRow[i + 1];
			auto ib = ir * m_numPlayers + ic - (1 + ir) * (2 + ir) / 2;
			ASSERT(ib / 8 >= m_lenMask);
			SET_MASK_BIT(bm, ib);
		}
	}
	CC void initMaskStorage(uint numObjects) {
		m_pMaskStorage = new CStorage<tchar>(numObjects, (((m_numPlayers * (m_numPlayers - 1) / 2) + 63) / 64) * 8);
		memset(m_pRowSolutionCntr, 0, m_numPlayers * sizeof(m_pRowSolutionCntr[0]));
		m_numObjects = 0;
	}

	const int m_numPreconstructedRows;     // Number of preconstructed matrix rows
	const int m_numPlayers;
	const kSysParam* m_pSysParam;
	CStorage<tchar>* m_pMaskStorage = NULL;
#if 0
	long long cnt1 = 0;
	long long cnt2[18] = { 0 };
	long long cnt3[18] = { 0 };
	long long cnt4[18] = { 0 };
	long long cnt5[18] = { 0 };
#endif
	uint m_numObjects;
	uint m_numObjectsMax;
	int m_lenMask;
	uint* m_pRowSolutionCntr = NULL;
	uint m_numSolutionTotalB;
	tchar* m_fullExcludeTable = NULL;
	// For each row of the matrix, we define two masks, each containing an interval of consecutive bits set to 1
	// These intervals represent the row's first and last sets of solutions that lie outside the separately tested 64-bit intervals.
	tmask* m_pRowSolutionMasks = NULL;
	//  ... and the the set of indices of the long long elements which corresponds to two mask's sets.                           
	uint* m_pRowSolutionMasksIdx = NULL; 
	ctchar* m_u1fCycles = NULL;
	uint *m_pNumLongs2Skip = NULL; // Pointer to the number of long long's that we don't need to copy for each row.
	int m_useCliquesAfterRow;
	int m_numRecAdj = 0;
};

class CRowUsage : public CompSolStorage {
public:
	CC CRowUsage(const CRowStorage* const pRowStorage) : CompSolStorage(pRowStorage) {
		const auto numPlayers = pRowStorage->numPlayers();
		m_pRowSolutionIdx = new uint[numPlayers + 1];
		memset(m_pRowSolutionIdx, 0, numPlayers * sizeof(m_pRowSolutionIdx[0]));
	}

	CC ~CRowUsage() {
		delete[] m_pRowSolutionIdx;
		delete[] m_pCompatibleSolutions;
	}
	CC void init(int iThread = 0, int numThreads = 1) {
		m_numSolutionTotalB = m_pRowStorage->numSolutionTotalB();
		const auto len = (m_pRowStorage->numPlayers() - m_pRowStorage->numPreconstructedRows() - 1) * m_numSolutionTotalB;
		m_pCompatibleSolutions = new tchar[len];
		m_pRowSolutionIdx[m_pRowStorage->numPreconstructedRows()] = iThread;
		m_step = numThreads;
	}
	CC bool getMatrix2(tchar* row, tchar* neighbors, int nRows, int iRow) {
		getMatrix(row, neighbors, ++iRow);
		return completeMatrix(row, neighbors, nRows, iRow);
	}
	CC void getMatrix(tchar * row, tchar * neighbors, int nRows) {
		m_pRowStorage->getMatrix(row, neighbors, nRows, m_step, m_pRowSolutionIdx);
	}
	CC int getRow(int iRow, int ipx);
private:
	uint m_numSolutionTotalB;
	uint* m_pRowSolutionIdx = NULL;
	tchar* m_pCompatibleSolutions = NULL;
	int m_step;
#if NEW_GET_ROW == 0
	tchar* m_excludeForRow[MAX_PLAYER_NUMBER];
#endif
};
