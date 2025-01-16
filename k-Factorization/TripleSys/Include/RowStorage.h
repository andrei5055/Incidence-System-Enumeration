#pragma once
#include <vector>
#include "k-SysSupport.h"
#include "Storage.h"

#define USE_64_BIT_MASK		!USE_CUDA
#define UseSolutionCliques	!USE_CUDA	
										// The graph whose vertices are the remaining solutions must have a maximum 
										// clique whose size is equal to the number of unconstructed rows of the matrix.
typedef long long ll;

#if USE_64_BIT_MASK
typedef ll tmask;
#define SHIFT						6
#else
typedef tchar tmask;
#define SHIFT						3
#endif

#define MASK_BIT(idx)				(tmask)1 << ((idx) & ((1<<SHIFT) - 1))	
#define IDX(n)						(n + (1<<SHIFT) - 1) >> SHIFT
#define REM(n)						(n % ((tmask)1<<SHIFT))			// remainder from division
#define SET_MASK_BIT(mask, idx)		mask[(idx) >> SHIFT] |= MASK_BIT(idx)
#define RESET_MASK_BIT(mask, idx)	mask[(idx) >> SHIFT] ^= MASK_BIT(idx)
#define CHECK_MASK_BIT(mask, idx)	(mask[(idx) >> SHIFT] & MASK_BIT(idx))


#include "CompSolGraph.h"

class alldata;

class CRowStorage : public CStorage<tchar> {
	typedef void(CRowStorage::* rowToBitmask)(ctchar* pRow, tmask *pMask) const;
public:
	CC CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects = 1000, alldata* pAllData = NULL);
	CC ~CRowStorage();
	CC void generateCompatibilityMasks(tmask* pMask, uint solIdx, uint idx) const;
	CC bool maskForCombinedSolutions(tmask* pMaskOut, uint& solIdx) const;
	CC void init();
	CC inline void reset()								{ m_numObjects = 0; }
	CC inline auto numPlayers() const					{ return m_numPlayers; }
	CC inline auto numDaysResult() const				{ return m_numDaysResult; }
	CC void addRow(ctchar* pRow, ctchar* pNeighbors);
	CC void initCompatibilityMasks(ctchar* u1fCycles = NULL);
	CC int initRowUsage(tchar** ppCompatibleSolutions, uint* pRowSolutionLastIdx) const;
	CC inline auto numPreconstructedRows() const		{ return m_numPreconstructedRows; }
	CC inline auto numSolutionTotalB() const			{ return m_numSolutionTotalB; }
	CC inline auto numRowSolutionsPtr() const			{ return m_pRowSolutionCntr; }
	CC inline auto getSolutionMask(uint solNumb) const  { return m_pRowsCompatMasks[1] + (solNumb + m_solAdj) * m_lenSolutionMask; }
	CC inline auto numLongs2Skip(int iRow) const		{ return m_pNumLongs2Skip[iRow]; }
	CC inline const auto rowSolutionMasksIdx() const	{ return m_pRowSolutionMasksIdx; }
	CC inline const auto rowSolutionMasks() const		{ return m_pRowSolutionMasks; }
	CC inline auto useCliquesAfterRow() const			{ return m_useCliquesAfterRow; }
	CC inline auto useCliques(int iRow) const			{ return iRow > m_useCliquesAfterRow; }
	CC inline const kSysParam* sysParam() const			{ return m_pSysParam; }
	CC inline const auto numRecAdj() const				{ return m_numRecAdj; }
	CC inline const auto numRec(int idx) const			{ return m_numRec[1]; }
	CC inline const auto lastInFirstSet() const			{ return m_lastInFirstSet; }
	CC void getMatrix(tchar* row, tchar* neighbors, int nRows, uint* pRowSolutionIdx) const {
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
#if !USE_64_BIT_MASK
	CC inline auto firstOnePosition(tchar byte) const	{ return m_FirstOnePosition[byte]; }
private:
	tchar m_FirstOnePosition[256]; // Table for fast determination of the first 1's position in byte.
#endif


#define SetMask(bm, ir, ic) \
	const int ib = ir * m_numPlayers + ic - (1 + ir) * (2 + ir) / 2; \
	ASSERT(ib / 8 >= m_lenMask); \
	SET_MASK_BIT(bm, ib)

#define SetMaskWithOffset(bm, iOffset, ic) \
	ib = iOffset + ic; \
	ASSERT(ib / 8 >= m_lenMask); \
	SET_MASK_BIT(bm, ib)

private:
	CC void rowToBitmask2(ctchar* pRow, tmask* bm) const
	{
		memset(bm, 0, m_lenMask);
		for (int i = 0; i < m_numPlayers; i += 2) {
			SetMask(bm, pRow[i], pRow[i + 1]);
		}
	}
	CC void rowToBitmask3(ctchar* pRow, tmask* bm) const
	{
		int ib;
		memset(bm, 0, m_lenMask);
		for (int i = 0; i < m_numPlayers; i += 3, pRow += 3) {
			const int iOffset = pRow[0] * m_numPlayers - (1 + pRow[0]) * (2 + pRow[0]) / 2;
			SetMaskWithOffset(bm, iOffset, pRow[1]);
			SetMaskWithOffset(bm, iOffset, pRow[2]);
			SetMask(bm, pRow[1], pRow[2]);
		}
	}
	CC void initMaskStorage(uint numObjects) {
		m_pMaskStorage = new CStorage<tchar>(numObjects, (((m_numPlayers * (m_numPlayers - 1) / 2) + 63) / 64) * 8);
		memset(m_pRowSolutionCntr, 0, m_numPlayers * sizeof(m_pRowSolutionCntr[0]));
		reset();
	}
	CC bool p1fCheck2(ctchar* neighborsi, ctchar* neighborsj) const;
	CC bool checkCompatibility(ctchar* neighborsi, const long long* rm, uint idx) const;

	const kSysParam* m_pSysParam;
	const int m_numPlayers;
	const int m_numPreconstructedRows;     // Number of preconstructed matrix rows
	const int m_numDaysResult;
	const alldata* m_pAllData;

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
	rowToBitmask m_pRowToBitmask;
	uint* m_pRowSolutionCntr = NULL;
	uint m_numSolutionTotal;
	uint m_numSolutionTotalB;
	uint m_lenSolutionMask;
	tmask* m_pRowsCompatMasks[2];
	// For each row of the matrix, we define two masks, each containing an interval of consecutive bits set to 1
	// These intervals represent the row's first and last sets of solutions that lie outside the separately tested 64-bit intervals.
	tmask* m_pRowSolutionMasks = NULL;
	//  ... and the the set of indices of the long long elements which corresponds to two mask's sets.                           
	uint* m_pRowSolutionMasksIdx = NULL; 
	ctchar* m_u1fCycles = NULL;
	uint *m_pNumLongs2Skip = NULL; // Pointer to the number of long long's that we don't need to copy for each row.
	int m_useCliquesAfterRow;
	uint m_numRecAdj = 0;
	uint m_numRec[2];				   // Number of solutions for first and second nonfixed rows.
	uint m_numRecAdj2;
	uint m_lastInFirstSet;
	const bool m_bUseCombinedSolutions;
	const int m_step;
	int m_solAdj = 0;
	long long m_playersMask = 0;       // Mask with bits corresponding to players from first group of predefined rows equal to zeros.
};

class CRowUsage : public CompSolStorage {
public:
	CC CRowUsage(const CRowStorage* const pRowStorage) : CompSolStorage(pRowStorage) {
		const auto len = 2 * m_pRowStorage->numDaysResult();
		m_pRowSolutionIdx = new uint[len];
		memset(m_pRowSolutionIdx, 0, len * sizeof(m_pRowSolutionIdx[0]));
		m_pRowSolutionLastIdx = m_pRowSolutionIdx + (len >> 1);
		m_bUseCombinedSolutions = pRowStorage->sysParam()->val[t_useCombinedSolutions];
	}
	CC ~CRowUsage() {
		delete[] m_pRowSolutionIdx;
		delete[] m_pCompatibleSolutions;
	}
	CC void init(int iThread = 0, int numThreads = 1);
	CC int getRow(int iRow, int ipx);
	CC inline bool getMatrix2(tchar* row, tchar* neighbors, int nRows, int iRow) {
		getMatrix(row, neighbors, ++iRow);
		return completeMatrix(row, neighbors, nRows, iRow);
	}
	CC inline void getMatrix(tchar* row, tchar* neighbors, int nRows) {
		m_pRowStorage->getMatrix(row, neighbors, nRows, m_pRowSolutionIdx);
	}
private:
	uint m_numSolutionTotalB;
	uint* m_pRowSolutionIdx = NULL;
	uint* m_pRowSolutionLastIdx = NULL;
	tchar* m_pCompatibleSolutions = NULL;
	int m_step;
	bool m_bUseCombinedSolutions;
	bool m_bSolutionReady;   // true, when solution was prepared ar a part of combined solution.  
};
