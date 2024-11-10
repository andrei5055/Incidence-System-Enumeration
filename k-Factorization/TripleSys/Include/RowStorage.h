#pragma once
#include "k-SysSupport.h"
#include "Storage.h"

#define NEW_GET_ROW 1
#define UseSolutionMasks 1
#define USE_64_BIT_MASK		!USE_CUDA

#if USE_64_BIT_MASK
typedef long long tmask;
#define IDX(n)					(n + 63) >> 6
#define REM(n)					(n % 64)			// remainder from division
#define SET_MASK_BIT(mask, idx)	mask[idx >> 6] |= (tmask)1 << (idx & 0x3f)
#else
typedef tchar tmask;
#define IDX(n)					(n + 7) >> 3
#define REM(n)					(n % 8)				// remainder from division
#define SET_MASK_BIT(mask, idx)	mask[idx >> 3] |= (tmask)1 << (idx & 0x07)
#endif

typedef unsigned int uint;

class CRowStorage : public CStorage<tchar> {
public:
	CC CRowStorage(int numPreconstructedRows, int numPlayers, int numObjects = 1000) : m_numPlayers(numPlayers), 
		m_numPreconstructedRows(numPreconstructedRows), CStorage<tchar>(numObjects, 2 * numPlayers)	
	{
		m_numObjectsMax = numObjects;
		m_pRowSolutionCntr = new uint[2*numPlayers];
		m_pNumLongs2Skip = m_pRowSolutionCntr + numPlayers;
		initMaskStorage(numObjects);
		m_lenMask = m_pMaskStorage->lenObject();
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
	CC void initCompatibilityMasks(tchar* u1fCycles = NULL)
	{
		m_u1fCycles = u1fCycles;
		for (int i = 1; i < m_numPlayers; i++)
			m_pRowSolutionCntr[i] += m_pRowSolutionCntr[i - 1];

		const auto& numSolutionTotal = m_pRowSolutionCntr[m_numPlayers - 1];
		delete[] m_fullExcludeTable;
		m_numSolutionTotalB = ((numSolutionTotal + 7) / 8 + 7) / 8 * 8;
		auto len = numSolutionTotal * m_numSolutionTotalB;
		m_fullExcludeTable = new tchar[len];
		memset(m_fullExcludeTable, 0, len);


		delete[] m_pRowSolutionMasksIdx;
		delete[] m_pRowSolutionMasks;
		len = numPlayers();
		m_pRowSolutionMasksIdx = new uint[len];

		m_pRowSolutionMasks = new tmask[len];
		memset(m_pRowSolutionMasks, 0, len * sizeof(tmask));

#if !USE_64_BIT_MASK || !NEW_GET_ROW
		// Filling the lookup table m_FirstOnePosition
		memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
		for (int i = 2; i < 256; i += 2)
			m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

		// Define the number of first long long's we don't need to copy to the next row.
		memset(m_pNumLongs2Skip, 0, m_numPlayers * sizeof(m_pNumLongs2Skip[0]));
		for (int i = m_numPreconstructedRows; i < m_numPlayers - 1; i++)
			m_pNumLongs2Skip[i+1] = m_pRowSolutionCntr[i] >> 6;

		const auto jMax = m_lenMask >> 3;
		auto* pFullIncludeTable = (tmask *)m_fullExcludeTable;
		unsigned int last = 0;
		m_pRowSolutionMasksIdx[0] = 0;
		int i = m_numPreconstructedRows;
		while (i < m_numPlayers) {
			auto first = last;
			last = m_pRowSolutionCntr[i++];
			m_pRowSolutionMasksIdx[i - 1] = IDX(last);
			if (REM(last))
				m_pRowSolutionMasks[i-1] = (tmask)(-1) << REM(last);

			while (first < last) {
				auto* rm = (const long long*)m_pMaskStorage->getObject(first);
				const auto pRow = getObject(first++);
				ASSERT(pRow[1] != i);
				const auto pNeighbors = pRow + m_numPlayers;
				auto idx = last - 1;
				while (++idx < numSolutionTotal) {
					// Let's check if the masks are mutually compatible
					auto* pMask = (const long long*)(m_pMaskStorage->getObject(idx));
					int j = jMax;
					while (j-- && !(rm[j] & pMask[j]));

					if (j < 0 && p1fCheck2(m_u1fCycles, pNeighbors, getObject(idx) + m_numPlayers, m_numPlayers)) {
						// The masks are compatible and the length of the cycle is equal to m_numPlayers
						SET_MASK_BIT(pFullIncludeTable, idx);     // 1 - means OK
					}
				}

				pFullIncludeTable += m_numSolutionTotalB / sizeof(tmask);
			}
		}
		delete m_pMaskStorage;
		m_pMaskStorage = NULL;
	}

	CC inline auto numPlayers() const					{ return m_numPlayers; }
	CC inline auto numPreconstructedRows() const		{ return m_numPreconstructedRows; }
	CC inline auto numSolutionTotalB() const			{ return m_numSolutionTotalB; }
	CC inline auto numRowSolutions(int nRow) const		{ return m_pRowSolutionCntr[nRow]; }
	CC inline auto getSolutionMask(uint solNumb) const	{ return (const long long *)(m_fullExcludeTable + solNumb * m_numSolutionTotalB); }
	CC inline auto numLongs2Skip(int iRow) const		{ return m_pNumLongs2Skip[iRow]; }
	CC inline auto rowSolutionMasksIdx() const			{ return m_pRowSolutionMasksIdx; }
	CC inline auto rowSolutionMasks() const				{ return m_pRowSolutionMasks; }
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
	tchar* m_u1fCycles = NULL;
	int m_lenType[4];              // The number of long long's, int's short's and bytes  in the solution mask that need to be copied to the next row.
	uint *m_pNumLongs2Skip = NULL; // Pointer to the number of long long's that we don't need to copy for each row.
};

class CRowUsage {
public:
	CC CRowUsage(const CRowStorage* pRowStorage) : m_pRowStorage(pRowStorage), m_nRowMax(pRowStorage->numPlayers() - 2) {
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
	CC void getMatrix(tchar* row, tchar* neighbors, int nRows) {
		const auto numPlayers = m_pRowStorage->numPlayers();
		auto iStep = m_step;
		for (int iRow = m_pRowStorage->numPreconstructedRows(); iRow < nRows; iRow++) {
			auto& first = m_pRowSolutionIdx[iRow];
			auto* pObj = m_pRowStorage->getObject(first - iStep);
			iStep = 1;
			memcpy(row + iRow * numPlayers, pObj, numPlayers);
			memcpy(neighbors + iRow * numPlayers, pObj + numPlayers, numPlayers);
		}

		FILE* f;
		extern int cntr;
		static int nMatr= 0;
		fopen_s(&f, "aaa.txt", nMatr++? "a" : "w");
		fprintf(f, "===========  nMatr = %2d  cntr = %d \n", nRows, cntr);
		fclose(f);
	}
	CC bool getRow(int iRow, int ipx) const;
private:
	const CRowStorage* m_pRowStorage;
	const int m_nRowMax;				// Maximum value of iRow
	uint m_numSolutionTotalB;
	uint* m_pRowSolutionIdx = NULL;
	tchar* m_pCompatibleSolutions = NULL;
	int m_step;
#if NEW_GET_ROW == 0
	tchar* m_excludeForRow[MAX_PLAYER_NUMBER];
#endif
};