#pragma once
#include <vector>
#include "k-SysSupport.h"
#include "Storage.h"

#define SAME_MASK_IDX		0			// We allow the same mask index to be used for tree consecutive row. 
                                        // If 0, we will not apply the acceleration method that analyzes valid solutions for the remaining rows in such situations. 
#define USE_64_BIT_MASK		!USE_CUDA
#define UseSolutionCliques	!USE_CUDA	
										// The graph whose vertices are the remaining solutions must have a maximum 
										// clique whose size is equal to the number of unconstructed rows of the matrix.
#define COUNT_GET_ROW_CALLS  0          // Trace and print the number of calls of CRowUsage::getRow() method 

typedef long long ll;

#if USE_64_BIT_MASK
typedef ll tmask;
#define SHIFT						6
#else
typedef tchar tmask;
#define SHIFT						3
#endif

#define MASK_BIT(idx)				((tmask)1 << ((idx) & ((1<<SHIFT) - 1)))	
#define IDX(n)						(n + (1<<SHIFT) - 1) >> SHIFT
#define REM(n)						(n % ((tmask)1<<SHIFT))			// remainder from division
#define SET_MASK_BIT(mask, idx)		mask[(idx) >> SHIFT] |= MASK_BIT(idx)
#define RESET_MASK_BIT(mask, idx)	mask[(idx) >> SHIFT] ^= MASK_BIT(idx)
#define CHECK_MASK_BIT(mask, idx)	(mask[(idx) >> SHIFT] & MASK_BIT(idx))


#include "CompSolGraph.h"

class alldata;

class CRowStorage : public CStorage<tchar> {
	typedef void (CRowStorage::* rowToBitmask)(ctchar* pRow, tmask *pMask) const;
	typedef uint& (CRowStorage::* solutionInterval)(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const;
public:
	CC CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects = 1000, const alldata* pAllData = NULL);
	CC ~CRowStorage();
	CC inline void init()								{ initMaskStorage(m_numObjectsMax); }
	CC void initPlayerMask(ctchar* pFirstMatr = NULL);
	CC bool maskForCombinedSolutions(tmask* pMaskOut, uint& solIdx) const;
	CC inline uint& getSolutionInterval(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const {
		return (this->*m_fSolutionInterval)(pRowSolutionIdx, pLast, availablePlayers);
	}
	CC void passCompatibilityMask(tmask *pCompatibleSolutions, uint first, uint last) const;
	CC bool checkSolutionByMask(int iRow, const tmask* pToASol) const;
	CC inline void reset()								{ m_numObjects = 0; }
	CC inline auto numPlayers() const					{ return m_numPlayers; }
	CC inline auto numDaysResult() const				{ return m_numDaysResult; }
	CC bool addRow(ctchar* pRow, ctchar* pNeighbors);
	CC bool initCompatibilityMasks();
	CC int initRowUsage(tmask** ppCompatibleSolutions, bool *pUsePlayersMask) const;
	CC inline auto numPreconstructedRows() const		{ return m_numPreconstructedRows; }
	CC inline auto numSolutionTotalB() const			{ return m_numSolutionTotalB; }
	CC inline auto numPlayerSolutionsPtr() const		{ return m_pPlayerSolutionCntr; }
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
	CC inline const auto getPlayersMask() const			{ return m_playersMask[0]; }
	CC inline auto maskTestingCompleted() const			{ return m_pMaskTestingCompleted; }
	CC inline auto groupSize2() const					{ return m_bGroupSize2; }
	CC bool initRowSolution(uint **ppRowSolutionIdx) const {
		*ppRowSolutionIdx = new uint[m_lenDayResults * (m_bGroupSize2 ? 1 : 2)];
		(*ppRowSolutionIdx)[numPreconstructedRows()] = 0;
		return sysParam()->val[t_useCombinedSolutions];
	}
	CC void releaseSolMaskInfo() {
		delete[] m_pRowSolutionMasks;
		delete[] m_pRowSolutionMasksIdx;
		delete[] m_pMaskTestingCompleted;
	}

	CC void getMatrix(tchar* row, tchar* neighbors, int nRows, uint* pRowSolutionIdx) const;
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
	CC uint& solInterval(uint* pRowSolutionIdx, int iRow, uint* pLast, ll availablePlayers) const;
	CC bool generateCompatibilityMasks(tmask* pMask, uint solIdx, uint idx, ll* pUsedPlayers = NULL) const;
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
	CC void initMaskStorage(uint numObjects);
	CC bool p1fCheck2(ctchar* neighborsi, ctchar* neighborsj) const;
	CC bool checkCompatibility(ctchar* neighborsi, const ll* rm, uint idx) const;
	CC uint& solutionInterval2(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const;
	CC uint& solutionInterval3(uint* pRowSolutionIdx, uint* pLast, ll availablePlayers) const;
	CC void updateMasksByAut(const CGroupInfo* pGroupInfo) const;
	CC void updateMasksByAutForSolution(ctchar* pSolution, const CGroupInfo* pGroupInfo, tmask* pMask, uint solIdx, uint last, uint idxMin = 0) const;
	CC uint getSolutionRange(uint& last, ll& availablePlayers, int i) const;
	CC inline unsigned long minPlayer(ll availablePlayers) const {
#if USE_64_BIT_MASK
		unsigned long iBit;
		_BitScanForward64(&iBit, availablePlayers);
		return iBit;
#else
		#pragma message("A GPU-equivalent function similar to `_BitScanForward64` needs to be implemented.")
		return 0;
#endif
	}
	CC uint getTransformerSolIndex(ctchar* pSol, ctchar* pPerm, uint last, uint first = 0) const;

	const kSysParam* m_pSysParam;
	const int m_numPlayers;
	const int m_numPreconstructedRows;     // Number of preconstructed matrix rows
	const int m_numDaysResult;
	const alldata* m_pAllData;
	const bool m_bGroupSize2;
	const bool m_bUseCombinedSolutions;
	const int m_step;

	ctchar* m_pFirstMatr = NULL;                         // Pointer to the array of initial matrices
	CStorage<tchar>* m_pMaskStorage = NULL;
	CRepository<uint>* m_pTRTSN_Storage = NULL;   // Storage for T(hird) R(ow) T(ransformed) S(olution) N(numbers)
	                                              // Each solution number is stored alongside its corresponding permutation number
#if 0
	ll cnt1 = 0;
	ll cnt2[18] = { 0 };
	ll cnt3[18] = { 0 };
	ll cnt4[18] = { 0 };
	ll cnt5[18] = { 0 };
#endif
	uint m_numObjects;
	uint m_numObjectsMax;
	int m_lenDayResults;

	int m_lenMask;
	rowToBitmask m_fRowToBitmask;
	solutionInterval m_fSolutionInterval;

	uint* m_pPlayerSolutionCntr = NULL;
	uint m_numSolutionTotal;
	uint m_numSolutionTotalB;
	uint m_lenSolutionMask;
	tmask* m_pRowsCompatMasks[2];
	// For each row of the matrix, we define two masks, each containing an interval of consecutive bits set to 1
	// These intervals represent the row's first and last sets of solutions that lie outside the separately tested 64-bit intervals.
	tmask* m_pRowSolutionMasks = NULL;
	//  ... and the the set of indices of the long long elements which corresponds to two mask's sets.                           
	uint* m_pRowSolutionMasksIdx = NULL;
	bool* m_pMaskTestingCompleted = NULL;
	uint *m_pNumLongs2Skip = NULL; // Pointer to the number of long long's that we don't need to copy for each row.
	int m_useCliquesAfterRow;
	uint m_numRecAdj = 0;
	uint m_numRec[2];			   // Number of solutions for first and second nonfixed rows.
	uint m_numRecAdj2;
	uint m_lastInFirstSet;
	int m_solAdj = 0;
	ll m_playersMask[2] = { 0, 0 };// Mask with bits corresponding to players from first group of predefined rows equal to zeros.
	tchar* m_pSolMemory = NULL;    // Memory allocated to support the use of the automorphism group of the matrix with with pre-constructed rows.
	bool m_bUseAut;
};

class CRowUsage : public CompSolStorage {
public:
	CC CRowUsage(const CRowStorage* const pRowStorage) : CompSolStorage(pRowStorage) {
		m_bUseCombinedSolutions = pRowStorage->initRowSolution(&m_pRowSolutionIdx);
		m_bGroupSize2 = pRowStorage->groupSize2();
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
	uint m_lenMask = 0;				// Count of tmask elements, each encoding mutually compatible solutions
	uint* m_pRowSolutionIdx = NULL;
	tmask* m_pCompatibleSolutions = NULL;
	int m_step = 0;
	bool m_bUseCombinedSolutions;
	bool m_bSolutionReady = false;   // true, when solution was prepared ar a part of combined solution. 
	bool m_bUsePlayersMask = false;
	bool m_bGroupSize2 = false;
};

#define PERMUTATION_OF_PLAYERS(numPlayers, pLayersIn, permut, pLayersOut)	for (auto j = numPlayers; j--;) \
																				pLayersOut[j] = permut[pLayersIn[j]];
