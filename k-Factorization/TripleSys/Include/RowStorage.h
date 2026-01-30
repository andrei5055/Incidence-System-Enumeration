#pragma once
#include "CompatMasks.h"

class CRowStorage : public CStorage<tchar>, public CCompatMasks {
	typedef void (CRowStorage::* rowToBitmask)(ctchar* pRow, tmask *pMask) const;
public:
	CC CRowStorage(const kSysParam* pSysParam, int numPlayers, int numObjects = 1000, const alldata* pAllData = NULL);
	CC ~CRowStorage();
	CC inline void init()								{ initMaskStorage(m_numObjectsMax); }
	CC void initPlayerMask(ctchar* pFirstMatr = NULL, ctchar lastNeighborOfPlayer0 = 0);
	CC bool maskForCombinedSolutions(tmask* pMaskOut, uint& solIdx) const;
	CC void passCompatibilityMask(tmask *pCompatibleSolutions, uint first, uint last, const alldata* pAllData, CCompatMasks** pCompMaskHandle = NULL) const;
	CC bool checkSolutionByMask(int iRow, const tmask* pToASol) const;
	CC inline auto useFeature(featureFlags mask) const	{ return m_pSysParam->useFeature(mask); }
	CC bool addRow(ctchar* pRow, ctchar* pNeighbors, ctchar* pNeighbors2);
	CC int initCompatibilityMasks(CStorageIdx<tchar>** ppSolRecast = NULL);
	CC inline auto useCliquesAfterRow() const			{ return m_useCliquesAfterRow; }
	CC inline auto useCliques(int iRow) const			{ return iRow > m_useCliquesAfterRow; }
	CC inline const kSysParam* sysParam() const			{ return m_pSysParam; }
	CC inline const auto numRec(int idx) const			{ return m_numRec[1]; }
	CC inline auto groupSize2() const					{ return m_bGroupSize2; }
	CC bool initRowSolution(uint **ppRowSolutionIdx) const {
		*ppRowSolutionIdx = new uint[lenDayResults() * (selectPlayerByMask() ? 2 : 1)];
		(*ppRowSolutionIdx)[numPreconstructedRows()] = 0;
		return sysParam()->val[t_useCombinedSolutions];
	}
	const alldata* allData() const						{ return m_pAllData; }

	CC void getMatrix(tchar* row, tchar* neighbors, int nRows, uint* pRowSolutionIdx, const CCompatMasks* pCompMask) const;
	CC void outSelectedSolution(int iRow, uint first, uint last, int threadID = 0) const;
#if !USE_64_BIT_MASK
	CC inline auto firstOnePosition(tchar byte) const	{ return m_FirstOnePosition[byte]; }
private:
	tchar m_FirstOnePosition[256]; // Table for fast determination of the first 1's position in byte.
#endif

#define SetMask(bm, ir, ic) \
	const int ib = ir * numPlayers() + ic - (1 + ir) * (2 + ir) / 2; \
	ASSERT_IF(ib / 8 >= m_lenMask); \
	SET_MASK_BIT(bm, ib)

#define SetMaskWithOffset(bm, iOffset, ic) \
	ib = iOffset + ic; \
	ASSERT_IF(ib / 8 >= m_lenMask); \
	SET_MASK_BIT(bm, ib)

private:
	CC uint& solInterval(uint* pRowSolutionIdx, int iRow, uint* pLast, ll availablePlayers) const;
	CC bool generateCompatibilityMasks(tmask* pMask, uint solIdx, uint idx, const alldata* pAllData = NULL, ll* pUsedPlayers = NULL) const;
	CC void rowToBitmask2(ctchar* pRow, tmask* bm) const
	{
		memset(bm, 0, m_lenMask);
		for (int i = 0; i < numPlayers(); i += 2) {
			SetMask(bm, pRow[i], pRow[i + 1]);
		}
	}
	CC void rowToBitmask3(ctchar* pRow, tmask* bm) const
	{
		int ib;
		memset(bm, 0, m_lenMask);
		for (int i = 0; i < numPlayers(); i += 3, pRow += 3) {
			const int iOffset = pRow[0] * numPlayers() - (1 + pRow[0]) * (2 + pRow[0]) / 2;
			SetMaskWithOffset(bm, iOffset, pRow[1]);
			SetMaskWithOffset(bm, iOffset, pRow[2]);
			SetMask(bm, pRow[1], pRow[2]);
		}
	}
	CC void initMaskStorage(uint numObjects);
	CC bool p1fCheck2P1F(ctchar* neighborsi, ctchar* neighborsj) const;
	CC bool checkCompatibility(ctchar* neighborsi, const ll* rm, uint idx, const alldata* pAllData, int& sameP1) const;
	CC void updateMasksByAut(const CGroupInfo* pGroupInfo) const;
	CC void updateMasksByAutForSolution(ctchar* pSolution, const CGroupInfo* pGroupInfo, tmask* pMask, uint solIdx, uint last, uint idxMin = 0) const;
	CC uint getTransformerSolIndex(ctchar* pSol, ctchar* pPerm, uint last, uint first = 0) const;
	CC void modifyMask(CStorageIdx<tchar>** ppSolRecast);
	CC int findIndexInRange(int left, int right, ctchar* pSol) const;

	const kSysParam* m_pSysParam;
	const alldata* m_pAllData;
	const bool m_bGroupSize2;
	const bool m_bUseCombinedSolutions;
	const int m_step;
	const bool m_use3RowCheck;
	int m_stepCombSolution;

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
	uint m_numObjectsMax;

	int m_lenMask;
	rowToBitmask m_fRowToBitmask;

	int m_useCliquesAfterRow;
	uint m_numRec[2];			   // Number of solutions for first and second nonfixed rows.
	uint m_numRecAdj2 = 0;         // used only in combined solutions mode (it seems, that it never been properly initialized)
	tchar* m_pSolMemory = NULL;    // Memory allocated to support the use of the automorphism group of the matrix with with pre-constructed rows.
	bool m_bUseAut;
};

class CRowUsage : public CompSolStorage {
public:
	CC CRowUsage(const CRowStorage* const pRowStorage) : CompSolStorage(pRowStorage) {
		m_bUseCombinedSolutions = pRowStorage->initRowSolution(&m_pRowSolutionIdx);
		m_bSelectPlayerByMask = pRowStorage->selectPlayerByMask();
		m_pCompatMasks = pRowStorage;
	}
	CC ~CRowUsage() {
		delete[] m_pRowSolutionIdx;
		delete[] m_pCompatibleSolutions;
		if (m_pCompatMasks != m_pRowStorage)
			delete m_pCompatMasks;
	}
	CC void init(int iThread = 0, int numThreads = 1);
	CC int getRow(int iRow, int ipx, const alldata* pAllData);
	CC inline bool getMatrix2(tchar* row, tchar* neighbors, int nRows, int iRow) {
		getMatrix(row, neighbors, ++iRow);
		return completeMatrix(row, neighbors, nRows, iRow);
	}
	CC inline void getMatrix(tchar* row, tchar* neighbors, int nRows) {
		m_pRowStorage->getMatrix(row, neighbors, nRows, m_pRowSolutionIdx, m_pCompatMasks);
	}
private:
	CC inline auto selectPlayerByMask() const				{ return m_bSelectPlayerByMask; }
	uint* m_pRowSolutionIdx = NULL;
	tmask* m_pCompatibleSolutions = NULL;
	int m_step = 0;
	bool m_bUseCombinedSolutions;
	bool m_bSolutionReady = false;   // true, when solution was prepared ar a part of combined solution. 
	bool m_bSelectPlayerByMask = false;
	bool m_bPlayerByMask = false;
	int m_threadID = 0;
	const CCompatMasks *m_pCompatMasks = NULL;
};

#define PERMUTATION_OF_PLAYERS(numPlayers, pLayersIn, permut, pLayersOut)	for (auto j = numPlayers; j--;) \
																				pLayersOut[j] = permut[pLayersIn[j]];
