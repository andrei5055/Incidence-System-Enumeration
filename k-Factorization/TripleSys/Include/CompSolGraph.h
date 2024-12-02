#pragma once

class CRowStorage;

#if !USE_CUDA

class CompSol {
public:
	CC CompSol()					{ m_ppCompatibleSol = NULL; }
	CC ~CompSol()					{ delete [] m_ppCompatibleSol; }
	CC void Init(int numPrevGroups) {
		if (numPrevGroups)
			m_ppCompatibleSol = new std::vector<CompSol *>[numPrevGroups];
	}
	CompSol* Initialize(uint idx, int rowIdx) {
		initialized() ? Reset(rowIdx) : Init(rowIdx);
		setSolIndex(idx);
		return this;
	}
	CC void Reset(int numPrevGroups) {
		setUnused();
		for (int i = numPrevGroups; i--;)
			m_ppCompatibleSol[i].clear();
	}
	CC inline bool initialized() const				{ return m_ppCompatibleSol != NULL; }
	CC void setSolIndex(uint val)					{ m_solIdx = val; }
	inline auto solIdx() const						{ return m_solIdx; }
	inline void addNeighbor(CompSol *sol, int lev)	{ m_ppCompatibleSol[lev].push_back(sol); }
	inline bool noNeighborsOnLevel(int lev) const	{ return m_ppCompatibleSol[lev].empty(); }
	inline auto compatibleSolutions(int lev) const	{ return m_ppCompatibleSol[lev]; }
	inline void setUsed()							{ m_bUsed = true; }
	inline void setUnused()							{ m_bUsed = false; }
	inline auto isUsed() const						{ return m_bUsed; }
private:
	uint m_solIdx;
	bool m_bUsed = false;						// true, if the solution is compatible with at least one solution of next level. 
	std::vector<CompSol *>* m_ppCompatibleSol;	// Solutions from previous groups, that are compatible with the solution indexed by m_solIdx
};

class CompSolSet : public CStorageSet<CompSol> {
public:
	CC CompSolSet(int numObjects) : CStorageSet<CompSol>(numObjects) {}
private:
};

class CompSolStorage {
public:
	CC CompSolStorage(const CRowStorage* const pRowStorage, int lenGroup = 100);
	CC ~CompSolStorage();
	bool ConstructCompatibleSolutionGraph(long long* pToA, int iRow);
protected:
	CC bool completeMatrix(tchar* row, tchar* neighbors, int nRows, int iRow);
private:
	CC void addCompatibleSolutions(uint jBase, tmask& mask, int kMax);
	CC void releaseSolDB() {
		for (int i = m_nGroups; i--;) {
			m_solDB[i].clear();
			m_idxUsedSol[i] = -1;
			m_ppCompSol[i]->releaseAllObjects();
		}
	}
	bool removeUnreachableVertices(int rowIdx);

	CC CompSol* compatibleSolutions(uint idx, int rowIdx) const {
		return m_ppCompSol[rowIdx]->getNextObject()->Initialize(idx, rowIdx);
	}
	CC inline void releaseCompatibleSolutions(int rowIdx) {
		m_ppCompSol[rowIdx]->releaseObject();
	}
	CC inline auto currentSolution(int idx) const {
		return m_solDB[idx][m_idxUsedSol[idx]];
	}
protected:
	const CRowStorage* const m_pRowStorage;
	int m_nRowMax;				// Maximum value of iRow
private:
	int m_nGroups;
	CompSolSet** m_ppCompSol;
	std::vector<CompSol *> *m_solDB;
	int* m_idxUsedSol = NULL;  // Indices of the solutions from m_solDB used for each row.  
};
#else
// Dummy class just to make program compiled for GPU
class CompSolStorage {
public:
	CC CompSolStorage(const CRowStorage* const pRowStorage) : m_pRowStorage(pRowStorage) { 
		m_nRowMax = 0;//  pRowStorage->numPlayers() - 2;
	}
	CC bool completeMatrix(tchar* row, tchar* neighbors, int nRows, int iRow) { return false; }
protected:
	const CRowStorage* const m_pRowStorage;
	int m_nRowMax;
};
#endif