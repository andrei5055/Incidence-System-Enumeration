#pragma once

#if !USE_CUDA
class CompSol {
public:
	CC CompSol()					{ m_ppCompatibleSol = NULL; }
	CC void Init(int numPrevGroups) {
		if (numPrevGroups) {
			m_ppCompatibleSol = new std::vector<uint>*[numPrevGroups];
			for (int i = numPrevGroups; i--;)
				m_ppCompatibleSol[i] = new std::vector<uint>;
		}
	}
	CompSol* Initialize(uint idx, int rowIdx) {
		initialized() ? Reset(rowIdx) : Init(rowIdx);
		setSolIndex(idx);
		return this;
	}
	CC void Release(int numPrevGroups) {
		for (int i = numPrevGroups; i--;)
			delete m_ppCompatibleSol[i];
	}
	CC void Reset(int numPrevGroups) {
		setUnused();
		for (int i = numPrevGroups; i--;)
			m_ppCompatibleSol[i]->clear();
	}
	CC inline bool initialized() const				{ return m_ppCompatibleSol != NULL; }
	CC void setSolIndex(uint val)					{ m_solIdx = val; }
	inline auto solIdx() const						{ return m_solIdx; }
	inline void addNeighbor(uint idx, int lev)		{ m_ppCompatibleSol[lev]->push_back(idx); }
	inline bool noNeighborsOnLevel(int lev) const	{ return m_ppCompatibleSol[lev]->empty(); }
	inline auto compatibleSolutions(int lev) const	{ return m_ppCompatibleSol[lev]; }
	inline void setUsed()							{ m_bUsed = true; }
	inline void setUnused()							{ m_bUsed = false; }
	inline auto isUsed() const						{ return m_bUsed; }
private:
	uint m_solIdx;
	bool m_bUsed = false;					// true, if the solution is compatible with at least one solution of next level. 
	std::vector<uint>** m_ppCompatibleSol;  // Indices of solutions from previous groups, that are
											// compatible with the solution indexed by m_solIdx.
};

class CompSolSet : public CStorageSet<CompSol> {
public:
	CC CompSolSet(int numObjects) : CStorageSet<CompSol>(numObjects) {}
private:
};

class CRowStorage;

class CompSolStorage {
public:
	CC CompSolStorage(int nGroups, int lenGroup) : m_nGroups(nGroups) {
		m_ppCompSol = new CompSolSet* [nGroups];
		m_ppCompSol[0] = new CompSolSet(lenGroup);
		for (int i = 1; i < nGroups; i++) {
			m_ppCompSol[i] = new CompSolSet(lenGroup);
			for (int j = 0; j < lenGroup; j++)
				m_ppCompSol[i]->getObject(j)->Init(i);
		}
		m_solDB = new std::vector<CompSol *>[nGroups];
	}
	CC ~CompSolStorage() {
		for (int i = 0; i < m_nGroups; i++) {
			for (int j = i; j--;)
				m_ppCompSol[i]->getObject(j)->Release(i);

			delete m_ppCompSol[i];
		}

		delete[] m_ppCompSol;
		delete[] m_solDB;
		delete[] getBuffer();

	}
	CC void allocateBuffer(size_t nLongs)		{ m_pBuffer = new long long[nLongs]; }
	CC auto getBuffer() const					{ return m_pBuffer; }
	CC void addCompatibleSolutions(uint jBase, tmask mask, int kMax, const CRowStorage* pRowStorage, long long* pSolMask);
	CC void releaseSolDB() {
		for (int i = m_nGroups; i--;)
			m_solDB[i].clear();
	}
	void removeUnreachableVertices(int rowIdx);
private:
	CC CompSol* compatibleSolutions(uint idx, int rowIdx) const {
		return m_ppCompSol[rowIdx]->getNextObject()->Initialize(idx, rowIdx);
	}
	CC inline void releaseCompatibleSolutions(int rowIdx) {
		m_ppCompSol[rowIdx]->releaseObject();
	}

	CompSolSet** m_ppCompSol;
	const int m_nGroups;
	std::vector<CompSol *> *m_solDB;
	long long* m_pBuffer = NULL;
};
#else
// Dummy class just to make program compiled for GPU
class CompSolStorage {
public: 
	CC void allocateBuffer(size_t nLongs) {}
};
#endif