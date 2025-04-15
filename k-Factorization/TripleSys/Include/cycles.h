#pragma once

class Cycles {
public:
	bool defined = false;
	tchar rowInd[2];
	tchar cycles[MAX_3PF_SETS * MAX_CYCLES_PER_SET];
	int nCycleSets = 0;
	tchar* cycleStart(int ic) {
		return (ic >= nCycleSets) ? NULL : cycles + ic * MAX_CYCLES_PER_SET + 1;
	}
};

class CycleMapping {
protected:
	CC CycleMapping(int numPlayers, bool bBipartiteGraph = false) : m_bBipartite(bBipartiteGraph) {
		m_pNumCycles = new tchar[3 * numPlayers];
	}
	CC ~CycleMapping()				{ delete [] m_pNumCycles; }
	CC ctchar* InitCycleMapping(ctchar* const pLenCycles, ctchar* pStartCycles, int nCycles, int lenGroup, ctchar** const pDirection, ctchar** const pStartCycleOut);
	CC bool ProceedToNextMapping();
	CC auto permCycles() const		{ return m_permCycles; }
	CC bool nextPerm(tchar* pPerm, int lenPerm) const;
private:
	const bool m_bBipartite;
	tchar* m_pNumCycles = NULL;		// number of cycles in each group of cycles
	tchar* m_permCycles;			// permutations of cycles of the same length
	tchar* m_pIdxElem;				// indices of current elements in each cycle
	tchar* m_pDirection;
	tchar* m_pLenCycles;
	tchar* m_pStartCyclesIn;
	tchar* m_pStartCyclesOut;
	int m_nLenGroup;
	int m_nCycles;
	int m_numCycleGroups;
	size_t m_totalVariants;
	size_t m_currVariant;

};

class CycleSupport : protected CycleMapping {
protected:
	CC CycleSupport(int numPlayers, bool bBipartiteGraph = false) : CycleMapping(numPlayers, bBipartiteGraph) {}
	CC ~CycleSupport() { delete[] m_pV0; }
    CC void InitCycleSupport(int groupNumber) {
		auto const len = groupNumber * MAX_3PF_SETS;
		m_pV0 = new tchar[2 * len];
		m_pV1 = m_pV0 + len;
	}
	CC inline auto getV0() const { return m_pV0; }
	CC inline auto getV1() const { return m_pV1; }

	//CC void combineCycles(ctchar* pLenCycles, int nCycles);
private:
	tchar* m_pV0 = NULL;
	tchar* m_pV1 = NULL;
};