#pragma once

class CycleMapping {
protected:
	CC CycleMapping(int numPlayers, bool cmpGraph = false) : m_bCMP_Graph(cmpGraph) {
		m_pNumCycles = new tchar[3 * numPlayers];
	}
	CC ~CycleMapping()					{ delete [] m_pNumCycles; }
	CC ctchar* InitCycleMapping(ctchar* const pLenCycles, ctchar* pStartCycles, int nCycles, int lenGroup, ctchar** const pDirection, ctchar** const pStartCycleOut);
	CC bool ProceedToNextMapping();
	CC inline auto permCycles() const	{ return m_permCycles; }
	CC bool nextPerm(tchar* pPerm, int lenPerm) const;
private:
	const bool m_bCMP_Graph;
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
	CC CycleSupport(int numPlayers, bool cmpGraph = false) : CycleMapping(numPlayers, cmpGraph) {}
	CC ~CycleSupport() { delete[] m_pV0; }
    CC void InitCycleSupport(int groupNumber, int maxSets) {
		auto const len = groupNumber * maxSets;
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