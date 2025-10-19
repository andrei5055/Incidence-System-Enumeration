#pragma once
#include "Table.h"

#define DEBUGING_CANON		0
#define PRINT_MATRICES		0
#define PRINT_NUM_CUR_GRAPH TRACE_GROUP_ORDER

#if PRINT_MATRICES || PRINT_NUM_CUR_GRAPH
extern int numCurrGraph;
#endif

#if PRINT_MATRICES || DEBUGING_CANON
#define PRINT_ADJ_MATRIX(...) printAdjMatrix(__VA_ARGS__)
#else
#define PRINT_ADJ_MATRIX(...)
#endif

#if PRINT_MATRICES
#define DO_PRINT(nIter)   (printFlag && FFF != -1 && nIter >= FFF)
// Parameters specific to the bug we are trying to fix.
#define N_MATR 3     // Number of matrix to activate the output
#define NUM_GENERATOR	1  // Generator number that is presumably missing from the matrix with the smaller group
#define FFF 22  // should be set equal to index NM of last bbb_NM.txt file
#define FF_ 1 // 2   // 4 - for 26. 2 for 22
#endif

class CGraphCanonizer : public Generators<ushort>
{
public:
	CGraphCanonizer(int nVer = 0);
	~CGraphCanonizer() { releaseCanonizerMemory(); }

	ctchar* canonize_graph(ctchar* pGraph = NULL, int* pCanonIdx = NULL);
	inline auto* graphPntr(int idx = 0) const { ASSERT_IF(idx < 0 || idx > 1); return m_pGraph[idx]; }
protected:
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVert = 0);
	inline auto lenGraphMatr() const { return m_lenGraphMatr; }
	inline auto* groupOrbits() const { return m_pGroupOrbits; }
	void outAdjMatrix(ctchar* pGraphOut, FILE* f, int endVertex = 0) const;
	void printAdjMatrix(ctchar* pGraphOut, int idx = 0, int endVertex = 0, ushort* pOrb = NULL, const char* fName = "aaa", cchar *mode = "w", cchar* frmt = "%s_%02d.txt") const;
	void printAdjMatrix(ctchar* pGraphOut, const char* fileName, int endVertex, int idx = 0, ushort* pOrb = NULL, cchar * mode = "w") const;
private:
	int Init(int nVer);
	void initCanonizer();
	void releaseCanonizerMemory();
	void initVertexGroupOrbits();
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx,
		ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex);
	tchar* createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex = 0, int endVertex = 0, const ushort* pOrb = NULL) const;
	ushort* restoreParam(int& i, int iStart, ushort* pLenOrbits);
	inline auto* vertOrbits() const { return m_pOrbits; }

	int m_v;
	size_t m_lenGraphMatr;
	tchar* m_pGraph[2] = { nullptr };
	ushort* m_pOrbits = nullptr;
	ushort* m_pNumOrbits = nullptr;     // Number of orbits for each vertex
	ushort* m_pGroupOrbits = nullptr;
	tchar* m_bUsedFlags = nullptr;
	ushort* m_pLenOrbits = nullptr;
	ushort* m_pSavedOrbits = nullptr;
	ushort* m_pSavedOrbIdx = nullptr;
};


