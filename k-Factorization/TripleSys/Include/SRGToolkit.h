#pragma once
#include "k-SysSupport.h"
#include "GroupOrder.h"

typedef struct {
	unsigned int m_cntr[4];
	int k;
	int λ;
	int μ;
	int α;
	int β;
} SRGParam;

typedef unsigned short ushort;

class SRGToolkit : public CGroupOrder<ushort>
{
public:
	SRGToolkit(int nCols, int nRows, int groupSize);
	~SRGToolkit();
	void exploreMatrix(ctchar* pMatr);
	void printStat();
private:
	void exploreMatrixOfType(int typeIdx, ctchar* pMatr);
	bool checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	void initCanonizer();
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVertex = 0);
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx, 
		ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex);
	ushort* restoreParam(int& i, int iStart, ushort* pLenOrbits);
	void printAdjMatrix(ctchar* pGraph, tchar* pGraphOut, int idx = 0, int startVertex = 0, int endVertex = 0) const;
	tchar* createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex = 0, int endVertex = 0) const;
	void initVertexGroupOrbits();

	const int m_nCols;
	const int m_nRows; 
	const int m_groupSize;
	const int m_v;
	
	int m_len;

	SRGParam m_graphParam[2];
	tchar* m_pGraph[2] = { nullptr };
	ushort* m_pOrbits = nullptr;
	ushort* m_pSavedOrbits = nullptr;
	ushort* m_pGroupOrbits = nullptr;
	ushort* m_subgraphVertex = nullptr;
	ushort* m_pNumOrbits = nullptr;     // Number of orbits for each vertex
	ushort* m_pLenOrbits = nullptr;
	ushort* m_pSavedOrbIdx = nullptr;
};

