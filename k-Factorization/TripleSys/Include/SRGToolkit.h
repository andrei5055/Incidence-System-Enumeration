#pragma once
#include "k-SysSupport.h"

typedef struct {
	unsigned int m_cntr[4];
	int k;
	int λ;
	int μ;
	int α;
	int β;
} SRGParam;

class SRGToolkit
{
public:
	SRGToolkit(int nCols, int nRows, int groupSize);
	~SRGToolkit();
	void exploreMatrix(ctchar* pMatr);
	void printStat();
private:
	void exploreMatrixOfType(int typeIdx, ctchar* pMatr);
	bool checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVertex = 0);
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx, 
		tchar** ppLenOrbits, int& idxRight, bool flag, int& lastUnfixedVertexIndex);
	void printAdjMatrix(ctchar* pGraph, tchar* pGraphOut, int idx = 0);

	const int m_nCols;
	const int m_nRows; 
	const int m_groupSize;
	const int m_v;

	SRGParam m_graphParam[2];
	tchar* m_pGraph[2] = { nullptr };
	unsigned short* m_pOrbits = nullptr;
	unsigned short* m_subgraphVertex = nullptr;
	unsigned short* m_pNumOrbits = nullptr;     // Number of orbits for each vertex
	tchar* m_pLenOrbits = nullptr;
};

