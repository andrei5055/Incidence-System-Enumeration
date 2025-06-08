#pragma once
#include "TripleSys.h"

typedef struct SRGParam {
	unsigned int m_cntr[4];
	int k;
	int λ;
	int μ;
	int α;
	int β;
	SRGParam() { memset(this, 0, sizeof(*this)); }
};

typedef enum {
	t_nonregular	= 0,
	t_regular		= 1 << 0,
	t_srg			= t_regular + (1 << 1),
	t_4_vert		= t_srg + (1 << 2),
	t_complete		= -1
} t_graphType;

class SRGToolkit : public CGroupOrder<ushort>
{
public:
	SRGToolkit(int nCols, int nRows, int groupSize);
	~SRGToolkit();
	bool exploreMatrix(ctchar* pMatr);
	void printStat();
private:
	bool exploreMatrixOfType(int typeIdx, ctchar* pMatr);
	t_graphType checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	void initCanonizer();
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVertex = 0);
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx, 
		ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex);
	ushort* restoreParam(int& i, int iStart, ushort* pLenOrbits);
	void printAdjMatrix(ctchar* pGraph, tchar* pGraphOut = NULL, int idx = 0, int startVertex = 0, int endVertex = 0) const;
	tchar* createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex = 0, int endVertex = 0) const;
	void initVertexGroupOrbits();

	const int m_nCols;
	const int m_nRows; 
	const int m_groupSize;
	const int m_v;
	
	int m_len;
	bool m_bChekMatr[2];
	SRGParam *m_pGraphParam[2] = { nullptr };
	tchar* m_pGraph[2] = { nullptr };
	ushort* m_pOrbits = nullptr;
	ushort* m_pSavedOrbits = nullptr;
	ushort* m_pGroupOrbits = nullptr;
	ushort* m_subgraphVertex = nullptr;
	ushort* m_pNumOrbits = nullptr;     // Number of orbits for each vertex
	ushort* m_pLenOrbits = nullptr;
	ushort* m_pSavedOrbIdx = nullptr;
	CBinaryMatrixStorage* m_pMarixStorage[2] = { nullptr };
};

