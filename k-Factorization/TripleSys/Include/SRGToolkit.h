#pragma once
#include "Table.h"

#define CHECK_NON_SRG				1   // Set this value to 1, if you want to see the graph which are not strongly-regular
#define OUT_SRG_TO_SEPARATE_FILE	0	// Set this value to 1 if you want to see the constructed SRG in a separate file.

typedef struct SRGParam {
	unsigned int m_cntr[5];   // 0 - total, 1 - regular, 2 - SRG, 3 - 4-vert cond; 4 - rank 3 
	int k;
	int λ;
	int μ;
	int α;
	int β;
	SRGParam() { memset(this, 0, sizeof(*this)); }
	t_graphType updateParam(int* pCommon, bool flag_4_ver);
} SRGParam;

class SRGToolkit : public  Generators<ushort>
{
public:
	SRGToolkit(int nCols, int nRows, int groupSize, const std::string& resFileName, bool semiSymmetric, int exploreMatrices);
	~SRGToolkit();
	bool exploreMatrix(ctchar* pMatr, GraphDB *ppGraphDB, uint sourceMatrID, uint srcGroupOrder);
	void printStat();
private:
	bool exploreMatrixOfType(int typeIdx, ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID, uint srcGroupOrder);
	t_graphType checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	void initCanonizer();
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVertex = 0);
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx, 
		ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex);
	ushort* restoreParam(int& i, int iStart, ushort* pLenOrbits);
	void printAdjMatrix(ctchar* pGraphOut, int idx = 0, int endVertex = 0, ushort* pOrb = NULL, const char* fName = "aaa") const;
	void printAdjMatrix(ctchar* pGraphOut, const char* fileName, int endVertex, int idx = 0, ushort* pOrb = NULL) const;
	void outAdjMatrix(ctchar* pGraphOut, FILE* f, int endVertex = 0) const;
	tchar* createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex = 0, int endVertex = 0, const ushort* pOrb = NULL) const;
	void initVertexGroupOrbits();

	const int m_nCols;
	const int m_nRows; 
	const int m_groupSize;
	const int m_v;
	const std::string m_resFileName;
	const int m_nExploreMatrices;
	int m_nPrevMatrNumb = 0;
	
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
	tchar* m_bConnectFlags = nullptr;
};

