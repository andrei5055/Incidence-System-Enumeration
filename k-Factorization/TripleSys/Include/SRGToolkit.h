﻿#pragma once
#include "Table.h"

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

class CGraphCanonizer : public Generators<ushort>
{
public:
	CGraphCanonizer(int nVer = 0);
	~CGraphCanonizer()					{ releaseCanonizerMemory(); }

	ctchar* canonize_graph(ctchar* pGraph = NULL, int* pCanonIdx = NULL);
	inline auto* graphPntr(int idx = 0) const { ASSERT_IF(idx < 0 || idx > 1); return m_pGraph[idx]; }
protected:
	int canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVert = 0);
	inline auto lenGraphMatr() const			{ return m_lenGraphMatr; }
	inline auto* groupOrbits() const			{ return m_pGroupOrbits; }
private:
	int Init(int nVer);
	void initCanonizer();
	void releaseCanonizerMemory();
	void initVertexGroupOrbits();
	int canonizeMatrixRow(ctchar* pGraph, tchar* pVertOut, int vertIdx,
		ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex);
	tchar* createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex = 0, int endVertex = 0, const ushort* pOrb = NULL) const;
	ushort* restoreParam(int& i, int iStart, ushort* pLenOrbits);

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

class SRGToolkit : public CGraphCanonizer
{
public:
	SRGToolkit(ctchar* pGraph, int nVert);
	SRGToolkit(const kSysParam* pParam, int nRows, const std::string& resFileName, int exploreMatrices);
	~SRGToolkit();
	bool exploreMatrix(ctchar* pMatr, GraphDB *ppGraphDB, uint sourceMatrID, uint srcGroupOrder);
	void printStat();
private:
	bool exploreMatrixOfType(int typeIdx, ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID, uint srcGroupOrder);
	t_graphType checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	t_graphType checkSRG(const tchar *pGraph, int graphDegree, int* nCommon, size_t lenCommon, bool& flag) const;
	void printAdjMatrix(ctchar* pGraphOut, int idx = 0, int endVertex = 0, ushort* pOrb = NULL, const char* fName = "aaa") const;
	void printAdjMatrix(ctchar* pGraphOut, const char* fileName, int endVertex, int idx = 0, ushort* pOrb = NULL) const;
	void outAdjMatrix(ctchar* pGraphOut, FILE* f, int endVertex = 0) const;
	inline int param(paramID id) const { return m_pParam->val[id]; }

	const int m_nRows; 
	const std::string m_resFileName;
	const int m_nExploreMatrices;

	int m_nPrevMatrNumb = 0;
	bool m_bChekMatr[2];
	ushort* m_subgraphVertex = nullptr;
	SRGParam *m_pGraphParam[2] = { nullptr };
	CBinaryMatrixStorage* m_pMarixStorage[2] = { nullptr };
	const kSysParam* m_pParam;
};

