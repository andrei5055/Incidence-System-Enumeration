#pragma once
#include "GraphCanonizer.h"

#define OUT_SRG_TO_SEPARATE_FILE	0	// Set this value to 1 if you want to see the constructed SRG in a separate file.

typedef struct SRGParam {
	unsigned int m_cntr[5];   // 0 - total, 1 - regular, 2 - SRG, 3 - 4-vert cond; 4 - rank 3 
	int k;
#ifndef USE_CUDA
	int λ;
	int μ;
	int α;
	int β;
#endif
	SRGParam() { memset(this, 0, sizeof(*this)); }
	t_graphType updateParam(int* pCommon, bool flag_4_ver);
} SRGParam;

class SRGToolkit : public CGraphCanonizer
{
public:
	SRGToolkit(const kSysParam* pParam, int nRows, const std::string& resFileName, int exploreMatrices);
	~SRGToolkit();
	bool exploreMatrix(ctchar* pMatr, GraphDB *ppGraphDB, uint sourceMatrID, uint srcGroupOrder);
	void printStat();
private:
	bool exploreMatrixOfType(int typeIdx, ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID, uint srcGroupOrder);
	t_graphType checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	t_graphType checkSRG(const tchar *pGraph, int graphDegree, int* nCommon, size_t lenCommon, bool& flag) const;
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

