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
	SRGToolkit(const kSysParam* pParam, int nRows, const std::string& resFileName, int exploreMatrices, SRGToolkit* pMaster = NULL);
	~SRGToolkit();
	bool exploreMatrix(ctchar* pMatr, uint sourceMatrID, CBinaryMatrixStorage** ppMarixStorage);
	void printStat();
	inline SRGParam* graphParam(int i) const	{ return m_pGraphParam[i]; }
	inline void reportOnScreen(bool val)		{ m_reportOnScreen = val; }
	inline bool reportOnScreen() const			{ return m_reportOnScreen; }
	inline auto srcGroupOrderPntr()		        { return &m_srcGroupOrder; }
	inline void setMaster(SRGToolkit* pMaster)	{ m_pMaster = pMaster; }
	bool outputGraph(int typeIdx, t_graphType graphType, uint sourceMatrID, CBinaryMatrixStorage* pMarixStorage, 
		bool rank3, ctchar* pResGraph, ctchar* pUpperDiag, SRGToolkit* pSlaveToolKit);
private:
	bool exploreMatrixOfType(int typeIdx, ctchar* pMatr, uint sourceMatrID, CBinaryMatrixStorage* pMarixStorage);
	t_graphType checkSRG(tchar* pGraph, SRGParam* pGraphParam = nullptr);
	t_graphType checkSRG(const tchar *pGraph, int graphDegree, int* nCommon, size_t lenCommon, bool& flag) const;
	inline int param(paramID id) const { return m_pParam->val[id]; }
	inline auto getMaster() const				{ return m_pMaster; }
	void outputGraph(int typeIdx, uint prevMatrNumb, t_graphType graphType, bool rank3, ctchar *pResGraph, SRGToolkit* pSlaveToolKit);
	
	const int m_nRows; 
	const std::string m_resFileName;
	const int m_nExploreMatrices;

	uint m_srcGroupOrder;		// The order of the group of the source matrix (the matrix from which the SRG is constructed).
	int m_nPrevMatrNumb = 0;
	bool m_bChekMatr[2];
	ushort* m_subgraphVertex = nullptr;
	SRGParam *m_pGraphParam[2] = { nullptr };
	const kSysParam* m_pParam;
	bool m_reportOnScreen = false;
	SRGToolkit *m_pMaster = nullptr;
};

