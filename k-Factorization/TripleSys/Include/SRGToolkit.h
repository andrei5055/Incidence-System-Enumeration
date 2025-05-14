#pragma once

typedef struct {
	unsigned int m_cntr[3];
	int k;
	int λ;
	int μ;
} SRGParam;

class SRGToolkit
{
public:
	SRGToolkit(int nCols, int nRows, int groupSize);
	~SRGToolkit()			{ delete[] m_pAdjacencyMatrix; }
	void exploreMatrix(const unsigned char* pMatr);
	void printStat();
	void SetSRG_params(int k, int λ, int μ) {
		m_graphParam[0].k = k;
		m_graphParam[0].λ = λ;
		m_graphParam[0].μ = μ;
		m_graphParam[1].k = m_v - k - 1;
		m_graphParam[1].λ = m_v - 2 * k + μ - 2;
		m_graphParam[1].μ = m_v - 2 * k + λ;
	}
private:
	void exploreMatrixOfType(int typeIdx, const unsigned char* pMatr);
	const int m_nCols;
	const int m_nRows; 
	const int m_groupSize;
	const int m_v;

	SRGParam m_graphParam[2];
	unsigned char* m_pAdjacencyMatrix = nullptr;
};

