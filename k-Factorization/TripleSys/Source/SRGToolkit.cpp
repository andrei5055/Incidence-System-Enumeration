#include "SRGToolkit.h"
#include <cstring>
#include "k-SysSupport.h"

typedef unsigned char tchar;
typedef const tchar ctchar;

SRGToolkit::SRGToolkit(int nCols, int nRows, int groupSize) : 
	m_nCols(nCols), m_nRows(nRows), m_groupSize(groupSize), m_v(nRows * nCols/groupSize) {
	m_pAdjacencyMatrix = new tchar[m_v * m_v];
	memset(m_graphParam, 0, sizeof(m_graphParam));
}

static bool one_common_element(ctchar* pArray1, ctchar* pArray2, int len) {
	int i = 0, j = 0;
	while (i < len && j < len) {
		if (pArray1[i] == pArray2[j])
			return true;

		if (pArray1[i] < pArray2[j])
			++i;
		else
			++j;
	}
	return false;
}

void SRGToolkit::exploreMatrix(const unsigned char* pMatr) {
	for (int i = 0; i < 2; i++)
		exploreMatrixOfType(i, pMatr);
}

void SRGToolkit::exploreMatrixOfType(int typeIdx, const unsigned char* pMatr) {
	const auto numGroups = m_nCols / m_groupSize;
	const auto pVertexLast = pMatr + m_nRows * m_nCols;

	auto& graphParam = m_graphParam[typeIdx];
	const auto cond = typeIdx != 0;
	// Create graph
	memset(m_pAdjacencyMatrix, 0, m_v * m_v * sizeof(m_pAdjacencyMatrix[0]));
	auto pVertex = pMatr;
	int numVertex = 0;
	for (int i = 0; i < m_nRows; i++) {
		for (int j = 0; j < numGroups; j++, pVertex += m_groupSize, numVertex++) {
			// Add edges to the graph for the vertices in the same matrix row
			auto numNextVertex = numVertex;
			for (int k = j; ++k < numGroups;) {
				++numNextVertex;
				m_pAdjacencyMatrix[numVertex * m_v + numNextVertex] =
					m_pAdjacencyMatrix[numNextVertex * m_v + numVertex] = 1;
			}

			// Add edges to the graph for the vertices intersecting by one element
			auto* pNextVertex = pMatr + (i + 1) * m_nCols;
			while (pNextVertex < pVertexLast) {
				numNextVertex++;
				// We do have two pointers (pVertex and pNextVertex).
				// If one out of groupSize elements is the same, then we have an edge between two vertices.
				if (one_common_element(pVertex, pNextVertex, m_groupSize) == cond) {
					// Add edges to the graph for the vertices intersecting by one element
					m_pAdjacencyMatrix[numVertex * m_v + numNextVertex] =
						m_pAdjacencyMatrix[numNextVertex * m_v + numVertex] = 1;
				}

				// Move to the next vertex
				pNextVertex += m_groupSize;
			}
		}
	}

	graphParam.m_cntr[0]++;

	// Check if constructed graph is regular
	int graphDegree = 0;
	pVertex = m_pAdjacencyMatrix;
	for (int i = 0; i < m_v; i++, pVertex += m_v) {
		int vertexDegree = 0;
		for (int j = 0; j < m_v; j++)
			if (pVertex[j])
				vertexDegree++;

		if (graphDegree) {
			if (graphDegree != vertexDegree) {
				printfRed("Graph is not regular\n");
				return;
			}
		}
		else
			graphDegree = vertexDegree;
	}

	if (!graphParam.m_cntr[1]++)
		graphParam.k = graphDegree;

	// Check if constructed graph is strongly regular
	int nCommon[2] = { 0, 0 };
	auto pFirstVertex = m_pAdjacencyMatrix;
	for (int i = 0; i < m_v; i++, pFirstVertex += m_v) {
		auto pSecondVertex = pFirstVertex;
		for (int j = i; ++j < m_v;) {
			pSecondVertex += m_v;
			const auto idx = pFirstVertex[j] ? 0 : 1;
			int nCommonCurr[2] = { 0, 0 };
			for (int k = 0; k < m_v; k++) {
				if (pFirstVertex[k] && pSecondVertex[k])
					nCommonCurr[idx]++;
			}
			if (nCommon[idx]) {
				if (nCommon[idx] != nCommonCurr[idx]) {
					printfRed("Graph is not strongly regular\n");
					return;
				}
			}
			else
				nCommon[idx] = nCommonCurr[idx];
		}
	}

	if (!graphParam.m_cntr[2]++) {
		graphParam.λ = nCommon[0];
		graphParam.μ = nCommon[1];
	}
}

void SRGToolkit::printStat() {
	for (int i = 0; i < 2; i++) {
		auto& graphParam = m_graphParam[i];
		const bool plural = graphParam.m_cntr[0] > 1;
		const char* pntr0 = plural ? "s" : "";
		const char* pntr1 = plural ? "are" : "is";
		const char* pntr2 = plural ? "  " : "";

		printfYellow("Out of %d graph%s with %d vertices\n", graphParam.m_cntr[0], pntr0, m_v);
		printfYellow("       %d %s regular of degree %d\n", graphParam.m_cntr[1], pntr1, graphParam.k);

		printfYellow("Graph%s %s strongly regular with parameters: (%d,%2d,%2d,%2d)\n", pntr0, pntr1, m_v, graphParam.k, graphParam.λ, graphParam.μ);
		const auto v_2k = m_v - 2 * graphParam.k;
		graphParam.k = m_v - graphParam.k - 1;
		const auto λ = graphParam.λ;
		graphParam.λ = v_2k + graphParam.μ - 2;
		graphParam.μ = v_2k + λ;
		printfYellow("%s        Parameters of complimentary graph: (%d,%2d,%2d,%2d)\n", pntr2, m_v, graphParam.k, graphParam.λ, graphParam.μ);
	}
}
