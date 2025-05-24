#include "SRGToolkit.h"
#include "TripleSys.h"
#include <cstring>

SRGToolkit::SRGToolkit(int nCols, int nRows, int groupSize) : 
	m_nCols(nCols), m_nRows(nRows), m_groupSize(groupSize), m_v(nRows * nCols/groupSize) {
	m_pGraph[0] = new tchar[2 * m_v * m_v];
	m_pGraph[1] = m_pGraph[0] + m_v * m_v;
	memset(m_graphParam, 0, sizeof(m_graphParam));
	m_subgraphVertex = new unsigned short[m_v];
	m_pOrbits = new unsigned short[m_v * m_v];
	m_pNumOrbits = new unsigned short[m_v];
	m_pLenOrbits = new tchar[m_v * (m_v - 1) / 2];
}

SRGToolkit::~SRGToolkit() { 
	delete[] m_pGraph[0];
	delete[] m_subgraphVertex;
	delete[] m_pNumOrbits;
	delete[] m_pOrbits;
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
	auto pAdjacencyMatrix = m_pGraph[0];
	memset(pAdjacencyMatrix, 0, m_v * m_v * sizeof(pAdjacencyMatrix[0]));
	auto pVertex = pMatr;
	int numVertex = 0;
	for (int i = 0; i < m_nRows; i++) {
		for (int j = 0; j < numGroups; j++, pVertex += m_groupSize, numVertex++) {
			// Add edges to the graph for the vertices in the same matrix row
			auto numNextVertex = numVertex;
			for (int k = j; ++k < numGroups;) {
				++numNextVertex;
				pAdjacencyMatrix[numVertex * m_v + numNextVertex] =
					pAdjacencyMatrix[numNextVertex * m_v + numVertex] = 1;
			}

			// Add edges between vertex pairs whose corresponding tuples of size `m_groupSize`
			// either share exactly one common element (if `cont == true`)
			// or have no elements in common (if `cont == false`).
			auto* pNextVertex = pMatr + (i + 1) * m_nCols;
			while (pNextVertex < pVertexLast) {
				numNextVertex++;
				// We do have two pointers (pVertex and pNextVertex).
				// If one out of groupSize elements is the same, then we have an edge between two vertices.
				if (one_common_element(pVertex, pNextVertex, m_groupSize) == cond) {
					// Add edges to the graph for the vertices intersecting by one element
					pAdjacencyMatrix[numVertex * m_v + numNextVertex] =
						pAdjacencyMatrix[numNextVertex * m_v + numVertex] = 1;
				}

				// Move to the next vertex
				pNextVertex += m_groupSize;
			}
		}
	}

	checkSRG(pAdjacencyMatrix, &graphParam);

	int i, idx, firstVert = 0;
	i = idx = 0;
	while (firstVert = canonizeGraph(m_pGraph[i], m_pGraph[1 - i], firstVert)) {
		printAdjMatrix(m_pGraph[i], m_pGraph[1 - i], idx++);
		checkSRG(m_pGraph[1 - i]);
		i = 1 - i;
	}
}

bool SRGToolkit::checkSRG(tchar* pGraph, SRGParam* pGraphParam) {
	if (pGraphParam)
		pGraphParam->m_cntr[0]++;

	// Check if constructed graph is regular
	int graphDegree = 0;
	auto pVertex = pGraph;
	for (int i = 0; i < m_v; i++, pVertex += m_v) {
		int vertexDegree = 0;
		for (int j = 0; j < m_v; j++)
			if (pVertex[j])
				vertexDegree++;

		if (graphDegree) {
			if (graphDegree != vertexDegree) {
				printfRed("Graph is not regular\n");
				return false;
			}
		}
		else
			graphDegree = vertexDegree;
	}

	if (pGraphParam && !pGraphParam->m_cntr[1]++)
		pGraphParam->k = graphDegree;

	// Check if constructed graph is strongly regular
	int nCommon[4] = { 0 };
	bool flag = true;
	auto pFirstVertex = pGraph;
	for (int i = 0; i < m_v; i++, pFirstVertex += m_v) {
		auto pSecondVertex = pFirstVertex;
		for (int j = i; ++j < m_v;) {
			pSecondVertex += m_v;
			int nCommonCurr, alpha;
			nCommonCurr = alpha = 0;
			for (int k = 0; k < m_v; k++) {
				if (pFirstVertex[k] && pSecondVertex[k])
					m_subgraphVertex[nCommonCurr++] = k;
			}

			if (flag) {
				// Count the number of edges in the subgraph of two vertices common neighbors 
				for (int k = 0; k < nCommonCurr; k++) {
					for (int l = k + 1; l < nCommonCurr; l++) {
						if (pGraph[m_subgraphVertex[k] * m_v + m_subgraphVertex[l]])
							alpha++;
					}
				}
			}

			const auto idx = pFirstVertex[j] ? 0 : 1;
			if (nCommon[idx]) {
				if (nCommon[idx] != nCommonCurr) {
					printfRed("Graph is not strongly regular\n");
					return false;
				}
				if (flag) {
					//printfRed("Graph does not satisfy 4-vertex condition\n");
					flag = (nCommon[idx + 2] == alpha);
				}
			}
			else {
				nCommon[idx] = nCommonCurr;
				nCommon[idx + 2] = alpha;
			}
		}
	}

	if (pGraphParam) {
		if (!pGraphParam->m_cntr[2]++) {
			pGraphParam->λ = nCommon[0];
			pGraphParam->μ = nCommon[1];
			pGraphParam->α = nCommon[2];
			pGraphParam->β = nCommon[3];
		}
		if (!flag)
			pGraphParam->m_cntr[3]++;
	}
	return true;
}

void SRGToolkit::printStat() {
	for (int i = 0; i < 2; i++) {
		auto& graphParam = m_graphParam[i];
		const bool plural = graphParam.m_cntr[0] > 1;
		const char* pntr0 = plural ? "s" : "";
		const char* pntr1 = plural ? "are" : "is";
		const char* pntr2 = plural ? "    " : "";

		printfYellow("Out of %d graph%s with %d vertices\n", graphParam.m_cntr[0], pntr0, m_v);
		printfYellow("       %d %s regular of degree %d\n", graphParam.m_cntr[1], pntr1, graphParam.k);
		printfYellow("       %d %s strongly regular with parameters: (%d,%2d,%2d,%2d)\n", 
			graphParam.m_cntr[2], pntr1, m_v, graphParam.k, graphParam.λ, graphParam.μ);
		const auto n4VertCond = graphParam.m_cntr[2] - graphParam.m_cntr[3];
		if (n4VertCond)
			printfYellow("       %d of them satisfy 4-vertex conditions: (%d, %d)\n", n4VertCond, graphParam.α, graphParam.β);

		const auto v_2k = m_v - 2 * graphParam.k;
		graphParam.k = m_v - graphParam.k - 1;
		const auto λ = graphParam.λ;
		const auto μ = graphParam.μ;
		graphParam.λ = v_2k + graphParam.μ - 2;
		graphParam.μ = v_2k + λ;
		graphParam.α = λ * (λ - 1) / 2 - graphParam.α;
		graphParam.β = μ * (μ - 1) / 2 - graphParam.β;
		printfYellow("%s           complimentary graph parameters: (%d,%2d,%2d,%2d)\n", pntr2, m_v, graphParam.k, graphParam.λ, graphParam.μ);
//		if (n4VertCond)
//			printfYellow("%s            4-vertex condition parameters: (%d, %d)\n", pntr2, graphParam.α, graphParam.β);		
	}
}

int SRGToolkit::canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVert) {
	// Initialize orbits
	for (int i = m_v; i--;)
		m_pOrbits[i] = i;

	int flag = 0;
	// Loop over all vertices
	m_pNumOrbits[0] = 1;
	m_pLenOrbits[0] = m_v;
	auto pLenOrbitsNext = m_pLenOrbits;
	auto pOrbits = m_pOrbits;
	auto pCurVertex = pGraph;
	int len, idxOrb, lastUnfixedVertexIndex = 0;
	int i = 0;// firstVert; // Smalest vertex number that will 
	int idxRight = 1;  // Continue, until at least one orbit length > 1
	while (idxRight && i < m_v) {
		// Loop over all orbits generated by first i vertices
		// Copying current orbits and orbit's lengths to the next level
		auto pLenOrbits = pLenOrbitsNext;
		// Decrease the first orbit's length after we starting use the current vertex
		pLenOrbitsNext += (len = m_v - i - 1);
		if (!--pLenOrbits[idxOrb = 0]) {
			// The current vertex is the sole member of its orbit; 
			// removing the orbit entirely
			memcpy(pLenOrbitsNext, pLenOrbits + 1, len);
			memcpy(pLenOrbits, pLenOrbitsNext, len);
			m_pNumOrbits[i]--;
		}
		else {
			memcpy(pLenOrbitsNext, pLenOrbits, len);
			lastUnfixedVertexIndex = i;
		}

		auto idxLast = i;
		const auto idxOrbMax = m_pNumOrbits[i] = m_pNumOrbits[i++];
//		if (pLenOrbits[idxOrb] > 1)
//			lastUnfixedVertexIndex = idxLast;

		idxRight = 0;
		int idxOrbNext = -1;
		auto pCurVert = pGraph + m_v * m_pOrbits[idxLast];
		while (idxOrb < idxOrbMax) {
			idxOrbNext++;
			auto idxLeft = idxLast;
			const auto lenOrb = pLenOrbits[idxOrb++];
			ASSERT(!lenOrb);
			idxLast += lenOrb;
			if (lenOrb == 1)
				continue;

			// Orbit's length > 1, loop over the vertices from the same orbit
			idxRight = idxLast;
			const auto idxLeftStart = idxLeft;
			int splitPos = 0;
			while (true) {
				while (++idxLeft < idxRight && pCurVert[m_pOrbits[idxLeft]]) splitPos = idxLeft;
				if (idxLeft == idxRight) {
					if (pCurVert[m_pOrbits[idxLeft]])
						splitPos = idxRight;

					break;
				}

				while (idxLeft < idxRight && !pCurVert[m_pOrbits[idxRight]]) idxRight--;

				if (idxLeft >= idxRight)
					break;

				// Mark this swap in orbits
				const auto tmp = m_pOrbits[splitPos = idxLeft];
				m_pOrbits[idxLeft] = m_pOrbits[idxRight];
				m_pOrbits[idxRight] = tmp;
			}
#if 1
			idxLeft = splitPos;
			if (idxLeft && idxLeft < idxLast) {
				// Current orbits was split in two
				m_pNumOrbits[i]++;
				// Moving lengths of the orbits to reseve a place for a new one
				auto j = idxOrbMax + idxOrbNext - idxOrb;
				while (j > idxOrbNext)
					pLenOrbitsNext[j-- + 1] = pLenOrbitsNext[j];

				pLenOrbitsNext[idxOrbNext++] = (idxLeft -= idxLeftStart);
				pLenOrbitsNext[idxOrbNext] = lenOrb - idxLeft;
				ASSERT(!idxLeft || lenOrb == idxLeft);
			}
#else
			if (splitPos) {
				// Current orbits was split in two
				m_pNumOrbits[i]++;
				// Moving lengths of the orbits to reseve a place for a new one
				auto j = idxOrbMax + idxOrbNext - idxOrb;
				while (j > idxOrbNext)
					pLenOrbitsNext[j-- + 1] = pLenOrbitsNext[j];

				pLenOrbitsNext[idxOrbNext++] = (splitPos -= idxLeftStart);
				pLenOrbitsNext[idxOrbNext] = lenOrb - splitPos;
				ASSERT(lenOrb == splitPos);
			}
#endif
		}


		if (!flag) {
			const auto pVertexIn = pGraph + m_pOrbits[i-1] * m_v;
			for (int j = 0; j < m_v; j++)
				pGraphOut[j] = pVertexIn[m_pOrbits[j]];

			flag = memcmp(pGraphOut, pCurVertex, m_v);
		}

		pCurVertex += m_v;
	}

	return flag > 0 ? lastUnfixedVertexIndex : 0;
}

void SRGToolkit::printAdjMatrix(ctchar* pGraph, tchar* pGraphOut, int idx) {
	for (int i = 0; i < m_v; i++) {
		auto pVertexOut = pGraphOut + i * m_v;
		const auto pVertexIn = pGraph + m_pOrbits[i] * m_v;
		for (int j = 0; j < m_v; j++) {
			pVertexOut[j] = pVertexIn[m_pOrbits[j]];
		}
	}

	char buf[256], * pBuf;
	snprintf(buf, sizeof(buf), "aaa_%02d.txt", idx);
	FOPEN_F(f, buf, "w");
	fprintf(f, "Adjacency matrix for graph% d\n", idx);
	
	pBuf = buf;
	int iMax = 20;
	int k = 0;
	for (int j = 0; j < 3; j++) {
		if (j == 2)
			iMax = m_v % iMax;
		for (int i = 0; i < iMax; i++) {
			SPRINTFD(pBuf, buf, "%3d", m_pOrbits[k++]);
		}
		SPRINTFD(pBuf, buf, "\n");
	}
	fprintf(f, "%s\n", buf);


	pBuf = buf;
	for (int i = 0; i < m_v; i += 10) {
		SPRINTFD(pBuf, buf, "%1d", i / 10);
		for (int j = 0; j < 9; j++)
			SPRINTFD(pBuf, buf, " ");
	}
	fprintf(f, "%s\n", buf);

	pBuf = buf;
	for (int i = 0; i < m_v; i++)
		SPRINTFD(pBuf, buf, "%1d", i%10);

	fprintf(f, "%s\n", buf);

	auto pVertexOut = pGraphOut;
	for (int i = 0; i < m_v; i++, pVertexOut += m_v) {
		pBuf = buf;
		for (int j = 0; j < m_v; j++)
			SPRINTFD(pBuf, buf, "%1d", pVertexOut[j]);

		*(pBuf - m_v + i) = '*';		// *'s on diagonal
		SPRINTFD(pBuf, buf, " %2d", i); // row #'s from left
		fprintf(f, "%s\n", buf);
	}
	FCLOSE_F(f);
}