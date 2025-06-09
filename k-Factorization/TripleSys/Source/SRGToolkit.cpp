#include "SRGToolkit.h"
#include <cstring>

#define PRINT_MATRICES 0
#if PRINT_MATRICES
#define PRINT_ADJ_MATRIX(...) printAdjMatrix(__VA_ARGS__)
#else
#define PRINT_ADJ_MATRIX(...)
#endif

SRGToolkit::SRGToolkit(int nCols, int nRows, int groupSize) : 
	m_nCols(nCols), m_nRows(nRows), m_groupSize(groupSize), m_v(nRows * nCols/groupSize) {
	m_pGraph[0] = new tchar[2 * m_v * m_v];
	m_pGraph[1] = m_pGraph[0] + m_v * m_v;
	m_subgraphVertex = new ushort[m_v];
	m_pOrbits = new ushort[m_v];
	m_pNumOrbits = new ushort[2 * m_v];
	m_pGroupOrbits = m_pNumOrbits + m_v;
	m_len = m_v * (m_v - 1) / 2;
	m_pLenOrbits = new ushort [3 * m_len];
	m_pSavedOrbits = m_pLenOrbits + m_len;
	m_pSavedOrbIdx = m_pSavedOrbits + m_len;
	for (int i = 0; i < 2; i++) {
		m_bChekMatr[i] = true;
		m_pGraphParam[i] = new SRGParam();
		m_pMarixStorage[i] = new CBinaryMatrixStorage(m_len, 50 * m_len);
	}
}

SRGToolkit::~SRGToolkit() { 
	delete[] m_pGraph[0];
	delete[] m_subgraphVertex;
	delete[] m_pNumOrbits;
	delete[] m_pOrbits;
	for (int i = 0; i < 2; i++) {
		delete m_pMarixStorage[i];
		delete m_pGraphParam[i];
	}
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

bool SRGToolkit::exploreMatrix(ctchar* pMatr) {
	int counter = 0;
	for (int i = 0; i < 2; i++) {
		if (m_bChekMatr[i] && (m_bChekMatr[i] = exploreMatrixOfType(i, pMatr)))
			counter++;
	}

	return counter > 0;
}

bool SRGToolkit::exploreMatrixOfType(int typeIdx, ctchar* pMatr) {
	const auto numGroups = m_nCols / m_groupSize;
	const auto pVertexLast = pMatr + m_nRows * m_nCols;

	auto graphParam = m_pGraphParam[typeIdx];
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

	const auto graphType = checkSRG(pAdjacencyMatrix, graphParam);
	switch (graphType) {
	case t_nonregular:
	case t_regular:
	case t_complete:
		delete graphParam;
		m_pGraphParam[typeIdx] = NULL;
		return false;
	}

	initCanonizer();
	int i, firstVert = 0;
	i = 0;
	while (firstVert = canonizeGraph(m_pGraph[i], m_pGraph[1 - i], firstVert)) {
		createGraphOut(m_pGraph[i], m_pGraph[1 - i]);
		i = 1 - i;
	}

	// Copy elements above the main diagonal into the array.
	auto pFrom = m_pGraph[i];
	auto pTo = m_pGraph[1 - i];
	auto pGraph = pTo;
	for (int j = m_v, i = 0; --j; pTo += j, pFrom += m_v)
		memcpy(pTo, pFrom + ++i, j);

	m_pMarixStorage[typeIdx]->updateRepo(pGraph);
	PRINT_ADJ_MATRIX(m_pGraph[i], m_pMarixStorage[typeIdx]->numObjects(), m_v);
	return true;
}

void SRGToolkit::initCanonizer() {
	// Initialize orbits
	m_pNumOrbits[0] = 1;
	m_pLenOrbits[0] = m_v;
}

void SRGToolkit::initVertexGroupOrbits() {
	setGroupOrder(1);
	for (int i = m_v; i--;)
		m_pGroupOrbits[i] = i;
}

tchar *SRGToolkit::createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex, int endVertex) const {
	if (!endVertex)
		endVertex = m_v;

	int i = startVertex;
	auto pVertexOutRet = pGraphOut + i * m_v;
	auto pVertexOut = pVertexOutRet;
	for (; i < endVertex; i++, pVertexOut += m_v) {
		const auto pVertexIn = pGraph + m_pOrbits[i] * m_v;
		for (int j = 0; j < m_v; j++) {
			pVertexOut[j] = pVertexIn[m_pOrbits[j]];
		}
	}

	return pVertexOutRet;
}

ushort* SRGToolkit::restoreParam(int &i, int iStart, ushort * pLenOrbits) {
	do {
		// Restore previous parameters;
		while (i > iStart && !*(pLenOrbits -= m_v - i--));

		if (i == iStart)
			return NULL;

	} while (m_pSavedOrbIdx[i]++ >= *pLenOrbits);

	return pLenOrbits;
}

int SRGToolkit::canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVert) {
#if PRINT_MATRICES
	if (firstVert) {
		// Copy the top part of the adjacency matrix
		memcpy(pGraphOut, pGraph, m_v * (m_v - firstVert));
	}
#endif
	const auto defineAut = firstVert > 0;
	int lastUnfixedVertexIndex = 0;
	// Vertex index which will will be swapped with the vertex firstVert
	int indVertMax = 1;
	int indVert = 0;
	ushort* pLenOrbitsPrev = m_pLenOrbits;
	if (defineAut) {
		pLenOrbitsPrev += (2 * m_v - firstVert - 1) * firstVert / 2;
		initVertexGroupOrbits();
	}

	while (true) {
		// Loop over all vertices
		// Initialize orbits
		for (int i = m_v; i--;)
			m_pOrbits[i] = i;

		ushort* pLenOrbits = pLenOrbitsPrev;
		int i = firstVert;
		while (defineAut) {
			if (indVert >= indVertMax) {
				while (--i >= 0) {
					// The smallest index of the vertex that will be considered first
					pLenOrbits = m_pLenOrbits + (2 * m_v - i - 1) * i / 2;
					if (*pLenOrbits)
						break;	
				}

				if (i < 0) {
					updateGroupOrder(m_v, m_pGroupOrbits);
					return 0;
				}
				
				firstVert = i;
				indVert = 0;
				indVertMax = *(pLenOrbitsPrev = pLenOrbits);
			}

			const auto pGroupOrbs = m_pGroupOrbits + i;
			while (++indVert <= indVertMax && pGroupOrbs[indVert] != i + indVert);
#if 0
			static int ff;
			FOPEN_F(f, "bbb.txt", ff++ ? "a" : "w");
			fprintf(f, "ff = %2d  i = %2d  indVert = %2d\n", ff, i, indVert);
			FCLOSE_F(f);
#endif
			if (indVert <= indVertMax) {
				// swapping of i-th and (i+indVert)-th vertices
				m_pOrbits[m_pOrbits[i] = i + indVert] = i;
				pLenOrbits[0]++;
				break;
			}
		}

		int flag = 0;
		int idxRight = 1;  // Continue, until at least one orbit length > 1
		if (defineAut) {
			const auto iStart = i;
			while (idxRight && i < m_v) {
				if (*pLenOrbits > 1)
					m_pSavedOrbIdx[i] = 0;
canonRow:
				flag = canonizeMatrixRow(pGraph, pGraphOut, i++, &pLenOrbits, idxRight, flag, lastUnfixedVertexIndex);
				if (flag < 0) {
					if (i == iStart)
						break;

					pLenOrbits = restoreParam(i, iStart, pLenOrbits);
					if (!pLenOrbits)
						break;
						
					idxRight = 1;
					flag = 0;
					const auto j = m_pSavedOrbIdx[i];
					const auto tmp = m_pOrbits[i];
					m_pOrbits[i] = m_pOrbits[i + j];
					m_pOrbits[i + j] = tmp;
					++*pLenOrbits;
					goto canonRow;
				}
			}

			if (!flag) {
#if 0
				static int fff; fff = 0;
				FOPEN_F(f, "aaa.txt", fff++ ? "a" : "w");
				char buf[256], * pBuf = buf;
				SPRINTFD(pBuf, buf, " %2d (%2d): ", fff, idx);
				for (int k = 0; k < m_v; k++)
					SPRINTFD(pBuf, buf, " %2d", m_pGroupOrbits[k]);
				fprintf(f, "%s\n", buf);
				FCLOSE_F(f);

				FOPEN_F(f1, "ccc.txt", fff>1 ? "a" : "w");
				pBuf = buf;
				SPRINTFD(pBuf, buf, " %2d (%2d): ", fff, idx);
				for (int k = 0; k < m_v; k++)
					SPRINTFD(pBuf, buf, " %2d", m_pOrbits[k]);
				fprintf(f1, "%s\n", buf);
				FCLOSE_F(f1);
#endif
				// Build the matrix of the graph to the end
				const auto pGraphLast = createGraphOut(pGraph, pGraphOut, i);
				flag = memcmp(pGraphLast, pGraph + i * m_v, m_v * (m_v - i));
				if (!flag)
					addAutomorphism(m_v, m_pOrbits, m_pGroupOrbits, true);
			}
		}
		else {
			while (idxRight && i < m_v) {
				flag = canonizeMatrixRow(pGraph, pGraphOut, i++, &pLenOrbits, idxRight, flag, lastUnfixedVertexIndex);
			}
		}

		if (flag > 0 || !defineAut)
			return lastUnfixedVertexIndex;
	}

	return 0;
}

int SRGToolkit::canonizeMatrixRow(ctchar* pGraph, tchar* pGraphOut, int vertIdx,
	ushort ** ppLenOrbits, int& idxRight, int flag, int &lastUnfixedVertexIndex) {

	auto idxLast = vertIdx;
	const auto vertIdxNext = vertIdx + 1;
	auto len = m_v - vertIdxNext;
	auto pLenOrbits = *ppLenOrbits;
	auto pLenOrbitsNext = pLenOrbits + len;
	len *= sizeof(pLenOrbits[0]);

	// Loop over all orbits generated by first vertIdx vertices
	// Copying current orbits and orbit's lengths to the next level
	// Decrease the first orbit's length after we starting use the current vertex
	if (!--pLenOrbits[0]) {
		// The current vertex is the sole member of its orbit; 
		// removing the orbit entirely
		memcpy(pLenOrbitsNext, pLenOrbits + 1, len);
		pLenOrbits++;  // compensation for removed orbits
		m_pNumOrbits[vertIdx]--;
	}
	else {
		memcpy(pLenOrbitsNext, pLenOrbits, len);
		lastUnfixedVertexIndex = vertIdx;
	}

	const auto idxOrbMax = m_pNumOrbits[vertIdxNext] = m_pNumOrbits[vertIdx];

	int idxOrb = idxRight = 0;
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

		if (splitPos && splitPos < idxLast) {
			// Current orbits was split in two
			m_pNumOrbits[vertIdxNext]++;
			// Moving lengths of the orbits to reseve a place for a new one
			auto j = idxOrbMax + idxOrbNext - idxOrb;
			while (j > idxOrbNext)
				pLenOrbitsNext[j-- + 1] = pLenOrbitsNext[j];

			pLenOrbitsNext[idxOrbNext++] = (splitPos -= idxLeftStart);
			pLenOrbitsNext[idxOrbNext] = lenOrb - splitPos;
			ASSERT(!splitPos || lenOrb == splitPos);
		}
	}

	if (!flag) {
		const auto pVertexIn = pGraph + m_pOrbits[vertIdx] * m_v;
		auto pVertOut = pGraphOut + vertIdx * m_v;
		for (int j = 0; j < m_v; j++)
			pVertOut[j] = pVertexIn[m_pOrbits[j]];

		flag = memcmp(pVertOut, pGraph + vertIdx * m_v, m_v);
	}

	*ppLenOrbits = pLenOrbitsNext;
	return flag;
}


t_graphType SRGToolkit::checkSRG(tchar* pGraph, SRGParam* pGraphParam) {
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
				return t_nonregular;
			}
		}
		else
			graphDegree = vertexDegree;
	}

	if (pGraphParam && !pGraphParam->m_cntr[1]++)
		pGraphParam->k = graphDegree;

	// Check if constructed graph is strongly regular
	int nCommon[10] = { 0 };
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
				if (nCommonCurr == graphDegree - 1)
					return t_complete;   // Complete graph, we don't need them

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
#if 0
					printfRed("Graph is not strongly regular:\n"
						"For(%d, %d) %s is %d and not %d as it was for (%d, %d)\n", 
						i, j, idx? "mu" : "lambda", nCommonCurr, nCommon[idx], nCommon[4 * idx + 4], nCommon[4 * idx + 5]);
					printAdjMatrix(pGraph);
#endif
					return t_regular;
				}
				if (flag) {
					flag = (nCommon[idx + 2] == alpha);
					if (flag)
						continue;
#if 0
					printfRed("Graph does not satisfy 4-vertex condition\n"
						"For (%d, %d) %s is %d and not %d as it was for (%d, %d)\n",
						i, j, idx ? "beta" : "alpha", alpha, nCommon[idx+2], nCommon[4 * idx + 4], nCommon[4 * idx + 5]);
					printAdjMatrix(pGraph);
#endif
				}
			}
			else {
				nCommon[idx] = nCommonCurr;
				nCommon[idx + 2] = alpha;
				nCommon[4 * idx + 4] = i;
				nCommon[4 * idx + 5] = j;
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
	return t_srg;
}

void SRGToolkit::printStat() {
	for (int i = 0; i < 2; i++) {
		if (!m_pGraphParam[i])
			continue;

		auto& graphParam = *m_pGraphParam[i];
		if (!m_bChekMatr[i]) {
			printfRed("At least one out of %d graphs of type %d is not SRG.\n", graphParam.m_cntr[0], i);
			continue;
		}

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
		printfYellow("       %d of these graphs %s non-isomorphic\n", m_pMarixStorage[i]->numObjects(), pntr1);
	}
}

#if PRINT_MATRICES
void SRGToolkit::printAdjMatrix(tchar* pGraphOut, int idx, int endVertex) const {
	char buf[256], * pBuf;
	snprintf(buf, sizeof(buf), "aaa_%02d.txt", idx);
	FOPEN_F(f, buf, "w");
	fprintf(f, "Adjacency matrix for graph% d\n", idx);
	
	pBuf = buf;
	int iMax = m_v <= 30? m_v : 20;
	int k = 0;
	const auto jMax = (m_v + iMax - 1) / iMax;
	for (int j = 0; j++ < jMax;) {
		if (j == jMax && jMax > 1)
			iMax = m_v % iMax;

		for (int i = 0; i < iMax; i++)
			SPRINTFD(pBuf, buf, "%3d", m_pOrbits[k++]);

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

	if (!endVertex)
		endVertex = m_v;

	auto pVertexOut = pGraphOut;
	for (int i = 0; i < endVertex; i++, pVertexOut += m_v) {
		pBuf = buf;
		for (int j = 0; j < m_v; j++)
			SPRINTFD(pBuf, buf, "%1d", pVertexOut[j]);

		*(pBuf - m_v + i) = '*';		// *'s on diagonal
		SPRINTFD(pBuf, buf, " %2d", i); // row #'s from left
		fprintf(f, "%s\n", buf);
	}
	FCLOSE_F(f);
}
#endif