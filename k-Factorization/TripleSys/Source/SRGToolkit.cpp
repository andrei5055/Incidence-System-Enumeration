#include "SRGToolkit.h"
#include <cstring>

#pragma execution_character_set("utf-8")

#define PRINT_MATRICES	0
#if PRINT_MATRICES
#define PRINT_ADJ_MATRIX(...) printAdjMatrix(__VA_ARGS__)
#define FFF 8    // 6 - for 26, 9 - for 22
#define FF_ 2    // 4 - for 26. 2 for 22
tchar* pGraph[2] = { NULL };
int nIter = 0;
#else
#define PRINT_ADJ_MATRIX(...)
#endif

SRGToolkit::SRGToolkit(int nCols, int nRows, int groupSize, const std::string& resFileName) :
	m_nCols(nCols), m_nRows(nRows), m_groupSize(groupSize), m_v(nRows * nCols/groupSize),
	m_resFileName(resFileName), Generators<ushort>(0, "\nVertex orbits and group generators of graph", nRows * nCols / groupSize) {
	setOutFileName(NULL);
	const int coeff = PRINT_MATRICES? 2 : 0;
	m_pGraph[0] = new tchar[(2 + coeff) * m_v * m_v];
	m_pGraph[1] = m_pGraph[0] + m_v * m_v;
	m_subgraphVertex = new ushort[m_v];
	m_pOrbits = new ushort[(1 + coeff) * m_v];
	m_pNumOrbits = new ushort[2 * m_v];
	m_pGroupOrbits = m_pNumOrbits + m_v;
	m_len = m_v * (m_v - 1) / 2;
	m_pLenOrbits = new ushort [3 * m_len];
	m_pSavedOrbits = m_pLenOrbits + m_len;
	m_pSavedOrbIdx = m_pSavedOrbits + m_len;
	for (int i = 0; i < 2; i++) {
		m_bChekMatr[i] = true;
		m_pGraphParam[i] = new SRGParam();
		m_pMarixStorage[i] = new CBinaryMatrixStorage(m_len, 50);
	}
}

SRGToolkit::~SRGToolkit() { 
	delete[] m_pGraph[0];
	delete[] m_subgraphVertex;
	delete[] m_pNumOrbits;
	delete[] m_pOrbits;
	delete[] outFileName();
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

bool SRGToolkit::exploreMatrix(ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID) {
	int counter = 0;
	for (int i = 0; i < 2; i++) {
		if (m_bChekMatr[i])
			if (exploreMatrixOfType(i, pMatr, pGraphDB+i, sourceMatrID))
				counter++;
#if !CHECK_NON_SRG			
			else 
				m_bChekMatr[i] = false;
#endif
	}

	return counter > 0;
}

bool SRGToolkit::exploreMatrixOfType(int typeIdx, ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID) {
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
	pGraphDB->setGraphType(graphType);
	switch (graphType) {
	case t_nonregular:
	case t_complete:
		delete graphParam;
		m_pGraphParam[typeIdx] = NULL;
		return false;
	case t_regular:
		if (!CHECK_NON_SRG) {
			delete graphParam;
			m_pGraphParam[typeIdx] = NULL;
			return false;
		}
	}

	initCanonizer();
	int i, firstVert = 0;
	i = 0;
#if PRINT_MATRICES
	new ushort[m_v];
	for (int j = m_v; j--;)
		m_pOrbits[j] = j;

	auto pInitOrbits = m_pOrbits + m_v;
	auto pResOrbits = pInitOrbits + m_v;
	memcpy(pGraph[0] = m_pGraph[i] + 2 * m_v * m_v, m_pGraph[0], m_v * m_v * sizeof(m_pGraph[0]));
	pGraph[1] = pGraph[0] + m_v * m_v;
	memcpy(pInitOrbits, m_pOrbits, m_v * sizeof(m_pOrbits[0]));
	PRINT_ADJ_MATRIX(m_pGraph[i], -1, m_v);
#endif
	while (firstVert = canonizeGraph(m_pGraph[i], m_pGraph[1 - i], firstVert)) {
		createGraphOut(m_pGraph[i], m_pGraph[1 - i]);
#if PRINT_MATRICES
		for (int j = m_v; j--;)
			pResOrbits[j] = pInitOrbits[m_pOrbits[j]];

		memcpy(pInitOrbits, pResOrbits, m_v * sizeof(m_pOrbits[0]));
		createGraphOut(pGraph[0], pGraph[1], 0, 0, pInitOrbits);
		PRINT_ADJ_MATRIX(pGraph[1], nIter, m_v, pInitOrbits, "bbb");
		PRINT_ADJ_MATRIX(m_pGraph[1-i], nIter++, m_v);
#endif
		i = 1 - i;
	}

	// Copy elements above the main diagonal into the array.
	const auto pResGraph = m_pGraph[i];
	auto pFrom = m_pGraph[i];
	auto pTo = m_pGraph[1 - i];
	auto pGraph = pTo;
	for (int j = m_v, i = 0; --j; pTo += j, pFrom += m_v)
		memcpy(pTo, pFrom + ++i, j);

	// Copying vertex orbits and trivial permutation
	auto pntr = this->getObject(0);
	memcpy(pntr, m_pGroupOrbits, m_v * sizeof(*pntr));
	pntr += 2 * m_v;
	bool rank3 = true;
	for (i = 0; i < m_v; i++) {
		pntr[i] = i;
		rank3 &= !m_pGroupOrbits[i];
	}

	// Analyze the stabilizer of first vertex
	pntr -= m_v;
	const auto j = graphParam->k + 1;
	for (i = 1; i < j; i++)
		rank3 &= pntr[i] == 1;

	while (i < m_v)
		rank3 &= pntr[i++] == j;

	const char* pGraphDescr = "";
	if (rank3)
		pGraphDescr = "rank 3 graph";
	else
	if (graphType == t_4_vert)
		pGraphDescr = "4-vertex condition";

	char buf[512], *pBuf = buf;
	if (graphType != t_regular)
		SPRINTFD(pBuf, buf, "Strongly regular graphs with parameters: (v,k,λ μ) = (%d,%2d,%d,%d)",
			m_v, graphParam->k, graphParam->λ, graphParam->μ);
	else
		SPRINTFD(pBuf, buf, "Regular graphs with parameters: (v,k) = (%d,%2d)", m_v, graphParam->k);

	pGraphDB->setTableTitle(buf);
	if (rank3)
		graphParam->m_cntr[4]++;

	const auto prevMatrNumb = m_pMarixStorage[typeIdx]->numObjects();
	m_pMarixStorage[typeIdx]->updateRepo(pGraph);
	const auto newGraph = prevMatrNumb < m_pMarixStorage[typeIdx]->numObjects() ? 1 : 0;
	if (newGraph) {
		pBuf = buf;
		// New SRG constructed
		if (!prevMatrNumb) {
			std::string fileName;
#if OUT_SRG_TO_SEPARATE_FILE
			SPRINTFD(pBuf, buf, "%d_matrices.txt", typeIdx + 1);
			fileName = m_resFileName + buf;
#else
			fileName = m_resFileName;
#endif			
			const auto len = fileName.length() + 1;
			delete[] outFileName();
			memcpy(pBuf = new char[len], fileName.c_str(), len);
			setOutFileName(pBuf, false);
		}

		FOPEN_F(f, outFileName(), prevMatrNumb || !OUT_SRG_TO_SEPARATE_FILE ? "a" : "w");

		pBuf = buf;
#if OUT_SRG_TO_SEPARATE_FILE
		if (!prevMatrNumb) {
			if (graphType != t_regular)
				fprintf(f, "List of SRGs of type %d with parameters (v,k,λ μ) = (%d,%2d,%d,%d):\n", typeIdx + 1, m_v, graphParam->k, graphParam->λ, graphParam->μ);
			else
				fprintf(f, "List of regular graphs of type %d with parameters (v,k) = (%d,%2d):\n", typeIdx + 1, m_v, graphParam->k);
		}

		SPRINTFD(pBuf, buf, "\nGraph #%d:  |Aut(G)| = %zd", prevMatrNumb + 1, groupOrder());
#else
		if (graphType != t_regular)
			SPRINTFD(pBuf, buf, "\nSRG #%d of type %d with parameters (v,k,λ μ) = (%d,%2d,%d,%d): |Aut(G)| = %zd", 
				prevMatrNumb + 1, typeIdx + 1, m_v, graphParam->k, graphParam->λ, graphParam->μ, groupOrder());
		else
			SPRINTFD(pBuf, buf, "\nRegular graph #%d of type %d with parameters (v,k) = (%d,%2d): |Aut(G)| = %zd",
				prevMatrNumb + 1, typeIdx + 1, m_v, graphParam->k, groupOrder());
#endif
		if (rank3)
			SPRINTFD(pBuf, buf, "\nIt's a rank 3 graph with");
		else
		if (graphType == t_4_vert)
			SPRINTFD(pBuf, buf, "\n4-vertex condition satisfied");


		if (rank3 || graphType == t_4_vert)
			SPRINTFD(pBuf, buf, " (α, β) = (%d, %d)", graphParam->α, graphParam->β);

		fprintf(f, "%s\n", buf);
		outAdjMatrix(pResGraph, f);
		FCLOSE_F(f);

		makeGroupOutput(NULL, false, false);
	}
	else {
		// Graph is isomorphic to a previously constructed one — adjust statistics to avoid duplicate counting.
		--graphParam->m_cntr[0];
		--graphParam->m_cntr[1];
		if (graphType != t_regular) 
			--graphParam->m_cntr[2];

		if (graphType != t_4_vert)
			--graphParam->m_cntr[3];

		if (rank3)
			--graphParam->m_cntr[4];
	}

	pGraphDB->addObjDescriptor(groupOrder(), pGraphDescr, newGraph, sourceMatrID);

//	PRINT_ADJ_MATRIX(pResGraph, m_pMarixStorage[typeIdx]->numObjects(), m_v);
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

tchar *SRGToolkit::createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex, int endVertex, const ushort* pOrb) const {
	if (!endVertex)
		endVertex = m_v;

	if (!pOrb)
		pOrb = m_pOrbits;

	int i = startVertex;
	auto pVertexOutRet = pGraphOut + i * m_v;
	auto pVertexOut = pVertexOutRet;
	for (; i < endVertex; i++, pVertexOut += m_v) {
		const auto pVertexIn = pGraph + pOrb[i] * m_v;
		for (int j = 0; j < m_v; j++)
			pVertexOut[j] = pVertexIn[pOrb[j]];
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
	resetGroup();
	releaseAllObjects();
	reserveObjects(3);
	bool stabFlag = true;

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
#if FFF
		static int ff;
		char buffer[32];
		if (nIter >= FFF) {
			FOPEN_F(f, "bbb.txt", ff++ ? "a" : "w");
			fprintf(f, "==> ff = %2d  firstVert = %2d\n", ff, firstVert);
			FCLOSE_F(f);
		}
#endif
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
#if FFF
			if (nIter >= FFF) {
				FOPEN_F(f, "bbb.txt", ff++ ? "a" : "w");
				fprintf(f, "ff = %2d  i = %2d  indVert = %2d\n", ff, i, indVert);
				FCLOSE_F(f);
			}
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
#if FFF
				if (nIter >= FFF) {
					static int hhh;
					const bool flg = i == FF_ && idxRight == 1;
					if (flg) {
						hhh++;
						FOPEN_F(f, "bbb.txt", ff++ ? "a" : "w");
						fprintf(f, "ff = %2d  canonizeMatrixRow: hhh = %3d\n", ff, hhh);
						FCLOSE_F(f);
						sprintf_s(buffer, sizeof(buffer), "ccc_%04d.txt", hhh);
						printAdjMatrix(pGraphOut, buffer, i+1);
					}
					sprintf_s(buffer, sizeof(buffer), "ccc_%04d.txt", hhh);
					FOPEN_F(f, buffer, /*flg ? "w" :*/ "a");
					fprintf(f, "canonizeMatrixRow: i = %2d  flag = %d idxRight = %2d\n", i, flag, idxRight);
					FCLOSE_F(f);
				}
#endif
				flag = canonizeMatrixRow(pGraph, pGraphOut, i++, &pLenOrbits, idxRight, flag, lastUnfixedVertexIndex);
				if (!flag && !idxRight) {
					// All remaining orbits consist of single elements
					const auto pGraphLast = createGraphOut(pGraph, pGraphOut, i, i + 1);
					flag = memcmp(pGraphLast, pGraph + i * m_v, m_v);
				}

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
				const auto pGraphLast = createGraphOut(pGraph, pGraphOut, ++i);
#if FFF
				if (nIter >= FFF)
					printAdjMatrix(pGraphOut, buffer, m_v);
#endif
				flag = memcmp(pGraphLast, pGraph + i * m_v, m_v * (m_v - i));
				if (!flag) {
					if (m_pOrbits[0] && stabFlag) {
						stabFlag = false;
						memcpy(this->getObject(1), m_pGroupOrbits, sizeof(m_pGroupOrbits[0]) * m_v);
					}

					addAutomorphism(m_v, m_pOrbits, m_pGroupOrbits, true);
				}
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
			// Moving lengths of the orbits to reserve a place for a new one
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
						i, j, idx ? "β" : "α", alpha, nCommon[idx+2], nCommon[4 * idx + 4], nCommon[4 * idx + 5]);
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

	return pGraphParam->updateParam(nCommon, flag);
}

t_graphType SRGParam::updateParam(int* pCommon, bool flag_4_ver) {
	if (this) {
		if (flag_4_ver) {
			if (α && α != pCommon[2] || β && β != pCommon[3])
				printfRed("Found graph with 4-vertex condition for different (α, β) = (%d, %d) != (%d, %d)\n",
					pCommon[2], pCommon[3], α, β);
			α = pCommon[2];
			β = pCommon[3];

		}
		else
			m_cntr[3]++;  // # of graphs NOT satisfying 4-vertex condition

		if (!m_cntr[2]++) {
			λ = pCommon[0];
			μ = pCommon[1];
		}
	}
	return flag_4_ver ? t_4_vert : t_srg;
}

void SRGToolkit::printStat() {
	for (int i = 0; i < 2; i++) {
		if (!m_pGraphParam[i])
			continue;

		auto& graphParam = *m_pGraphParam[i];
		if (!m_bChekMatr[i]) {
			printfRed("At least one out of %d graphs of type %d is not SRG.\n", graphParam.m_cntr[0], i);
			//continue;
		}

		const bool plural = graphParam.m_cntr[0] > 1;
		const char* pntr0 = plural ? "s" : "";
		const char* pntr1 = plural ? "are" : "is";

		printfYellow("\nConstructed %d graph%s of type %d with %d vertices\n", graphParam.m_cntr[0], pntr0, i + 1, m_v);
		if (!graphParam.m_cntr[2]) {
			printfYellow(" • %d %s regular of degree %d\n", graphParam.m_cntr[1], pntr1, graphParam.k);
			continue;
		}

		printfYellow(" • %d %s strongly regular with parameters: (v, k, λ, μ) = (%d,%2d,%d,%d)\n",
			graphParam.m_cntr[2], pntr1, m_v, graphParam.k, graphParam.λ, graphParam.μ);

		if (graphParam.m_cntr[4])
			printfYellow(" • %d - rank 3 graph%s\n", graphParam.m_cntr[4], (graphParam.m_cntr[4] > 1? "s" : ""));

		const auto n4VertCond = graphParam.m_cntr[2] - graphParam.m_cntr[3] - graphParam.m_cntr[4];
		if (n4VertCond)
			printfRed(" • %d graph%s satisf%s 4-vertex conditions: (α = %d, β = %d)\n", 
				n4VertCond, (n4VertCond > 1? "s" : ""), (n4VertCond > 1 ? "y" : "ies"), graphParam.α, graphParam.β);

		const auto v_2k = m_v - 2 * graphParam.k;
		const auto k = m_v - graphParam.k - 1;
		const auto λ = v_2k + graphParam.μ - 2;
		const auto μ = v_2k + graphParam.λ;
		//graphParam.α = λ * (λ - 1) / 2 - graphParam.α;
		//graphParam.β = μ * (μ - 1) / 2 - graphParam.β;
		printfYellow("Complementary graph parameters: (v, k, λ, μ) = (%d,%2d,%2d,%2d)\n", m_v, k, λ, μ);
		const auto α = graphParam.k - graphParam.β;
		const auto β = v_2k + graphParam.μ - graphParam.α - 2;
		//if (n4VertCond)
		//	printfYellow("            4-vertex condition parameters: (%d, %d)\n",  α, β);
	}
}

#if PRINT_MATRICES
void SRGToolkit::printAdjMatrix(ctchar* pGraphOut, int idx, int endVertex, ushort* pOrb, const char *fName) const {
	char buf[1024];
	snprintf(buf, sizeof(buf), "%s_%02d.txt", fName, idx);
	printAdjMatrix(pGraphOut, buf, endVertex, idx, pOrb);
}

void SRGToolkit::printAdjMatrix(ctchar* pGraphOut, const char* fileName, int endVertex, int idx, ushort *pOrb) const {
	FOPEN_F(f, fileName, "w");
	fprintf(f, "Adjacency matrix for graph %d\n", idx);
	
	char buf[1024], * pBuf;
	pBuf = buf;
	int iMax = m_v <= 30? m_v : 20;
	int k = 0;
	if (!pOrb)
		pOrb = m_pOrbits;

	const auto jMax = (m_v + iMax - 1) / iMax;
	for (int j = 0; j++ < jMax;) {
		if (j == jMax && jMax > 1)
			iMax = m_v % iMax;

		for (int i = 0; i < iMax; i++)
			SPRINTFD(pBuf, buf, "%3d", pOrb[k++]);

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

	outAdjMatrix(pGraphOut, f, endVertex);
	FCLOSE_F(f);
}
#endif

void SRGToolkit::outAdjMatrix(ctchar* pGraphOut, FILE *f, int endVertex) const {
	if (!endVertex)
		endVertex = m_v;

	char buf[1024], * pBuf;
	for (int i = 0; i < endVertex; i++, pGraphOut += m_v) {
		pBuf = buf;
		for (int j = 0; j < m_v; j++)
			SPRINTFD(pBuf, buf, "%1d", pGraphOut[j]);

		*(pBuf - m_v + i) = '*';		// *'s on diagonal
		SPRINTFD(pBuf, buf, " %2d", i); // row #'s from left
		fprintf(f, "%s\n", buf);
	}
}

template<>
void Generators<ushort>::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	this->printTable(this->getObject(0), false, false, this->numObjects(), "");
}

