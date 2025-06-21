#include "SRGToolkit.h"
#include <cstring>

#pragma execution_character_set("utf-8")

#define PRINT_MATRICES	0
#if PRINT_MATRICES
#define PRINT_ADJ_MATRIX(...) printAdjMatrix(__VA_ARGS__)
#define FFF 8    // 6 - for 26, 9 - for 22
#define FF_ 2    // 4 - for 26. 2 for 22
tchar* pGraph[2] = { NULL };
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

bool SRGToolkit::exploreMatrix(ctchar* pMatr) {
	int counter = 0;
	for (int i = 0; i < 2; i++) {
		if (m_bChekMatr[i])
			if (exploreMatrixOfType(i, pMatr))
				counter++;
#if !CHECK_NON_SRG			
			else 
				m_bChekMatr[i] = false;
#endif
	}

	return counter > 0;
}
int nIter = 0;
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
	case t_complete:
		delete graphParam;
		m_pGraphParam[typeIdx] = NULL;
		return false;
	case t_regular:
		if (!CHECK_NON_SRG) {
			delete graphParam;
			m_pGraphParam[typeIdx] = NULL;
		}
		return false;
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
//	PRINT_ADJ_MATRIX(m_pGraph[i], -1, m_v);
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

#if 0
	ushort perm[] = { 0, 1, 5, 4, 6, 7, 3, 2,18,16,17,19,23,21,20,22,58,56,57,59,63,61,60,
		62,26,24,25,29,27,28,30,31,50,48,49,53,51,52,55,54,34,32,33,37,35,36,38,39,42,40,41,
		45,43,44,47,46,10,8, 9,11,15,13,12,14 };
	ushort permA[64];
	for (int i = 0; i < 64; i++)
		permA[i] = perm[pInitOrbits[i]];

	createGraphOut(pGraph[0], pGraph[1], 0, 0, permA);
	PRINT_ADJ_MATRIX(pGraph[1], nIter, m_v, permA, "ccc");
#endif

#if 0
	extern ushort kSetPerm[];
	const auto* pGr = m_pGraph[i];
	const auto len = m_v * m_v * sizeof(*pGr);
	ushort* perm;
	ushort permA[64], permB[64];
	perm = kSetPerm;
	int jInd = 0;
	for (int j = 0; j < 1344; j++, perm += 64) {
		for (int i = 0; i < 64; i++)
			permA[i] = perm[pInitOrbits[i]];

		createGraphOut(pGraph[0], pGraph[1], 0, 0, permA);
		const auto cmp = memcmp(pGraph[1], pGr, len);
		if (!cmp) {
			const auto permX = permA; // perm;
			if (j && memcmp(permB, permX, m_v * sizeof(*permA)) < 0)
				continue;
			
			jInd = j;
			memcpy(permB, permX, m_v*sizeof(*permA));
			continue;
		}

		jInd += 0;
	}

	Table<ushort> BestPerm(NULL, 1, 64, 0, 0);
	BestPerm.setOutFileName("bestPerm.txt");
	BestPerm.printTable(permB, false, false, 1);// , NULL, NULL);
#endif

	// Copy elements above the main diagonal into the array.
	const auto pResGraph = m_pGraph[i];
	auto pFrom = m_pGraph[i];
	auto pTo = m_pGraph[1 - i];
	auto pGraph = pTo;
	for (int j = m_v, i = 0; --j; pTo += j, pFrom += m_v)
		memcpy(pTo, pFrom + ++i, j);

	const auto prevMatrNumb = m_pMarixStorage[typeIdx]->numObjects();
	m_pMarixStorage[typeIdx]->updateRepo(pGraph);
	if (prevMatrNumb < m_pMarixStorage[typeIdx]->numObjects()) {
		char buf[512], *pBuf = buf;
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

		if (rank3)
			graphParam->m_cntr[4]++;

		FOPEN_F(f, outFileName(), prevMatrNumb || !OUT_SRG_TO_SEPARATE_FILE ? "a" : "w");

		pBuf = buf;
#if OUT_SRG_TO_SEPARATE_FILE
		if (!prevMatrNumb)
			fprintf(f, "List of SRGs of type %d with parameters (%d, %d, %d, %d):\n", typeIdx + 1, m_v, graphParam->k, graphParam->λ, graphParam->μ);

		SPRINTFD(pBuf, buf, "\nGraph #%d:  |Aut(G)| = %zd", prevMatrNumb + 1, groupOrder());
#else
		SPRINTFD(pBuf, buf, "\nSRG #%d of type %d with parameters (%d, %d, %d, %d): |Aut(G)| = %zd", 
			prevMatrNumb + 1, typeIdx + 1, m_v, graphParam->k, graphParam->λ, graphParam->μ, groupOrder());
#endif
		if (rank3)
			SPRINTFD(pBuf, buf, "\nIt's a rank 3 graph with");
		else
		if (graphType == t_4_vert)
			SPRINTFD(pBuf, buf, "\n4-vertex condition satisfied");

		if (rank3 || graphType == t_4_vert)
			SPRINTFD(pBuf, buf, " (alpha = % d, beta = % d)", graphParam->α, graphParam->β);

		fprintf(f, "%s\n", buf);
		outAdjMatrix(pResGraph, f);
		FCLOSE_F(f);

		makeGroupOutput(NULL, false, false);
	}
	else {
		// Graph is isomorphic to a previously constructed one — adjust statistics to avoid duplicate counting.
		--graphParam->m_cntr[0];
		--graphParam->m_cntr[1];
		--graphParam->m_cntr[2];
		if (graphType != t_4_vert)
			--graphParam->m_cntr[3];
	}

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
				const auto pGraphLast = createGraphOut(pGraph, pGraphOut, i);
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
		if (flag) {
			if (pGraphParam->m_cntr[2] == pGraphParam->m_cntr[3]) {
				pGraphParam->α = nCommon[2];
				pGraphParam->β = nCommon[3];
			}
			else {
				if (pGraphParam->α != nCommon[2] || pGraphParam->β != nCommon[3])
					printfRed("Found graph with 4-vertex condition for different alpha an beta");
			}
		} else
			pGraphParam->m_cntr[3]++;  // # of graphs not satisfiing 4-vertex condition

		if (!pGraphParam->m_cntr[2]++) {
			pGraphParam->λ = nCommon[0];
			pGraphParam->μ = nCommon[1];
		}
	}
	return flag? t_4_vert : t_srg;
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
			return;
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

ushort kSetPerm[] =
{
	0, 1, 2,3, 4, 5, 6, 7, 8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
		0, 1, 5, 4, 6, 7, 3, 2,18,16,17,19,23,21,20,22,58,56,57,59,63,61,60,62,26,24,25,29,27,28,30,31,50,48,49,53,51,52,55,54,34,32,33,37,35,36,38,39,42,40,41,45,43,44,47,46,10, 8,9,11,15,13,12,14,
		0, 1, 7, 6, 3, 2,4, 5,57,58,56,59,62,61,63,60, 9,10, 8,11,14,13,15,12,25,26,24,28,29,27,30,31,41,42,40,44,45,43,46,47,49,50,48,52,53,51,55,54,33,34,32,36,37,35,39,38,17,18,16,19,22,21,23,20,
		0, 2,3, 1, 4, 6, 7, 5,10, 8,9,11,14,12,13,15,50,48,49,54,52,55,53,51,42,40,41,43,45,46,47,44,18,16,17,22,20,21,23,19,58,56,57,62,60,61,63,59,34,32,33,38,36,35,39,37,26,24,25,30,31,28,27,29,
		0, 2,5, 7, 1, 3, 4, 6,24,25,26,30,27,28,29,31,8,9,10,11,13,12,15,14,40,41,42,45,46,43,47,44,56,57,58,60,61,62,63,59,32,33,34,36,35,38,37,39,16,17,18,20,21,22,19,23,48,49,50,54,53,55,51,52,
		0, 2,6, 4, 7, 5, 1, 3,49,50,48,54,51,55,52,53,25,26,24,30,29,28,31,27,41,42,40,46,43,45,47,44,33,34,32,35,38,36,37,39,17,18,16,21,22,20,23,19,57,58,56,61,62,60,59,63, 9,10, 8,11,15,12,14,13,
		0, 3, 1, 2,4, 7, 5, 6, 9,10, 8,11,13,14,12,15,33,34,32,39,36,37,35,38,57,58,56,62,61,63,59,60,49,50,48,53,52,55,51,54,25,26,24,27,31,28,29,30,17,18,16,23,20,22,19,21,41,42,40,47,44,45,43,46,
		0, 3, 6, 5, 2,1, 4, 7,42,40,41,47,43,45,46,44,10, 8,9,11,12,14,15,13,58,56,57,61,63,62,59,60,26,24,25,31,28,27,29,30,18,16,17,20,22,23,21,19,50,48,49,52,55,53,54,51,34,32,33,39,35,37,38,36,
		0, 3, 7, 4, 5, 6, 2,1,32,33,34,39,38,37,36,35,40,41,42,47,46,45,44,43,56,57,58,63,62,61,59,60,16,17,18,22,23,20,21,19,48,49,50,55,53,52,51,54,24,25,26,28,27,31,30,29, 8,9,10,11,15,14,13,12,
		0, 4, 1, 5, 6, 2,7, 3,16,17,18,19,21,20,23,22,48,49,50,54,51,52,53,55, 8,9,10,12,13,14,11,15,40,41,42,44,43,46,45,47,24,25,26,29,31,27,28,30,56,57,58,62,63,60,59,61,32,33,34,39,35,36,37,38,
		0, 4, 2,6, 7, 3, 5, 1,50,48,49,54,55,52,51,53,34,32,33,39,38,36,35,37,10, 8,9,14,12,13,11,15,58,56,57,60,62,63,61,59,42,40,41,46,44,43,45,47,26,24,25,27,29,31,30,28,18,16,17,19,22,20,21,23,
		0, 4, 3, 7, 5, 1, 6, 2,33,34,32,39,37,36,38,35,17,18,16,19,23,20,22,21,9,10, 8,13,14,12,11,15,25,26,24,31,27,29,28,30,57,58,56,63,60,62,61,59,41,42,40,43,46,44,47,45,49,50,48,54,53,52,55,51,
		0, 5, 3, 6, 2,7, 1, 4,40,41,42,47,45,46,43,44,24,25,26,30,28,27,31,29,32,33,34,38,37,36,39,35,48,49,50,53,55,51,52,54,56,57,58,61,60,63,62,59, 8,9,10,13,12,15,11,14,16,17,18,19,22,23,20,21,
		0, 5, 4, 1, 6, 3, 2,7,17,18,16,19,20,23,21,22,41,42,40,47,43,46,44,45,33,34,32,37,36,38,39,35,57,58,56,60,63,61,62,59, 9,10, 8,12,15,13,14,11,49,50,48,55,51,53,54,52,25,26,24,30,31,27,29,28,
		0, 5, 7, 2,1, 4, 6, 3,26,24,25,30,29,27,28,31,18,16,17,19,21,23,22,20,34,32,33,36,38,37,39,35,10, 8,9,15,13,12,14,11,50,48,49,51,53,55,52,54,58,56,57,63,61,60,59,62,42,40,41,47,44,46,45,43,
		0, 6, 1, 7, 3, 5, 2,4,58,56,57,59,61,63,62,60,42,40,41,47,45,43,44,46,18,16,17,23,21,20,19,22,34,32,33,35,37,38,36,39,26,24,25,28,31,29,27,30,10, 8,9,12,14,15,11,13,50,48,49,54,53,51,52,55,
		0, 6, 4, 2,7, 1, 3, 5,48,49,50,54,52,51,55,53,56,57,58,59,62,63,60,61,16,17,18,21,20,23,19,22,24,25,26,31,29,28,27,30, 8,9,10,14,15,12,13,11,32,33,34,37,38,35,39,36,40,41,42,47,44,43,46,45,
		0, 6, 5, 3, 2,4, 7, 1,41,42,40,47,46,43,45,44,49,50,48,54,55,51,53,52,17,18,16,20,23,21,19,22, 9,10, 8,15,12,14,13,11,33,34,32,38,35,37,36,39,25,26,24,29,28,31,30,27,57,58,56,59,60,63,61,62,
		0, 7, 2,5, 1, 6, 3, 4,25,26,24,30,28,29,27,31,57,58,56,59,61,62,60,63,49,50,48,51,55,52,54,53,17,18,16,22,21,23,20,19,41,42,40,45,44,46,43,47, 9,10, 8,14,13,15,11,12,33,34,32,39,35,38,36,37,
		0, 7, 4, 3, 5, 2,1, 6,34,32,33,39,36,38,37,35,26,24,25,30,27,29,31,28,50,48,49,55,52,51,54,53,42,40,41,44,46,45,43,47,10, 8,9,13,15,14,12,11,18,16,17,21,23,22,19,20,58,56,57,59,60,62,63,61,
		0, 7, 6, 1, 3, 4, 5, 2,56,57,58,59,63,62,61,60,32,33,34,39,37,38,35,36,48,49,50,52,51,55,54,53, 8,9,10,15,14,13,12,11,16,17,18,23,22,21,20,19,40,41,42,46,45,44,47,43,24,25,26,30,31,29,28,27,
		1, 0, 3, 2,5, 4, 7, 6, 9, 8,11,10,13,12,15,14,18,19,16,17,22,23,20,21,27,28,29,24,25,26,31,30,35,36,37,32,33,34,39,38,43,44,45,40,41,42,47,46,51,52,53,48,49,50,55,54,59,58,57,56,63,62,61,60,
		1, 0, 4, 5, 7, 6, 2,3,16,18,19,17,21,23,22,20,57,59,58,56,60,62,63,61,29,27,28,26,24,25,31,30,53,51,52,50,48,49,54,55,37,35,36,34,32,33,39,38,45,43,44,42,40,41,46,47,11,9, 8,10,14,12,13,15,
		1, 0, 6, 7, 2,3, 5, 4,58,57,59,56,61,62,60,63, 8,11,9,10,15,12,14,13,28,29,27,25,26,24,31,30,44,45,43,41,42,40,47,46,52,53,51,49,50,48,54,55,36,37,35,33,34,32,38,39,19,16,18,17,20,23,21,22,
		1, 2,0, 3, 5, 6, 4, 7, 8,11,9,10,12,15,13,14,36,37,35,38,33,34,32,39,58,57,59,61,62,60,56,63,52,53,51,50,49,54,48,55,28,29,27,24,30,25,26,31,19,16,18,21,22,20,17,23,44,45,43,46,41,42,40,47,
		1, 2,6, 5, 4, 7, 3, 0,35,36,37,38,39,34,33,32,43,44,45,46,47,42,41,40,59,58,57,60,61,62,56,63,18,19,16,20,21,22,23,17,51,52,53,54,50,49,48,55,27,28,29,25,24,30,31,26, 9, 8,11,10,14,15,12,13,
		1, 2,7, 4, 3, 0, 5, 6,45,43,44,46,40,42,47,41,11,9, 8,10,13,15,14,12,57,59,58,62,60,61,56,63,29,27,28,30,25,24,26,31,16,18,19,22,20,21,23,17,53,51,52,49,54,50,55,48,37,35,36,38,32,34,39,33,
		1, 3, 2,0, 5, 7, 6, 4,11,9, 8,10,15,13,12,14,53,51,52,55,49,54,50,48,45,43,44,40,42,47,46,41,16,18,19,20,22,23,21,17,57,59,58,61,63,62,60,56,37,35,36,39,33,32,38,34,29,27,28,31,30,25,24,26,
		1, 3, 4, 6, 0, 2,5, 7,27,28,29,31,24,25,26,30, 9, 8,11,10,12,13,14,15,43,44,45,42,47,40,46,41,59,58,57,63,62,61,60,56,35,36,37,33,32,39,34,38,18,19,16,22,23,20,17,21,51,52,53,55,50,54,48,49,
		1, 3, 7, 5, 6, 4, 0, 2,52,53,51,55,48,54,49,50,28,29,27,31,26,25,30,24,44,45,43,47,40,42,46,41,36,37,35,32,39,33,34,38,19,16,18,23,20,22,21,17,58,57,59,62,61,63,56,60, 8,11,9,10,14,13,15,12,
		1, 4, 2,7, 3, 6, 0, 5,43,44,45,46,42,47,40,41,27,28,29,31,25,24,30,26,35,36,37,39,34,33,38,32,51,52,53,50,54,48,49,55,59,58,57,62,63,60,61,56, 9, 8,11,12,13,14,10,15,18,19,16,17,20,21,22,23,
		1, 4, 5, 0, 7, 2,3, 6,19,16,18,17,22,21,23,20,44,45,43,46,40,47,41,42,36,37,35,34,33,39,38,32,58,57,59,63,60,62,61,56, 8,11,9,13,14,12,15,10,52,53,51,54,48,50,55,49,28,29,27,31,30,24,26,25,
		1, 4, 6, 3, 0, 5, 7, 2,29,27,28,31,26,24,25,30,16,18,19,17,23,21,20,22,37,35,36,33,39,34,38,32,11,9, 8,14,12,13,15,10,53,51,52,48,50,54,49,55,57,59,58,60,62,63,56,61,45,43,44,46,41,47,42,40,
		1, 5, 0, 4, 7, 3, 6, 2,18,19,16,17,23,22,21,20,51,52,53,55,48,49,50,54, 9, 8,11,13,12,15,10,14,43,44,45,41,40,47,42,46,27,28,29,26,30,24,25,31,59,58,57,61,60,63,56,62,35,36,37,38,32,33,34,39,
		1, 5, 2,6, 4, 0, 7, 3,36,37,35,38,34,33,39,32,19,16,18,17,21,22,20,23, 8,11,9,12,15,13,10,14,28,29,27,30,24,26,25,31,58,57,59,60,63,61,62,56,44,45,43,40,47,41,46,42,52,53,51,55,50,49,54,48,
		1, 5, 3, 7, 6, 2,4, 0,53,51,52,55,54,49,48,50,37,35,36,38,39,33,32,34,11,9, 8,15,13,12,10,14,57,59,58,63,61,60,62,56,45,43,44,47,41,40,42,46,29,27,28,24,26,30,31,25,16,18,19,17,20,22,23,21,
		1, 6, 3, 4, 0, 7, 2,5,28,29,27,31,25,26,24,30,58,57,59,56,62,61,63,60,52,53,51,48,54,49,55,50,19,16,18,20,23,21,22,17,44,45,43,42,41,47,40,46, 8,11,9,15,12,14,10,13,36,37,35,38,32,39,33,34,
		1, 6, 5, 2,4, 3, 0, 7,37,35,36,38,33,39,34,32,29,27,28,31,24,26,30,25,53,51,52,54,49,48,55,50,45,43,44,41,47,42,40,46,11,9, 8,12,14,15,13,10,16,18,19,23,21,20,17,22,57,59,58,56,63,61,60,62,
		1, 6, 7, 0, 2,5, 4, 3,59,58,57,56,60,61,62,63,35,36,37,38,34,39,32,33,51,52,53,49,48,54,55,50, 9, 8,11,14,15,12,13,10,18,19,16,21,20,23,22,17,43,44,45,47,42,41,46,40,27,28,29,31,30,26,25,24,
		1, 7, 0, 6, 2,4, 3, 5,57,59,58,56,62,60,61,63,45,43,44,46,42,40,41,47,16,18,19,21,23,22,17,20,37,35,36,32,34,39,33,38,29,27,28,25,30,26,24,31,11,9, 8,13,15,14,10,12,53,51,52,55,50,48,49,54,
		1, 7, 4, 2,3, 5, 6, 0,44,45,43,46,47,40,42,41,52,53,51,55,54,48,50,49,19,16,18,22,21,23,17,20, 8,11,9,14,13,15,12,10,36,37,35,39,32,34,33,38,28,29,27,26,25,30,31,24,58,57,59,56,63,60,62,61,
		1, 7, 5, 3, 6, 0, 2,4,51,52,53,55,49,48,54,50,59,58,57,56,61,60,63,62,18,19,16,23,22,21,17,20,27,28,29,30,26,25,24,31,9, 8,11,15,14,13,12,10,35,36,37,34,39,32,38,33,43,44,45,46,41,40,47,42,
		2,0, 1, 3, 6, 4, 5, 7, 8,10,11,9,12,14,15,13,49,54,50,48,53,51,52,55,43,45,46,42,40,41,44,47,22,20,21,18,16,17,19,23,62,60,61,58,56,57,59,63,38,36,35,34,32,33,37,39,30,25,24,26,29,27,28,31,
		2,0, 4, 6, 5, 7, 3, 1,50,49,54,48,55,51,53,52,24,30,25,26,31,27,29,28,46,43,45,41,42,40,44,47,35,38,36,33,34,32,39,37,21,22,20,17,18,16,19,23,61,62,60,57,58,56,63,59,11,8,10, 9,13,14,12,15,
		2,0, 7, 5, 3, 1, 6, 4,25,24,30,26,28,27,31,29,10,11,8,9,15,14,13,12,45,46,43,40,41,42,44,47,60,61,62,56,57,58,59,63,36,35,38,32,33,34,39,37,20,21,22,16,17,18,23,19,54,50,49,48,52,51,55,53,
		2,1, 3, 0, 6, 5, 7, 4,11,8,10, 9,15,12,14,13,35,38,36,37,32,39,33,34,61,62,60,58,57,59,63,56,50,49,54,52,53,51,55,48,24,30,25,28,29,27,31,26,21,22,20,19,16,18,23,17,46,43,45,44,47,40,42,41,
		2,1, 4, 7, 0, 3, 6, 5,43,45,46,44,42,40,41,47, 8,10,11,9,14,12,13,15,62,60,61,57,59,58,63,56,30,25,24,29,27,28,31,26,22,20,21,16,18,19,17,23,49,54,50,53,51,52,48,55,38,36,35,37,33,39,34,32,
		2,1, 5, 6, 7, 4, 0, 3,36,35,38,37,34,39,32,33,45,46,43,44,41,40,47,42,60,61,62,59,58,57,63,56,20,21,22,18,19,16,17,23,54,50,49,51,52,53,55,48,25,24,30,27,28,29,26,31,10,11,8,9,13,12,15,14,
		2,3, 0, 1, 6, 7, 4, 5,10,11,8,9,14,15,12,13,20,21,22,23,16,17,18,19,25,24,30,28,27,31,26,29,36,35,38,33,32,39,34,37,45,46,43,42,47,40,41,44,54,50,49,55,53,52,48,51,60,61,62,63,56,57,58,59,
		2,3, 5, 4, 1, 0, 6, 7,61,62,60,63,58,57,59,56,11,8,10, 9,12,15,13,14,24,30,25,27,31,28,26,29,46,43,45,47,40,42,41,44,50,49,54,53,52,55,51,48,35,38,36,32,39,33,37,34,21,22,20,23,18,17,19,16,
		2,3, 7, 6, 4, 5, 1, 0,22,20,21,23,19,17,16,18,62,60,61,63,59,57,56,58,30,25,24,31,28,27,26,29,49,54,50,52,55,53,51,48,38,36,35,39,33,32,34,37,43,45,46,40,42,47,44,41,8,10,11,9,13,15,14,12,
		2,4, 3, 5, 1, 7, 0, 6,62,60,61,63,57,59,58,56,43,45,46,44,40,42,47,41,22,20,21,19,17,16,23,18,38,36,35,33,39,34,32,37,30,25,24,27,29,31,28,26, 8,10,11,14,12,13, 9,15,49,54,50,48,52,55,53,51,
		2,4, 6, 0, 5, 3, 1, 7,54,50,49,48,53,55,51,52,60,61,62,63,58,59,56,57,20,21,22,17,16,19,23,18,25,24,30,29,31,27,28,26,10,11,8,12,13,14,15, 9,36,35,38,39,34,33,37,32,45,46,43,44,47,42,41,40,
		2,4, 7, 1, 0, 6, 5, 3,46,43,45,44,41,42,40,47,50,49,54,48,51,55,52,53,21,22,20,16,19,17,23,18,11,8,10,13,14,12,15, 9,35,38,36,34,33,39,32,37,24,30,25,31,27,29,26,28,61,62,60,63,56,59,57,58,
		2,5, 0, 7, 3, 4, 1, 6,24,30,25,26,27,31,28,29,61,62,60,63,57,58,56,59,50,49,54,55,51,53,48,52,21,22,20,18,17,19,16,23,46,43,45,40,47,41,42,44,11,8,10,12,15,13, 9,14,35,38,36,37,33,34,32,39,
		2,5, 4, 3, 1, 6, 7, 0,60,61,62,63,59,58,57,56,36,35,38,37,39,34,33,32,54,50,49,53,55,51,48,52,10,11,8,13,12,15,14, 9,20,21,22,19,18,17,16,23,45,46,43,41,40,47,44,42,25,24,30,26,29,31,27,28,
		2,5, 6, 1, 7, 0, 3, 4,38,36,35,37,32,34,39,33,30,25,24,26,28,31,29,27,49,54,50,51,53,55,48,52,43,45,46,47,41,40,42,44, 8,10,11,15,13,12,14, 9,22,20,21,17,19,18,23,16,62,60,61,63,56,58,59,57,
		2,6, 0, 4, 5, 1, 7, 3,49,54,50,48,51,53,55,52,38,36,35,37,34,32,33,39, 8,10,11,12,14,15, 9,13,62,60,61,56,58,59,57,63,43,45,46,41,47,42,40,44,30,25,24,28,31,29,26,27,22,20,21,23,18,16,17,19,
		2,6, 1, 5, 7, 3, 4, 0,35,38,36,37,39,32,34,33,21,22,20,23,19,16,18,17,11,8,10,15,12,14, 9,13,24,30,25,29,28,31,27,26,61,62,60,59,56,58,57,63,46,43,45,42,41,47,44,40,50,49,54,48,52,53,51,55,
		2,6, 3, 7, 4, 0, 5, 1,20,21,22,23,17,16,19,18,54,50,49,48,55,53,52,51,10,11,8,14,15,12, 9,13,45,46,43,47,42,41,40,44,25,24,30,31,29,28,27,26,60,61,62,58,59,56,63,57,36,35,38,37,33,32,39,34,
		2,7, 1, 4, 0, 5, 3, 6,45,46,43,44,40,41,42,47,25,24,30,26,27,28,29,31,36,35,38,34,39,32,37,33,54,50,49,52,51,55,53,48,60,61,62,57,56,59,58,63,10,11,8,15,14,13, 9,12,20,21,22,23,18,19,16,17,
		2,7, 5, 0, 3, 6, 4, 1,30,25,24,26,31,28,27,29,22,20,21,23,17,19,18,16,38,36,35,32,34,39,37,33, 8,10,11,13,15,14,12, 9,49,54,50,55,52,51,53,48,62,60,61,59,57,56,63,58,43,45,46,44,47,41,40,42,
		2,7, 6, 3, 4, 1, 0, 5,21,22,20,23,16,19,17,18,46,43,45,44,42,41,47,40,35,38,36,39,32,34,37,33,61,62,60,56,59,57,58,63,11,8,10,14,13,15,12, 9,50,49,54,51,55,52,48,53,24,30,25,26,29,28,31,27,
		3, 0, 2,1, 7, 4, 6, 5,10, 9,11,8,14,13,15,12,32,39,33,34,35,38,36,37,62,61,63,57,58,56,60,59,53,52,55,49,50,48,54,51,27,31,28,25,26,24,30,29,23,20,22,17,18,16,21,19,47,40,42,41,46,43,45,44,
		3, 0, 4, 7, 6, 5, 1, 2,33,32,39,34,37,38,35,36,42,47,40,41,44,43,46,45,63,62,61,56,57,58,60,59,22,23,20,16,17,18,19,21,55,53,52,48,49,50,54,51,28,27,31,24,25,26,29,30,11,10, 9, 8,12,13,14,15,
		3, 0, 5, 6, 1, 2,7, 4,40,42,47,41,45,43,44,46, 9,11,10, 8,15,13,12,14,61,63,62,58,56,57,60,59,31,28,27,26,24,25,30,29,20,22,23,18,16,17,19,21,52,55,53,50,48,49,51,54,39,33,32,34,36,38,37,35,
		3, 1, 0, 2,7, 5, 4, 6, 9,11,10, 8,13,15,14,12,52,55,53,51,50,48,49,54,40,42,47,45,43,44,41,46,20,22,23,16,18,19,17,21,61,63,62,57,59,58,56,60,39,33,32,37,35,36,34,38,31,28,27,29,26,24,25,30,
		3, 1, 5, 7, 4, 6, 2,0,53,52,55,51,54,48,50,49,27,31,28,29,30,24,26,25,47,40,42,44,45,43,41,46,32,39,33,36,37,35,38,34,23,20,22,19,16,18,17,21,62,61,63,58,57,59,60,56,10, 9,11,8,12,15,13,14,
		3, 1, 6, 4, 2,0, 7, 5,28,27,31,29,25,24,30,26,11,10, 9, 8,14,15,12,13,42,47,40,43,44,45,41,46,63,62,61,59,58,57,56,60,33,32,39,35,36,37,38,34,22,23,20,18,19,16,21,17,55,53,52,51,49,48,54,50,
		3, 2,1, 0, 7, 6, 5, 4,11,10, 9, 8,15,14,13,12,22,23,20,21,18,19,16,17,28,27,31,25,24,30,29,26,33,32,39,36,35,38,37,34,42,47,40,45,46,43,44,41,55,53,52,54,50,49,51,48,63,62,61,60,59,58,57,56,
		3, 2,4, 5, 0, 1, 7, 6,62,61,63,60,57,58,56,59,10, 9,11,8,13,14,12,15,27,31,28,24,30,25,29,26,47,40,42,46,43,45,44,41,53,52,55,50,49,54,48,51,32,39,33,35,38,36,34,37,23,20,22,21,16,19,17,18,
		3, 2,6, 7, 5, 4, 0, 1,20,22,23,21,17,19,18,16,61,63,62,60,56,58,59,57,31,28,27,30,25,24,29,26,52,55,53,49,54,50,48,51,39,33,32,38,36,35,37,34,40,42,47,43,45,46,41,44, 9,11,10, 8,12,14,15,13,
		3, 4, 1, 6, 2,5, 0, 7,27,31,28,29,24,30,25,26,62,61,63,60,58,57,59,56,53,52,55,54,48,50,51,49,23,20,22,16,19,17,18,21,47,40,42,43,46,44,45,41,10, 9,11,13,14,12, 8,15,32,39,33,34,36,37,35,38,
		3, 4, 5, 2,0, 7, 6, 1,63,62,61,60,56,57,58,59,33,32,39,34,38,37,36,35,55,53,52,50,54,48,51,49,11,10, 9,12,13,14,15, 8,22,23,20,17,16,19,18,21,42,47,40,44,43,46,41,45,28,27,31,29,26,30,24,25,
		3, 4, 7, 0, 6, 1, 2,5,39,33,32,34,35,37,38,36,31,28,27,29,25,30,26,24,52,55,53,48,50,54,51,49,40,42,47,46,44,43,45,41,9,11,10,14,12,13,15, 8,20,22,23,19,17,16,21,18,61,63,62,60,59,57,56,58,
		3, 5, 2,4, 0, 6, 1, 7,61,63,62,60,58,56,57,59,40,42,47,41,43,45,46,44,20,22,23,17,19,18,21,16,39,33,32,36,38,37,35,34,31,28,27,24,26,30,25,29, 9,11,10,15,13,12, 8,14,52,55,53,51,49,54,50,48,
		3, 5, 6, 0, 1, 7, 4, 2,47,40,42,41,44,45,43,46,53,52,55,51,48,54,49,50,23,20,22,18,17,19,21,16,10, 9,11,12,15,13,14, 8,32,39,33,37,36,38,35,34,27,31,28,30,24,26,29,25,62,61,63,60,59,56,58,57,
		3, 5, 7, 1, 4, 2,0, 6,55,53,52,51,50,54,48,49,63,62,61,60,57,56,59,58,22,23,20,19,18,17,21,16,28,27,31,26,30,24,25,29,11,10, 9,13,12,15,14, 8,33,32,39,38,37,36,34,35,42,47,40,41,46,45,44,43,
		3, 6, 0, 5, 1, 4, 2,7,42,47,40,41,43,44,45,46,28,27,31,29,24,25,26,30,33,32,39,37,38,35,34,36,55,53,52,49,48,54,50,51,63,62,61,58,59,56,57,60,11,10, 9,14,15,12, 8,13,22,23,20,21,16,17,18,19,
		3, 6, 4, 1, 2,7, 5, 0,31,28,27,29,30,25,24,26,20,22,23,21,19,17,16,18,39,33,32,35,37,38,34,36, 9,11,10,12,14,15,13, 8,52,55,53,54,49,48,50,51,61,63,62,56,58,59,60,57,40,42,47,41,46,44,43,45,
		3, 6, 7, 2,5, 0, 1, 4,23,20,22,21,18,17,19,16,47,40,42,41,45,44,46,43,32,39,33,38,35,37,34,36,62,61,63,59,56,58,57,60,10, 9,11,15,12,14,13, 8,53,52,55,48,54,49,51,50,27,31,28,29,26,25,30,24,
		3, 7, 0, 4, 6, 2,5, 1,32,39,33,34,38,35,37,36,23,20,22,21,17,18,16,19,10, 9,11,14,13,15, 8,12,27,31,28,26,25,30,24,29,62,61,63,56,59,57,58,60,47,40,42,45,44,46,41,43,53,52,55,51,49,50,48,54,
		3, 7, 1, 5, 4, 0, 6, 2,52,55,53,51,48,50,54,49,39,33,32,34,37,35,36,38,9,11,10,13,15,14, 8,12,61,63,62,59,57,56,58,60,40,42,47,44,46,45,43,41,31,28,27,25,30,26,29,24,20,22,23,21,16,18,19,17,
		3, 7, 2,6, 5, 1, 4, 0,22,23,20,21,19,18,17,16,55,53,52,51,54,50,49,48,11,10, 9,15,14,13, 8,12,42,47,40,46,45,44,43,41,28,27,31,30,26,25,24,29,63,62,61,57,56,59,60,58,33,32,39,34,36,35,38,37,
		4, 0, 5, 1, 2,6, 3, 7,17,16,19,18,20,21,22,23,50,54,48,49,53,55,51,52,12,13,14, 8,9,10,15,11,44,43,46,40,41,42,47,45,29,31,27,24,25,26,30,28,62,63,60,56,57,58,61,59,39,34,33,32,38,37,36,35,
		4, 0, 6, 2,3, 7, 1, 5,48,50,54,49,52,55,53,51,33,39,34,32,35,37,38,36,14,12,13,10, 8,9,15,11,60,62,63,58,56,57,59,61,46,44,43,42,40,41,47,45,27,29,31,26,24,25,28,30,19,17,16,18,23,21,20,22,
		4, 0, 7, 3, 1, 5, 2,6,34,33,39,32,36,37,35,38,16,19,17,18,22,21,23,20,13,14,12, 9,10, 8,15,11,31,27,29,25,26,24,30,28,63,60,62,57,58,56,59,61,43,46,44,41,42,40,45,47,54,48,50,49,51,55,52,53,
		4, 1, 0, 5, 2,7, 6, 3,16,19,17,18,21,22,20,23,43,46,44,45,41,42,40,47,34,33,39,36,37,35,32,38,63,60,62,58,57,59,56,61,13,14,12, 8,11,9,10,15,54,48,50,52,53,51,49,55,31,27,29,28,25,26,24,30,
		4, 1, 3, 6, 5, 0, 2,7,27,29,31,28,24,26,30,25,19,17,16,18,20,22,23,21,33,39,34,37,35,36,32,38,14,12,13,11,9, 8,10,15,48,50,54,53,51,52,55,49,60,62,63,57,59,58,61,56,46,44,43,45,40,42,47,41,
		4, 1, 7, 2,6, 3, 5, 0,44,43,46,45,47,42,41,40,29,31,27,28,30,26,25,24,39,34,33,35,36,37,32,38,50,54,48,51,52,53,55,49,62,63,60,59,58,57,56,61,12,13,14, 9, 8,11,15,10,17,16,19,18,23,22,21,20,
		4, 2,0, 6, 3, 5, 7, 1,50,54,48,49,55,53,52,51,62,63,60,61,56,57,58,59,17,16,19,20,21,22,18,23,29,31,27,25,24,30,26,28,12,13,14,10,11,8,9,15,39,34,33,36,35,38,32,37,44,43,46,45,40,41,42,47,
		4, 2,1, 7, 6, 0, 3, 5,43,46,44,45,42,41,47,40,54,48,50,49,52,53,51,55,16,19,17,21,22,20,18,23,13,14,12,11,8,10, 9,15,34,33,39,35,38,36,37,32,31,27,29,24,30,25,28,26,63,60,62,61,58,57,59,56,
		4, 2,5, 3, 7, 1, 6, 0,60,62,63,61,59,57,56,58,46,44,43,45,47,41,40,42,19,17,16,22,20,21,18,23,33,39,34,38,36,35,37,32,27,29,31,30,25,24,26,28,14,12,13, 8,10,11,15, 9,48,50,54,49,51,53,55,52,
		4, 3, 0, 7, 1, 6, 5, 2,33,39,34,32,37,35,36,38,27,29,31,28,26,24,25,30,48,50,54,52,55,53,49,51,46,44,43,40,42,47,41,45,14,12,13, 9,11,10, 8,15,19,17,16,20,22,23,18,21,60,62,63,61,58,56,57,59,
		4, 3, 2,5, 7, 0, 1, 6,62,63,60,61,57,56,59,58,39,34,33,32,36,35,38,37,50,54,48,55,53,52,49,51,12,13,14,11,10, 9, 8,15,17,16,19,22,23,20,21,18,44,43,46,42,47,40,45,41,29,31,27,28,25,24,30,26,
		4, 3, 6, 1, 5, 2,7, 0,31,27,29,28,30,24,26,25,63,60,62,61,59,56,58,57,54,48,50,53,52,55,49,51,16,19,17,23,20,22,21,18,43,46,44,47,40,42,41,45,13,14,12,10, 9,11,15, 8,34,33,39,32,38,35,37,36,
		4, 5, 1, 0, 2,3, 7, 6,19,17,16,18,22,20,21,23,60,62,63,61,57,59,58,56,27,29,31,24,26,30,28,25,48,50,54,51,53,55,52,49,33,39,34,36,38,37,35,32,46,44,43,47,41,40,45,42,14,12,13,15,11,9, 8,10,
		4, 5, 3, 2,7, 6, 0, 1,63,60,62,61,56,59,57,58,13,14,12,15,10, 9,11,8,31,27,29,30,24,26,28,25,43,46,44,40,47,41,42,45,54,48,50,55,51,53,52,49,34,33,39,37,36,38,32,35,16,19,17,18,23,20,22,21,
		4, 5, 6, 7, 0, 1, 2,3,12,13,14,15, 8,9,10,11,17,16,19,18,21,20,23,22,29,31,27,26,30,24,28,25,39,34,33,38,37,36,35,32,44,43,46,41,40,47,42,45,50,54,48,53,55,51,49,52,62,63,60,61,58,59,56,57,
		4, 6, 1, 3, 5, 7, 0, 2,29,31,27,28,26,30,24,25,12,13,14,15, 9, 8,11,10,44,43,46,47,42,41,45,40,62,63,60,58,59,56,57,61,39,34,33,37,38,35,36,32,17,16,19,21,20,23,18,22,50,54,48,49,51,52,53,55,
		4, 6, 2,0, 3, 1, 5, 7,54,48,50,49,53,52,55,51,31,27,29,28,24,30,25,26,43,46,44,42,41,47,45,40,34,33,39,38,35,37,36,32,16,19,17,20,23,21,22,18,63,60,62,59,56,58,61,57,13,14,12,15,11,8,10, 9,
		4, 6, 7, 5, 0, 2,3, 1,14,12,13,15,10, 8,9,11,48,50,54,49,55,52,51,53,46,44,43,41,47,42,45,40,19,17,16,23,21,20,22,18,60,62,63,56,58,59,57,61,33,39,34,35,37,38,32,36,27,29,31,28,25,30,26,24,
		4, 7, 2,1, 6, 5, 0, 3,46,44,43,45,41,47,42,40,14,12,13,15, 8,10,11,9,60,62,63,59,57,56,61,58,27,29,31,25,30,26,24,28,19,17,16,21,23,22,20,18,48,50,54,55,52,51,49,53,33,39,34,32,38,36,35,37,
		4, 7, 3, 0, 1, 2,6, 5,39,34,33,32,35,36,37,38,44,43,46,45,42,47,40,41,62,63,60,57,56,59,61,58,17,16,19,23,22,21,20,18,50,54,48,52,51,55,53,49,29,31,27,30,26,25,28,24,12,13,14,15,11,10, 9, 8,
		4, 7, 5, 6, 0, 3, 1, 2,13,14,12,15, 9,10, 8,11,34,33,39,32,37,36,38,35,63,60,62,56,59,57,61,58,54,48,50,51,55,52,53,49,31,27,29,26,25,30,24,28,16,19,17,22,21,23,18,20,43,46,44,45,40,47,41,42,
		5, 0, 1, 4, 3, 6, 7, 2,18,17,19,16,23,20,22,21,40,47,41,42,44,45,43,46,37,36,38,33,34,32,35,39,60,63,61,57,58,56,59,62,12,15,13, 9,10, 8,11,14,55,51,53,49,50,48,52,54,30,24,26,25,28,29,27,31,
		5, 0, 2,7, 4, 1, 3, 6,24,26,30,25,27,29,31,28,17,19,18,16,22,20,21,23,36,38,37,34,32,33,35,39,15,13,12,10, 8,9,11,14,51,53,55,50,48,49,54,52,63,61,60,58,56,57,62,59,47,41,40,42,43,45,46,44,
		5, 0, 6, 3, 7, 2,4, 1,41,40,47,42,46,45,44,43,26,30,24,25,31,29,28,27,38,37,36,32,33,34,35,39,53,55,51,48,49,50,54,52,61,60,63,56,57,58,59,62,13,12,15, 8,9,10,14,11,19,18,17,16,21,20,23,22,
		5, 1, 4, 0, 3, 7, 2,6,19,18,17,16,22,23,20,21,53,55,51,52,50,54,48,49,13,12,15, 9, 8,11,14,10,41,40,47,43,44,45,46,42,26,30,24,27,28,29,31,25,61,60,63,59,58,57,62,56,38,37,36,35,39,34,33,32,
		5, 1, 6, 2,0, 4, 3, 7,37,36,38,35,33,34,32,39,18,17,19,16,20,23,21,22,12,15,13, 8,11,9,14,10,30,24,26,28,29,27,31,25,60,63,61,58,57,59,56,62,40,47,41,44,45,43,42,46,55,51,53,52,48,54,49,50,
		5, 1, 7, 3, 2,6, 0, 4,51,53,55,52,49,54,50,48,36,38,37,35,32,34,39,33,15,13,12,11,9, 8,14,10,63,61,60,57,59,58,56,62,47,41,40,45,43,44,46,42,24,26,30,29,27,28,25,31,17,19,18,16,21,23,22,20,
		5, 2,1, 6, 0, 7, 4, 3,36,38,37,35,34,32,33,39,24,26,30,25,29,27,28,31,51,53,55,49,54,50,52,48,47,41,40,43,45,46,44,42,15,13,12, 8,10,11,9,14,17,19,18,22,20,21,16,23,63,61,60,62,57,59,58,56,
		5, 2,3, 4, 6, 1, 0, 7,61,60,63,62,58,59,56,57,38,37,36,35,33,32,39,34,53,55,51,54,50,49,52,48,13,12,15,10,11,8,9,14,19,18,17,20,21,22,23,16,41,40,47,45,46,43,42,44,26,30,24,25,28,27,31,29,
		5, 2,7, 0, 4, 3, 6, 1,30,24,26,25,31,27,29,28,60,63,61,62,56,59,57,58,55,51,53,50,49,54,52,48,18,17,19,21,22,20,23,16,40,47,41,46,43,45,44,42,12,15,13,11,8,10,14, 9,37,36,38,35,39,32,34,33,
		5, 3, 0, 6, 7, 1, 2,4,40,47,41,42,45,44,46,43,55,51,53,52,49,50,48,54,18,17,19,23,20,22,16,21,12,15,13,10, 9,11,8,14,37,36,38,32,39,33,34,35,30,24,26,27,31,28,25,29,60,63,61,62,57,58,56,59,
		5, 3, 1, 7, 2,4, 6, 0,53,55,51,52,54,50,49,48,61,60,63,62,59,58,57,56,19,18,17,22,23,20,16,21,26,30,24,28,27,31,29,25,13,12,15,11,10, 9, 8,14,38,37,36,33,32,39,35,34,41,40,47,42,43,44,45,46,
		5, 3, 4, 2,6, 0, 7, 1,63,61,60,62,56,58,59,57,47,41,40,42,46,44,43,45,17,19,18,20,22,23,16,21,36,38,37,39,33,32,34,35,24,26,30,31,28,27,29,25,15,13,12, 9,11,10,14, 8,51,53,55,52,48,50,54,49,
		5, 4, 0, 1, 3, 2,6, 7,17,19,18,16,20,22,23,21,63,61,60,62,58,56,57,59,24,26,30,27,29,31,25,28,51,53,55,48,50,54,49,52,36,38,37,33,39,34,32,35,47,41,40,46,44,43,42,45,15,13,12,14,10, 8,9,11,
		5, 4, 2,3, 6, 7, 1, 0,60,63,61,62,59,56,58,57,12,15,13,14,11,8,10, 9,30,24,26,31,27,29,25,28,40,47,41,43,46,44,45,42,55,51,53,54,48,50,49,52,37,36,38,34,33,39,35,32,18,17,19,16,21,22,20,23,
		5, 4, 7, 6, 1, 0, 3, 2,13,12,15,14, 9, 8,11,10,19,18,17,16,23,22,21,20,26,30,24,29,31,27,25,28,38,37,36,39,34,33,32,35,41,40,47,44,43,46,45,42,53,55,51,50,54,48,52,49,61,60,63,62,57,56,59,58,
		5, 6, 2,1, 0, 3, 7, 4,38,37,36,35,32,33,34,39,41,40,47,42,45,46,43,44,61,60,63,58,59,56,62,57,19,18,17,21,20,23,22,16,53,55,51,49,48,54,50,52,26,30,24,31,29,28,25,27,13,12,15,14,10,11,8,9,
		5, 6, 3, 0, 7, 4, 1, 2,47,41,40,42,44,46,45,43,15,13,12,14, 9,11,10, 8,63,61,60,56,58,59,62,57,24,26,30,28,31,29,27,25,17,19,18,23,21,20,22,16,51,53,55,54,49,48,52,50,36,38,37,35,39,33,32,34,
		5, 6, 4, 7, 1, 2,0, 3,12,15,13,14, 8,11,9,10,37,36,38,35,34,33,39,32,60,63,61,59,56,58,62,57,55,51,53,48,54,49,50,52,30,24,26,29,28,31,27,25,18,17,19,20,23,21,16,22,40,47,41,42,43,46,44,45,
		5, 7, 0, 2,4, 6, 1, 3,26,30,24,25,29,31,27,28,13,12,15,14, 8,9,10,11,41,40,47,46,45,44,42,43,61,60,63,57,56,59,58,62,38,37,36,34,39,32,33,35,19,18,17,23,22,21,16,20,53,55,51,52,48,49,50,54,
		5, 7, 3, 1, 2,0, 4, 6,55,51,53,52,50,49,54,48,30,24,26,25,27,31,28,29,40,47,41,45,44,46,42,43,37,36,38,39,32,34,33,35,18,17,19,22,21,23,20,16,60,63,61,56,59,57,62,58,12,15,13,14,10, 9,11,8,
		5, 7, 6, 4, 1, 3, 2,0,15,13,12,14,11,9, 8,10,51,53,55,52,54,49,48,50,47,41,40,44,46,45,42,43,17,19,18,21,23,22,20,16,63,61,60,59,57,56,58,62,36,38,37,32,34,39,35,33,24,26,30,25,28,31,29,27,
		6, 0, 2,4, 1, 7, 5, 3,49,48,54,50,51,52,53,55,58,59,56,57,60,61,62,63,21,20,23,16,17,18,22,19,31,29,28,24,25,26,30,27,14,15,12, 8,9,10,11,13,37,38,35,32,33,34,36,39,47,42,41,40,45,46,43,44,
		6, 0, 3, 5, 4, 2,1, 7,42,41,47,40,43,46,44,45,48,54,49,50,53,52,55,51,20,23,21,17,18,16,22,19,15,12,14, 9,10, 8,11,13,38,35,37,33,34,32,39,36,29,28,31,25,26,24,27,30,59,56,58,57,62,61,63,60,
		6, 0, 7, 1, 5, 3, 4, 2,56,58,59,57,63,61,60,62,41,47,42,40,44,46,45,43,23,21,20,18,16,17,22,19,35,37,38,34,32,33,39,36,28,31,29,26,24,25,30,27,12,14,15,10, 8,9,13,11,54,49,48,50,55,52,51,53,
		6, 1, 0, 7, 5, 2,3, 4,58,59,56,57,61,60,63,62,37,38,35,36,32,33,34,39,49,48,54,51,52,53,50,55,14,15,12, 9, 8,11,10,13,21,20,23,18,19,16,17,22,47,42,41,43,44,45,40,46,31,29,28,27,24,25,26,30,
		6, 1, 2,5, 3, 4, 7, 0,35,37,38,36,39,33,32,34,28,31,29,27,30,25,24,26,54,49,48,53,51,52,50,55,41,47,42,45,43,44,46,40,12,14,15,11,9, 8,10,13,23,21,20,16,18,19,22,17,56,58,59,57,62,60,61,63,
		6, 1, 4, 3, 7, 0, 5, 2,29,28,31,27,26,25,30,24,59,56,58,57,63,60,62,61,48,54,49,52,53,51,50,55,20,23,21,19,16,18,17,22,42,41,47,44,45,43,46,40,15,12,14, 8,11,9,13,10,38,35,37,36,34,33,39,32,
		6, 2,4, 0, 1, 5, 3, 7,54,49,48,50,53,51,52,55,35,37,38,36,33,39,34,32,12,14,15, 8,10,11,13, 9,56,58,59,62,60,61,63,57,41,47,42,43,45,46,44,40,28,31,29,30,25,24,27,26,23,21,20,22,19,17,16,18,
		6, 2,5, 1, 3, 7, 0, 4,38,35,37,36,32,39,33,34,20,23,21,22,18,17,19,16,15,12,14,11,8,10,13, 9,29,28,31,24,30,25,26,27,59,56,58,61,62,60,63,57,42,41,47,46,43,45,40,44,48,54,49,50,55,51,53,52,
		6, 2,7, 3, 0, 4, 1, 5,21,20,23,22,16,17,18,19,49,48,54,50,52,51,55,53,14,15,12,10,11,8,13, 9,47,42,41,45,46,43,44,40,31,29,28,25,24,30,26,27,58,59,56,60,61,62,57,63,37,38,35,36,34,39,32,33,
		6, 3, 1, 4, 7, 2,0, 5,28,31,29,27,25,30,26,24,23,21,20,22,16,18,19,17,35,37,38,39,33,32,36,34,12,14,15, 9,11,10, 8,13,54,49,48,52,55,53,51,50,56,58,59,61,63,62,57,60,41,47,42,40,45,43,44,46,
		6, 3, 2,7, 0, 5, 4, 1,20,23,21,22,17,18,16,19,42,41,47,40,46,43,45,44,38,35,37,32,39,33,36,34,59,56,58,62,61,63,60,57,15,12,14,10, 9,11,8,13,48,54,49,53,52,55,50,51,29,28,31,27,24,30,25,26,
		6, 3, 5, 0, 4, 1, 7, 2,47,42,41,40,44,43,46,45,31,29,28,27,26,30,24,25,37,38,35,33,32,39,36,34,49,48,54,55,53,52,51,50,58,59,56,63,62,61,60,57,14,15,12,11,10, 9,13, 8,21,20,23,22,19,18,17,16,
		6, 4, 0, 2,1, 3, 7, 5,48,54,49,50,52,53,51,55,29,28,31,27,25,26,24,30,42,41,47,43,46,44,40,45,38,35,37,34,33,39,32,36,20,23,21,16,19,17,18,22,59,56,58,63,60,62,57,61,15,12,14,13, 9,10, 8,11,
		6, 4, 3, 1, 7, 5, 2,0,31,29,28,27,30,26,25,24,14,15,12,13,11,10, 9, 8,47,42,41,44,43,46,40,45,58,59,56,62,63,60,61,57,37,38,35,39,34,33,32,36,21,20,23,17,16,19,22,18,49,48,54,50,55,53,52,51,
		6, 4, 5, 7, 2,0, 1, 3,12,14,15,13, 8,10,11,9,54,49,48,50,51,53,55,52,41,47,42,46,44,43,40,45,23,21,20,19,17,16,18,22,56,58,59,60,62,63,61,57,35,37,38,33,39,34,36,32,28,31,29,27,24,26,30,25,
		6, 5, 0, 3, 4, 7, 2,1,41,47,42,40,46,44,43,45,12,14,15,13,10, 8,9,11,56,58,59,63,61,60,57,62,28,31,29,24,26,30,25,27,23,21,20,17,19,18,16,22,54,49,48,51,53,55,50,52,35,37,38,36,34,32,33,39,
		6, 5, 1, 2,3, 0, 4, 7,37,38,35,36,33,32,39,34,47,42,41,40,43,44,45,46,58,59,56,61,60,63,57,62,21,20,23,19,18,17,16,22,49,48,54,53,55,51,52,50,31,29,28,26,30,24,27,25,14,15,12,13, 9, 8,11,10,
		6, 5, 7, 4, 2,1, 3, 0,15,12,14,13,11,8,10, 9,38,35,37,36,39,32,34,33,59,56,58,60,63,61,57,62,48,54,49,55,51,53,52,50,29,28,31,30,24,26,25,27,20,23,21,18,17,19,22,16,42,41,47,40,45,44,46,43,
		6, 7, 1, 0, 5, 4, 2,3,59,56,58,57,60,63,61,62,15,12,14,13, 8,11,9,10,29,28,31,26,25,30,27,24,42,41,47,45,44,46,43,40,48,54,49,51,55,52,53,50,38,35,37,39,32,34,36,33,20,23,21,22,19,16,18,17,
		6, 7, 3, 2,0, 1, 5, 4,23,21,20,22,18,16,17,19,56,58,59,57,61,63,62,60,28,31,29,25,30,26,27,24,54,49,48,55,52,51,53,50,35,37,38,32,34,39,33,36,41,47,42,44,46,45,40,43,12,14,15,13, 9,11,10, 8,
		6, 7, 4, 5, 2,3, 0, 1,14,15,12,13,10,11,8,9,21,20,23,22,17,16,19,18,31,29,28,30,26,25,27,24,37,38,35,34,39,32,33,36,47,42,41,46,45,44,43,40,49,48,54,52,51,55,50,53,58,59,56,57,62,63,60,61,
		7, 0, 1, 6, 4, 3, 2,5,57,56,59,58,62,63,60,61,34,39,32,33,35,36,37,38,52,51,55,48,49,50,53,54,15,14,13, 8,9,10,11,12,23,22,21,16,17,18,19,20,46,45,44,40,41,42,43,47,30,26,25,24,27,28,29,31,
		7, 0, 3, 4, 2,5, 6, 1,32,34,39,33,38,36,35,37,25,30,26,24,31,28,27,29,55,52,51,50,48,49,53,54,44,46,45,42,40,41,47,43,13,15,14,10, 8,9,11,12,21,23,22,18,16,17,20,19,59,57,56,58,61,63,62,60,
		7, 0, 5, 2,6, 1, 4, 3,26,25,30,24,29,28,31,27,56,59,57,58,60,63,61,62,51,55,52,49,50,48,53,54,22,21,23,17,18,16,19,20,45,44,46,41,42,40,47,43,14,13,15, 9,10, 8,12,11,39,32,34,33,37,36,38,35,
		7, 1, 2,4, 5, 3, 0, 6,45,44,46,43,40,47,41,42,51,55,52,53,50,49,54,48,22,21,23,19,16,18,20,17,14,13,15, 8,11,9,10,12,39,32,34,36,37,35,38,33,26,25,30,28,29,27,24,31,56,59,57,58,61,62,60,63,
		7, 1, 3, 5, 0, 6, 4, 2,52,51,55,53,48,49,50,54,57,56,59,58,63,62,61,60,23,22,21,18,19,16,20,17,30,26,25,27,28,29,31,24,15,14,13, 9, 8,11,10,12,34,39,32,35,36,37,33,38,46,45,44,43,42,47,40,41,
		7, 1, 6, 0, 4, 2,5, 3,59,57,56,58,60,62,63,61,44,46,45,43,41,47,42,40,21,23,22,16,18,19,20,17,32,34,39,37,35,36,38,33,25,30,26,29,27,28,31,24,13,15,14,11,9, 8,12,10,55,52,51,53,54,49,48,50,
		7, 2,0, 5, 6, 3, 1, 4,25,30,26,24,28,31,29,27,21,23,22,20,18,16,17,19,32,34,39,38,36,35,33,37,13,15,14, 8,10,11,9,12,55,52,51,49,54,50,48,53,59,57,56,62,60,61,58,63,44,46,45,43,42,40,41,47,
		7, 2,3, 6, 1, 4, 5, 0,22,21,23,20,19,16,18,17,45,44,46,43,47,40,42,41,39,32,34,35,38,36,33,37,56,59,57,61,62,60,63,58,14,13,15,11,8,10, 9,12,51,55,52,50,49,54,53,48,26,25,30,24,27,31,28,29,
		7, 2,4, 1, 5, 0, 6, 3,46,45,44,43,41,40,47,42,30,26,25,24,29,31,27,28,34,39,32,36,35,38,33,37,52,51,55,54,50,49,48,53,57,56,59,60,61,62,63,58,15,14,13,10,11,8,12, 9,23,22,21,20,17,16,19,18,
		7, 3, 4, 0, 2,6, 1, 5,39,32,34,33,35,38,36,37,22,21,23,20,16,19,17,18,14,13,15,10, 9,11,12, 8,26,25,30,27,31,28,29,24,56,59,57,62,61,63,60,58,45,44,46,47,40,42,43,41,51,55,52,53,54,48,50,49,
		7, 3, 5, 1, 0, 4, 2,6,55,52,51,53,50,48,49,54,32,34,39,33,36,38,37,35,13,15,14, 9,11,10,12, 8,59,57,56,61,63,62,60,58,44,46,45,40,42,47,41,43,25,30,26,31,28,27,24,29,21,23,22,20,17,19,18,16,
		7, 3, 6, 2,1, 5, 0, 4,23,22,21,20,18,19,16,17,52,51,55,53,49,48,54,50,15,14,13,11,10, 9,12, 8,46,45,44,42,47,40,41,43,30,26,25,28,27,31,29,24,57,56,59,63,62,61,58,60,34,39,32,33,37,38,35,36,
		7, 4, 0, 3, 2,1, 5, 6,34,39,32,33,36,35,38,37,46,45,44,43,40,41,42,47,57,56,59,62,63,60,58,61,23,22,21,17,16,19,18,20,52,51,55,50,54,48,49,53,30,26,25,29,31,27,24,28,15,14,13,12, 8,9,10,11,
		7, 4, 1, 2,5, 6, 3, 0,44,46,45,43,47,41,40,42,13,15,14,12,11,9, 8,10,59,57,56,60,62,63,58,61,25,30,26,27,29,31,28,24,21,23,22,19,17,16,18,20,55,52,51,48,50,54,53,49,32,34,39,33,37,35,36,38,
		7, 4, 6, 5, 3, 0, 2,1,14,13,15,12,10, 9,11,8,39,32,34,33,38,35,37,36,56,59,57,63,60,62,58,61,51,55,52,54,48,50,49,53,26,25,30,31,27,29,28,24,22,21,23,16,19,17,20,18,45,44,46,43,42,41,47,40,
		7, 5, 1, 3, 0, 2,6, 4,51,55,52,53,49,50,48,54,26,25,30,24,28,29,27,31,45,44,46,40,47,41,43,42,39,32,34,37,36,38,35,33,22,21,23,18,17,19,16,20,56,59,57,60,63,61,58,62,14,13,15,12, 8,11,9,10,
		7, 5, 2,0, 6, 4, 3, 1,30,26,25,24,31,29,28,27,15,14,13,12,10,11,8,9,46,45,44,41,40,47,43,42,57,56,59,61,60,63,62,58,34,39,32,38,37,36,35,33,23,22,21,19,18,17,20,16,52,51,55,53,54,50,49,48,
		7, 5, 4, 6, 3, 1, 0, 2,13,15,14,12, 9,11,10, 8,55,52,51,53,48,50,54,49,44,46,45,47,41,40,43,42,21,23,22,17,19,18,16,20,59,57,56,63,61,60,62,58,32,34,39,36,38,37,33,35,25,30,26,24,27,29,31,28,
		7, 6, 0, 1, 4, 5, 3, 2,56,59,57,58,63,60,62,61,14,13,15,12, 9,10, 8,11,26,25,30,29,28,31,24,27,45,44,46,42,41,47,40,43,51,55,52,48,54,49,50,53,39,32,34,38,35,37,33,36,22,21,23,20,17,18,16,19,
		7, 6, 2,3, 1, 0, 4, 5,21,23,22,20,16,18,19,17,59,57,56,58,62,60,61,63,25,30,26,28,31,29,24,27,55,52,51,54,49,48,50,53,32,34,39,35,37,38,36,33,44,46,45,41,47,42,43,40,13,15,14,12, 8,10,11,9,
		7, 6, 5, 4, 3, 2,1, 0,15,14,13,12,11,10, 9, 8,23,22,21,20,19,18,17,16,30,26,25,31,29,28,24,27,34,39,32,37,38,35,36,33,46,45,44,47,42,41,40,43,52,51,55,49,48,54,53,50,57,56,59,58,61,60,63,62,
		8,9,10,11,12,13,14,15, 0, 1, 2,3, 4, 5, 6, 7,40,44,43,41,45,47,42,46,32,36,38,35,33,39,34,37,24,28,30,27,25,31,26,29,16,19,22,18,17,20,23,21,56,58,62,59,57,61,60,63,48,52,49,51,54,53,50,55,
		8,9,13,12,14,15,11,10,43,40,44,41,46,47,45,42,49,48,52,51,55,53,54,50,38,32,36,39,35,33,34,37,62,56,58,61,59,57,63,60,30,24,28,31,27,25,26,29,22,16,19,20,18,17,21,23, 2,0, 1, 3, 7, 5, 4, 6,
		8,9,15,14,11,10,12,13,52,49,48,51,50,53,55,54, 1, 2,0, 3, 6, 5, 7, 4,36,38,32,33,39,35,34,37,19,22,16,17,20,18,23,21,58,62,56,57,61,59,63,60,28,30,24,25,31,27,29,26,44,43,40,41,42,47,46,45,
		8,10,11,9,12,14,15,13, 2,0, 1, 3, 6, 4, 5, 7,62,56,58,60,57,63,61,59,22,16,19,18,20,23,21,17,43,40,44,42,45,47,46,41,49,48,52,50,54,53,55,51,30,24,28,26,25,27,29,31,38,32,36,34,37,33,35,39,
		8,10,13,15, 9,11,12,14,32,36,38,34,35,33,39,37, 0, 1, 2,3, 5, 4, 7, 6,16,19,22,20,23,18,21,17,48,52,49,54,53,50,55,51,24,28,30,25,27,26,31,29,40,44,43,45,47,42,41,46,56,58,62,60,61,63,59,57,
		8,10,14,12,15,13, 9,11,58,62,56,60,59,63,57,61,36,38,32,34,39,33,37,35,19,22,16,23,18,20,21,17,28,30,24,27,26,25,31,29,44,43,40,47,42,45,46,41,52,49,48,53,50,54,51,55, 1, 2,0, 3, 7, 4, 6, 5,
		8,11,9,10,12,15,13,14, 1, 2,0, 3, 5, 6, 4, 7,28,30,24,29,25,31,27,26,52,49,48,50,53,55,51,54,58,62,56,61,57,63,59,60,36,38,32,35,37,33,39,34,44,43,40,46,45,42,41,47,19,22,16,21,17,20,18,23,
		8,11,14,13,10, 9,12,15,22,16,19,21,18,20,23,17, 2,0, 1, 3, 4, 6, 7, 5,49,48,52,53,55,50,51,54,38,32,36,37,33,35,39,34,43,40,44,45,42,46,47,41,62,56,58,57,63,61,60,59,30,24,28,29,27,31,26,25,
		8,11,15,12,13,14,10, 9,24,28,30,29,26,31,25,27,16,19,22,21,23,20,17,18,48,52,49,55,50,53,51,54,40,44,43,42,46,45,47,41,56,58,62,63,61,57,59,60,32,36,38,33,35,37,34,39, 0, 1, 2,3, 7, 6, 5, 4,
		8,12, 9,13,14,10,15,11,40,44,43,41,47,45,46,42,56,58,62,60,59,57,61,63, 0, 1, 2,4, 5, 6, 3, 7,16,19,22,17,18,23,20,21,32,36,38,39,37,35,33,34,48,52,49,50,55,54,51,53,24,28,30,29,27,25,31,26,
		8,12,10,14,15,11,13, 9,62,56,58,60,63,57,59,61,30,24,28,29,26,25,27,31,2,0, 1, 6, 4, 5, 3, 7,49,48,52,54,50,55,53,51,22,16,19,23,17,18,20,21,38,32,36,35,39,37,34,33,43,40,44,41,42,45,47,46,
		8,12,11,15,13, 9,14,10,28,30,24,29,31,25,26,27,44,43,40,41,46,45,42,47, 1, 2,0, 5, 6, 4, 3, 7,36,38,32,37,35,39,33,34,52,49,48,55,54,50,53,51,19,22,16,18,23,17,21,20,58,62,56,60,61,57,63,59,
		8,13,11,14,10,15, 9,12,16,19,22,21,20,23,18,17,32,36,38,34,33,35,37,39,24,28,30,26,31,25,29,27,56,58,62,61,63,59,57,60,48,52,49,53,54,55,50,51,0, 1, 2,5, 4, 7, 3, 6,40,44,43,41,42,46,45,47,
		8,13,12, 9,14,11,10,15,44,43,40,41,45,46,47,42,19,22,16,21,18,23,17,20,28,30,24,31,25,26,29,27,52,49,48,54,55,53,50,51,1, 2,0, 4, 7, 5, 6, 3,58,62,56,63,59,61,60,57,36,38,32,34,37,35,39,33,
		8,13,15,10, 9,12,14,11,38,32,36,34,39,35,33,37,43,40,44,41,47,46,42,45,30,24,28,25,26,31,29,27, 2,0, 1, 7, 5, 4, 6, 3,62,56,58,59,61,63,57,60,49,48,52,55,53,54,51,50,22,16,19,21,17,23,20,18,
		8,14, 9,15,11,13,10,12,49,48,52,51,53,55,50,54,22,16,19,21,20,18,17,23,43,40,44,46,47,45,41,42,30,24,28,27,31,26,25,29,38,32,36,33,37,39,35,34, 2,0, 1, 4, 6, 7, 3, 5,62,56,58,60,61,59,57,63,
		8,14,12,10,15, 9,11,13,56,58,62,60,57,59,63,61,48,52,49,51,50,55,54,53,40,44,43,47,45,46,41,42,32,36,38,37,39,33,35,34, 0, 1, 2,6, 7, 4, 5, 3,24,28,30,31,26,27,29,25,16,19,22,21,17,18,23,20,
		8,14,13,11,10,12,15, 9,19,22,16,21,23,18,20,17,58,62,56,60,63,59,61,57,44,43,40,45,46,47,41,42, 1, 2,0, 7, 4, 6, 5, 3,28,30,24,26,27,31,25,29,36,38,32,39,33,37,34,35,52,49,48,51,54,55,53,50,
		8,15,10,13, 9,14,11,12,36,38,32,34,33,39,35,37,52,49,48,51,53,50,54,55,58,62,56,59,63,57,60,61,44,43,40,42,47,46,45,41,19,22,16,20,17,23,18,21,1, 2,0, 6, 5, 7, 3, 4,28,30,24,29,27,26,25,31,
		8,15,12,11,13,10, 9,14,30,24,28,29,25,26,31,27,38,32,36,34,35,39,37,33,62,56,58,63,57,59,60,61,22,16,19,17,23,20,18,21,2,0, 1, 5, 7, 6, 4, 3,43,40,44,47,46,42,41,45,49,48,52,51,54,50,55,53,
		8,15,14, 9,11,12,13,10,48,52,49,51,55,50,53,54,24,28,30,29,31,26,27,25,56,58,62,57,59,63,60,61,0, 1, 2,7, 6, 5, 4, 3,40,44,43,46,42,47,45,41,16,19,22,23,20,17,21,18,32,36,38,34,37,39,33,35,
		9, 8,11,10,13,12,15,14, 1, 0, 3, 2,5, 4, 7, 6,43,41,40,44,42,46,45,47,35,33,39,32,36,38,37,34,27,25,31,24,28,30,29,26,18,17,20,16,19,22,21,23,59,57,61,56,58,62,63,60,51,49,52,48,55,50,53,54,
		9, 8,12,13,15,14,10,11,40,43,41,44,47,46,42,45,52,51,49,48,54,50,55,53,39,35,33,38,32,36,37,34,61,59,57,62,56,58,60,63,31,27,25,30,24,28,29,26,20,18,17,22,16,19,23,21,3, 1, 0, 2,6, 4, 5, 7,
		9, 8,14,15,10,11,13,12,49,52,51,48,53,50,54,55, 0, 3, 1, 2,7, 4, 6, 5,33,39,35,36,38,32,37,34,17,20,18,19,22,16,21,23,57,61,59,58,62,56,60,63,25,31,27,28,30,24,26,29,41,40,43,44,45,46,47,42,
		9,10, 8,11,13,14,12,15, 0, 3, 1, 2,4, 7, 5, 6,25,31,27,26,28,30,24,29,49,52,51,53,50,54,48,55,57,61,59,62,58,60,56,63,33,39,35,32,34,36,38,37,41,40,43,47,42,45,44,46,17,20,18,23,19,22,16,21,
		9,10,14,13,12,15,11,8,27,25,31,26,29,30,28,24,18,17,20,23,21,22,19,16,51,49,52,54,53,50,48,55,43,41,40,45,47,42,46,44,59,57,61,60,62,58,56,63,35,33,39,36,32,34,37,38,1, 0, 3, 2,6, 7, 4, 5,
		9,10,15,12,11,8,13,14,20,18,17,23,16,22,21,19, 3, 1, 0, 2,5, 7, 6, 4,52,51,49,50,54,53,48,55,39,35,33,34,36,32,38,37,40,43,41,42,45,47,46,44,61,59,57,58,60,62,63,56,31,27,25,26,24,30,29,28,
		9,11,10, 8,13,15,14,12, 3, 1, 0, 2,7, 5, 4, 6,61,59,57,63,58,60,62,56,20,18,17,16,22,21,23,19,40,43,41,45,42,46,47,44,52,51,49,53,55,50,54,48,31,27,25,29,28,24,26,30,39,35,33,37,34,36,32,38,
		9,11,12,14, 8,10,13,15,35,33,39,37,32,36,38,34, 1, 0, 3, 2,4, 5, 6, 7,18,17,20,22,21,16,23,19,51,49,52,55,50,53,54,48,27,25,31,28,24,29,30,26,43,41,40,42,46,45,44,47,59,57,61,63,62,60,56,58,
		9,11,15,13,14,12, 8,10,57,61,59,63,56,60,58,62,33,39,35,37,38,36,34,32,17,20,18,21,16,22,23,19,25,31,27,24,29,28,30,26,41,40,43,46,45,42,47,44,49,52,51,50,53,55,48,54, 0, 3, 1, 2,6, 5, 7, 4,
		9,12,10,15,11,14, 8,13,18,17,20,23,22,21,16,19,35,33,39,37,36,32,34,38,27,25,31,29,30,28,26,24,59,57,61,62,60,56,58,63,51,49,52,50,55,54,53,48,1, 0, 3, 4, 5, 6, 2,7,43,41,40,44,45,47,42,46,
		9,12,13, 8,15,10,11,14,41,40,43,44,42,47,46,45,17,20,18,23,16,21,19,22,25,31,27,30,28,29,26,24,49,52,51,55,54,50,53,48,0, 3, 1, 5, 6, 4, 7, 2,57,61,59,60,56,62,63,58,33,39,35,37,34,32,38,36,
		9,12,14,11,8,13,15,10,39,35,33,37,38,32,36,34,40,43,41,44,46,47,45,42,31,27,25,28,29,30,26,24, 3, 1, 0, 6, 4, 5, 7, 2,61,59,57,56,62,60,58,63,52,51,49,54,50,55,48,53,20,18,17,23,19,21,22,16,
		9,13, 8,12,15,11,14,10,43,41,40,44,46,42,47,45,59,57,61,63,56,58,62,60, 1, 0, 3, 5, 4, 7, 2,6,18,17,20,19,16,21,22,23,35,33,39,38,34,32,36,37,51,49,52,53,54,55,48,50,27,25,31,26,24,28,30,29,
		9,13,10,14,12, 8,15,11,25,31,27,26,30,28,29,24,41,40,43,44,47,42,45,46, 0, 3, 1, 4, 7, 5, 2,6,33,39,35,34,32,38,36,37,49,52,51,54,55,53,50,48,17,20,18,16,21,19,23,22,57,61,59,63,62,58,60,56,
		9,13,11,15,14,10,12, 8,61,59,57,63,60,58,56,62,31,27,25,26,29,28,24,30, 3, 1, 0, 7, 5, 4, 2,6,52,51,49,55,53,54,50,48,20,18,17,21,19,16,22,23,39,35,33,32,38,34,37,36,40,43,41,44,45,42,46,47,
		9,14,11,12, 8,15,10,13,33,39,35,37,36,38,32,34,49,52,51,48,50,53,55,54,57,61,59,56,60,58,63,62,41,40,43,45,46,47,42,44,17,20,18,22,19,21,16,23, 0, 3, 1, 7, 4, 6, 2,5,25,31,27,26,24,29,28,30,
		9,14,13,10,12,11,8,15,31,27,25,26,28,29,30,24,39,35,33,37,32,38,34,36,61,59,57,60,58,56,63,62,20,18,17,19,21,22,16,23, 3, 1, 0, 4, 6, 7, 5, 2,40,43,41,46,47,45,44,42,52,51,49,48,55,53,54,50,
		9,14,15, 8,10,13,12,11,51,49,52,48,54,53,50,55,27,25,31,26,30,29,24,28,59,57,61,58,56,60,63,62, 1, 0, 3, 6, 7, 4, 5, 2,43,41,40,47,45,46,42,44,18,17,20,21,22,19,23,16,35,33,39,37,34,38,36,32,
		9,15, 8,14,10,12,11,13,52,51,49,48,50,54,53,55,20,18,17,23,22,16,19,21,40,43,41,47,46,42,44,45,31,27,25,24,30,29,28,26,39,35,33,36,34,38,32,37, 3, 1, 0, 5, 7, 6, 2,4,61,59,57,63,62,56,58,60,
		9,15,12,10,11,13,14, 8,17,20,18,23,21,16,22,19,57,61,59,63,60,56,62,58,41,40,43,42,47,46,44,45, 0, 3, 1, 6, 5, 7, 4, 2,25,31,27,29,24,30,28,26,33,39,35,38,36,34,37,32,49,52,51,48,55,54,50,53,
		9,15,13,11,14, 8,10,12,59,57,61,63,58,56,60,62,51,49,52,48,53,54,55,50,43,41,40,46,42,47,44,45,35,33,39,34,38,36,32,37, 1, 0, 3, 7, 6, 5, 4, 2,27,25,31,30,29,24,26,28,18,17,20,23,19,16,21,22,
		10, 8,9,11,14,12,13,15, 0, 2,3, 1, 4, 6, 7, 5,58,60,62,56,61,59,57,63,18,20,23,22,16,19,17,21,42,45,47,43,40,44,41,46,50,54,53,49,48,52,51,55,26,25,27,30,24,28,31,29,34,36,32,38,39,35,33,37,
		10, 8,12,14,13,15,11,9,62,58,60,56,63,59,61,57,32,34,36,38,37,35,39,33,23,18,20,19,22,16,17,21,27,26,25,28,30,24,29,31,47,42,45,44,43,40,41,46,53,50,54,52,49,48,55,51,3, 0, 2,1, 5, 6, 4, 7,
		10, 8,15,13,11,9,14,12,36,32,34,38,33,35,37,39, 2,3, 0, 1, 7, 6, 5, 4,20,23,18,16,19,22,17,21,54,53,50,48,52,49,51,55,25,27,26,24,28,30,29,31,45,47,42,40,44,43,46,41,60,62,58,56,57,59,63,61,
		10, 9,11,8,14,13,15,12, 3, 0, 2,1, 7, 4, 6, 5,27,26,25,31,24,29,28,30,53,50,54,49,52,51,55,48,62,58,60,57,61,59,63,56,32,34,36,33,39,35,37,38,47,42,45,41,40,43,46,44,23,18,20,17,21,16,22,19,
		10, 9,12,15, 8,11,14,13,18,20,23,17,22,16,19,21,0, 2,3, 1, 6, 4, 5, 7,50,54,53,52,51,49,55,48,34,36,32,39,35,33,37,38,42,45,47,40,43,41,44,46,58,60,62,61,59,57,56,63,26,25,27,31,28,29,30,24,
		10, 9,13,14,15,12, 8,11,25,27,26,31,30,29,24,28,20,23,18,17,19,16,21,22,54,53,50,51,49,52,55,48,45,47,42,43,41,40,44,46,60,62,58,59,57,61,63,56,36,32,34,35,33,39,38,37, 2,3, 0, 1, 5, 4, 7, 6,
		10,11,8,9,14,15,12,13, 2,3, 0, 1, 6, 7, 4, 5,45,47,42,46,40,44,43,41,36,32,34,33,35,37,38,39,25,27,26,28,24,29,30,31,20,23,18,22,21,16,19,17,60,62,58,63,61,57,56,59,54,53,50,55,48,52,49,51,
		10,11,13,12, 9, 8,14,15,53,50,54,55,49,52,51,48,3, 0, 2,1, 4, 7, 5, 6,32,34,36,35,37,33,38,39,23,18,20,21,16,22,19,17,62,58,60,61,57,63,59,56,27,26,25,24,29,28,31,30,47,42,45,46,43,44,41,40,
		10,11,15,14,12,13, 9, 8,42,45,47,46,41,44,40,43,50,54,53,55,51,52,48,49,34,36,32,37,33,35,38,39,58,60,62,57,63,61,59,56,26,25,27,29,28,24,30,31,18,20,23,16,22,21,17,19, 0, 2,3, 1, 5, 7, 6, 4,
		10,12,11,13, 9,15, 8,14,50,54,53,55,52,51,49,48,18,20,23,17,16,22,21,19,42,45,47,41,44,40,46,43,26,25,27,28,29,30,24,31,34,36,32,35,39,37,33,38,0, 2,3, 6, 4, 5, 1, 7,58,60,62,56,57,63,61,59,
		10,12,14, 8,13,11,9,15,60,62,58,56,61,63,59,57,54,53,50,55,49,51,48,52,45,47,42,44,40,41,46,43,36,32,34,39,37,35,33,38,2,3, 0, 4, 5, 6, 7, 1,25,27,26,29,30,28,31,24,20,23,18,17,21,22,19,16,
		10,12,15, 9, 8,14,13,11,23,18,20,17,19,22,16,21,62,58,60,56,59,63,57,61,47,42,45,40,41,44,46,43, 3, 0, 2,5, 6, 4, 7, 1,27,26,25,30,28,29,24,31,32,34,36,37,35,39,38,33,53,50,54,55,48,51,52,49,
		10,13, 8,15,11,12, 9,14,32,34,36,38,35,37,33,39,53,50,54,55,52,49,48,51,62,58,60,63,59,61,56,57,47,42,45,43,44,41,40,46,23,18,20,16,21,19,22,17, 3, 0, 2,4, 7, 5, 1, 6,27,26,25,31,28,30,24,29,
		10,13,12,11,9,14,15, 8,54,53,50,55,51,49,52,48,25,27,26,31,29,30,28,24,60,62,58,61,63,59,56,57, 2,3, 0, 5, 4, 7, 6, 1,45,47,42,41,43,44,40,46,20,23,18,19,16,21,17,22,36,32,34,38,39,37,35,33,
		10,13,14, 9,15, 8,11,12,26,25,27,31,24,30,29,28,34,36,32,38,33,37,39,35,58,60,62,59,61,63,56,57,18,20,23,21,19,16,22,17, 0, 2,3, 7, 5, 4, 6, 1,42,45,47,44,41,43,46,40,50,54,53,55,48,49,51,52,
		10,14, 8,12,13, 9,15,11,58,60,62,56,59,61,63,57,26,25,27,31,30,24,28,29, 0, 2,3, 4, 6, 7, 1, 5,50,54,53,48,49,51,52,55,18,20,23,19,21,22,16,17,34,36,32,33,37,39,38,35,42,45,47,46,43,40,44,41,
		10,14, 9,13,15,11,12, 8,27,26,25,31,29,24,30,28,47,42,45,46,41,40,43,44, 3, 0, 2,7, 4, 6, 1, 5,32,34,36,39,33,37,35,38,53,50,54,51,48,49,52,55,23,18,20,22,19,21,17,16,62,58,60,56,57,61,59,63,
		10,14,11,15,12, 8,13, 9,45,47,42,46,44,40,41,43,60,62,58,56,63,61,57,59, 2,3, 0, 6, 7, 4, 1, 5,20,23,18,21,22,19,16,17,36,32,34,37,39,33,35,38,54,53,50,49,51,48,55,52,25,27,26,31,28,24,29,30,
		10,15, 9,12, 8,13,11,14,20,23,18,17,16,19,22,21,36,32,34,38,35,33,39,37,25,27,26,30,29,24,31,28,60,62,58,57,59,63,61,56,54,53,50,52,48,51,49,55, 2,3, 0, 7, 6, 5, 1, 4,45,47,42,46,43,41,40,44,
		10,15,13, 8,11,14,12, 9,34,36,32,38,37,33,35,39,42,45,47,46,44,41,43,40,26,25,27,24,30,29,31,28,0, 2,3, 5, 7, 6, 4, 1,58,60,62,63,57,59,61,56,50,54,53,51,52,48,55,49,18,20,23,17,21,19,16,22,
		10,15,14,11,12, 9, 8,13,47,42,45,46,40,41,44,43,23,18,20,17,22,19,21,16,27,26,25,29,24,30,31,28,53,50,54,48,51,52,49,55, 3, 0, 2,6, 5, 7, 4, 1,62,58,60,59,63,57,56,61,32,34,36,38,39,33,37,35,
		11,8,10, 9,15,12,14,13, 2,1, 3, 0, 6, 5, 7, 4,24,29,28,30,27,26,25,31,50,53,55,52,49,48,54,51,61,57,63,58,62,56,60,59,35,37,33,36,38,32,34,39,46,45,42,44,43,40,47,41,21,16,22,19,23,18,20,17,
		11,8,12,15,14,13, 9,10,28,24,29,30,31,26,27,25,22,21,16,19,17,18,23,20,55,50,53,48,52,49,54,51,42,46,45,40,44,43,41,47,63,61,57,56,58,62,60,59,33,35,37,32,36,38,39,34, 3, 2,1, 0, 4, 5, 6, 7,
		11,8,13,14, 9,10,15,12,16,22,21,19,20,18,17,23, 1, 3, 2,0, 7, 5, 4, 6,53,55,50,49,48,52,54,51,37,33,35,38,32,36,34,39,45,42,46,43,40,44,41,47,57,63,61,62,56,58,59,60,29,28,24,30,25,26,31,27,
		11,9, 8,10,15,13,12,14, 1, 3, 2,0, 5, 7, 6, 4,57,63,61,59,62,56,58,60,16,22,21,20,18,17,19,23,45,42,46,40,43,41,44,47,53,55,50,52,51,49,48,54,29,28,24,31,27,25,30,26,37,33,35,39,38,32,36,34,
		11,9,13,15,12,14,10, 8,61,57,63,59,60,56,62,58,35,37,33,39,34,32,38,36,21,16,22,17,20,18,19,23,24,29,28,25,31,27,26,30,46,45,42,41,40,43,44,47,50,53,55,49,52,51,54,48,2,1, 3, 0, 4, 7, 5, 6,
		11,9,14,12,10, 8,15,13,33,35,37,39,36,32,34,38,3, 2,1, 0, 6, 7, 4, 5,22,21,16,18,17,20,19,23,55,50,53,51,49,52,48,54,28,24,29,27,25,31,26,30,42,46,45,43,41,40,47,44,63,61,57,59,58,56,60,62,
		11,10, 9, 8,15,14,13,12, 3, 2,1, 0, 7, 6, 5, 4,42,46,45,47,43,41,40,44,33,35,37,36,32,34,39,38,28,24,29,25,27,26,31,30,22,21,16,20,23,18,17,19,63,61,57,60,62,58,59,56,55,50,53,54,51,49,52,48,
		11,10,12,13, 8,9,15,14,50,53,55,54,52,49,48,51,2,1, 3, 0, 5, 6, 4, 7,35,37,33,32,34,36,39,38,21,16,22,23,18,20,17,19,61,57,63,62,58,60,56,59,24,29,28,27,26,25,30,31,46,45,42,47,40,41,44,43,
		11,10,14,15,13,12, 8,9,45,42,46,47,44,41,43,40,53,55,50,54,48,49,51,52,37,33,35,34,36,32,39,38,57,63,61,58,60,62,56,59,29,28,24,26,25,27,31,30,16,22,21,18,20,23,19,17, 1, 3, 2,0, 4, 6, 7, 5,
		11,12, 9,14,10,13, 8,15,35,37,33,39,32,34,36,38,50,53,55,54,49,52,51,48,61,57,63,60,56,62,59,58,46,45,42,40,41,44,43,47,21,16,22,18,23,17,20,19, 2,1, 3, 5, 6, 4, 0, 7,24,29,28,30,25,31,27,26,
		11,12,13,10, 8,15,14, 9,55,50,53,54,48,52,49,51,28,24,29,30,26,31,25,27,63,61,57,62,60,56,59,58,3, 2,1, 4, 5, 6, 7, 0,42,46,45,44,40,41,43,47,22,21,16,17,18,23,19,20,33,35,37,39,38,34,32,36,
		11,12,15, 8,14, 9,10,13,29,28,24,30,27,31,26,25,37,33,35,39,36,34,38,32,57,63,61,56,62,60,59,58,16,22,21,23,17,18,20,19, 1, 3, 2,6, 4, 5, 7, 0,45,42,46,41,44,40,47,43,53,55,50,54,51,52,48,49,
		11,13,10,12, 8,14, 9,15,53,55,50,54,49,48,52,51,16,22,21,19,18,20,23,17,45,42,46,44,41,43,47,40,29,28,24,25,26,31,27,30,37,33,35,32,38,34,36,39, 1, 3, 2,7, 5, 4, 0, 6,57,63,61,59,58,60,62,56,
		11,13,14, 8,9,15,12,10,21,16,22,19,17,20,18,23,61,57,63,59,56,60,58,62,46,45,42,43,44,41,47,40, 2,1, 3, 4, 7, 5, 6, 0,24,29,28,31,25,26,27,30,35,37,33,34,32,38,39,36,50,53,55,54,51,48,49,52,
		11,13,15, 9,12,10, 8,14,63,61,57,59,62,60,56,58,55,50,53,54,52,48,51,49,42,46,45,41,43,44,47,40,33,35,37,38,34,32,36,39, 3, 2,1, 5, 4, 7, 6, 0,28,24,29,26,31,25,30,27,22,21,16,19,23,20,17,18,
		11,14, 8,13, 9,12,10,15,22,21,16,19,18,17,20,23,33,35,37,39,32,36,38,34,28,24,29,31,26,27,30,25,63,61,57,58,56,60,62,59,55,50,53,49,51,48,52,54, 3, 2,1, 6, 7, 4, 0, 5,42,46,45,47,40,44,43,41,
		11,14,12, 9,10,15,13, 8,37,33,35,39,34,36,32,38,45,42,46,47,41,44,40,43,29,28,24,27,31,26,30,25, 1, 3, 2,4, 6, 7, 5, 0,57,63,61,60,58,56,62,59,53,55,50,48,49,51,54,52,16,22,21,19,23,17,18,20,
		11,14,15,10,13, 8,9,12,46,45,42,47,43,44,41,40,21,16,22,19,20,17,23,18,24,29,28,26,27,31,30,25,50,53,55,51,48,49,52,54, 2,1, 3, 7, 4, 6, 5, 0,61,57,63,56,60,58,59,62,35,37,33,39,38,36,34,32,
		11,15, 8,12,14,10,13, 9,24,29,28,30,26,27,31,25,46,45,42,47,44,43,40,41,2,1, 3, 6, 5, 7, 0, 4,35,37,33,38,36,34,32,39,50,53,55,48,51,52,49,54,21,16,22,20,17,23,19,18,61,57,63,59,58,62,56,60,
		11,15, 9,13,12, 8,14,10,57,63,61,59,56,62,60,58,29,28,24,30,31,27,25,26, 1, 3, 2,5, 7, 6, 0, 4,53,55,50,51,52,48,49,54,16,22,21,17,23,20,18,19,37,33,35,36,34,38,39,32,45,42,46,47,40,43,41,44,
		11,15,10,14,13, 9,12, 8,42,46,45,47,41,43,44,40,63,61,57,59,60,62,58,56, 3, 2,1, 7, 6, 5, 0, 4,22,21,16,23,20,17,18,19,33,35,37,34,38,36,32,39,55,50,53,52,48,51,54,49,28,24,29,30,25,27,26,31,
		12, 8,13, 9,10,14,11,15,44,40,41,43,45,47,42,46,62,60,56,58,61,63,59,57, 4, 5, 6, 0, 1, 2,7, 3,17,18,23,16,19,22,21,20,39,37,35,32,36,38,34,33,50,55,54,48,52,49,53,51,29,30,28,24,26,31,25,27,
		12, 8,14,10,11,15, 9,13,56,62,60,58,57,63,61,59,28,29,30,24,27,31,26,25, 6, 4, 5, 2,0, 1, 7, 3,54,50,55,49,48,52,51,53,23,17,18,22,16,19,21,20,35,39,37,38,32,36,33,34,41,44,40,43,46,47,45,42,
		12, 8,15,11,9,13,10,14,30,28,29,24,25,31,27,26,40,41,44,43,42,47,46,45, 5, 6, 4, 1, 2,0, 7, 3,37,35,39,36,38,32,34,33,55,54,50,52,49,48,51,53,18,23,17,19,22,16,20,21,60,56,62,58,59,63,57,61,
		12, 9, 8,13,10,15,14,11,40,41,44,43,47,42,45,46,18,23,17,20,19,22,16,21,30,28,29,25,31,27,24,26,55,54,50,49,52,51,48,53, 5, 6, 4, 0, 3, 1, 2,7,60,56,62,57,61,59,58,63,37,35,39,33,36,38,32,34,
		12, 9,11,14,13, 8,10,15,35,39,37,33,32,38,34,36,41,44,40,43,45,42,46,47,28,29,30,31,27,25,24,26, 6, 4, 5, 3, 1, 0, 2,7,56,62,60,61,59,57,63,58,54,50,55,52,51,49,53,48,23,17,18,20,16,22,21,19,
		12, 9,15,10,14,11,13, 8,17,18,23,20,21,22,19,16,39,37,35,33,34,38,36,32,29,30,28,27,25,31,24,26,62,60,56,59,57,61,63,58,50,55,54,51,49,52,48,53, 4, 5, 6, 1, 0, 3, 7, 2,44,40,41,43,46,42,47,45,
		12,10, 8,14,11,13,15, 9,62,60,56,58,63,61,57,59,50,55,54,53,48,52,49,51,44,40,41,45,47,42,43,46,39,37,35,36,32,34,38,33, 4, 5, 6, 2,3, 0, 1, 7,29,30,28,25,27,26,24,31,17,18,23,20,16,19,22,21,
		12,10, 9,15,14, 8,11,13,18,23,17,20,22,19,21,16,60,56,62,58,57,61,59,63,40,41,44,47,42,45,43,46, 5, 6, 4, 3, 0, 2,1, 7,30,28,29,27,26,25,31,24,37,35,39,32,34,36,33,38,55,54,50,53,49,52,51,48,
		12,10,13,11,15, 9,14, 8,54,50,55,53,51,52,48,49,23,17,18,20,21,19,16,22,41,44,40,42,45,47,43,46,28,29,30,26,25,27,31,24,35,39,37,34,36,32,38,33, 6, 4, 5, 0, 2,3, 7, 1,56,62,60,58,59,61,63,57,
		12,11,8,15, 9,14,13,10,28,29,30,24,31,27,25,26,35,39,37,33,38,32,36,34,56,62,60,57,63,61,58,59,23,17,18,16,22,21,19,20, 6, 4, 5, 1, 3, 2,0, 7,41,44,40,45,42,46,43,47,54,50,55,53,49,48,52,51,
		12,11,10,13,15, 8,9,14,50,55,54,53,52,48,51,49,29,30,28,24,25,27,26,31,62,60,56,63,61,57,58,59, 4, 5, 6, 3, 2,1, 0, 7,44,40,41,42,46,45,47,43,17,18,23,22,21,16,20,19,39,37,35,33,36,32,34,38,
		12,11,14, 9,13,10,15, 8,37,35,39,33,34,32,38,36,55,54,50,53,51,48,49,52,60,56,62,61,57,63,58,59,40,41,44,46,45,42,47,43,18,23,17,21,16,22,19,20, 5, 6, 4, 2,1, 3, 7, 0,30,28,29,24,26,27,31,25,
		12,13, 9, 8,10,11,15,14,41,44,40,43,42,45,47,46,54,50,55,53,52,51,49,48,35,39,37,32,38,34,33,36,56,62,60,59,61,63,57,58,28,29,30,25,26,31,27,24,23,17,18,21,19,16,20,22, 6, 4, 5, 7, 3, 1, 0, 2,
		12,13,11,10,15,14, 8,9,55,54,50,53,48,51,52,49, 5, 6, 4, 7, 2,1, 3, 0,37,35,39,34,32,38,33,36,18,23,17,16,21,19,22,20,60,56,62,63,59,61,57,58,30,28,29,31,25,26,24,27,40,41,44,43,46,45,42,47,
		12,13,14,15, 8,9,10,11,4, 5, 6, 7, 0, 1, 2,3,44,40,41,43,47,45,46,42,39,37,35,38,34,32,33,36,29,30,28,26,31,25,27,24,17,18,23,19,16,21,22,20,62,60,56,61,63,59,58,57,50,55,54,53,49,51,48,52,
		12,14, 9,11,13,15, 8,10,39,37,35,33,38,34,32,36, 4, 5, 6, 7, 1, 0, 3, 2,17,18,23,21,22,19,20,16,50,55,54,49,51,48,52,53,29,30,28,31,26,27,25,24,44,40,41,47,45,46,43,42,62,60,56,58,59,57,61,63,
		12,14,10, 8,11,9,13,15,60,56,62,58,61,57,63,59,37,35,39,33,32,34,36,38,18,23,17,22,19,21,20,16,30,28,29,26,27,31,25,24,40,41,44,45,46,47,42,43,55,54,50,51,48,49,53,52, 5, 6, 4, 7, 3, 0, 2,1,
		12,14,15,13, 8,10,11,9, 6, 4, 5, 7, 2,0, 1, 3,56,62,60,58,63,57,59,61,23,17,18,19,21,22,20,16,41,44,40,46,47,45,42,43,54,50,55,48,49,51,52,53,28,29,30,27,31,26,24,25,35,39,37,33,36,34,38,32,
		12,15,10, 9,14,13, 8,11,23,17,18,20,19,21,22,16, 6, 4, 5, 7, 0, 2,3, 1,54,50,55,51,52,48,53,49,35,39,37,36,34,38,32,33,41,44,40,47,46,42,45,43,56,62,60,63,57,59,58,61,28,29,30,24,26,25,27,31,
		12,15,11,8,9,10,14,13,29,30,28,24,27,25,31,26,17,18,23,20,22,21,16,19,50,55,54,52,48,51,53,49,44,40,41,46,42,47,45,43,62,60,56,57,59,63,61,58,39,37,35,34,38,36,33,32, 4, 5, 6, 7, 3, 2,1, 0,
		12,15,13,14, 8,11,9,10, 5, 6, 4, 7, 1, 2,0, 3,30,28,29,24,31,25,26,27,55,54,50,48,51,52,53,49,60,56,62,59,63,57,61,58,37,35,39,38,36,34,32,33,40,41,44,42,47,46,43,45,18,23,17,20,16,21,19,22,
		13, 8,9,12,11,14,15,10,43,44,41,40,46,45,42,47,16,21,19,22,17,20,18,23,31,25,26,28,30,24,27,29,54,55,53,52,49,48,51,50, 4, 7, 5, 1, 2,0, 3, 6,63,59,61,58,62,56,57,60,34,32,38,36,33,39,35,37,
		13, 8,10,15,12, 9,11,14,32,38,34,36,35,39,37,33,44,41,43,40,42,45,47,46,25,26,31,30,24,28,27,29, 7, 5, 4, 2,0, 1, 3, 6,59,61,63,62,56,58,60,57,55,53,54,49,48,52,50,51,21,19,16,22,18,20,23,17,
		13, 8,14,11,15,10,12, 9,19,16,21,22,23,20,17,18,38,34,32,36,37,39,33,35,26,31,25,24,28,30,27,29,61,63,59,56,58,62,60,57,53,54,55,48,52,49,51,50, 5, 4, 7, 0, 1, 2,6, 3,41,43,44,40,47,45,46,42,
		13, 9,12, 8,11,15,10,14,41,43,44,40,42,46,45,47,61,63,59,57,62,60,56,58,5, 4, 7, 1, 0, 3, 6, 2,19,16,21,18,17,20,23,22,38,34,32,35,33,39,37,36,53,54,55,51,49,52,50,48,26,31,25,27,29,30,28,24,
		13, 9,14,10, 8,12,11,15,31,25,26,27,28,30,24,29,43,44,41,40,45,46,47,42, 4, 7, 5, 0, 3, 1, 6, 2,34,32,38,33,39,35,37,36,54,55,53,49,52,51,48,50,16,21,19,17,20,18,22,23,63,59,61,57,56,60,58,62,
		13, 9,15,11,10,14, 8,12,59,61,63,57,58,60,62,56,25,26,31,27,24,30,29,28,7, 5, 4, 3, 1, 0, 6, 2,55,53,54,52,51,49,48,50,21,19,16,20,18,17,23,22,32,38,34,39,35,33,36,37,44,41,43,40,47,46,42,45,
		13,10, 9,14, 8,15,12,11,25,26,31,27,30,24,28,29,32,38,34,36,39,35,33,37,59,61,63,58,60,62,57,56,21,19,16,18,20,23,17,22, 7, 5, 4, 0, 2,3, 1, 6,44,41,43,42,45,47,40,46,55,53,54,50,52,51,49,48,
		13,10,11,12,14, 9, 8,15,53,54,55,50,49,51,48,52,26,31,25,27,28,24,29,30,61,63,59,60,62,58,57,56, 5, 4, 7, 2,3, 0, 1, 6,41,43,44,45,47,42,46,40,19,16,21,20,23,18,22,17,38,34,32,36,33,35,37,39,
		13,10,15, 8,12,11,14, 9,34,32,38,36,37,35,39,33,54,55,53,50,48,51,52,49,63,59,61,62,58,60,57,56,43,44,41,47,42,45,46,40,16,21,19,23,18,20,17,22, 4, 7, 5, 3, 0, 2,6, 1,31,25,26,27,29,24,30,28,
		13,11,8,14,15, 9,10,12,16,21,19,22,20,17,23,18,63,59,61,57,58,62,56,60,43,44,41,46,45,42,40,47, 4, 7, 5, 2,1, 3, 0, 6,31,25,26,24,29,28,30,27,34,32,38,35,37,33,36,39,54,55,53,50,52,49,48,51,
		13,11,9,15,10,12,14, 8,61,63,59,57,60,62,58,56,53,54,55,50,51,49,52,48,41,43,44,42,46,45,40,47,38,34,32,33,35,37,39,36, 5, 4, 7, 3, 2,1, 0, 6,26,31,25,28,24,29,27,30,19,16,21,22,18,17,20,23,
		13,11,12,10,14, 8,15, 9,55,53,54,50,48,49,51,52,21,19,16,22,23,17,18,20,44,41,43,45,42,46,40,47,25,26,31,29,28,24,30,27,32,38,34,37,33,35,39,36, 7, 5, 4, 1, 3, 2,6, 0,59,61,63,57,56,62,60,58,
		13,12, 8,9,11,10,14,15,44,41,43,40,45,42,46,47,55,53,54,50,49,48,52,51,32,38,34,35,39,37,36,33,59,61,63,56,62,60,58,57,25,26,31,28,29,30,24,27,21,19,16,23,17,18,22,20, 7, 5, 4, 6, 2,0, 1, 3,
		13,12,10,11,14,15, 9, 8,54,55,53,50,51,48,49,52, 4, 7, 5, 6, 3, 0, 2,1,34,32,38,37,35,39,36,33,16,21,19,18,23,17,20,22,63,59,61,60,56,62,58,57,31,25,26,30,28,29,27,24,43,44,41,40,47,42,45,46,
		13,12,15,14, 9, 8,11,10, 5, 4, 7, 6, 1, 0, 3, 2,41,43,44,40,46,42,47,45,38,34,32,39,37,35,36,33,26,31,25,29,30,28,24,27,19,16,21,17,18,23,20,22,61,63,59,62,60,56,57,58,53,54,55,50,52,48,51,49,
		13,14,10, 9, 8,11,15,12,26,31,25,27,24,28,30,29,19,16,21,22,20,23,18,17,53,54,55,49,51,48,50,52,41,43,44,47,45,46,42,40,61,63,59,58,56,60,62,57,38,34,32,37,39,33,36,35, 5, 4, 7, 6, 2,3, 0, 1,
		13,14,11,8,15,12, 9,10,21,19,16,22,17,23,20,18,7, 5, 4, 6, 1, 3, 2,0,55,53,54,48,49,51,50,52,32,38,34,33,37,39,35,36,44,41,43,46,47,45,42,40,59,61,63,60,58,56,57,62,25,26,31,27,29,28,24,30,
		13,14,12,15, 9,10, 8,11,4, 7, 5, 6, 0, 3, 1, 2,31,25,26,27,30,28,29,24,54,55,53,51,48,49,50,52,63,59,61,56,60,58,62,57,34,32,38,39,33,37,35,36,43,44,41,45,46,47,40,42,16,21,19,22,18,23,17,20,
		13,15, 8,10,12,14, 9,11,38,34,32,36,39,37,35,33, 5, 4, 7, 6, 0, 1, 2,3,19,16,21,23,20,17,22,18,53,54,55,52,48,51,49,50,26,31,25,30,29,24,28,27,41,43,44,46,42,47,40,45,61,63,59,57,56,58,62,60,
		13,15,11,9,10, 8,12,14,63,59,61,57,62,58,60,56,34,32,38,36,35,37,33,39,16,21,19,20,17,23,22,18,31,25,26,29,24,30,28,27,43,44,41,42,47,46,45,40,54,55,53,48,51,52,50,49, 4, 7, 5, 6, 2,1, 3, 0,
		13,15,14,12, 9,11,10, 8,7, 5, 4, 6, 3, 1, 0, 2,59,61,63,57,60,58,56,62,21,19,16,17,23,20,22,18,44,41,43,47,46,42,45,40,55,53,54,51,52,48,49,50,25,26,31,24,30,29,27,28,32,38,34,36,33,37,39,35,
		14, 8,10,12, 9,15,13,11,58,56,60,62,59,57,61,63,49,51,48,52,54,53,50,55,47,45,46,40,44,43,42,41,37,39,33,32,36,38,34,35, 6, 7, 4, 0, 1, 2,3, 5,31,26,27,24,28,30,25,29,21,22,19,16,20,23,18,17,
		14, 8,11,13,12,10, 9,15,22,19,21,16,18,23,17,20,56,60,58,62,61,57,63,59,45,46,47,44,43,40,42,41,7, 4, 6, 1, 2,0, 3, 5,26,27,31,28,30,24,29,25,39,33,37,36,38,32,35,34,51,48,49,52,50,53,55,54,
		14, 8,15, 9,13,11,12,10,48,49,51,52,55,53,54,50,19,21,22,16,17,23,20,18,46,47,45,43,40,44,42,41,27,31,26,30,24,28,29,25,33,37,39,38,32,36,34,35, 4, 6, 7, 2,0, 1, 5, 3,60,58,56,62,63,57,59,61,
		14, 9, 8,15,13,10,11,12,49,51,48,52,53,54,55,50,31,26,27,25,24,28,30,29,58,56,60,59,57,61,62,63, 6, 7, 4, 1, 0, 3, 2,5,47,45,46,43,41,40,44,42,21,22,19,18,17,20,16,23,37,39,33,35,32,36,38,34,
		14, 9,10,13,11,12,15, 8,27,31,26,25,29,28,24,30,33,37,39,35,34,36,32,38,60,58,56,61,59,57,62,63,19,21,22,20,18,17,23,16, 4, 6, 7, 3, 1, 0, 2,5,46,47,45,40,43,41,42,44,48,49,51,52,50,54,53,55,
		14, 9,12,11,15, 8,13,10,39,33,37,35,38,36,34,32,51,48,49,52,55,54,50,53,56,60,58,57,61,59,62,63,45,46,47,41,40,43,44,42,22,19,21,17,20,18,23,16, 7, 4, 6, 0, 3, 1, 5, 2,26,27,31,25,30,28,29,24,
		14,10,12, 8,9,13,11,15,60,58,56,62,61,59,57,63,27,31,26,25,28,29,30,24, 4, 6, 7, 0, 2,3, 5, 1,48,49,51,50,54,53,55,52,19,21,22,18,20,23,17,16,33,37,39,34,36,32,35,38,46,47,45,42,41,44,40,43,
		14,10,13, 9,11,15, 8,12,26,27,31,25,24,29,28,30,45,46,47,42,43,44,41,40, 7, 4, 6, 3, 0, 2,5, 1,39,33,37,32,34,36,38,35,51,48,49,53,50,54,55,52,22,19,21,23,18,20,16,17,56,60,58,62,63,59,61,57,
		14,10,15,11,8,12, 9,13,47,45,46,42,40,44,43,41,58,56,60,62,57,59,63,61,6, 7, 4, 2,3, 0, 5, 1,21,22,19,20,23,18,17,16,37,39,33,36,32,34,38,35,49,51,48,54,53,50,52,55,31,26,27,25,30,29,24,28,
		14,11,9,12,15,10, 8,13,33,37,39,35,36,34,38,32,46,47,45,42,40,43,41,44,27,31,26,29,28,24,25,30, 4, 6, 7, 1, 3, 2,0, 5,60,58,56,57,63,61,59,62,48,49,51,53,55,50,52,54,19,21,22,16,20,18,17,23,
		14,11,10,15, 8,13,12, 9,45,46,47,42,44,43,40,41,22,19,21,16,23,18,20,17,26,27,31,24,29,28,25,30,51,48,49,50,53,55,54,52, 7, 4, 6, 2,1, 3, 0, 5,56,60,58,61,57,63,62,59,39,33,37,35,32,34,36,38,
		14,11,13, 8,12, 9,15,10,21,22,19,16,17,18,23,20,37,39,33,35,38,34,32,36,31,26,27,28,24,29,25,30,58,56,60,63,61,57,59,62,49,51,48,55,50,53,54,52, 6, 7, 4, 3, 2,1, 5, 0,47,45,46,42,41,43,44,40,
		14,12, 8,10, 9,11,15,13,56,60,58,62,57,61,59,63,39,33,37,35,36,38,32,34,22,19,21,18,23,17,16,20,26,27,31,30,28,29,24,25,45,46,47,40,41,44,43,42,51,48,49,55,54,50,52,53, 7, 4, 6, 5, 1, 2,0, 3,
		14,12,11,9,15,13,10, 8,37,39,33,35,34,38,36,32, 6, 7, 4, 5, 3, 2,1, 0,21,22,19,17,18,23,16,20,49,51,48,50,55,54,53,52,31,26,27,29,30,28,24,25,47,45,46,44,40,41,42,43,58,56,60,62,63,61,57,59,
		14,12,13,15,10, 8,9,11,4, 6, 7, 5, 0, 2,3, 1,60,58,56,62,59,61,63,57,19,21,22,23,17,18,16,20,46,47,45,41,44,40,43,42,48,49,51,54,50,55,53,52,27,31,26,28,29,30,25,24,33,37,39,35,32,38,34,36,
		14,13, 8,11,12,15,10, 9,19,21,22,16,23,17,18,20, 4, 6, 7, 5, 2,0, 1, 3,48,49,51,55,53,54,52,50,33,37,39,32,38,34,36,35,46,47,45,44,41,43,40,42,60,58,56,59,61,63,62,57,27,31,26,25,30,24,28,29,
		14,13, 9,10,11,8,12,15,31,26,27,25,28,24,29,30,21,22,19,16,18,17,20,23,49,51,48,53,54,55,52,50,47,45,46,41,43,44,40,42,58,56,60,61,63,59,57,62,37,39,33,38,34,32,35,36, 6, 7, 4, 5, 1, 0, 3, 2,
		14,13,15,12,10, 9,11,8,7, 4, 6, 5, 3, 0, 2,1,26,27,31,25,29,24,30,28,51,48,49,54,55,53,52,50,56,60,58,63,59,61,57,62,39,33,37,34,32,38,36,35,45,46,47,43,44,41,42,40,22,19,21,16,20,17,23,18,
		14,15, 9, 8,13,12,10,11,51,48,49,52,54,55,53,50, 7, 4, 6, 5, 0, 3, 1, 2,39,33,37,38,36,34,35,32,22,19,21,20,17,23,18,16,56,60,58,59,63,57,61,62,26,27,31,29,24,30,25,28,45,46,47,42,41,40,43,44,
		14,15,11,10, 8,9,13,12,46,47,45,42,43,40,44,41,48,49,51,52,53,55,50,54,33,37,39,36,34,38,35,32,60,58,56,63,57,59,61,62,27,31,26,24,30,29,28,25,19,21,22,17,23,20,16,18,4, 6, 7, 5, 1, 3, 2,0,
		14,15,12,13,10,11,8,9, 6, 7, 4, 5, 2,3, 0, 1,47,45,46,42,44,40,41,43,37,39,33,34,38,36,35,32,31,26,27,30,29,24,28,25,21,22,19,23,20,17,18,16,58,56,60,57,59,63,62,61,49,51,48,52,50,55,54,53,
		15, 8,9,14,12,11,10,13,52,48,51,49,50,55,54,53,30,29,24,28,27,25,31,26,57,59,63,56,58,62,61,60, 7, 6, 5, 0, 1, 2,3, 4,46,42,47,40,44,43,41,45,23,20,17,16,19,22,18,21,34,38,36,32,35,33,39,37,
		15, 8,11,12,10,13,14, 9,24,30,29,28,26,25,27,31,36,34,38,32,37,33,35,39,63,57,59,62,56,58,61,60,17,23,20,22,16,19,21,18,5, 7, 6, 2,0, 1, 3, 4,47,46,42,43,40,44,45,41,51,52,48,49,53,55,50,54,
		15, 8,13,10,14, 9,12,11,38,36,34,32,39,33,37,35,48,51,52,49,54,55,53,50,59,63,57,58,62,56,61,60,42,47,46,44,43,40,41,45,20,17,23,19,22,16,21,18,6, 5, 7, 1, 2,0, 4, 3,29,24,30,28,31,25,26,27,
		15, 9,10,12,13,11,8,14,20,17,23,18,16,21,19,22,59,63,57,61,62,58,60,56,42,47,46,41,40,43,45,44, 6, 5, 7, 0, 3, 1, 2,4,29,24,30,25,31,27,26,28,38,36,34,33,39,35,32,37,48,51,52,49,53,50,54,55,
		15, 9,11,13, 8,14,12,10,57,59,63,61,56,58,62,60,52,48,51,49,55,50,53,54,46,42,47,43,41,40,45,44,34,38,36,35,33,39,37,32, 7, 6, 5, 1, 0, 3, 2,4,30,29,24,27,25,31,28,26,23,20,17,18,22,21,16,19,
		15, 9,14, 8,12,10,13,11,51,52,48,49,54,50,55,53,17,23,20,18,19,21,22,16,47,46,42,40,43,41,45,44,24,30,29,31,27,25,26,28,36,34,38,39,35,33,37,32, 5, 7, 6, 3, 1, 0, 4, 2,63,57,59,61,60,58,56,62,
		15,10, 8,13,14,11,9,12,36,34,38,32,33,37,39,35,47,46,42,45,43,40,44,41,24,30,29,26,25,27,28,31,5, 7, 6, 0, 2,3, 1, 4,63,57,59,58,60,62,56,61,51,52,48,50,54,53,49,55,17,23,20,18,22,16,19,21,
		15,10,11,14, 9,12,13, 8,42,47,46,45,41,40,43,44,20,17,23,18,21,16,22,19,29,24,30,27,26,25,28,31,48,51,52,53,50,54,55,49, 6, 5, 7, 3, 0, 2,1, 4,59,63,57,62,58,60,61,56,38,36,34,32,35,37,33,39,
		15,10,12, 9,13, 8,14,11,23,20,17,18,19,16,21,22,34,38,36,32,39,37,35,33,30,29,24,25,27,26,28,31,57,59,63,60,62,58,56,61,52,48,51,54,53,50,55,49, 7, 6, 5, 2,3, 0, 4, 1,46,42,47,45,44,40,41,43,
		15,11,12, 8,10,14, 9,13,29,24,30,28,27,26,25,31,42,47,46,45,40,41,44,43, 6, 5, 7, 2,1, 3, 4, 0,38,36,34,35,37,33,39,32,48,51,52,50,53,55,54,49,20,17,23,21,16,22,18,19,59,63,57,61,60,56,62,58,
		15,11,13, 9, 8,12,10,14,63,57,59,61,62,56,58,60,24,30,29,28,25,26,31,27, 5, 7, 6, 1, 3, 2,4, 0,51,52,48,53,55,50,54,49,17,23,20,16,22,21,19,18,36,34,38,37,33,35,32,39,47,46,42,45,44,41,43,40,
		15,11,14,10, 9,13, 8,12,46,42,47,45,43,41,40,44,57,59,63,61,58,56,60,62, 7, 6, 5, 3, 2,1, 4, 0,23,20,17,22,21,16,19,18,34,38,36,33,35,37,39,32,52,48,51,55,50,53,49,54,30,29,24,28,31,26,27,25,
		15,12, 8,11,10, 9,13,14,30,29,24,28,25,27,26,31,23,20,17,18,16,19,22,21,52,48,51,50,55,54,49,53,46,42,47,44,40,41,43,45,57,59,63,62,60,56,58,61,34,38,36,39,37,35,32,33, 7, 6, 5, 4, 0, 1, 2,3,
		15,12, 9,10,13,14,11,8,17,23,20,18,21,19,16,22, 5, 7, 6, 4, 3, 1, 0, 2,51,52,48,54,50,55,49,53,36,34,38,35,39,37,33,32,47,46,42,41,44,40,43,45,63,57,59,56,62,60,61,58,24,30,29,28,31,27,25,26,
		15,12,14,13,11,8,10, 9, 6, 5, 7, 4, 2,1, 3, 0,29,24,30,28,26,27,31,25,48,51,52,55,54,50,49,53,59,63,57,60,56,62,58,61,38,36,34,37,35,39,33,32,42,47,46,40,41,44,45,43,20,17,23,18,22,19,21,16,
		15,13, 9,11,8,10,14,12,59,63,57,61,58,62,56,60,38,36,34,32,33,39,35,37,20,17,23,16,21,19,18,22,29,24,30,31,25,26,27,28,42,47,46,43,44,41,40,45,48,51,52,54,55,53,49,50, 6, 5, 7, 4, 0, 3, 1, 2,
		15,13,10, 8,14,12,11,9,34,38,36,32,37,39,33,35, 7, 6, 5, 4, 2,3, 0, 1,23,20,17,19,16,21,18,22,52,48,51,53,54,55,50,49,30,29,24,26,31,25,27,28,46,42,47,41,43,44,45,40,57,59,63,61,60,62,58,56,
		15,13,12,14,11,9, 8,10, 5, 7, 6, 4, 1, 3, 2,0,63,57,59,61,56,62,60,58,17,23,20,21,19,16,18,22,47,46,42,44,41,43,40,45,51,52,48,55,53,54,50,49,24,30,29,25,26,31,28,27,36,34,38,32,35,39,37,33,
		15,14, 8,9,12,13,11,10,48,51,52,49,55,54,50,53, 6, 5, 7, 4, 1, 2,0, 3,38,36,34,39,33,37,32,35,20,17,23,22,19,21,16,18,59,63,57,56,60,58,62,61,29,24,30,26,27,31,28,25,42,47,46,45,44,43,40,41,
		15,14,10,11,9, 8,12,13,47,46,42,45,40,43,41,44,51,52,48,49,50,54,53,55,36,34,38,33,37,39,32,35,63,57,59,60,58,56,62,61,24,30,29,27,31,26,25,28,17,23,20,19,21,22,18,16, 5, 7, 6, 4, 0, 2,3, 1,
		15,14,13,12,11,10, 9, 8,7, 6, 5, 4, 3, 2,1, 0,46,42,47,45,41,43,44,40,34,38,36,37,39,33,32,35,30,29,24,31,26,27,25,28,23,20,17,21,22,19,16,18,57,59,63,58,56,60,61,62,52,48,51,49,53,54,55,50,
		16,17,18,19,21,20,23,22, 0, 4, 1, 5, 6, 2,7, 3,24,31,29,25,27,30,26,28,40,43,45,44,41,47,42,46, 8,13,11,12, 9,15,10,14,48,54,53,50,49,51,55,52,32,34,37,39,33,36,35,38,56,63,57,62,59,60,58,61,
		16,17,20,21,23,22,19,18,29,24,31,25,28,30,27,26,57,56,63,62,61,60,59,58,45,40,43,47,44,41,42,46,37,32,34,36,39,33,38,35,11,8,13,15,12, 9,10,14,53,48,54,51,50,49,52,55, 1, 0, 4, 5, 3, 2,6, 7,
		16,17,22,23,19,18,21,20,63,57,56,62,58,60,61,59, 4, 1, 0, 5, 7, 2,3, 6,43,45,40,41,47,44,42,46,54,53,48,49,51,50,55,52,34,37,32,33,36,39,38,35,13,11,8,9,15,12,14,10,31,29,24,25,26,30,28,27,
		16,18,19,17,21,23,22,20, 1, 0, 4, 5, 7, 6, 2,3,37,32,34,35,33,38,36,39,53,48,54,50,51,55,52,49,29,24,31,26,27,30,28,25,57,56,63,58,59,60,61,62,11,8,13,10, 9,12,14,15,45,40,43,42,46,41,44,47,
		16,18,20,22,17,19,21,23,40,43,45,42,44,41,47,46, 0, 4, 1, 5, 2,6, 3, 7,48,54,53,51,55,50,52,49,56,63,57,59,60,58,61,62, 8,13,11,9,12,10,15,14,24,31,29,27,30,26,25,28,32,34,37,35,36,38,39,33,
		16,18,23,21,22,20,17,19,34,37,32,35,39,38,33,36,43,45,40,42,47,41,46,44,54,53,48,55,50,51,52,49,13,11,8,12,10, 9,15,14,31,29,24,30,26,27,28,25,63,57,56,60,58,59,62,61,4, 1, 0, 5, 3, 6, 7, 2,
		16,19,17,18,21,22,20,23, 4, 1, 0, 5, 2,7, 6, 3,13,11,8,14, 9,15,12,10,63,57,56,58,60,61,62,59,34,37,32,36,33,38,39,35,43,45,40,44,46,41,47,42,31,29,24,28,27,26,25,30,54,53,48,52,49,51,50,55,
		16,19,22,21,20,23,18,17, 8,13,11,14,10,15, 9,12,48,54,53,52,55,51,49,50,56,63,57,61,58,60,62,59,24,31,29,26,28,27,30,25,32,34,37,38,36,33,39,35,40,43,45,41,44,46,42,47, 0, 4, 1, 5, 3, 7, 2,6,
		16,19,23,20,18,17,21,22,53,48,54,52,50,51,55,49, 1, 0, 4, 5, 6, 7, 3, 2,57,56,63,60,61,58,62,59,45,40,43,46,41,44,47,42,29,24,31,27,26,28,30,25,37,32,34,33,38,36,35,39,11,8,13,14,12,15,10, 9,
		16,20,19,23,18,22,17,21,48,54,53,52,51,55,50,49,40,43,45,42,41,44,46,47, 8,13,11,10,15, 9,14,12,32,34,37,36,38,39,33,35,56,63,57,60,59,61,58,62, 0, 4, 1, 2,6, 3, 5, 7,24,31,29,25,26,28,27,30,
		16,20,21,17,23,19,18,22,31,29,24,25,27,28,30,26,54,53,48,52,50,55,49,51,13,11,8,15, 9,10,14,12,63,57,56,59,61,60,58,62, 4, 1, 0, 6, 3, 2,7, 5,34,37,32,38,39,36,35,33,43,45,40,42,46,44,47,41,
		16,20,22,18,17,21,23,19,45,40,43,42,47,44,41,46,29,24,31,25,30,28,26,27,11,8,13, 9,10,15,14,12, 1, 0, 4, 3, 2,6, 7, 5,37,32,34,39,36,38,33,35,57,56,63,61,60,59,62,58,53,48,54,52,49,55,51,50,
		16,21,17,20,23,18,22,19,24,31,29,25,30,27,28,26,32,34,37,35,39,33,36,38,0, 4, 1, 6, 2,7, 5, 3,48,54,53,49,50,55,51,52,40,43,45,47,46,44,41,42,56,63,57,58,61,59,62,60, 8,13,11,14,12, 9,15,10,
		16,21,18,23,22,19,20,17,37,32,34,35,38,33,39,36,11,8,13,14,10, 9,12,15, 1, 0, 4, 7, 6, 2,5, 3,57,56,63,59,58,61,60,62,53,48,54,55,49,50,51,52,45,40,43,44,47,46,42,41,29,24,31,25,26,27,30,28,
		16,21,19,22,20,17,23,18,13,11,8,14,15, 9,10,12,31,29,24,25,28,27,26,30, 4, 1, 0, 2,7, 6, 5, 3,43,45,40,46,44,47,41,42,63,57,56,61,59,58,60,62,54,53,48,50,55,49,52,51,34,37,32,35,36,33,38,39,
		16,22,18,20,17,23,19,21,43,45,40,42,41,47,44,46,63,57,56,62,60,58,59,61,34,37,32,39,38,33,35,36,31,29,24,26,30,28,27,25,54,53,48,51,49,55,50,52, 4, 1, 0, 7, 2,3, 5, 6,13,11,8,14,12,10, 9,15,
		16,22,21,19,20,18,17,23,11,8,13,14, 9,10,15,12,45,40,43,42,44,47,46,41,37,32,34,38,33,39,35,36,53,48,54,49,55,51,50,52, 1, 0, 4, 2,3, 7, 6, 5,29,24,31,30,28,26,25,27,57,56,63,62,59,58,61,60,
		16,22,23,17,19,21,20,18,56,63,57,62,61,58,60,59, 8,13,11,14,15,10,12, 9,32,34,37,33,39,38,35,36, 0, 4, 1, 3, 7, 2,6, 5,24,31,29,28,26,30,27,25,48,54,53,55,51,49,52,50,40,43,45,42,46,47,41,44,
		16,23,17,22,19,20,18,21,57,56,63,62,60,61,58,59,53,48,54,52,51,50,49,55,29,24,31,28,30,27,25,26,11,8,13,12,15,10, 9,14,45,40,43,41,46,47,44,42, 1, 0, 4, 6, 7, 3, 5, 2,37,32,34,35,36,39,33,38,
		16,23,20,19,18,21,22,17,54,53,48,52,55,50,51,49,34,37,32,35,38,39,36,33,31,29,24,27,28,30,25,26, 4, 1, 0, 3, 6, 7, 2,5,13,11,8,10,12,15, 9,14,43,45,40,47,41,46,42,44,63,57,56,62,59,61,60,58,
		16,23,21,18,22,17,19,20,32,34,37,35,33,39,38,36,56,63,57,62,58,61,59,60,24,31,29,30,27,28,25,26,40,43,45,46,47,41,44,42, 0, 4, 1, 7, 3, 6, 2,5, 8,13,11,15,10,12,14, 9,48,54,53,52,49,50,55,51,
		17,16,19,18,20,21,22,23, 4, 0, 5, 1, 2,6, 3, 7,29,25,24,31,26,28,27,30,44,41,47,40,43,45,46,42,12, 9,15, 8,13,11,14,10,50,49,51,48,54,53,52,55,39,33,36,32,34,37,38,35,62,57,63,56,61,58,60,59,
		17,16,21,20,22,23,18,19,24,29,25,31,30,28,26,27,63,62,57,56,59,58,61,60,47,44,41,45,40,43,46,42,36,39,33,37,32,34,35,38,15,12, 9,11,8,13,14,10,51,50,49,53,48,54,55,52, 5, 4, 0, 1, 7, 6, 2,3,
		17,16,23,22,18,19,20,21,57,63,62,56,60,58,59,61,0, 5, 4, 1, 3, 6, 7, 2,41,47,44,43,45,40,46,42,49,51,50,54,53,48,52,55,33,36,39,34,37,32,35,38,9,15,12,13,11,8,10,14,25,24,29,31,27,28,30,26,
		17,18,16,19,20,23,21,22, 0, 5, 4, 1, 6, 3, 2,7, 9,15,12,10,13,11,8,14,57,63,62,60,58,59,56,61,33,36,39,37,34,35,32,38,41,47,44,40,42,43,45,46,25,24,29,30,26,27,31,28,49,51,50,55,54,53,48,52,
		17,18,22,21,19,16,20,23,51,50,49,55,48,53,52,54, 5, 4, 0, 1, 2,3, 7, 6,63,62,57,58,59,60,56,61,47,44,41,42,43,40,45,46,24,29,25,26,27,30,28,31,36,39,33,34,35,37,38,32,15,12, 9,10, 8,11,14,13,
		17,18,23,20,21,22,19,16,12, 9,15,10,14,11,13, 8,50,49,51,55,52,53,54,48,62,57,63,59,60,58,56,61,29,25,24,27,30,26,28,31,39,33,36,35,37,34,32,38,44,41,47,43,40,42,46,45, 4, 0, 5, 1, 7, 3, 6, 2,
		17,19,18,16,20,22,23,21,5, 4, 0, 1, 3, 2,6, 7,36,39,33,38,34,35,37,32,51,50,49,48,53,52,55,54,24,29,25,27,26,28,30,31,63,62,57,60,61,58,59,56,15,12, 9,14,13, 8,10,11,47,44,41,46,42,43,40,45,
		17,19,21,23,16,18,20,22,44,41,47,46,40,43,45,42, 4, 0, 5, 1, 6, 2,7, 3,50,49,51,53,52,48,55,54,62,57,63,61,58,60,59,56,12, 9,15,13, 8,14,11,10,29,25,24,26,28,27,31,30,39,33,36,38,37,35,32,34,
		17,19,22,20,23,21,16,18,33,36,39,38,32,35,34,37,41,47,44,46,45,43,42,40,49,51,50,52,48,53,55,54, 9,15,12, 8,14,13,11,10,25,24,29,28,27,26,30,31,57,63,62,58,60,61,56,59, 0, 5, 4, 1, 7, 2,3, 6,
		17,20,16,21,22,19,23,18,29,25,24,31,28,26,30,27,39,33,36,38,32,34,37,35, 4, 0, 5, 2,6, 3, 1, 7,50,49,51,54,48,52,53,55,44,41,47,45,42,40,43,46,62,57,63,60,59,61,56,58,12, 9,15,10, 8,13,11,14,
		17,20,18,23,21,16,22,19, 9,15,12,10,11,13,14, 8,25,24,29,31,30,26,27,28,0, 5, 4, 6, 3, 2,1, 7,41,47,44,42,40,45,43,46,57,63,62,59,61,60,58,56,49,51,50,48,52,54,55,53,33,36,39,38,37,34,35,32,
		17,20,19,22,23,18,21,16,36,39,33,38,35,34,32,37,15,12, 9,10,14,13, 8,11,5, 4, 0, 3, 2,6, 1, 7,63,62,57,61,60,59,58,56,51,50,49,52,54,48,53,55,47,44,41,40,45,42,46,43,24,29,25,31,27,26,28,30,
		17,21,18,22,19,23,16,20,50,49,51,55,53,52,48,54,44,41,47,46,43,40,42,45,12, 9,15,14,11,13,10, 8,39,33,36,37,35,32,34,38,62,57,63,58,61,59,60,56, 4, 0, 5, 6, 2,7, 1, 3,29,25,24,31,27,30,26,28,
		17,21,20,16,22,18,19,23,25,24,29,31,26,30,28,27,49,51,50,55,48,52,54,53, 9,15,12,11,13,14,10, 8,57,63,62,61,59,58,60,56, 0, 5, 4, 2,7, 6, 3, 1,33,36,39,35,32,37,38,34,41,47,44,46,42,40,45,43,
		17,21,23,19,16,20,22,18,47,44,41,46,45,40,43,42,24,29,25,31,28,30,27,26,15,12, 9,13,14,11,10, 8,5, 4, 0, 7, 6, 2,3, 1,36,39,33,32,37,35,34,38,63,62,57,59,58,61,56,60,51,50,49,55,54,52,53,48,
		17,22,16,23,18,21,19,20,63,62,57,56,58,59,60,61,51,50,49,55,53,48,54,52,24,29,25,30,28,26,31,27,15,12, 9, 8,11,14,13,10,47,44,41,43,42,45,40,46, 5, 4, 0, 2,3, 7, 1, 6,36,39,33,38,37,32,34,35,
		17,22,20,19,23,16,18,21,39,33,36,38,34,32,35,37,62,57,63,56,60,59,61,58,29,25,24,28,26,30,31,27,44,41,47,42,45,43,40,46, 4, 0, 5, 3, 7, 2,6, 1,12, 9,15,11,14, 8,10,13,50,49,51,55,54,48,52,53,
		17,22,21,18,19,20,23,16,49,51,50,55,52,48,53,54,33,36,39,38,35,32,37,34,25,24,29,26,30,28,31,27, 0, 5, 4, 7, 2,3, 6, 1, 9,15,12,14, 8,11,13,10,41,47,44,45,43,42,46,40,57,63,62,56,61,59,58,60,
		17,23,19,21,16,22,18,20,41,47,44,46,43,45,40,42,57,63,62,56,58,60,61,59,33,36,39,32,35,34,38,37,25,24,29,27,28,30,26,31,49,51,50,53,54,52,48,55, 0, 5, 4, 3, 6, 7, 1, 2,9,15,12,10, 8,14,13,11,
		17,23,20,18,21,19,16,22,15,12, 9,10,13,14,11,8,47,44,41,46,40,45,42,43,36,39,33,35,34,32,38,37,51,50,49,54,52,53,48,55, 5, 4, 0, 6, 7, 3, 2,1,24,29,25,28,30,27,31,26,63,62,57,56,61,60,59,58,
		17,23,22,16,18,20,21,19,62,57,63,56,59,60,58,61,12, 9,15,10,11,14, 8,13,39,33,36,34,32,35,38,37, 4, 0, 5, 7, 3, 6, 2,1,29,25,24,30,27,28,26,31,50,49,51,52,53,54,55,48,44,41,47,46,42,45,43,40,
		18,16,17,19,23,21,20,22, 0, 1, 5, 4, 6, 7, 3, 2,34,35,37,32,36,39,33,38,50,51,55,53,48,54,49,52,26,27,30,29,24,31,25,28,58,59,60,57,56,63,62,61,10, 9,12,11,8,13,15,14,42,43,40,45,47,44,41,46,
		18,16,21,23,20,22,19,17,37,34,35,32,38,39,36,33,40,42,43,45,46,44,47,41,55,50,51,54,53,48,49,52,12,10, 9,13,11,8,14,15,30,26,27,31,29,24,25,28,60,58,59,63,57,56,61,62, 5, 0, 1, 4, 2,7, 6, 3,
		18,16,22,20,19,17,23,21,43,40,42,45,41,44,46,47, 1, 5, 0, 4, 3, 7, 2,6,51,55,50,48,54,53,49,52,59,60,58,56,63,57,62,61,9,12,10, 8,13,11,14,15,27,30,26,24,31,29,28,25,35,37,34,32,33,39,38,36,
		18,17,19,16,23,20,22,21,5, 0, 1, 4, 3, 6, 7, 2,12,10, 9,15, 8,14,13,11,60,58,59,57,63,62,61,56,37,34,35,33,36,39,38,32,40,42,43,41,47,44,46,45,30,26,27,25,24,29,28,31,55,50,51,49,52,48,53,54,
		18,17,20,23,22,21,16,19, 9,12,10,15,11,14, 8,13,51,55,50,49,54,48,52,53,59,60,58,62,57,63,61,56,27,30,26,29,25,24,31,28,35,37,34,39,33,36,38,32,43,40,42,44,41,47,45,46, 1, 5, 0, 4, 2,6, 3, 7,
		18,17,21,22,16,19,23,20,50,51,55,49,53,48,54,52, 0, 1, 5, 4, 7, 6, 2,3,58,59,60,63,62,57,61,56,42,43,40,47,44,41,46,45,26,27,30,24,29,25,31,28,34,35,37,36,39,33,32,38,10, 9,12,15,13,14,11,8,
		18,19,16,17,23,22,21,20, 1, 5, 0, 4, 7, 3, 6, 2,27,30,26,28,24,31,29,25,43,40,42,41,44,46,45,47, 9,12,10,13, 8,14,11,15,51,55,50,53,52,48,54,49,35,37,34,38,36,33,32,39,59,60,58,61,56,63,57,62,
		18,19,20,21,17,16,23,22,60,58,59,61,57,63,62,56, 5, 0, 1, 4, 6, 3, 2,7,40,42,43,44,46,41,45,47,55,50,51,52,48,53,54,49,37,34,35,36,33,38,39,32,12,10, 9, 8,14,13,15,11,30,26,27,28,29,31,25,24,
		18,19,22,23,21,20,17,16,26,27,30,28,25,31,24,29,58,59,60,61,62,63,56,57,42,43,40,46,41,44,45,47,34,35,37,33,38,36,39,32,10, 9,12,14,13, 8,11,15,50,51,55,48,53,52,49,54, 0, 1, 5, 4, 2,3, 7, 6,
		18,20,16,22,19,21,17,23,40,42,43,45,44,46,41,47,60,58,59,61,63,57,56,62,37,34,35,38,39,36,32,33,30,26,27,29,31,25,24,28,55,50,51,48,52,54,53,49, 5, 0, 1, 6, 3, 2,4, 7,12,10, 9,15,13,11,8,14,
		18,20,21,19,17,23,22,16,59,60,58,61,62,57,63,56, 9,12,10,15,14,11,13, 8,35,37,34,36,38,39,32,33, 1, 5, 0, 2,6, 3, 7, 4,27,30,26,25,29,31,24,28,51,55,50,54,48,52,49,53,43,40,42,45,47,46,44,41,
		18,20,23,17,22,16,19,21,10, 9,12,15, 8,11,14,13,42,43,40,45,41,46,47,44,34,35,37,39,36,38,32,33,50,51,55,52,54,48,53,49, 0, 1, 5, 3, 2,6, 7, 4,26,27,30,31,25,29,28,24,58,59,60,61,56,57,62,63,
		18,21,19,20,17,22,16,23,58,59,60,61,63,62,57,56,50,51,55,49,48,53,52,54,26,27,30,25,31,24,28,29,10, 9,12,13,14,11,8,15,42,43,40,44,47,46,41,45, 0, 1, 5, 7, 6, 2,4, 3,34,35,37,32,33,38,36,39,
		18,21,22,17,16,23,20,19,55,50,51,49,54,53,48,52,37,34,35,32,39,38,33,36,30,26,27,24,25,31,28,29, 5, 0, 1, 2,7, 6, 3, 4,12,10, 9,11,13,14, 8,15,40,42,43,46,44,47,45,41,60,58,59,61,56,62,63,57,
		18,21,23,16,20,19,17,22,35,37,34,32,36,38,39,33,59,60,58,61,57,62,56,63,27,30,26,31,24,25,28,29,43,40,42,47,46,44,41,45, 1, 5, 0, 6, 2,7, 3, 4, 9,12,10,14,11,13,15, 8,51,55,50,49,52,53,54,48,
		18,22,17,21,16,20,19,23,51,55,50,49,48,54,53,52,43,40,42,45,44,41,47,46, 9,12,10,11,14, 8,15,13,35,37,34,33,39,38,36,32,59,60,58,63,56,62,57,61,1, 5, 0, 3, 7, 2,4, 6,27,30,26,28,29,25,24,31,
		18,22,20,16,19,23,21,17,42,43,40,45,46,41,44,47,26,27,30,28,31,25,29,24,10, 9,12, 8,11,14,15,13, 0, 1, 5, 2,3, 7, 6, 4,34,35,37,38,33,39,36,32,58,59,60,62,63,56,61,57,50,51,55,49,52,54,48,53,
		18,22,23,19,21,17,16,20,30,26,27,28,24,25,31,29,55,50,51,49,53,54,52,48,12,10, 9,14, 8,11,15,13,60,58,59,56,62,63,57,61,5, 0, 1, 7, 2,3, 6, 4,37,34,35,39,38,33,32,36,40,42,43,45,47,41,46,44,
		18,23,16,21,20,17,22,19,34,35,37,32,39,36,38,33,10, 9,12,15,11,8,13,14, 0, 1, 5, 6, 7, 3, 4, 2,58,59,60,56,57,62,63,61,50,51,55,54,52,53,48,49,42,43,40,41,46,47,45,44,26,27,30,28,29,24,31,25,
		18,23,17,20,22,19,21,16,12,10, 9,15,14, 8,11,13,30,26,27,28,25,24,29,31,5, 0, 1, 3, 6, 7, 4, 2,40,42,43,47,41,46,44,45,60,58,59,62,56,57,63,61,55,50,51,53,54,52,49,48,37,34,35,32,33,36,39,38,
		18,23,19,22,21,16,20,17,27,30,26,28,31,24,25,29,35,37,34,32,38,36,33,39, 1, 5, 0, 7, 3, 6, 4, 2,51,55,50,52,53,54,48,49,43,40,42,46,47,41,44,45,59,60,58,57,62,56,61,63, 9,12,10,15,13, 8,14,11,
		19,16,18,17,22,21,23,20, 1, 4, 5, 0, 7, 2,3, 6, 8,14,13,11,12,10, 9,15,58,60,61,63,57,56,59,62,36,33,38,34,37,32,35,39,44,46,41,43,45,40,42,47,28,27,26,31,29,24,30,25,52,48,53,54,55,50,51,49,
		19,16,20,23,17,18,22,21,48,53,52,54,51,50,49,55, 4, 5, 1, 0, 3, 2,6, 7,60,61,58,57,56,63,59,62,46,41,44,45,40,43,42,47,27,26,28,29,24,31,25,30,33,38,36,37,32,34,39,35,14,13, 8,11,9,10,15,12,
		19,16,21,22,23,20,17,18,13, 8,14,11,15,10,12, 9,53,52,48,54,49,50,55,51,61,58,60,56,63,57,59,62,26,28,27,24,31,29,25,30,38,36,33,32,34,37,35,39,41,44,46,40,43,45,47,42, 5, 1, 4, 0, 6, 2,7, 3,
		19,17,16,18,22,20,21,23, 4, 5, 1, 0, 2,3, 7, 6,33,38,36,39,37,32,34,35,48,53,52,51,50,49,54,55,27,26,28,24,29,25,31,30,60,61,58,63,62,57,56,59,14,13, 8,15,12, 9,11,10,46,41,44,47,45,40,43,42,
		19,17,20,22,21,23,18,16,36,33,38,39,35,32,37,34,44,46,41,47,42,40,45,43,52,48,53,49,51,50,54,55, 8,14,13, 9,15,12,10,11,28,27,26,25,24,29,31,30,58,60,61,57,63,62,59,56, 1, 4, 5, 0, 6, 3, 2,7,
		19,17,23,21,18,16,22,20,41,44,46,47,43,40,42,45, 5, 1, 4, 0, 7, 3, 6, 2,53,52,48,50,49,51,54,55,61,58,60,62,57,63,56,59,13, 8,14,12, 9,15,10,11,26,28,27,29,25,24,30,31,38,36,33,39,34,32,35,37,
		19,18,17,16,22,23,20,21,5, 1, 4, 0, 3, 7, 2,6,26,28,27,30,29,25,24,31,41,44,46,43,40,42,47,45,13, 8,14, 9,12,10,15,11,53,52,48,51,55,50,49,54,38,36,33,35,37,34,39,32,61,58,60,59,62,57,63,56,
		19,18,21,20,16,17,22,23,58,60,61,59,63,57,56,62, 1, 4, 5, 0, 2,7, 6, 3,44,46,41,40,42,43,47,45,52,48,53,55,50,51,49,54,36,33,38,37,34,35,32,39, 8,14,13,12,10, 9,11,15,28,27,26,30,24,25,31,29,
		19,18,23,22,20,21,16,17,27,26,28,30,31,25,29,24,60,61,58,59,56,57,62,63,46,41,44,42,43,40,47,45,33,38,36,34,35,37,32,39,14,13, 8,10, 9,12,15,11,48,53,52,50,51,55,54,49, 4, 5, 1, 0, 6, 7, 3, 2,
		19,20,18,21,16,23,17,22,60,61,58,59,57,56,63,62,48,53,52,54,50,51,55,49,27,26,28,31,25,29,30,24,14,13, 8,9,10,15,12,11,46,41,44,40,45,42,43,47, 4, 5, 1, 3, 2,6, 0, 7,33,38,36,39,34,35,37,32,
		19,20,22,17,21,18,16,23,38,36,33,39,37,35,32,34,61,58,60,59,63,56,62,57,26,28,27,25,29,31,30,24,41,44,46,45,42,40,43,47, 5, 1, 4, 2,6, 3, 7, 0,13, 8,14,10,15, 9,11,12,53,52,48,54,55,51,49,50,
		19,20,23,16,17,22,21,18,52,48,53,54,49,51,50,55,36,33,38,39,32,35,34,37,28,27,26,29,31,25,30,24, 1, 4, 5, 6, 3, 2,7, 0, 8,14,13,15, 9,10,12,11,44,46,41,42,40,45,47,43,58,60,61,59,62,56,57,63,
		19,21,17,23,18,20,16,22,44,46,41,47,40,42,43,45,58,60,61,59,57,63,62,56,36,33,38,35,32,37,39,34,28,27,26,24,25,31,29,30,52,48,53,50,55,49,51,54, 1, 4, 5, 2,7, 6, 0, 3, 8,14,13,11,9,15,12,10,
		19,21,20,18,16,22,23,17,61,58,60,59,56,63,57,62,13, 8,14,11,10,15, 9,12,38,36,33,37,35,32,39,34, 5, 1, 4, 6, 2,7, 3, 0,26,28,27,31,24,25,29,30,53,52,48,49,50,55,54,51,41,44,46,47,45,42,40,43,
		19,21,22,16,23,17,18,20,14,13, 8,11,12,15,10, 9,46,41,44,47,43,42,45,40,33,38,36,32,37,35,39,34,48,53,52,55,49,50,51,54, 4, 5, 1, 7, 6, 2,3, 0,27,26,28,25,31,24,30,29,60,61,58,59,62,63,56,57,
		19,22,16,21,23,18,20,17, 8,14,13,11,10,12,15, 9,28,27,26,30,31,29,24,25, 1, 4, 5, 7, 2,3, 0, 6,44,46,41,45,43,42,40,47,58,60,61,56,62,63,57,59,52,48,53,51,49,55,54,50,36,33,38,39,34,37,32,35,
		19,22,17,20,21,16,23,18,33,38,36,39,32,37,35,34,14,13, 8,11,15,12, 9,10, 4, 5, 1, 2,3, 7, 0, 6,60,61,58,62,63,56,57,59,48,53,52,49,55,51,50,54,46,41,44,43,42,45,47,40,27,26,28,30,24,29,25,31,
		19,22,18,23,20,17,21,16,26,28,27,30,25,29,31,24,38,36,33,39,35,37,34,32, 5, 1, 4, 3, 7, 2,0, 6,53,52,48,55,51,49,50,54,41,44,46,42,45,43,40,47,61,58,60,63,56,62,59,57,13, 8,14,11,9,12,10,15,
		19,23,16,20,17,21,18,22,53,52,48,54,50,49,51,55,41,44,46,47,40,43,45,42,13, 8,14,15,10,12,11,9,38,36,33,34,32,35,37,39,61,58,60,57,62,56,63,59, 5, 1, 4, 7, 3, 6, 0, 2,26,28,27,30,24,31,29,25,
		19,23,21,17,18,22,20,16,46,41,44,47,42,43,40,45,27,26,28,30,25,31,24,29,14,13, 8,12,15,10,11,9, 4, 5, 1, 6, 7, 3, 2,0,33,38,36,35,34,32,37,39,60,61,58,56,57,62,59,63,48,53,52,54,55,49,50,51,
		19,23,22,18,20,16,17,21,28,27,26,30,29,31,25,24,52,48,53,54,51,49,55,50, 8,14,13,10,12,15,11,9,58,60,61,62,56,57,63,59, 1, 4, 5, 3, 6, 7, 2,0,36,33,38,32,35,34,39,37,44,46,41,47,45,43,42,40,
		20,16,17,21,19,23,22,18,29,31,25,24,28,27,26,30,48,52,54,53,49,51,50,55,15, 9,10,13,11,8,12,14,59,61,60,63,57,56,62,58,6, 3, 2,4, 1, 0, 5, 7,38,39,36,34,37,32,33,35,42,40,45,43,41,47,44,46,
		20,16,18,22,21,17,19,23,40,45,42,43,44,47,46,41,31,25,29,24,26,27,30,28,9,10,15,11,8,13,12,14, 3, 2,6, 1, 0, 4, 5, 7,39,36,38,37,32,34,35,33,61,60,59,57,56,63,58,62,52,54,48,53,50,51,55,49,
		20,16,23,19,22,18,21,17,54,48,52,53,55,51,49,50,45,42,40,43,46,47,41,44,10,15, 9, 8,13,11,12,14,36,38,39,32,34,37,35,33,60,59,61,56,63,57,62,58,2,6, 3, 0, 4, 1, 7, 5,25,29,31,24,30,27,28,26,
		20,17,21,16,19,22,18,23,25,29,31,24,26,28,27,30,36,38,39,33,37,35,32,34, 2,6, 3, 4, 0, 5, 7, 1,54,48,52,50,49,51,55,53,45,42,40,44,41,47,46,43,60,59,61,62,57,63,58,56,10,15, 9,12,14,11,13, 8,
		20,17,22,19,18,23,16,21,39,36,38,33,34,35,37,32, 9,10,15,12, 8,11,14,13, 3, 2,6, 5, 4, 0, 7, 1,61,60,59,63,62,57,56,58,52,54,48,51,50,49,55,53,40,45,42,47,44,41,43,46,31,25,29,24,30,28,26,27,
		20,17,23,18,16,21,19,22,15, 9,10,12,13,11,8,14,29,31,25,24,27,28,30,26, 6, 3, 2,0, 5, 4, 7, 1,42,40,45,41,47,44,46,43,59,61,60,57,63,62,56,58,48,52,54,49,51,50,53,55,38,39,36,33,32,35,34,37,
		20,18,17,23,16,22,21,19, 9,10,15,12,11,8,13,14,40,45,42,43,47,44,41,46,39,36,38,34,35,37,33,32,52,54,48,50,51,55,49,53, 3, 2,6, 0, 1, 5, 4, 7,31,25,29,26,27,30,24,28,61,60,59,58,63,62,57,56,
		20,18,19,21,23,17,16,22,60,59,61,58,57,62,56,63,10,15, 9,12,13, 8,14,11,36,38,39,35,37,34,33,32, 2,6, 3, 1, 5, 0, 4, 7,25,29,31,27,30,26,28,24,54,48,52,51,55,50,53,49,45,42,40,43,41,44,46,47,
		20,18,22,16,21,19,23,17,42,40,45,43,46,44,47,41,59,61,60,58,56,62,63,57,38,39,36,37,34,35,33,32,29,31,25,30,26,27,28,24,48,52,54,55,50,51,49,53, 6, 3, 2,5, 0, 1, 7, 4,15, 9,10,12,14, 8,11,13,
		20,19,16,23,22,17,18,21,48,52,54,53,51,49,55,50,38,39,36,33,34,37,32,35,29,31,25,28,27,26,24,30, 6, 3, 2,1, 4, 5, 0, 7,15, 9,10, 8,14,13,11,12,42,40,45,44,46,41,43,47,59,61,60,58,63,57,56,62,
		20,19,17,22,18,21,23,16,36,38,39,33,35,37,34,32,60,59,61,58,62,57,63,56,25,29,31,26,28,27,24,30,45,42,40,41,44,46,47,43, 2,6, 3, 5, 1, 4, 0, 7,10,15, 9,13, 8,14,12,11,54,48,52,53,50,49,51,55,
		20,19,21,18,23,16,22,17,61,60,59,58,56,57,62,63,52,54,48,53,55,49,50,51,31,25,29,27,26,28,24,30, 9,10,15,14,13, 8,11,12,40,45,42,46,41,44,47,43, 3, 2,6, 4, 5, 1, 7, 0,39,36,38,33,32,37,35,34,
		20,21,16,17,19,18,23,22,31,25,29,24,27,26,28,30,61,60,59,58,57,56,63,62,40,45,42,44,47,46,43,41,39,36,38,32,37,35,34,33, 9,10,15,13,14,11,8,12,52,54,48,55,49,50,53,51,3, 2,6, 7, 1, 0, 4, 5,
		20,21,18,19,23,22,17,16,59,61,60,58,62,56,57,63, 6, 3, 2,7, 5, 0, 1, 4,42,40,45,46,44,47,43,41,48,52,54,50,55,49,51,53,38,39,36,35,32,37,34,33,15, 9,10,11,13,14,12, 8,29,31,25,24,30,26,27,28,
		20,21,22,23,17,16,19,18,2,6, 3, 7, 4, 0, 5, 1,25,29,31,24,28,26,30,27,45,42,40,47,46,44,43,41,10,15, 9,14,11,13, 8,12,54,48,52,49,50,55,51,53,36,38,39,37,35,32,33,34,60,59,61,58,63,56,62,57,
		20,22,16,18,21,23,17,19,45,42,40,43,47,46,44,41,2,6, 3, 7, 0, 4, 1, 5,54,48,52,55,51,49,53,50,60,59,61,63,56,62,57,58,10,15, 9,11,14, 8,13,12,25,29,31,28,26,30,24,27,36,38,39,33,32,34,37,35,
		20,22,19,17,18,16,21,23,38,39,36,33,37,34,35,32,42,40,45,43,44,46,41,47,48,52,54,51,49,55,53,50,15, 9,10,14, 8,11,13,12,29,31,25,26,30,28,27,24,59,61,60,56,62,63,58,57, 6, 3, 2,7, 1, 4, 5, 0,
		20,22,23,21,17,19,18,16, 3, 2,6, 7, 5, 4, 0, 1,39,36,38,33,35,34,32,37,52,54,48,49,55,51,53,50,31,25,29,30,28,26,27,24,61,60,59,62,63,56,57,58,9,10,15, 8,11,14,12,13,40,45,42,43,41,46,47,44,
		20,23,18,17,16,19,22,21,10,15, 9,12, 8,13,11,14,54,48,52,53,51,55,50,49,60,59,61,57,62,56,58,63,25,29,31,30,27,28,26,24,36,38,39,34,32,35,37,33,45,42,40,46,47,41,43,44, 2,6, 3, 7, 1, 5, 0, 4,
		20,23,19,16,22,21,17,18,52,54,48,53,49,55,51,50, 3, 2,6, 7, 4, 5, 1, 0,61,60,59,56,57,62,58,63,40,45,42,41,46,47,44,43,31,25,29,28,30,27,26,24,39,36,38,35,34,32,33,37, 9,10,15,12,14,13, 8,11,
		20,23,21,22,17,18,16,19, 6, 3, 2,7, 0, 5, 4, 1,15, 9,10,12,11,13,14, 8,59,61,60,62,56,57,58,63,38,39,36,32,35,34,37,33,42,40,45,47,41,46,44,43,29,31,25,27,28,30,24,26,48,52,54,53,50,55,49,51,
		21,16,20,17,18,23,19,22,31,24,25,29,27,30,26,28,37,35,32,34,36,38,39,33, 6, 2,7, 0, 4, 1, 3, 5,49,50,55,48,54,53,52,51,47,46,44,40,43,45,42,41,58,61,59,56,63,57,60,62,14,11,13, 8,10,15, 9,12,
		21,16,22,19,17,20,18,23,11,13,14, 8,9,15,12,10,24,25,31,29,26,30,28,27, 2,7, 6, 4, 1, 0, 3, 5,46,44,47,43,45,40,42,41,61,59,58,63,57,56,62,60,50,55,49,54,53,48,51,52,35,32,37,34,39,38,33,36,
		21,16,23,18,19,22,17,20,32,37,35,34,33,38,36,39,13,14,11,8,12,15,10, 9, 7, 6, 2,1, 0, 4, 3, 5,59,58,61,57,56,63,62,60,55,49,50,53,48,54,52,51,44,47,46,45,40,43,41,42,25,31,24,29,28,30,27,26,
		21,17,16,20,18,22,23,19,24,25,31,29,30,26,27,28,50,55,49,51,54,53,48,52,11,13,14, 9,15,12, 8,10,61,59,58,57,63,62,56,60, 2,7, 6, 0, 5, 4, 1, 3,35,32,37,33,36,39,34,38,46,44,47,41,43,45,40,42,
		21,17,19,23,20,16,18,22,44,47,46,41,40,45,42,43,25,31,24,29,27,26,28,30,13,14,11,15,12, 9, 8,10, 7, 6, 2,5, 4, 0, 1, 3,32,37,35,36,39,33,38,34,59,58,61,63,62,57,60,56,55,49,50,51,48,53,52,54,
		21,17,22,18,23,19,20,16,49,50,55,51,52,53,54,48,47,46,44,41,42,45,43,40,14,11,13,12, 9,15, 8,10,37,35,32,39,33,36,38,34,58,61,59,62,57,63,56,60, 6, 2,7, 4, 0, 5, 3, 1,31,24,25,29,28,26,30,27,
		21,18,16,23,19,20,22,17,37,35,32,34,38,36,33,39,58,61,59,60,56,63,57,62,31,24,25,27,30,26,29,28,47,46,44,43,40,42,45,41,6, 2,7, 1, 5, 0, 4, 3,14,11,13, 9,12,10, 8,15,49,50,55,51,48,54,53,52,
		21,18,17,22,23,16,19,20,50,55,49,51,53,54,52,48,35,32,37,34,33,36,39,38,24,25,31,30,26,27,29,28,2,7, 6, 5, 0, 1, 4, 3,11,13,14,12,10, 9,15, 8,46,44,47,40,42,43,41,45,61,59,58,60,57,63,62,56,
		21,18,20,19,22,17,23,16,59,58,61,60,62,63,56,57,55,49,50,51,52,54,48,53,25,31,24,26,27,30,29,28,13,14,11,10, 9,12,15, 8,44,47,46,42,43,40,45,41,7, 6, 2,0, 1, 5, 3, 4,32,37,35,34,39,36,38,33,
		21,19,16,22,17,23,20,18,13,14,11,8,15,12, 9,10,44,47,46,41,45,40,43,42,32,37,35,33,38,36,34,39,55,49,50,48,53,52,54,51,7, 6, 2,4, 5, 1, 0, 3,25,31,24,27,26,28,29,30,59,58,61,60,57,56,63,62,
		21,19,18,20,22,16,17,23,58,61,59,60,63,56,62,57,14,11,13, 8,9,12,10,15,37,35,32,38,36,33,34,39, 6, 2,7, 5, 1, 4, 0, 3,31,24,25,26,28,27,30,29,49,50,55,53,52,48,51,54,47,46,44,41,43,40,42,45,
		21,19,23,17,20,18,22,16,46,44,47,41,42,40,45,43,61,59,58,60,62,56,57,63,35,32,37,36,33,38,34,39,24,25,31,28,27,26,30,29,50,55,49,52,48,53,54,51,2,7, 6, 1, 4, 5, 3, 0,11,13,14, 8,10,12,15, 9,
		21,20,17,16,18,19,22,23,25,31,24,29,26,27,30,28,59,58,61,60,63,62,57,56,44,47,46,40,45,42,41,43,32,37,35,39,36,38,33,34,13,14,11,9,10,15,12, 8,55,49,50,52,54,48,51,53, 7, 6, 2,3, 5, 4, 0, 1,
		21,20,19,18,22,23,16,17,61,59,58,60,56,62,63,57, 2,7, 6, 3, 1, 4, 5, 0,46,44,47,42,40,45,41,43,50,55,49,48,52,54,53,51,35,32,37,38,39,36,33,34,11,13,14,15, 9,10, 8,12,24,25,31,29,28,27,26,30,
		21,20,23,22,16,17,18,19, 6, 2,7, 3, 0, 4, 1, 5,31,24,25,29,30,27,28,26,47,46,44,45,42,40,41,43,14,11,13,10,15, 9,12, 8,49,50,55,54,48,52,53,51,37,35,32,36,38,39,34,33,58,61,59,60,57,62,56,63,
		21,22,18,17,23,20,16,19,55,49,50,51,54,52,53,48,7, 6, 2,3, 0, 1, 5, 4,59,58,61,62,63,56,60,57,44,47,46,43,42,45,40,41,25,31,24,30,28,26,27,29,32,37,35,38,33,39,34,36,13,14,11,8,10, 9,12,15,
		21,22,19,16,17,18,23,20,14,11,13, 8,12, 9,15,10,49,50,55,51,53,52,48,54,58,61,59,63,56,62,60,57,31,24,25,28,26,30,27,29,37,35,32,33,39,38,36,34,47,46,44,42,45,43,41,40, 6, 2,7, 3, 5, 1, 4, 0,
		21,22,20,23,16,19,17,18,2,7, 6, 3, 4, 1, 0, 5,11,13,14, 8,15, 9,10,12,61,59,58,56,62,63,60,57,35,32,37,39,38,33,36,34,46,44,47,45,43,42,40,41,24,25,31,26,30,28,29,27,50,55,49,51,48,52,54,53,
		21,23,17,19,20,22,16,18,47,46,44,41,45,42,40,43, 6, 2,7, 3, 4, 0, 5, 1,49,50,55,52,53,54,51,48,58,61,59,57,62,56,63,60,14,11,13,15,10,12, 9, 8,31,24,25,30,27,28,29,26,37,35,32,34,39,33,36,38,
		21,23,18,16,19,17,20,22,35,32,37,34,36,33,38,39,46,44,47,41,40,42,43,45,50,55,49,53,54,52,51,48,11,13,14,10,12,15, 9, 8,24,25,31,27,28,30,26,29,61,59,58,62,56,57,60,63, 2,7, 6, 3, 5, 0, 1, 4,
		21,23,22,20,16,18,19,17, 7, 6, 2,3, 1, 0, 4, 5,32,37,35,34,38,33,39,36,55,49,50,54,52,53,51,48,25,31,24,28,30,27,26,29,59,58,61,56,57,62,63,60,13,14,11,12,15,10, 8,9,44,47,46,41,43,42,45,40,
		22,16,17,23,21,19,18,20,63,56,62,57,58,61,59,60,11,14, 8,13,12, 9,15,10,33,39,38,32,34,37,36,35, 3, 7, 2,0, 4, 1, 5, 6,28,26,30,24,31,29,25,27,55,51,49,48,54,53,50,52,42,45,43,40,44,41,47,46,
		22,16,19,21,18,20,23,17, 8,11,14,13,10, 9,12,15,43,42,45,40,46,41,44,47,38,33,39,37,32,34,36,35,49,55,51,53,48,54,52,50, 2,3, 7, 1, 0, 4, 5, 6,30,28,26,29,24,31,27,25,62,63,56,57,60,61,58,59,
		22,16,20,18,23,17,21,19,45,43,42,40,47,41,46,44,56,62,63,57,59,61,60,58,39,38,33,34,37,32,36,35,26,30,28,31,29,24,25,27,51,49,55,54,53,48,52,50, 7, 2,3, 4, 1, 0, 6, 5,14, 8,11,13,15, 9,10,12,
		22,17,18,21,20,19,16,23,51,49,55,50,48,52,54,53,39,38,33,36,37,34,35,32,26,30,28,25,24,29,27,31,7, 2,3, 0, 5, 4, 1, 6,14, 8,11,9,15,12,10,13,45,43,42,41,47,44,40,46,56,62,63,57,60,58,59,61,
		22,17,19,20,16,23,21,18,33,39,38,36,32,34,37,35,63,56,62,57,61,58,60,59,28,26,30,29,25,24,27,31,42,45,43,44,41,47,46,40, 3, 7, 2,4, 0, 5, 1, 6,11,14, 8,12, 9,15,13,10,55,51,49,50,53,52,48,54,
		22,17,23,16,21,18,20,19,62,63,56,57,59,58,61,60,49,55,51,50,54,52,53,48,30,28,26,24,29,25,27,31,8,11,14,15,12, 9,10,13,43,42,45,47,44,41,46,40, 2,3, 7, 5, 4, 0, 6, 1,38,33,39,36,35,34,32,37,
		22,18,16,20,23,19,17,21,43,42,45,40,41,46,47,44,30,28,26,27,29,24,31,25, 8,11,14,10, 9,12,13,15, 2,3, 7, 0, 1, 5, 4, 6,38,33,39,34,35,37,32,36,62,63,56,58,59,60,57,61,49,55,51,50,53,48,54,52,
		22,18,19,23,17,21,20,16,26,30,28,27,25,24,29,31,51,49,55,50,52,48,53,54,14, 8,11,12,10, 9,13,15,56,62,63,60,58,59,61,57, 7, 2,3, 5, 0, 1, 4, 6,39,38,33,37,34,35,36,32,45,43,42,40,44,46,41,47,
		22,18,21,17,20,16,23,19,55,51,49,50,54,48,52,53,42,45,43,40,47,46,44,41,11,14, 8,9,12,10,13,15,33,39,38,35,37,34,32,36,63,56,62,59,60,58,61,57, 3, 7, 2,1, 5, 0, 6, 4,28,26,30,27,31,24,25,29,
		22,19,20,17,16,21,18,23,38,33,39,36,37,32,34,35, 8,11,14,13, 9,10,15,12, 2,3, 7, 4, 5, 1, 6, 0,62,63,56,60,61,58,59,57,49,55,51,48,53,52,54,50,43,42,45,46,41,44,40,47,30,28,26,27,31,25,29,24,
		22,19,21,16,18,23,17,20,14, 8,11,13,12,10, 9,15,26,30,28,27,24,25,31,29, 7, 2,3, 1, 4, 5, 6, 0,45,43,42,44,46,41,47,40,56,62,63,58,60,61,59,57,51,49,55,52,48,53,50,54,39,38,33,36,35,32,37,34,
		22,19,23,18,17,20,16,21,28,26,30,27,29,25,24,31,33,39,38,36,34,32,35,37, 3, 7, 2,5, 1, 4, 6, 0,55,51,49,53,52,48,54,50,42,45,43,41,44,46,47,40,63,56,62,61,58,60,57,59,11,14, 8,13,15,10,12, 9,
		22,20,17,19,16,18,23,21,39,38,33,36,34,37,32,35,45,43,42,40,41,47,44,46,51,49,55,48,52,54,50,53,14, 8,11,15, 9,10,12,13,26,30,28,29,31,25,24,27,56,62,63,59,61,60,57,58,7, 2,3, 6, 0, 5, 4, 1,
		22,20,18,16,23,21,19,17,42,45,43,40,46,47,41,44, 3, 7, 2,6, 1, 5, 0, 4,55,51,49,54,48,52,50,53,63,56,62,60,59,61,58,57,11,14, 8,10,15, 9,12,13,28,26,30,25,29,31,27,24,33,39,38,36,35,37,34,32,
		22,20,21,23,19,17,16,18,2,3, 7, 6, 4, 5, 1, 0,38,33,39,36,32,37,35,34,49,55,51,52,54,48,50,53,30,28,26,31,25,29,24,27,62,63,56,61,60,59,58,57, 8,11,14, 9,10,15,13,12,43,42,45,40,44,47,46,41,
		22,21,16,19,18,17,20,23,11,14, 8,13, 9,12,10,15,55,51,49,50,48,54,53,52,63,56,62,58,61,59,57,60,28,26,30,31,24,25,29,27,33,39,38,37,35,32,34,36,42,45,43,47,46,44,40,41,3, 7, 2,6, 0, 4, 1, 5,
		22,21,17,18,20,23,19,16,49,55,51,50,52,54,48,53, 2,3, 7, 6, 5, 4, 0, 1,62,63,56,59,58,61,57,60,43,42,45,44,47,46,41,40,30,28,26,25,31,24,29,27,38,33,39,32,37,35,36,34, 8,11,14,13,15,12, 9,10,
		22,21,23,20,19,16,18,17, 7, 2,3, 6, 1, 4, 5, 0,14, 8,11,13,10,12,15, 9,56,62,63,61,59,58,57,60,39,38,33,35,32,37,34,36,45,43,42,46,44,47,41,40,26,30,28,24,25,31,27,29,51,49,55,50,53,54,52,48,
		22,23,16,17,21,20,19,18,56,62,63,57,61,59,58,60, 7, 2,3, 6, 4, 1, 0, 5,45,43,42,47,41,46,40,44,51,49,55,53,54,52,48,50,39,38,33,32,35,34,37,36,14, 8,11,10,12,15,13, 9,26,30,28,27,31,29,24,25,
		22,23,18,19,17,16,21,20,30,28,26,27,24,29,25,31,62,63,56,57,58,59,60,61,43,42,45,41,46,47,40,44,38,33,39,35,34,32,37,36, 8,11,14,12,15,10, 9,13,49,55,51,54,52,53,50,48,2,3, 7, 6, 0, 1, 5, 4,
		22,23,20,21,19,18,17,16, 3, 7, 2,6, 5, 1, 4, 0,28,26,30,27,25,29,31,24,42,45,43,46,47,41,40,44,11,14, 8,15,10,12, 9,13,55,51,49,52,53,54,48,50,33,39,38,34,32,35,36,37,63,56,62,57,60,59,61,58,
		23,16,18,21,17,22,20,19,34,32,35,37,39,33,36,38,57,62,56,63,59,60,58,61,30,27,28,24,31,29,26,25,46,47,41,40,43,45,42,44, 7, 3, 6, 0, 4, 1, 5, 2,15,10,12, 8,13,11,9,14,52,53,54,48,51,55,50,49,
		23,16,19,20,21,18,17,22,53,54,52,48,50,55,49,51,32,35,34,37,36,33,38,39,27,28,30,31,29,24,26,25, 3, 6, 7, 4, 1, 0, 5, 2,10,12,15,13,11,8,14, 9,47,41,46,43,45,40,44,42,62,56,57,63,58,60,61,59,
		23,16,22,17,20,19,21,18,56,57,62,63,61,60,59,58,54,52,53,48,49,55,51,50,28,30,27,29,24,31,26,25,12,15,10,11,8,13,14, 9,41,46,47,45,40,43,42,44, 6, 7, 3, 1, 0, 4, 2,5,35,34,32,37,38,33,39,36,
		23,17,16,22,20,18,19,21,57,62,56,63,60,59,61,58,15,10,12, 9, 8,13,11,14,34,32,35,39,33,36,37,38,7, 3, 6, 4, 0, 5, 1, 2,30,27,28,29,25,24,31,26,52,53,54,50,49,51,48,55,46,47,41,44,40,43,45,42,
		23,17,18,20,19,21,22,16,12,15,10, 9,14,13, 8,11,41,46,47,44,42,43,40,45,35,34,32,36,39,33,37,38,54,52,53,51,50,49,55,48,6, 7, 3, 5, 4, 0, 1, 2,28,30,27,24,29,25,26,31,56,57,62,63,58,59,60,61,
		23,17,21,19,22,16,20,18,47,41,46,44,45,43,42,40,62,56,57,63,61,59,58,60,32,35,34,33,36,39,37,38,27,28,30,25,24,29,31,26,53,54,52,49,51,50,55,48,3, 6, 7, 0, 5, 4, 2,1,10,12,15, 9,11,13,14, 8,
		23,18,20,17,19,22,16,21,10,12,15, 9, 8,14,13,11,27,28,30,26,29,31,25,24, 3, 6, 7, 5, 0, 1, 2,4,47,41,46,40,42,43,45,44,62,56,57,60,58,59,61,63,53,54,52,55,50,51,48,49,32,35,34,37,38,39,36,33,
		23,18,21,16,17,20,19,22,35,34,32,37,36,39,33,38,12,15,10, 9,13,14,11,8,6, 7, 3, 0, 1, 5, 2,4,56,57,62,58,59,60,61,63,54,52,53,50,51,55,49,48,41,46,47,42,43,40,44,45,28,30,27,26,25,31,24,29,
		23,18,22,19,16,21,17,20,30,27,28,26,24,31,29,25,34,32,35,37,33,39,38,36, 7, 3, 6, 1, 5, 0, 2,4,52,53,54,51,55,50,49,48,46,47,41,43,40,42,45,44,57,62,56,59,60,58,63,61,15,10,12, 9,11,14, 8,13,
		23,19,17,21,22,18,16,20,41,46,47,44,43,42,45,40,28,30,27,26,24,29,25,31,12,15,10,14,13, 8,9,11,6, 7, 3, 4, 5, 1, 0, 2,35,34,32,33,38,36,39,37,56,57,62,60,61,58,63,59,54,52,53,48,51,50,49,55,
		23,19,18,22,16,20,21,17,27,28,30,26,31,29,24,25,53,54,52,48,55,50,51,49,10,12,15, 8,14,13, 9,11,62,56,57,58,60,61,59,63, 3, 6, 7, 1, 4, 5, 0, 2,32,35,34,36,33,38,37,39,47,41,46,44,40,42,43,45,
		23,19,20,16,21,17,22,18,52,53,54,48,49,50,55,51,46,47,41,44,45,42,40,43,15,10,12,13, 8,14, 9,11,34,32,35,38,36,33,39,37,57,62,56,61,58,60,59,63, 7, 3, 6, 5, 1, 4, 2,0,30,27,28,26,25,29,31,24,
		23,20,16,19,21,22,18,17,54,52,53,48,55,49,50,51,6, 7, 3, 2,1, 0, 4, 5,56,57,62,61,60,59,63,58,41,46,47,40,45,42,43,44,28,30,27,31,25,29,24,26,35,34,32,39,36,38,37,33,12,15,10, 9,11,8,13,14,
		23,20,17,18,19,16,21,22,15,10,12, 9,13, 8,14,11,52,53,54,48,50,49,51,55,57,62,56,60,59,61,63,58,30,27,28,25,29,31,24,26,34,32,35,36,38,39,33,37,46,47,41,45,42,40,44,43, 7, 3, 6, 2,4, 0, 5, 1,
		23,20,22,21,18,17,19,16, 3, 6, 7, 2,5, 0, 1, 4,10,12,15, 9,14, 8,11,13,62,56,57,59,61,60,63,58,32,35,34,38,39,36,33,37,47,41,46,42,40,45,43,44,27,28,30,29,31,25,26,24,53,54,52,48,51,49,55,50,
		23,21,16,18,17,19,22,20,32,35,34,37,33,36,39,38,47,41,46,44,43,45,40,42,53,54,52,50,55,49,48,51,10,12,15,11,13,14, 8,9,27,28,30,24,25,31,29,26,62,56,57,61,59,58,63,60, 3, 6, 7, 2,4, 1, 0, 5,
		23,21,19,17,22,20,18,16,46,47,41,44,42,45,43,40, 7, 3, 6, 2,5, 1, 4, 0,52,53,54,49,50,55,48,51,57,62,56,58,61,59,60,63,15,10,12,14,11,13, 8,9,30,27,28,31,24,25,26,29,34,32,35,37,38,36,33,39,
		23,21,20,22,18,16,17,19, 6, 7, 3, 2,0, 1, 5, 4,35,34,32,37,39,36,38,33,54,52,53,55,49,50,48,51,28,30,27,25,31,24,29,26,56,57,62,59,58,61,60,63,12,15,10,13,14,11,9, 8,41,46,47,44,40,45,42,43,
		23,22,17,16,20,21,18,19,62,56,57,63,59,61,60,58,3, 6, 7, 2,0, 5, 4, 1,47,41,46,45,43,42,44,40,53,54,52,51,49,55,50,48,32,35,34,39,38,33,36,37,10,12,15,14, 8,11,9,13,27,28,30,26,25,24,29,31,
		23,22,19,18,16,17,20,21,28,30,27,26,29,24,31,25,56,57,62,63,60,61,58,59,41,46,47,43,42,45,44,40,35,34,32,38,33,39,36,37,12,15,10, 8,11,14,13, 9,54,52,53,49,55,51,48,50, 6, 7, 3, 2,4, 5, 1, 0,
		23,22,21,20,18,19,16,17, 7, 3, 6, 2,1, 5, 0, 4,30,27,28,26,31,24,25,29,46,47,41,42,45,43,44,40,15,10,12,11,14, 8,13, 9,52,53,54,55,51,49,50,48,34,32,35,33,39,38,37,36,57,62,56,63,58,61,59,60,
		24,25,26,30,27,28,29,31,0, 2,5, 7, 1, 3, 4, 6,32,35,36,33,38,39,34,37,56,61,63,60,57,59,58,62,40,46,47,45,41,44,42,43, 8,11,15,10, 9,13,14,12,48,50,51,54,49,55,53,52,16,21,17,20,19,22,18,23,
		24,25,28,27,29,31,30,26,36,32,35,33,37,39,38,34,17,16,21,20,23,22,19,18,63,56,61,59,60,57,58,62,51,48,50,55,54,49,52,53,47,40,46,44,45,41,42,43,15, 8,11,13,10, 9,12,14, 5, 0, 2,7, 6, 3, 1, 4,
		24,25,31,29,30,26,27,28,21,17,16,20,18,22,23,19, 2,5, 0, 7, 4, 3, 6, 1,61,63,56,57,59,60,58,62,11,15, 8,9,13,10,14,12,50,51,48,49,55,54,52,53,46,47,40,41,44,45,43,42,35,36,32,33,34,39,37,38,
		24,26,28,31,25,30,27,29,56,61,63,58,60,57,59,62, 0, 2,5, 7, 3, 1, 6, 4, 8,11,15,13,14,10,12, 9,16,21,17,19,22,18,23,20,40,46,47,41,45,42,44,43,32,35,36,38,39,34,33,37,48,50,51,53,55,52,54,49,
		24,26,29,27,31,28,25,30,50,51,48,53,54,52,49,55,61,63,56,58,59,57,62,60,11,15, 8,14,10,13,12, 9,46,47,40,45,42,41,44,43,35,36,32,39,34,38,37,33,21,17,16,22,18,19,20,23, 2,5, 0, 7, 6, 1, 4, 3,
		24,26,30,25,27,29,31,28,5, 0, 2,7, 4, 1, 3, 6,51,48,50,53,49,52,55,54,15, 8,11,10,13,14,12, 9,36,32,35,34,38,39,37,33,17,16,21,18,19,22,23,20,47,40,46,42,41,45,43,44,63,56,61,58,62,57,60,59,
		24,27,25,28,29,26,31,30,32,35,36,33,39,38,37,34,48,50,51,53,54,49,55,52, 0, 2,5, 1, 3, 4, 7, 6, 8,11,15, 9,10,14,13,12,56,61,63,59,62,60,57,58,16,21,17,18,23,19,20,22,40,46,47,43,45,41,44,42,
		24,27,26,29,31,30,28,25,51,48,50,53,52,49,54,55,47,40,46,43,42,41,45,44, 5, 0, 2,4, 1, 3, 7, 6,17,16,21,19,18,23,22,20,15, 8,11,14, 9,10,13,12,63,56,61,60,59,62,58,57,36,32,35,33,34,38,39,37,
		24,27,30,31,28,25,29,26,46,47,40,43,44,41,42,45,35,36,32,33,37,38,34,39, 2,5, 0, 3, 4, 1, 7, 6,61,63,56,62,60,59,57,58,21,17,16,23,19,18,22,20,11,15, 8,10,14, 9,12,13,50,51,48,53,55,49,52,54,
		24,28,27,25,29,30,26,31,35,36,32,33,38,37,39,34,11,15, 8,12,10,14, 9,13,46,47,40,44,41,42,43,45,21,17,16,19,23,22,18,20, 2,5, 0, 1, 6, 3, 4, 7,50,51,48,52,54,55,53,49,61,63,56,58,62,60,59,57,
		24,28,30,29,26,31,25,27, 8,11,15,12,13,14,10, 9,56,61,63,58,57,60,62,59,40,46,47,42,44,41,43,45,48,50,51,55,52,54,49,53,16,21,17,22,19,23,18,20, 0, 2,5, 3, 1, 6, 7, 4,32,35,36,33,34,37,38,39,
		24,28,31,26,25,27,29,30,63,56,61,58,59,60,57,62,36,32,35,33,39,37,34,38,47,40,46,41,42,44,43,45, 5, 0, 2,6, 3, 1, 4, 7,51,48,50,54,55,52,49,53,17,16,21,23,22,19,20,18,15, 8,11,12, 9,14,13,10,
		24,29,25,31,30,28,26,27,17,16,21,20,22,23,18,19,15, 8,11,12,13,10, 9,14,36,32,35,37,39,38,33,34,47,40,46,45,44,42,41,43,63,56,61,57,62,59,60,58,5, 0, 2,1, 4, 6, 7, 3,51,48,50,53,55,54,49,52,
		24,29,27,26,31,25,30,28,48,50,51,53,49,54,52,55,16,21,17,20,18,23,19,22,32,35,36,39,38,37,33,34,56,61,63,62,59,57,60,58,0, 2,5, 4, 6, 1, 3, 7,40,46,47,44,42,45,43,41,8,11,15,12, 9,10,14,13,
		24,29,28,30,26,27,31,25,11,15, 8,12,14,10,13, 9,50,51,48,53,52,54,55,49,35,36,32,38,37,39,33,34, 2,5, 0, 6, 1, 4, 3, 7,46,47,40,42,45,44,41,43,61,63,56,59,57,62,58,60,21,17,16,20,19,23,22,18,
		24,30,25,26,27,31,28,29, 2,5, 0, 7, 3, 4, 1, 6,46,47,40,43,41,44,45,42,21,17,16,18,22,23,20,19,50,51,48,55,49,52,54,53,61,63,56,60,62,57,59,58,35,36,32,37,38,34,33,39,11,15, 8,12, 9,13,10,14,
		24,30,29,28,26,25,27,31,15, 8,11,12,10,13,14, 9, 5, 0, 2,7, 1, 4, 6, 3,17,16,21,22,23,18,20,19,63,56,61,62,57,60,59,58,36,32,35,38,34,37,39,33,51,48,50,49,52,55,53,54,47,40,46,43,45,44,42,41,
		24,30,31,27,28,29,26,25,40,46,47,43,42,44,41,45, 8,11,15,12,14,13, 9,10,16,21,17,23,18,22,20,19,32,35,36,34,37,38,39,33,48,50,51,52,55,49,54,53,56,61,63,57,60,62,58,59, 0, 2,5, 7, 6, 4, 3, 1,
		24,31,26,28,25,29,30,27,61,63,56,58,57,59,60,62,21,17,16,20,22,18,19,23,50,51,48,54,52,49,53,55,35,36,32,34,39,37,38,33,11,15, 8,13, 9,14,10,12, 2,5, 0, 4, 3, 6, 7, 1,46,47,40,43,45,42,41,44,
		24,31,27,30,28,26,25,29,47,40,46,43,41,42,44,45,63,56,61,58,60,59,62,57,51,48,50,52,49,54,53,55,15, 8,11,9,14,13,10,12, 5, 0, 2,3, 6, 4, 1, 7,36,32,35,39,37,34,33,38,17,16,21,20,19,18,23,22,
		24,31,29,25,30,27,28,26,16,21,17,20,23,18,22,19,40,46,47,43,44,42,45,41,48,50,51,49,54,52,53,55, 0, 2,5, 6, 4, 3, 1, 7,32,35,36,37,34,39,38,33, 8,11,15,14,13, 9,12,10,56,61,63,58,62,59,57,60,
		25,24,27,28,31,29,26,30,32,36,33,35,39,37,34,38,21,20,17,16,19,18,23,22,59,60,57,63,56,61,62,58,55,54,49,51,48,50,53,52,44,45,41,47,40,46,43,42,13,10, 9,15, 8,11,14,12, 7, 2,0, 5, 4, 1, 3, 6,
		25,24,29,31,26,30,28,27,17,21,20,16,22,18,19,23, 0, 7, 2,5, 6, 1, 4, 3,57,59,60,61,63,56,62,58,9,13,10,11,15, 8,12,14,49,55,54,50,51,48,53,52,41,44,45,46,47,40,42,43,33,32,36,35,38,37,39,34,
		25,24,30,26,28,27,31,29, 2,0, 7, 5, 3, 1, 6, 4,36,33,32,35,34,37,38,39,60,57,59,56,61,63,62,58,45,41,44,40,46,47,43,42,10, 9,13, 8,11,15,12,14,54,49,55,48,50,51,52,53,20,17,21,16,23,18,22,19,
		25,26,24,30,28,29,27,31,0, 7, 2,5, 1, 6, 3, 4,41,44,45,42,46,47,40,43,17,21,20,22,18,19,16,23,49,55,54,51,50,53,48,52,57,59,60,56,58,61,63,62,33,32,36,39,34,38,35,37, 9,13,10,14,11,15, 8,12,
		25,26,29,28,27,31,30,24,45,41,44,42,43,47,46,40,10, 9,13,14,12,15,11,8,20,17,21,19,22,18,16,23,36,33,32,38,39,34,37,35,54,49,55,53,51,50,48,52,60,57,59,61,56,58,62,63, 2,0, 7, 5, 4, 6, 1, 3,
		25,26,31,27,30,24,28,29,13,10, 9,14, 8,15,12,11,7, 2,0, 5, 3, 6, 4, 1,21,20,17,18,19,22,16,23,59,60,57,58,61,56,63,62,32,36,33,34,38,39,37,35,55,54,49,50,53,51,52,48,44,45,41,42,40,47,43,46,
		25,27,26,31,30,29,24,28,10, 9,13,14,15,12, 8,11,60,57,59,62,61,56,58,63,45,41,44,43,47,46,42,40,54,49,55,51,53,48,50,52,20,17,21,18,23,19,22,16, 2,0, 7, 1, 3, 4, 5, 6,36,33,32,35,38,39,34,37,
		25,27,28,24,31,26,30,29,33,32,36,35,34,39,37,38,9,13,10,14, 8,12,11,15,41,44,45,47,46,43,42,40,17,21,20,23,19,18,22,16, 0, 7, 2,3, 4, 1, 6, 5,49,55,54,53,48,51,52,50,57,59,60,62,58,56,63,61,
		25,27,29,30,24,28,31,26,59,60,57,62,63,56,61,58,32,36,33,35,37,39,38,34,44,45,41,46,43,47,42,40, 7, 2,0, 4, 1, 3, 6, 5,55,54,49,48,51,53,50,52,21,20,17,19,18,23,16,22,13,10, 9,14,11,12,15, 8,
		25,28,24,27,31,30,29,26,36,33,32,35,37,34,39,38,54,49,55,52,48,50,51,53, 2,0, 7, 3, 1, 6, 5, 4,10, 9,13,11,8,12,15,14,60,57,59,63,58,56,61,62,20,17,21,22,19,23,16,18,45,41,44,42,40,46,47,43,
		25,28,26,29,27,24,31,30,41,44,45,42,47,46,43,40,33,32,36,35,39,34,38,37, 0, 7, 2,1, 6, 3, 5, 4,57,59,60,58,56,63,61,62,17,21,20,19,23,22,18,16, 9,13,10, 8,12,11,14,15,49,55,54,52,51,50,53,48,
		25,28,30,31,29,26,27,24,55,54,49,52,53,50,48,51,44,45,41,42,43,46,40,47, 7, 2,0, 6, 3, 1, 5, 4,21,20,17,23,22,19,18,16,13,10, 9,12,11,8,15,14,59,60,57,56,63,58,62,61,32,36,33,35,38,34,37,39,
		25,29,28,26,27,30,24,31,44,45,41,42,46,43,47,40,59,60,57,62,56,63,58,61,55,54,49,53,50,48,52,51,13,10, 9,11,12,15, 8,14, 7, 2,0, 1, 4, 6, 3, 5,32,36,33,37,39,38,35,34,21,20,17,16,23,22,19,18,
		25,29,30,27,24,31,26,28,57,59,60,62,61,63,56,58,17,21,20,16,18,22,23,19,49,55,54,48,53,50,52,51,33,32,36,38,37,39,34,35, 9,13,10,15,11,12, 8,14, 0, 7, 2,6, 1, 4, 5, 3,41,44,45,42,40,43,46,47,
		25,29,31,24,26,28,27,30,20,17,21,16,19,22,18,23,45,41,44,42,47,43,40,46,54,49,55,50,48,53,52,51,2,0, 7, 4, 6, 1, 3, 5,36,33,32,39,38,37,34,35,10, 9,13,12,15,11,14, 8,60,57,59,62,58,63,61,56,
		25,30,26,24,28,31,29,27, 7, 2,0, 5, 6, 3, 1, 4,55,54,49,52,50,53,51,48,13,10, 9, 8,15,12,14,11,32,36,33,38,34,37,39,35,21,20,17,22,23,18,19,16,44,45,41,43,46,40,42,47,59,60,57,62,58,61,56,63,
		25,30,27,29,24,26,28,31,60,57,59,62,56,61,63,58,2,0, 7, 5, 1, 3, 4, 6,10, 9,13,15,12, 8,14,11,20,17,21,23,18,22,19,16,45,41,44,46,40,43,47,42,36,33,32,34,37,38,35,39,54,49,55,52,51,53,48,50,
		25,30,31,28,29,27,24,26,49,55,54,52,48,53,50,51,57,59,60,62,63,61,58,56, 9,13,10,12, 8,15,14,11,41,44,45,40,43,46,47,42,33,32,36,37,38,34,39,35,17,21,20,18,22,23,16,19, 0, 7, 2,5, 4, 3, 6, 1,
		25,31,24,29,26,27,30,28,21,20,17,16,18,19,22,23,13,10, 9,14,15, 8,11,12,32,36,33,39,37,34,35,38,44,45,41,40,47,43,46,42,59,60,57,61,58,63,56,62, 7, 2,0, 3, 6, 4, 5, 1,55,54,49,52,51,48,50,53,
		25,31,27,26,30,28,29,24, 9,13,10,14,12, 8,15,11,49,55,54,52,53,48,51,50,33,32,36,34,39,37,35,38,0, 7, 2,4, 3, 6, 1, 5,41,44,45,43,40,47,46,42,57,59,60,63,61,58,62,56,17,21,20,16,23,19,18,22,
		25,31,28,30,29,24,26,27,54,49,55,52,50,48,53,51,20,17,21,16,22,19,23,18,36,33,32,37,34,39,35,38,60,57,59,58,63,61,56,62, 2,0, 7, 6, 4, 3, 1, 5,45,41,44,47,43,40,42,46,10, 9,13,14,11,8,12,15,
		26,24,25,30,29,27,28,31,0, 5, 7, 2,1, 4, 6, 3,50,53,51,48,55,54,49,52,10,13,14,15, 8,11,9,12,34,38,39,36,32,35,33,37,18,19,22,17,16,21,20,23,42,41,45,47,40,46,44,43,58,61,56,63,59,60,57,62,
		26,24,27,29,28,31,30,25,51,50,53,48,52,54,55,49,56,58,61,63,62,60,59,57,14,10,13,11,15, 8,9,12,45,42,41,46,47,40,43,44,39,34,38,35,36,32,33,37,22,18,19,21,17,16,23,20, 7, 0, 5, 2,3, 4, 1, 6,
		26,24,31,28,30,25,29,27,61,56,58,63,57,60,62,59, 5, 7, 0, 2,6, 4, 3, 1,13,14,10, 8,11,15, 9,12,19,22,18,16,21,17,20,23,41,45,42,40,46,47,43,44,38,39,34,32,35,36,37,33,53,51,50,48,49,54,52,55,
		26,25,27,31,24,30,29,28,10,13,14, 9,15, 8,11,12, 0, 5, 7, 2,4, 1, 3, 6,18,19,22,21,20,17,23,16,58,61,56,59,60,57,62,63,34,38,39,32,36,33,35,37,50,53,51,55,54,49,48,52,42,41,45,44,46,43,47,40,
		26,25,28,29,31,27,24,30,41,45,42,44,47,43,40,46,13,14,10, 9,11,8,12,15,19,22,18,20,17,21,23,16,38,39,34,36,33,32,35,37,53,51,50,54,49,55,52,48,61,56,58,60,57,59,63,62, 5, 7, 0, 2,3, 1, 6, 4,
		26,25,30,24,29,28,31,27, 7, 0, 5, 2,6, 1, 4, 3,45,42,41,44,40,43,46,47,22,18,19,17,21,20,23,16,51,50,53,49,55,54,52,48,56,58,61,57,59,60,62,63,39,34,38,33,32,36,37,35,14,10,13, 9,12, 8,15,11,
		26,27,29,24,28,30,25,31,53,51,50,48,55,52,54,49,19,22,18,23,17,20,16,21,38,39,34,35,32,33,37,36,61,56,58,59,62,60,57,63, 5, 7, 0, 1, 3, 4, 6, 2,41,45,42,43,47,46,44,40,13,14,10, 9,12,15,11,8,
		26,27,30,28,25,31,24,29,18,19,22,23,21,20,17,16,10,13,14, 9, 8,15,12,11,34,38,39,33,35,32,37,36,42,41,45,46,43,47,40,44,58,61,56,60,59,62,57,63, 0, 5, 7, 4, 1, 3, 2,6,50,53,51,48,49,52,55,54,
		26,27,31,25,24,29,28,30,14,10,13, 9,11,15, 8,12,51,50,53,48,54,52,49,55,39,34,38,32,33,35,37,36, 7, 0, 5, 3, 4, 1, 6, 2,45,42,41,47,46,43,40,44,56,58,61,62,60,59,63,57,22,18,19,23,16,20,21,17,
		26,28,24,31,30,27,25,29,56,58,61,63,60,62,57,59,22,18,19,23,21,17,16,20,51,50,53,52,54,55,48,49,39,34,38,36,35,33,32,37,14,10,13, 8,12,11,15, 9, 7, 0, 5, 1, 6, 3, 2,4,45,42,41,44,46,47,40,43,
		26,28,27,30,25,29,31,24,19,22,18,23,20,17,21,16,41,45,42,44,43,47,46,40,53,51,50,55,52,54,48,49, 5, 7, 0, 3, 1, 6, 4, 2,38,39,34,33,36,35,32,37,13,14,10,11,8,12, 9,15,61,56,58,63,59,62,60,57,
		26,28,29,25,31,24,30,27,42,41,45,44,40,47,43,46,58,61,56,63,57,62,59,60,50,53,51,54,55,52,48,49,10,13,14,12,11,8,15, 9, 0, 5, 7, 6, 3, 1, 4, 2,34,38,39,35,33,36,37,32,18,19,22,23,16,17,20,21,
		26,29,24,27,28,25,31,30,50,53,51,48,54,55,52,49,42,41,45,44,47,40,46,43, 0, 5, 7, 1, 4, 6, 2,3,18,19,22,16,17,20,21,23,10,13,14,11,12,15, 8,9,58,61,56,57,62,59,63,60,34,38,39,37,36,32,35,33,
		26,29,25,28,31,30,27,24,45,42,41,44,43,40,47,46,39,34,38,37,33,32,36,35, 7, 0, 5, 6, 1, 4, 2,3,56,58,61,59,57,62,60,63,22,18,19,20,16,17,21,23,14,10,13,15,11,12, 9, 8,51,50,53,48,49,55,54,52,
		26,29,30,31,27,24,28,25,38,39,34,37,35,32,33,36,53,51,50,48,52,55,49,54, 5, 7, 0, 4, 6, 1, 2,3,13,14,10,12,15,11,8,9,61,56,58,62,59,57,60,63,19,22,18,17,20,16,23,21,41,45,42,44,46,40,43,47,
		26,30,24,25,29,31,27,28,5, 7, 0, 2,4, 6, 1, 3,38,39,34,37,32,35,36,33,61,56,58,57,60,62,63,59,41,45,42,46,40,43,47,44,13,14,10,15,12, 8,11,9,53,51,50,52,55,49,48,54,19,22,18,23,16,21,17,20,
		26,30,28,27,25,24,29,31,22,18,19,23,17,21,20,16, 7, 0, 5, 2,1, 6, 3, 4,56,58,61,60,62,57,63,59,14,10,13,12, 8,15,11,9,51,50,53,55,49,52,54,48,45,42,41,40,43,46,44,47,39,34,38,37,36,35,33,32,
		26,30,31,29,27,28,25,24,34,38,39,37,33,35,32,36,18,19,22,23,20,21,16,17,58,61,56,62,57,60,63,59,50,53,51,49,52,55,54,48,42,41,45,43,46,40,47,44,10,13,14, 8,15,12, 9,11,0, 5, 7, 2,3, 6, 4, 1,
		26,31,25,27,24,28,30,29,13,14,10, 9, 8,11,15,12,61,56,58,63,60,57,59,62,41,45,42,47,43,40,44,46,53,51,50,49,54,52,55,48,19,22,18,21,16,20,17,23, 5, 7, 0, 6, 4, 3, 2,1,38,39,34,37,36,33,32,35,
		26,31,28,24,30,29,27,25,58,61,56,63,62,57,60,59,34,38,39,37,35,33,36,32,42,41,45,40,47,43,44,46, 0, 5, 7, 3, 6, 4, 1, 2,50,53,51,52,49,54,55,48,18,19,22,20,21,16,23,17,10,13,14, 9,12,11,8,15,
		26,31,29,30,27,25,24,28,39,34,38,37,32,33,35,36,14,10,13, 9,15,11,12, 8,45,42,41,43,40,47,44,46,22,18,19,16,20,21,17,23, 7, 0, 5, 4, 3, 6, 1, 2,51,50,53,54,52,49,48,55,56,58,61,63,59,57,62,60,
		27,24,28,25,26,29,30,31,35,32,33,36,38,39,34,37,51,53,48,50,55,52,54,49, 1, 3, 4, 0, 2,5, 6, 7, 9,10,14, 8,11,15,12,13,59,62,60,56,61,63,58,57,18,23,19,16,21,17,22,20,43,47,46,40,42,44,41,45,
		27,24,29,26,30,31,25,28,48,51,53,50,49,52,55,54,46,43,47,40,45,44,42,41,4, 1, 3, 5, 0, 2,6, 7,19,18,23,17,16,21,20,22,14, 9,10,15, 8,11,12,13,60,59,62,63,56,61,57,58,33,35,32,36,37,39,38,34,
		27,24,31,30,25,28,26,29,47,46,43,40,41,44,45,42,32,33,35,36,34,39,37,38,3, 4, 1, 2,5, 0, 6, 7,62,60,59,61,63,56,58,57,23,19,18,21,17,16,20,22,10,14, 9,11,15, 8,13,12,53,48,51,50,54,52,49,55,
		27,25,24,28,26,31,29,30,32,33,35,36,39,34,38,37,10,14, 9,13,11,15, 8,12,47,46,43,41,44,45,40,42,23,19,18,17,21,20,16,22, 3, 4, 1, 0, 7, 2,5, 6,53,48,51,49,55,54,50,52,62,60,59,57,61,63,56,58,
		27,25,30,29,28,24,26,31,60,59,62,57,56,63,58,61,33,35,32,36,38,34,37,39,46,43,47,44,45,41,40,42, 4, 1, 3, 7, 2,0, 5, 6,48,51,53,55,54,49,52,50,19,18,23,21,20,17,22,16,14, 9,10,13, 8,15,12,11,
		27,25,31,26,29,30,28,24, 9,10,14,13,12,15,11,8,59,62,60,57,58,63,61,56,43,47,46,45,41,44,40,42,51,53,48,54,49,55,52,50,18,23,19,20,17,21,16,22, 1, 3, 4, 2,0, 7, 6, 5,35,32,33,36,37,34,39,38,
		27,26,24,29,30,28,31,25,51,53,48,50,52,55,49,54,18,23,19,22,16,21,17,20,35,32,33,38,39,34,36,37,59,62,60,61,56,58,63,57, 1, 3, 4, 5, 7, 0, 2,6,43,47,46,41,45,42,40,44, 9,10,14,13, 8,11,15,12,
		27,26,25,31,29,24,30,28,10,14, 9,13,15,11,12, 8,53,48,51,50,49,55,54,52,32,33,35,39,34,38,36,37, 3, 4, 1, 7, 0, 5, 2,6,47,46,43,45,42,41,44,40,62,60,59,56,58,61,57,63,23,19,18,22,17,21,20,16,
		27,26,28,30,31,25,29,24,19,18,23,22,20,21,16,17,14, 9,10,13,12,11,8,15,33,35,32,34,38,39,36,37,46,43,47,42,41,45,44,40,60,59,62,58,61,56,63,57, 4, 1, 3, 0, 5, 7, 6, 2,48,51,53,50,54,55,52,49,
		27,28,25,24,26,30,31,29,33,35,32,36,34,38,39,37,19,18,23,22,21,20,17,16,60,59,62,56,63,58,57,61,48,51,53,54,55,52,49,50,46,43,47,41,42,44,45,40,14, 9,10,12,11,8,13,15, 4, 1, 3, 6, 7, 2,0, 5,
		27,28,29,31,24,25,26,30, 1, 3, 4, 6, 0, 2,5, 7,35,32,33,36,39,38,37,34,59,62,60,63,58,56,57,61,43,47,46,42,44,41,45,40, 9,10,14,11,8,12,15,13,51,53,48,55,52,54,50,49,18,23,19,22,17,20,16,21,
		27,28,30,26,31,29,24,25,23,19,18,22,16,20,21,17, 3, 4, 1, 6, 5, 2,7, 0,62,60,59,58,56,63,57,61,10,14, 9, 8,12,11,15,13,53,48,51,52,54,55,49,50,47,46,43,44,41,42,40,45,32,33,35,36,37,38,34,39,
		27,29,25,30,28,31,24,26,59,62,60,57,63,58,56,61,1, 3, 4, 6, 2,0, 7, 5, 9,10,14,12,15,11,13, 8,18,23,19,17,20,16,21,22,43,47,46,44,42,45,41,40,35,32,33,39,38,37,36,34,51,53,48,50,54,49,55,52,
		27,29,26,24,30,25,28,31,53,48,51,50,55,49,52,54,62,60,59,57,56,58,61,63,10,14, 9,15,11,12,13, 8,47,46,43,42,45,44,41,40,32,33,35,38,37,39,34,36,23,19,18,20,16,17,22,21,3, 4, 1, 6, 7, 0, 5, 2,
		27,29,31,28,24,26,30,25, 4, 1, 3, 6, 5, 0, 2,7,48,51,53,50,52,49,54,55,14, 9,10,11,12,15,13, 8,33,35,32,37,39,38,34,36,19,18,23,16,17,20,21,22,46,43,47,45,44,42,40,41,60,59,62,57,61,58,63,56,
		27,30,24,31,25,29,28,26,46,43,47,40,44,45,41,42,60,59,62,57,63,56,61,58,48,51,53,49,52,55,50,54,14, 9,10, 8,15,12,11,13, 4, 1, 3, 2,7, 5, 0, 6,33,35,32,38,34,37,36,39,19,18,23,22,17,16,21,20,
		27,30,26,28,31,24,25,29,18,23,19,22,21,16,20,17,43,47,46,40,41,45,42,44,51,53,48,52,55,49,50,54, 1, 3, 4, 7, 5, 2,0, 6,35,32,33,34,37,38,39,36, 9,10,14,15,12, 8,13,11,59,62,60,57,61,56,58,63,
		27,30,29,25,28,26,31,24,62,60,59,57,58,56,63,61,23,19,18,22,20,16,17,21,53,48,51,55,49,52,50,54,32,33,35,37,38,34,39,36,10,14, 9,12, 8,15,11,13, 3, 4, 1, 5, 2,7, 6, 0,47,46,43,40,42,45,44,41,
		27,31,26,25,29,28,24,30,14, 9,10,13,11,12,15, 8,4, 1, 3, 6, 0, 5, 7, 2,19,18,23,20,21,16,22,17,60,59,62,61,58,63,56,57,33,35,32,39,37,34,38,36,48,51,53,52,49,54,50,55,46,43,47,40,42,41,45,44,
		27,31,28,29,24,30,25,26, 3, 4, 1, 6, 2,5, 0, 7,47,46,43,40,44,41,42,45,23,19,18,16,20,21,22,17,53,48,51,54,52,49,55,50,62,60,59,63,61,58,56,57,32,33,35,34,39,37,36,38,10,14, 9,13, 8,12,11,15,
		27,31,30,24,25,26,29,28,43,47,46,40,45,41,44,42, 9,10,14,13,15,12, 8,11,18,23,19,21,16,20,22,17,35,32,33,37,34,39,38,36,51,53,48,49,54,52,55,50,59,62,60,58,63,61,57,56, 1, 3, 4, 6, 7, 5, 2,0,
		28,24,25,27,30,29,31,26,36,35,33,32,37,38,34,39, 8,12,11,15, 9,13,10,14,44,41,42,46,47,40,45,43,19,23,22,21,17,16,20,18,1, 6, 3, 2,5, 0, 7, 4,52,54,55,50,51,48,49,53,58,56,63,61,57,59,60,62,
		28,24,26,31,27,25,30,29,56,63,58,61,60,59,62,57,35,33,36,32,34,38,39,37,41,42,44,47,40,46,45,43, 6, 3, 1, 5, 0, 2,7, 4,54,55,52,51,48,50,53,49,23,22,19,17,16,21,18,20,12,11,8,15,10,13,14, 9,
		28,24,29,30,31,26,27,25,11,8,12,15,14,13, 9,10,63,58,56,61,62,59,57,60,42,44,41,40,46,47,45,43,55,52,54,48,50,51,53,49,22,19,23,16,21,17,20,18,3, 1, 6, 0, 2,5, 4, 7,33,36,35,32,39,38,37,34,
		28,25,27,24,30,31,26,29,33,36,35,32,34,37,38,39,55,52,54,49,51,53,48,50, 3, 1, 6, 2,0, 7, 4, 5,11,8,12,10, 9,13,14,15,63,58,56,60,57,59,62,61,22,19,23,20,17,21,18,16,42,44,41,45,43,47,46,40,
		28,25,29,26,24,27,30,31,44,41,42,45,46,47,40,43,36,35,33,32,38,37,39,34, 1, 6, 3, 0, 7, 2,4, 5,58,56,63,57,59,60,62,61,19,23,22,17,21,20,16,18,8,12,11,9,13,10,15,14,52,54,55,49,48,53,50,51,
		28,25,31,30,26,29,24,27,54,55,52,49,50,53,51,48,41,42,44,45,40,47,43,46, 6, 3, 1, 7, 2,0, 4, 5,23,22,19,21,20,17,16,18,12,11,8,13,10, 9,14,15,56,63,58,59,60,57,61,62,35,33,36,32,39,37,34,38,
		28,26,25,29,24,31,27,30,41,42,44,45,47,40,46,43,56,63,58,61,59,60,57,62,54,55,52,50,53,51,49,48,12,11,8,10,13,14, 9,15, 6, 3, 1, 0, 5, 7, 2,4,35,33,36,34,38,39,32,37,23,22,19,18,21,20,17,16,
		28,26,30,27,29,25,24,31,22,19,23,18,17,20,16,21,42,44,41,45,46,40,43,47,55,52,54,53,51,50,49,48,3, 1, 6, 5, 7, 0, 2,4,33,36,35,38,39,34,37,32,11,8,12,13,14,10,15, 9,63,58,56,61,57,60,62,59,
		28,26,31,24,27,30,29,25,58,56,63,61,62,60,59,57,19,23,22,18,16,20,21,17,52,54,55,51,50,53,49,48,36,35,33,39,34,38,37,32, 8,12,11,14,10,13, 9,15, 1, 6, 3, 7, 0, 5, 4, 2,44,41,42,45,43,40,47,46,
		28,27,24,25,30,26,29,31,35,33,36,32,38,34,37,39,23,22,19,18,17,16,21,20,56,63,58,60,59,62,61,57,54,55,52,48,51,53,50,49,41,42,44,46,43,47,40,45,12,11,8,14, 9,10,15,13, 6, 3, 1, 4, 5, 0, 2,7,
		28,27,26,30,29,31,25,24,19,23,22,18,20,16,17,21,1, 6, 3, 4, 7, 0, 5, 2,58,56,63,62,60,59,61,57, 8,12,11,10,14, 9,13,15,52,54,55,53,48,51,50,49,44,41,42,47,46,43,45,40,36,35,33,32,39,34,38,37,
		28,27,31,29,25,24,30,26, 3, 1, 6, 4, 2,0, 7, 5,33,36,35,32,37,34,39,38,63,58,56,59,62,60,61,57,42,44,41,43,47,46,40,45,11,8,12, 9,10,14,13,15,55,52,54,51,53,48,49,50,22,19,23,18,21,16,20,17,
		28,29,26,25,24,30,31,27,42,44,41,45,40,46,47,43,11,8,12,15,13,14,10, 9,22,19,23,17,20,16,18,21,33,36,35,39,38,37,34,32,55,52,54,50,48,53,51,49,63,58,56,62,59,57,61,60, 3, 1, 6, 4, 5, 7, 0, 2,
		28,29,27,31,25,26,24,30, 1, 6, 3, 4, 0, 7, 2,5,44,41,42,45,47,46,43,40,19,23,22,20,16,17,18,21,52,54,55,48,53,50,51,49,58,56,63,59,57,62,60,61,36,35,33,38,37,39,32,34, 8,12,11,15,10,14, 9,13,
		28,29,30,24,31,27,25,26,12,11,8,15, 9,14,13,10, 6, 3, 1, 4, 2,7, 5, 0,23,22,19,16,17,20,18,21,56,63,58,57,62,59,60,61,35,33,36,37,39,38,34,32,54,55,52,53,50,48,49,51,41,42,44,45,43,46,40,47,
		28,30,24,29,31,25,26,27, 8,12,11,15,13, 9,14,10,52,54,55,49,50,51,48,53,36,35,33,37,38,34,32,39, 1, 6, 3, 5, 2,7, 0, 4,44,41,42,40,43,46,47,45,58,56,63,60,62,57,61,59,19,23,22,18,21,17,16,20,
		28,30,25,31,26,27,29,24,55,52,54,49,53,51,50,48,22,19,23,18,20,17,21,16,33,36,35,34,37,38,32,39,63,58,56,57,60,62,59,61,3, 1, 6, 7, 5, 2,0, 4,42,44,41,46,40,43,45,47,11,8,12,15,10, 9,13,14,
		28,30,27,26,29,24,31,25,23,22,19,18,16,17,20,21,12,11,8,15,14, 9,10,13,35,33,36,38,34,37,32,39,41,42,44,43,46,40,47,45,56,63,58,62,57,60,59,61,6, 3, 1, 2,7, 5, 4, 0,54,55,52,49,48,51,53,50,
		28,31,24,26,27,29,25,30,63,58,56,61,59,62,60,57, 3, 1, 6, 4, 0, 2,5, 7,11,8,12,14,13, 9,15,10,22,19,23,21,16,20,17,18,42,44,41,47,43,40,46,45,33,36,35,37,34,39,32,38,55,52,54,49,48,50,51,53,
		28,31,29,27,25,30,26,24, 6, 3, 1, 4, 7, 2,0, 5,54,55,52,49,53,50,48,51,12,11,8,9,14,13,15,10,35,33,36,39,37,34,38,32,23,22,19,20,21,16,17,18,41,42,44,40,47,43,45,46,56,63,58,61,57,62,59,60,
		28,31,30,25,26,24,27,29,52,54,55,49,51,50,53,48,58,56,63,61,60,62,57,59, 8,12,11,13, 9,14,15,10,44,41,42,43,40,47,46,45,36,35,33,34,39,37,38,32,19,23,22,16,20,21,18,17, 1, 6, 3, 4, 5, 2,7, 0,
		29,24,26,27,25,31,28,30,50,48,53,51,54,49,55,52,17,20,16,21,19,22,18,23,39,38,37,32,35,36,34,33,62,59,57,56,61,63,58,60, 4, 6, 1, 0, 2,5, 7, 3,44,42,45,40,46,47,41,43,12,15,11,8,13,14,10, 9,
		29,24,30,28,27,26,25,31,15,11,12, 8,10,14, 9,13,48,53,50,51,55,49,52,54,38,37,39,35,36,32,34,33, 6, 1, 4, 2,5, 0, 7, 3,42,45,44,46,47,40,43,41,59,57,62,61,63,56,60,58,20,16,17,21,18,22,23,19,
		29,24,31,25,28,30,27,26,16,17,20,21,23,22,19,18,11,12,15, 8,9,14,13,10,37,39,38,36,32,35,34,33,45,44,42,47,40,46,43,41,57,62,59,63,56,61,58,60, 1, 4, 6, 5, 0, 2,3, 7,53,50,48,51,52,49,54,55,
		29,25,24,31,28,26,30,27,17,20,16,21,22,19,23,18,44,42,45,41,40,46,47,43,50,48,53,54,49,55,51,52, 4, 6, 1, 2,0, 7, 5, 3,39,38,37,36,33,32,35,34,12,15,11,10, 9,13, 8,14,62,59,57,60,56,61,63,58,
		29,25,26,28,30,27,31,24,45,44,42,41,43,46,40,47,57,62,59,60,58,61,56,63,53,50,48,55,54,49,51,52,11,12,15,13,10, 9,14, 8,1, 4, 6, 7, 2,0, 5, 3,37,39,38,32,36,33,34,35,16,17,20,21,18,19,22,23,
		29,25,27,30,31,24,28,26,59,57,62,60,63,61,58,56,20,16,17,21,23,19,18,22,48,53,50,49,55,54,51,52,38,37,39,33,32,36,35,34,15,11,12, 9,13,10,14, 8,6, 1, 4, 0, 7, 2,3, 5,42,45,44,41,47,46,43,40,
		29,26,27,24,25,28,30,31,53,50,48,51,55,54,49,52,45,44,42,41,46,43,47,40, 1, 4, 6, 0, 5, 7, 3, 2,16,17,20,18,19,22,23,21,11,12,15,10,13,14, 9, 8,57,62,59,58,61,56,60,63,37,39,38,34,33,35,32,36,
		29,26,28,25,30,31,24,27,42,45,44,41,40,43,46,47,38,37,39,34,36,35,33,32, 6, 1, 4, 7, 0, 5, 3, 2,59,57,62,56,58,61,63,60,20,16,17,22,18,19,23,21,15,11,12,14,10,13, 8,9,48,53,50,51,52,54,55,49,
		29,26,31,30,24,27,25,28,39,38,37,34,32,35,36,33,50,48,53,51,49,54,52,55, 4, 6, 1, 5, 7, 0, 3, 2,12,15,11,13,14,10, 9, 8,62,59,57,61,56,58,63,60,17,20,16,19,22,18,21,23,44,42,45,41,47,43,40,46,
		29,27,24,26,25,30,31,28,48,53,50,51,49,55,54,52,59,57,62,60,61,63,56,58,15,11,12,10,14, 9, 8,13,42,45,44,47,46,43,40,41,38,37,39,32,33,35,36,34,20,16,17,23,19,18,21,22, 6, 1, 4, 3, 2,5, 0, 7,
		29,27,28,31,26,24,25,30, 1, 4, 6, 3, 0, 5, 7, 2,53,50,48,51,54,55,52,49,11,12,15,14, 9,10, 8,13,37,39,38,33,35,32,36,34,16,17,20,19,18,23,22,21,45,44,42,46,43,47,41,40,57,62,59,60,56,63,58,61,
		29,27,30,25,31,28,26,24,62,59,57,60,58,63,61,56, 4, 6, 1, 3, 7, 5, 2,0,12,15,11,9,10,14, 8,13,17,20,16,18,23,19,22,21,44,42,45,43,47,46,40,41,39,38,37,35,32,33,34,36,50,48,53,51,52,55,49,54,
		29,28,24,30,27,31,26,25,11,12,15, 8,14, 9,10,13, 1, 4, 6, 3, 5, 0, 2,7,16,17,20,23,22,19,21,18,57,62,59,56,63,58,61,60,37,39,38,35,33,36,32,34,53,50,48,54,55,52,51,49,45,44,42,41,47,40,46,43,
		29,28,25,26,30,24,27,31,44,42,45,41,46,40,43,47,12,15,11,8,10, 9,13,14,17,20,16,22,19,23,21,18,39,38,37,33,36,35,32,34,50,48,53,55,52,54,49,51,62,59,57,63,58,56,60,61,4, 6, 1, 3, 2,0, 7, 5,
		29,28,31,27,26,25,30,24, 6, 1, 4, 3, 7, 0, 5, 2,42,45,44,41,43,40,47,46,20,16,17,19,23,22,21,18,48,53,50,52,54,55,49,51,59,57,62,58,56,63,61,60,38,37,39,36,35,33,34,32,15,11,12, 8,13, 9,14,10,
		29,30,25,27,31,26,24,28,57,62,59,60,61,58,63,56,37,39,38,34,32,36,33,35,45,44,42,43,46,40,41,47, 1, 4, 6, 2,7, 5, 0, 3,53,50,48,49,52,55,54,51,16,17,20,22,23,18,21,19,11,12,15, 8,13,10, 9,14,
		29,30,26,31,24,28,27,25,38,37,39,34,35,36,32,33,15,11,12, 8,14,10,13, 9,42,45,44,40,43,46,41,47,20,16,17,18,22,23,19,21,6, 1, 4, 5, 2,7, 0, 3,48,53,50,55,49,52,51,54,59,57,62,60,56,58,61,63,
		29,30,28,24,27,25,31,26,12,15,11,8,9,10,14,13,62,59,57,60,63,58,56,61,44,42,45,46,40,43,41,47,50,48,53,52,55,49,54,51,17,20,16,23,18,22,19,21,4, 6, 1, 7, 5, 2,3, 0,39,38,37,34,33,36,35,32,
		29,31,25,24,28,27,26,30,20,16,17,21,19,23,22,18,6, 1, 4, 3, 0, 7, 2,5,59,57,62,63,61,58,60,56,15,11,12,13, 9,14,10, 8,48,53,50,54,52,49,55,51,42,45,44,43,40,47,41,46,38,37,39,34,33,32,36,35,
		29,31,27,28,26,30,24,25, 4, 6, 1, 3, 5, 7, 0, 2,39,38,37,34,35,32,33,36,62,59,57,58,63,61,60,56,44,42,45,47,43,40,46,41,12,15,11,14,13, 9,10, 8,50,48,53,49,54,52,51,55,17,20,16,21,18,23,19,22,
		29,31,30,26,24,25,28,27,37,39,38,34,36,32,35,33,16,17,20,21,22,23,18,19,57,62,59,61,58,63,60,56,53,50,48,52,49,54,55,51,45,44,42,40,47,43,46,41,11,12,15, 9,14,13, 8,10, 1, 4, 6, 3, 2,7, 5, 0,
		30,24,26,25,31,27,29,28,5, 2,7, 0, 4, 3, 6, 1,40,43,46,47,45,42,41,44,18,22,23,21,17,16,19,20,55,49,52,50,51,48,53,54,60,62,57,61,63,56,58,59,37,38,34,35,36,32,39,33,12, 8,15,11,14,10,13, 9,
		30,24,27,31,29,28,25,26,46,40,43,47,44,42,45,41,15,12, 8,11,9,10,14,13,23,18,22,16,21,17,19,20,34,37,38,32,35,36,33,39,52,55,49,48,50,51,53,54,57,60,62,56,61,63,59,58,7, 5, 2,0, 1, 3, 4, 6,
		30,24,28,29,25,26,31,27, 8,15,12,11,13,10, 9,14, 2,7, 5, 0, 6, 3, 1, 4,22,23,18,17,16,21,19,20,62,57,60,63,56,61,58,59,38,34,37,36,32,35,33,39,49,52,55,51,48,50,54,53,43,46,40,47,41,42,44,45,
		30,25,24,26,31,28,27,29, 2,7, 5, 0, 3, 6, 4, 1,49,52,55,54,51,48,50,53, 8,15,12,13,10, 9,11,14,38,34,37,32,36,33,35,39,22,23,18,21,20,17,16,19,43,46,40,44,45,41,47,42,62,57,60,59,63,56,61,58,
		30,25,28,31,27,29,26,24,55,49,52,54,53,48,51,50,60,62,57,59,58,56,63,61,12, 8,15, 9,13,10,11,14,40,43,46,41,44,45,42,47,37,38,34,33,32,36,35,39,18,22,23,17,21,20,19,16, 5, 2,7, 0, 1, 6, 3, 4,
		30,25,29,27,26,24,31,28,57,60,62,59,61,56,58,63, 7, 5, 2,0, 4, 6, 1, 3,15,12, 8,10, 9,13,11,14,23,18,22,20,17,21,16,19,46,40,43,45,41,44,42,47,34,37,38,36,33,32,39,35,52,55,49,54,50,48,53,51,
		30,26,25,24,31,29,28,27, 7, 5, 2,0, 6, 4, 3, 1,34,37,38,39,36,33,32,35,57,60,62,61,56,58,59,63,46,40,43,41,45,42,44,47,15,12, 8,13,14,10, 9,11,52,55,49,53,51,50,54,48,23,18,22,19,20,17,21,16,
		30,26,27,28,24,25,31,29,18,22,23,19,21,17,16,20, 5, 2,7, 0, 3, 4, 1, 6,60,62,57,56,58,61,59,63,12, 8,15,14,10,13, 9,11,55,49,52,51,50,53,48,54,40,43,46,45,42,41,47,44,37,38,34,39,32,33,35,36,
		30,26,29,31,28,27,24,25,38,34,37,39,35,33,36,32,22,23,18,19,16,17,20,21,62,57,60,58,61,56,59,63,49,52,55,50,53,51,48,54,43,46,40,42,41,45,44,47, 8,15,12,10,13,14,11,9, 2,7, 5, 0, 1, 4, 6, 3,
		30,27,25,29,26,28,24,31,60,62,57,59,56,58,61,63,18,22,23,19,17,21,20,16,55,49,52,53,48,51,54,50,37,38,34,32,33,35,36,39,12, 8,15,10,14, 9,13,11,5, 2,7, 3, 4, 1, 0, 6,40,43,46,47,41,44,45,42,
		30,27,28,26,24,31,29,25,23,18,22,19,16,21,17,20,46,40,43,47,42,44,41,45,52,55,49,51,53,48,54,50, 7, 5, 2,1, 3, 4, 6, 0,34,37,38,35,32,33,36,39,15,12, 8,9,10,14,11,13,57,60,62,59,63,58,56,61,
		30,27,31,24,29,25,26,28,43,46,40,47,45,44,42,41,62,57,60,59,61,58,63,56,49,52,55,48,51,53,54,50, 8,15,12,14, 9,10,13,11,2,7, 5, 4, 1, 3, 6, 0,38,34,37,33,35,32,39,36,22,23,18,19,20,21,16,17,
		30,28,26,27,24,29,25,31,22,23,18,19,17,16,21,20, 8,15,12,11,10,13,14, 9,38,34,37,35,33,36,39,32,43,46,40,41,42,44,45,47,62,57,60,56,63,58,61,59, 2,7, 5, 6, 3, 1, 0, 4,49,52,55,54,50,53,51,48,
		30,28,29,24,25,31,27,26,12, 8,15,11,9,13,10,14,55,49,52,54,48,53,50,51,37,38,34,36,35,33,39,32, 5, 2,7, 1, 6, 3, 4, 0,40,43,46,44,41,42,45,47,60,62,57,58,56,63,59,61,18,22,23,19,20,16,17,21,
		30,28,31,25,27,26,24,29,52,55,49,54,51,53,48,50,23,18,22,19,21,16,20,17,34,37,38,33,36,35,39,32,57,60,62,63,58,56,61,59, 7, 5, 2,3, 1, 6, 4, 0,46,40,43,42,44,41,47,45,15,12, 8,11,14,13, 9,10,
		30,29,24,28,25,27,26,31,15,12, 8,11,10, 9,13,14,57,60,62,59,56,61,63,58,46,40,43,44,42,45,47,41,52,55,49,50,48,53,51,54,23,18,22,17,20,16,21,19, 7, 5, 2,4, 6, 1, 0, 3,34,37,38,39,32,35,36,33,
		30,29,27,25,26,31,28,24,62,57,60,59,58,61,56,63,38,34,37,39,33,35,32,36,43,46,40,45,44,42,47,41,2,7, 5, 1, 4, 6, 3, 0,49,52,55,53,50,48,51,54,22,23,18,16,17,20,19,21,8,15,12,11,14, 9,10,13,
		30,29,31,26,28,24,25,27,37,38,34,39,36,35,33,32,12, 8,15,11,13, 9,14,10,40,43,46,42,45,44,47,41,18,22,23,20,16,17,21,19, 5, 2,7, 6, 1, 4, 3, 0,55,49,52,48,53,50,54,51,60,62,57,59,63,61,58,56,
		30,31,24,27,29,26,28,25,40,43,46,47,42,45,44,41,37,38,34,39,35,36,32,33, 5, 2,7, 4, 3, 6, 0, 1,60,62,57,63,61,58,56,59,18,22,23,16,20,21,17,19,12, 8,15,13, 9,14,11,10,55,49,52,54,50,51,48,53,
		30,31,25,28,27,24,29,26,49,52,55,54,48,51,53,50,43,46,40,47,44,45,41,42, 2,7, 5, 3, 6, 4, 0, 1,22,23,18,20,21,16,17,19, 8,15,12, 9,14,13,10,11,62,57,60,61,58,63,59,56,38,34,37,39,32,36,33,35,
		30,31,26,29,28,25,27,24,34,37,38,39,33,36,35,32,52,55,49,54,53,51,50,48,7, 5, 2,6, 4, 3, 0, 1,15,12, 8,14,13, 9,10,11,57,60,62,58,63,61,56,59,23,18,22,21,16,20,19,17,46,40,43,47,41,45,42,44,
		31,24,25,29,27,30,26,28,21,16,20,17,18,23,19,22,47,43,40,46,45,41,44,42,49,54,52,48,50,51,55,53, 6, 4, 3, 0, 2,5, 7, 1,37,34,39,32,35,36,33,38,14,13, 9, 8,11,15,10,12,58,63,61,56,60,57,59,62,
		31,24,28,26,29,25,27,30,63,61,58,56,59,57,62,60,16,20,21,17,19,23,22,18,54,52,49,50,51,48,55,53,34,39,37,35,36,32,33,38,13, 9,14,11,15, 8,12,10, 4, 3, 6, 2,5, 0, 1, 7,43,40,47,46,44,41,42,45,
		31,24,30,27,26,28,29,25,40,47,43,46,42,41,45,44,61,58,63,56,62,57,60,59,52,49,54,51,48,50,55,53, 9,14,13,15, 8,11,12,10, 3, 6, 4, 5, 0, 2,7, 1,39,37,34,36,32,35,38,33,20,21,16,17,22,23,18,19,
		31,25,26,27,28,30,24,29,13, 9,14,10, 8,12,11,15,54,52,49,55,51,50,53,48,34,39,37,33,32,36,38,35, 4, 3, 6, 0, 7, 2,5, 1,43,40,47,41,44,45,42,46,63,61,58,57,59,60,56,62,16,20,21,17,22,18,19,23,
		31,25,29,24,27,26,28,30,20,21,16,17,19,18,23,22, 9,14,13,10,11,12,15, 8,39,37,34,32,36,33,38,35,40,47,43,44,45,41,42,46,61,58,63,59,60,57,62,56, 3, 6, 4, 7, 2,0, 1, 5,52,49,54,55,53,50,48,51,
		31,25,30,28,24,29,27,26,49,54,52,55,48,50,51,53,21,16,20,17,23,18,22,19,37,34,39,36,33,32,38,35,58,63,61,60,57,59,62,56, 6, 4, 3, 2,0, 7, 5, 1,47,43,40,45,41,44,46,42,14,13, 9,10,15,12, 8,11,
		31,26,24,28,29,30,25,27,61,58,63,56,57,62,59,60,39,37,34,38,36,32,35,33,40,47,43,42,41,45,46,44, 3, 6, 4, 0, 5, 7, 2,1,52,49,54,50,53,51,48,55,20,21,16,18,19,22,17,23, 9,14,13,10,15, 8,11,12,
		31,26,27,25,28,24,29,30,14,13, 9,10,11,8,12,15,58,63,61,56,59,62,60,57,47,43,40,41,45,42,46,44,49,54,52,53,51,50,48,55,21,16,20,19,22,18,23,17, 6, 4, 3, 5, 7, 0, 1, 2,37,34,39,38,35,32,33,36,
		31,26,30,29,25,27,28,24,34,39,37,38,33,32,36,35,13, 9,14,10,12, 8,15,11,43,40,47,45,42,41,46,44,16,20,21,22,18,19,23,17, 4, 3, 6, 7, 0, 5, 2,1,54,52,49,51,50,53,55,48,63,61,58,56,60,62,57,59,
		31,27,24,30,26,25,28,29,47,43,40,46,41,45,42,44,14,13, 9,10, 8,11,15,12,21,16,20,18,23,19,17,22,37,34,39,35,32,33,36,38,49,54,52,51,53,48,50,55,58,63,61,59,62,60,56,57, 6, 4, 3, 1, 0, 2,5, 7,
		31,27,25,26,28,29,30,24, 9,14,13,10,12,11,8,15, 3, 6, 4, 1, 7, 2,0, 5,20,21,16,19,18,23,17,22,61,58,63,60,59,62,57,56,39,37,34,33,35,32,36,38,52,49,54,48,51,53,55,50,40,47,43,46,44,45,41,42,
		31,27,29,28,30,24,26,25, 4, 3, 6, 1, 5, 2,7, 0,43,40,47,46,42,45,44,41,16,20,21,23,19,18,17,22,54,52,49,53,48,51,50,55,63,61,58,62,60,59,57,56,34,39,37,32,33,35,38,36,13, 9,14,10,15,11,12, 8,
		31,28,25,30,24,26,29,27,54,52,49,55,50,51,48,53,63,61,58,56,57,59,60,62,13, 9,14, 8,12,11,10,15,43,40,47,44,41,42,45,46,34,39,37,36,35,33,32,38,16,20,21,19,23,22,17,18,4, 3, 6, 1, 0, 7, 2,5,
		31,28,26,24,29,27,30,25,58,63,61,56,62,59,57,60, 6, 4, 3, 1, 5, 7, 0, 2,14,13, 9,11,8,12,10,15,21,16,20,22,19,23,18,17,47,43,40,42,44,41,45,46,37,34,39,33,36,35,38,32,49,54,52,55,53,51,50,48,
		31,28,27,29,30,25,24,26, 3, 6, 4, 1, 2,7, 5, 0,52,49,54,55,48,51,53,50, 9,14,13,12,11,8,10,15,39,37,34,35,33,36,32,38,20,21,16,23,22,19,18,17,40,47,43,41,42,44,46,45,61,58,63,56,60,59,62,57,
		31,29,24,25,27,28,30,26,16,20,21,17,23,19,18,22, 4, 3, 6, 1, 2,5, 0, 7,63,61,58,59,57,62,56,60,13, 9,14,15,11,12, 8,10,54,52,49,48,53,50,51,55,43,40,47,42,45,44,46,41,34,39,37,38,35,36,32,33,
		31,29,26,30,25,24,27,28,39,37,34,38,32,36,33,35,20,21,16,17,18,19,22,23,61,58,63,57,62,59,56,60,52,49,54,53,50,48,51,55,40,47,43,45,44,42,41,46, 9,14,13,11,12,15,10, 8,3, 6, 4, 1, 0, 5, 7, 2,
		31,29,28,27,30,26,25,24, 6, 4, 3, 1, 7, 5, 2,0,37,34,39,38,33,36,35,32,58,63,61,62,59,57,56,60,47,43,40,44,42,45,41,46,14,13, 9,12,15,11,8,10,49,54,52,50,48,53,55,51,21,16,20,17,22,19,23,18,
		31,30,27,24,26,29,25,28,43,40,47,46,45,42,41,44,34,39,37,38,32,33,35,36, 4, 3, 6, 5, 2,7, 1, 0,63,61,58,60,62,57,59,56,16,20,21,18,22,23,19,17,13, 9,14,12, 8,15,10,11,54,52,49,55,53,48,51,50,
		31,30,28,25,24,27,26,29,52,49,54,55,51,48,50,53,40,47,43,46,41,42,44,45, 3, 6, 4, 2,7, 5, 1, 0,20,21,16,22,23,18,19,17, 9,14,13, 8,15,12,11,10,61,58,63,62,57,60,56,59,39,37,34,38,35,33,36,32,
		31,30,29,26,25,28,24,27,37,34,39,38,36,33,32,35,49,54,52,55,50,48,53,51,6, 4, 3, 7, 5, 2,1, 0,14,13, 9,15,12, 8,11,10,58,63,61,57,60,62,59,56,21,16,20,23,18,22,17,19,47,43,40,46,44,42,45,41,
		32,33,34,39,38,37,36,35, 0, 3, 7, 4, 5, 6, 2,1,48,53,55,49,52,54,50,51,16,23,21,22,17,19,18,20,56,62,59,63,57,60,58,61,40,47,44,42,41,46,43,45, 8,10,13,11,9,14,15,12,24,27,25,28,30,31,26,29,
		32,33,35,36,39,34,38,37,27,25,24,28,26,31,29,30, 3, 7, 0, 4, 2,6, 1, 5,23,21,16,17,19,22,18,20,47,44,40,41,46,42,43,45,10,13, 8,9,14,11,12,15,62,59,56,57,60,63,61,58,53,55,48,49,50,54,51,52,
		32,33,37,38,36,35,39,34,55,48,53,49,51,54,52,50,25,24,27,28,29,31,30,26,21,16,23,19,22,17,18,20,13, 8,10,14,11,9,12,15,59,56,62,60,63,57,58,61,44,40,47,46,42,41,45,43, 7, 0, 3, 4, 1, 6, 5, 2,
		32,34,36,38,35,37,33,39,10,13, 8,15,11,12, 9,14,23,21,16,18,19,17,20,22,47,44,40,43,42,46,45,41,62,59,56,63,58,57,60,61,53,55,48,54,50,52,51,49,27,25,24,31,26,30,28,29, 3, 7, 0, 4, 1, 5, 2,6,
		32,34,37,35,33,39,38,36,16,23,21,18,22,17,19,20, 0, 3, 7, 4, 6, 5, 1, 2,40,47,44,46,43,42,45,41,24,27,25,30,31,26,29,28,56,62,59,57,63,58,60,61,48,53,55,52,54,50,49,51,8,10,13,15,14,12,11,9,
		32,34,39,33,38,36,35,37, 7, 0, 3, 4, 2,5, 6, 1,13, 8,10,15, 9,12,14,11,44,40,47,42,46,43,45,41,55,48,53,50,52,54,51,49,25,24,27,26,30,31,29,28,59,56,62,58,57,63,61,60,21,16,23,18,20,17,22,19,
		32,35,34,37,33,36,39,38,23,21,16,18,17,19,22,20,27,25,24,28,31,26,30,29,10,13, 8,11,12, 9,15,14,53,55,48,50,54,51,52,49,47,44,40,46,41,43,42,45, 3, 7, 0, 2,6, 1, 4, 5,62,59,56,61,63,58,57,60,
		32,35,36,33,39,38,37,34,24,27,25,28,29,26,31,30,56,62,59,61,60,58,63,57, 8,10,13, 9,11,12,15,14, 0, 3, 7, 1, 2,6, 5, 4,48,53,55,51,50,54,52,49,40,47,44,43,46,41,45,42,16,23,21,18,20,19,17,22,
		32,35,38,39,37,34,33,36,59,56,62,61,57,58,60,63,21,16,23,18,22,19,20,17,13, 8,10,12, 9,11,15,14,44,40,47,41,43,46,42,45, 7, 0, 3, 6, 1, 2,5, 4,55,48,53,54,51,50,49,52,25,24,27,28,30,26,29,31,
		32,36,33,35,39,37,34,38,25,24,27,28,31,29,26,30,44,40,47,45,46,42,41,43,55,48,53,51,54,52,49,50,59,56,62,63,60,58,57,61,21,16,23,17,20,19,22,18,7, 0, 3, 5, 2,1, 4, 6,13, 8,10,15,14,11,9,12,
		32,36,37,39,34,38,35,33,47,44,40,45,43,42,46,41,10,13, 8,15,12,11,14, 9,53,55,48,52,51,54,49,50, 3, 7, 0, 1, 5, 2,6, 4,62,59,56,58,63,60,57,61,23,21,16,19,17,20,18,22,27,25,24,28,30,29,31,26,
		32,36,38,34,35,33,39,37, 8,10,13,15, 9,11,12,14,24,27,25,28,26,29,30,31,48,53,55,54,52,51,49,50,16,23,21,20,19,17,22,18,0, 3, 7, 2,1, 5, 6, 4,56,62,59,60,58,63,61,57,40,47,44,45,41,42,43,46,
		32,37,35,34,33,38,36,39,21,16,23,18,19,22,17,20,55,48,53,49,54,51,50,52,59,56,62,57,58,60,61,63, 7, 0, 3, 1, 6, 5, 2,4,13, 8,10,11,14,12, 9,15,25,24,27,29,31,30,28,26,44,40,47,45,41,43,46,42,
		32,37,38,33,36,39,34,35,53,55,48,49,52,51,54,50,47,44,40,45,42,43,41,46,62,59,56,60,57,58,61,63,27,25,24,30,29,31,26,28,3, 7, 0, 5, 1, 6, 2,4,10,13, 8,12,11,14,15, 9,23,21,16,18,20,22,19,17,
		32,37,39,36,34,35,33,38,40,47,44,45,46,43,42,41,16,23,21,18,17,22,20,19,56,62,59,58,60,57,61,63, 8,10,13,14,12,11,9,15,24,27,25,31,30,29,26,28,0, 3, 7, 6, 5, 1, 4, 2,48,53,55,49,50,51,52,54,
		32,38,33,37,36,34,35,39,48,53,55,49,54,52,51,50, 8,10,13,15,11,9,14,12, 0, 3, 7, 5, 6, 2,4, 1,40,47,44,41,42,43,46,45,16,23,21,19,20,22,17,18,24,27,25,26,29,30,28,31,56,62,59,61,63,57,60,58,
		32,38,34,36,35,39,37,33,13, 8,10,15,12, 9,11,14,59,56,62,61,58,57,63,60, 7, 0, 3, 2,5, 6, 4, 1,25,24,27,30,26,29,31,28,44,40,47,43,41,42,46,45,21,16,23,22,19,20,18,17,55,48,53,49,50,52,54,51,
		32,38,39,35,37,33,36,34,62,59,56,61,60,57,58,63,53,55,48,49,51,52,50,54, 3, 7, 0, 6, 2,5, 4, 1,23,21,16,20,22,19,17,18,27,25,24,29,30,26,31,28,47,44,40,42,43,41,45,46,10,13, 8,15,14, 9,12,11,
		32,39,33,34,38,35,37,36, 3, 7, 0, 4, 6, 2,5, 1,62,59,56,61,57,60,63,58,27,25,24,26,31,29,28,30,10,13, 8,14, 9,12,11,15,23,21,16,22,20,17,19,18,53,55,48,51,52,50,49,54,47,44,40,45,41,46,42,43,
		32,39,35,38,37,36,34,33,56,62,59,61,58,60,57,63,40,47,44,45,43,46,41,42,24,27,25,29,26,31,28,30,48,53,55,50,51,52,54,49, 8,10,13,12,14, 9,11,15,16,23,21,17,22,20,18,19, 0, 3, 7, 4, 1, 2,6, 5,
		32,39,36,37,34,33,38,35,44,40,47,45,42,46,43,41,7, 0, 3, 4, 5, 2,1, 6,25,24,27,31,29,26,28,30,21,16,23,20,17,22,19,18,55,48,53,52,50,51,54,49,13, 8,10, 9,12,14,15,11,59,56,62,61,63,60,58,57,
		33,32,36,35,34,39,37,38,25,27,28,24,31,26,30,29, 0, 4, 3, 7, 1, 5, 2,6,17,19,22,23,21,16,20,18,41,46,42,47,44,40,45,43, 9,14,11,10,13, 8,15,12,57,60,63,62,59,56,58,61,49,48,55,53,52,51,54,50,
		33,32,38,37,35,36,34,39,48,55,49,53,54,51,50,52,27,28,25,24,30,26,29,31,19,22,17,21,16,23,20,18,14,11,9,13, 8,10,15,12,60,63,57,59,56,62,61,58,46,42,41,44,40,47,43,45, 4, 3, 0, 7, 2,5, 6, 1,
		33,32,39,34,37,38,35,36, 3, 0, 4, 7, 6, 5, 1, 2,55,49,48,53,50,51,52,54,22,17,19,16,23,21,20,18,63,57,60,56,62,59,61,58,42,41,46,40,47,44,45,43,11,9,14, 8,10,13,12,15,28,25,27,24,29,26,31,30,
		33,34,32,39,37,36,38,35, 0, 4, 3, 7, 5, 1, 6, 2,57,60,63,58,62,59,56,61,25,27,28,31,26,30,24,29, 9,14,11,13,10,15, 8,12,17,19,22,16,18,23,21,20,49,48,55,54,50,52,53,51,41,46,42,43,47,44,40,45,
		33,34,35,38,39,32,37,36,46,42,41,43,40,44,45,47, 4, 3, 0, 7, 6, 1, 2,5,27,28,25,26,30,31,24,29,19,22,17,18,23,16,21,20,48,55,49,50,52,54,51,53,14,11,9,10,15,13,12, 8,60,63,57,58,56,59,61,62,
		33,34,36,37,38,35,39,32,63,57,60,58,61,59,62,56,42,41,46,43,45,44,47,40,28,25,27,30,31,26,24,29,55,49,48,52,54,50,51,53,11,9,14,15,13,10, 8,12,22,17,19,23,16,18,20,21,3, 0, 4, 7, 2,1, 5, 6,
		33,35,32,36,34,38,39,37,27,28,25,24,26,30,31,29,46,42,41,43,44,40,47,45,48,55,49,54,51,50,53,52,60,63,57,56,59,61,62,58,19,22,17,23,18,21,16,20, 4, 3, 0, 6, 1, 2,7, 5,14,11,9,12,13, 8,10,15,
		33,35,37,39,36,32,34,38,11,9,14,12,10, 8,15,13,28,25,27,24,31,30,29,26,55,49,48,51,50,54,53,52,22,17,19,18,21,23,16,20, 3, 0, 4, 1, 2,6, 5, 7,63,57,60,59,61,56,58,62,42,41,46,43,47,40,45,44,
		33,35,38,34,39,37,36,32,41,46,42,43,45,40,44,47, 9,14,11,12,15, 8,13,10,49,48,55,50,54,51,53,52, 0, 4, 3, 2,6, 1, 5, 7,57,60,63,61,56,59,62,58,17,19,22,21,23,18,20,16,25,27,28,24,29,30,26,31,
		33,36,35,32,34,37,38,39,28,25,27,24,30,31,26,29,63,57,60,58,59,61,56,62,11,9,14,10, 8,15,12,13, 3, 0, 4, 2,1, 5, 6, 7,55,49,48,54,52,51,50,53,42,41,46,45,44,47,43,40,22,17,19,20,18,21,23,16,
		33,36,37,34,38,39,32,35,60,63,57,58,62,61,59,56,19,22,17,20,16,21,18,23,14,11,9,15,10, 8,12,13,46,42,41,47,45,44,40,43, 4, 3, 0, 5, 2,1, 6, 7,48,55,49,51,54,52,53,50,27,28,25,24,29,31,30,26,
		33,36,39,38,32,35,34,37,17,19,22,20,23,21,16,18,25,27,28,24,26,31,29,30, 9,14,11,8,15,10,12,13,49,48,55,52,51,54,50,53,41,46,42,44,47,45,40,43, 0, 4, 3, 1, 5, 2,7, 6,57,60,63,58,56,61,62,59,
		33,37,32,38,35,39,36,34,55,49,48,53,51,50,54,52,11,9,14,12, 8,10,13,15, 3, 0, 4, 6, 5, 1, 7, 2,42,41,46,47,40,45,44,43,22,17,19,21,18,16,23,20,28,25,27,31,30,29,24,26,63,57,60,58,56,62,59,61,
		33,37,34,36,38,32,35,39,57,60,63,58,59,62,61,56,49,48,55,53,54,50,52,51,0, 4, 3, 5, 1, 6, 7, 2,17,19,22,18,16,21,23,20,25,27,28,30,29,31,26,24,41,46,42,40,45,47,43,44, 9,14,11,12,13,10,15, 8,
		33,37,39,35,36,34,38,32,14,11,9,12,15,10, 8,13,60,63,57,58,61,62,56,59, 4, 3, 0, 1, 6, 5, 7, 2,27,28,25,29,31,30,26,24,46,42,41,45,47,40,44,43,19,22,17,16,21,18,20,23,48,55,49,53,52,50,51,54,
		33,38,34,35,39,36,32,37,42,41,46,43,44,45,40,47,22,17,19,20,23,16,18,21,63,57,60,61,59,62,58,56,11,9,14,13,15, 8,10,12,28,25,27,26,29,30,31,24, 3, 0, 4, 5, 6, 2,7, 1,55,49,48,53,52,54,50,51,
		33,38,36,39,32,37,35,34,19,22,17,20,21,16,23,18,48,55,49,53,51,54,52,50,60,63,57,62,61,59,58,56, 4, 3, 0, 2,5, 6, 1, 7,14,11,9, 8,13,15,10,12,27,28,25,30,26,29,24,31,46,42,41,43,47,45,44,40,
		33,38,37,32,35,34,39,36,49,48,55,53,50,54,51,52,41,46,42,43,40,45,47,44,57,60,63,59,62,61,58,56,25,27,28,29,30,26,31,24, 0, 4, 3, 6, 2,5, 1, 7, 9,14,11,15, 8,13,12,10,17,19,22,20,18,16,21,23,
		33,39,34,32,37,35,36,38,4, 3, 0, 7, 1, 6, 5, 2,14,11,9,12,10,15,13, 8,46,42,41,40,44,45,43,47,48,55,49,52,50,51,54,53,27,28,25,31,29,26,30,24,60,63,57,61,62,56,58,59,19,22,17,20,18,23,16,21,
		33,39,35,37,36,38,32,34, 9,14,11,12, 8,15,10,13,17,19,22,20,21,23,18,16,41,46,42,45,40,44,43,47,57,60,63,56,61,62,59,58,49,48,55,51,52,50,54,53,25,27,28,26,31,29,24,30, 0, 4, 3, 7, 2,6, 1, 5,
		33,39,38,36,32,34,37,35,22,17,19,20,16,23,21,18,3, 0, 4, 7, 5, 6, 2,1,42,41,46,44,45,40,43,47,28,25,27,29,26,31,30,24,63,57,60,62,56,61,59,58,55,49,48,50,51,52,53,54,11,9,14,12,13,15, 8,10,
		34,32,33,39,36,38,37,35, 0, 7, 4, 3, 5, 2,1, 6,10,15,13, 8,14,11,9,12,42,46,43,44,40,47,41,45,50,52,54,55,48,53,49,51,26,30,31,25,24,27,28,29,58,57,63,59,56,62,60,61,18,23,16,21,19,22,17,20,
		34,32,35,37,39,33,36,38,23,16,18,21,17,22,20,19, 7, 4, 0, 3, 1, 2,6, 5,46,43,42,40,47,44,41,45,30,31,26,24,27,25,28,29,57,63,58,56,62,59,61,60,52,54,50,48,53,55,51,49,15,13,10, 8,9,11,12,14,
		34,32,38,36,37,35,39,33,13,10,15, 8,12,11,14, 9,16,18,23,21,20,22,19,17,43,42,46,47,44,40,41,45,63,58,57,62,59,56,61,60,54,50,52,53,55,48,49,51,31,26,30,27,25,24,29,28,4, 0, 7, 3, 6, 2,5, 1,
		34,33,37,36,35,38,32,39,57,63,58,60,59,61,56,62,46,43,42,41,47,40,45,44,30,31,26,28,25,27,29,24,52,54,50,55,49,48,53,51,15,13,10,11,9,14,12, 8,23,16,18,22,17,19,21,20, 7, 4, 0, 3, 6, 5, 1, 2,
		34,33,38,35,32,39,36,37,42,46,43,41,44,40,47,45, 0, 7, 4, 3, 2,5, 6, 1,26,30,31,27,28,25,29,24,18,23,16,19,22,17,20,21,50,52,54,48,55,49,53,51,10,15,13,14,11,9, 8,12,58,57,63,60,62,61,59,56,
		34,33,39,32,36,37,35,38,4, 0, 7, 3, 1, 5, 2,6,63,58,57,60,56,61,62,59,31,26,30,25,27,28,29,24,13,10,15, 9,14,11,12, 8,16,18,23,17,19,22,20,21,54,50,52,49,48,55,51,53,43,42,46,41,45,40,44,47,
		34,35,33,38,32,37,39,36,46,43,42,41,40,47,44,45,23,16,18,21,22,17,19,20,57,63,58,59,61,56,60,62,15,13,10, 9,11,12,14, 8,30,31,26,27,24,28,25,29, 7, 4, 0, 1, 2,6, 3, 5,52,54,50,51,55,49,48,53,
		34,35,36,39,38,33,32,37,54,50,52,51,48,49,53,55,43,42,46,41,44,47,45,40,63,58,57,61,56,59,60,62,31,26,30,24,28,27,25,29, 4, 0, 7, 2,6, 1, 5, 3,13,10,15,11,12, 9, 8,14,16,18,23,21,19,17,20,22,
		34,35,37,32,39,36,38,33,18,23,16,21,20,17,22,19,50,52,54,51,53,49,55,48,58,57,63,56,59,61,60,62, 0, 7, 4, 6, 1, 2,5, 3,10,15,13,12, 9,11,14, 8,26,30,31,28,27,24,29,25,42,46,43,41,45,47,40,44,
		34,36,32,38,37,33,35,39,10,15,13, 8,11,14,12, 9,58,57,63,60,59,56,62,61,0, 7, 4, 5, 2,1, 3, 6,26,30,31,24,25,28,27,29,42,46,43,47,45,44,40,41,18,23,16,17,20,19,21,22,50,52,54,51,55,48,53,49,
		34,36,33,37,35,39,38,32,63,58,57,60,61,56,59,62,54,50,52,51,49,48,55,53, 4, 0, 7, 1, 5, 2,3, 6,16,18,23,19,17,20,22,21,31,26,30,28,24,25,27,29,43,42,46,44,47,45,41,40,13,10,15, 8,9,14,11,12,
		34,36,39,35,38,32,37,33,52,54,50,51,53,48,49,55,15,13,10, 8,12,14, 9,11,7, 4, 0, 2,1, 5, 3, 6,46,43,42,45,44,47,40,41,23,16,18,20,19,17,22,21,30,31,26,25,28,24,29,27,57,63,58,60,62,56,61,59,
		34,37,32,35,39,38,33,36,16,18,23,21,22,20,17,19,31,26,30,29,27,25,24,28,13,10,15,12,11,14, 8,9,54,50,52,55,53,49,48,51,43,42,46,40,45,47,44,41,4, 0, 7, 5, 1, 6, 3, 2,63,58,57,60,62,59,56,61,
		34,37,36,33,35,32,39,38,58,57,63,60,56,59,61,62,18,23,16,21,17,20,19,22,10,15,13,11,14,12, 8,9,42,46,43,45,47,40,44,41,0, 7, 4, 1, 6, 5, 2,3,50,52,54,53,49,55,51,48,26,30,31,29,24,25,28,27,
		34,37,38,39,33,36,35,32,30,31,26,29,28,25,27,24,57,63,58,60,61,59,62,56,15,13,10,14,12,11,8,9, 7, 4, 0, 6, 5, 1, 2,3,52,54,50,49,55,53,48,51,46,43,42,47,40,45,41,44,23,16,18,21,19,20,22,17,
		34,38,35,33,32,36,37,39,43,42,46,41,47,44,40,45,13,10,15, 8,11,12, 9,14,54,50,52,48,49,53,51,55, 4, 0, 7, 6, 2,5, 1, 3,63,58,57,59,62,61,56,60,16,18,23,20,22,19,21,17,31,26,30,29,24,28,27,25,
		34,38,36,32,37,39,33,35,15,13,10, 8,14,12,11,9,30,31,26,29,25,28,24,27,52,54,50,53,48,49,51,55,23,16,18,19,20,22,17,21,7, 4, 0, 5, 6, 2,1, 3,57,63,58,61,59,62,60,56,46,43,42,41,45,44,47,40,
		34,38,39,37,33,35,32,36,26,30,31,29,27,28,25,24,42,46,43,41,40,44,45,47,50,52,54,49,53,48,51,55,58,57,63,62,61,59,56,60,18,23,16,22,19,20,17,21,0, 7, 4, 2,5, 6, 3, 1,10,15,13, 8,9,12,14,11,
		34,39,32,33,36,35,38,37, 7, 4, 0, 3, 2,1, 5, 6,52,54,50,51,48,53,55,49,23,16,18,17,22,20,21,19,57,63,58,62,56,61,59,60,46,43,42,44,45,40,47,41,15,13,10,12,14, 9, 8,11,30,31,26,29,24,27,25,28,
		34,39,35,36,38,37,33,32,50,52,54,51,49,53,48,55,26,30,31,29,28,27,24,25,18,23,16,20,17,22,21,19,10,15,13, 9,12,14,11,8,58,57,63,61,62,56,59,60,42,46,43,40,44,45,41,47, 0, 7, 4, 3, 6, 1, 2,5,
		34,39,37,38,33,32,36,35,31,26,30,29,25,27,28,24, 4, 0, 7, 3, 5, 1, 6, 2,16,18,23,22,20,17,21,19,43,42,46,45,40,44,47,41,13,10,15,14, 9,12,11,8,63,58,57,56,61,62,60,59,54,50,52,51,55,53,49,48,
		35,32,33,36,38,39,34,37,27,24,28,25,26,29,30,31,59,61,56,62,63,57,60,58,9,11,12, 8,10,13,14,15, 1, 2,6, 0, 3, 7, 4, 5,51,50,54,48,53,55,49,52,43,46,41,40,47,44,42,45,18,21,23,16,22,17,19,20,
		35,32,37,34,36,33,38,39,21,23,18,16,19,17,20,22,24,28,27,25,30,29,31,26,11,12, 9,10,13, 8,14,15,50,54,51,53,55,48,49,52,46,41,43,47,44,40,45,42, 2,6, 1, 3, 7, 0, 5, 4,61,56,59,62,60,57,58,63,
		35,32,39,38,34,37,36,33,56,59,61,62,58,57,63,60,23,18,21,16,20,17,22,19,12, 9,11,13, 8,10,14,15,41,43,46,44,40,47,45,42, 6, 1, 2,7, 0, 3, 4, 5,54,51,50,55,48,53,52,49,28,27,24,25,31,29,26,30,
		35,33,34,38,37,39,32,36,46,41,43,42,40,45,47,44,11,12, 9,14,13,10,15, 8,50,54,51,49,48,55,52,53, 2,6, 1, 0, 4, 3, 7, 5,61,56,59,57,60,63,58,62,21,23,18,17,19,22,16,20,24,28,27,25,31,26,30,29,
		35,33,36,32,38,34,37,39,28,27,24,25,30,26,29,31,41,43,46,42,47,45,44,40,54,51,50,48,55,49,52,53,56,59,61,60,63,57,58,62,23,18,21,19,22,17,20,16, 6, 1, 2,4, 3, 0, 5, 7,12, 9,11,14,15,10, 8,13,
		35,33,39,37,32,36,38,34, 9,11,12,14, 8,10,13,15,27,24,28,25,29,26,31,30,51,50,54,55,49,48,52,53,18,21,23,22,17,19,20,16, 1, 2,6, 3, 0, 4, 7, 5,59,61,56,63,57,60,62,58,43,46,41,42,44,45,40,47,
		35,34,32,37,36,39,33,38,23,18,21,16,17,20,19,22,54,51,50,52,55,48,53,49,56,59,61,58,57,63,62,60, 6, 1, 2,0, 7, 4, 3, 5,12, 9,11,10,15,13, 8,14,28,27,24,26,30,31,25,29,41,43,46,42,44,40,47,45,
		35,34,38,33,37,32,36,39,43,46,41,42,47,40,45,44,18,21,23,16,19,20,22,17,59,61,56,57,63,58,62,60, 9,11,12,15,13,10, 8,14,27,24,28,30,31,26,29,25, 1, 2,6, 7, 4, 0, 5, 3,51,50,54,52,53,48,49,55,
		35,34,39,36,33,38,37,32,50,54,51,52,49,48,55,53,46,41,43,42,45,40,44,47,61,56,59,63,58,57,62,60,24,28,27,31,26,30,29,25, 2,6, 1, 4, 0, 7, 3, 5,11,12, 9,13,10,15,14, 8,21,23,18,16,22,20,17,19,
		35,36,32,33,38,37,39,34,24,28,27,25,29,30,26,31,2,6, 1, 5, 3, 7, 0, 4,21,23,18,19,17,20,16,22,46,41,43,44,47,45,40,42,11,12, 9, 8,15,10,13,14,61,56,59,58,63,60,62,57,50,54,51,52,53,55,48,49,
		35,36,34,39,33,32,38,37,54,51,50,52,48,55,49,53,28,27,24,25,26,30,31,29,23,18,21,17,20,19,16,22,12, 9,11,15,10, 8,13,14,56,59,61,63,60,58,57,62,41,43,46,47,45,44,42,40, 6, 1, 2,5, 0, 7, 4, 3,
		35,36,37,38,39,34,33,32, 1, 2,6, 5, 4, 7, 3, 0,51,50,54,52,49,55,53,48,18,21,23,20,19,17,16,22,59,61,56,60,58,63,57,62,43,46,41,45,44,47,40,42, 9,11,12,10, 8,15,14,13,27,24,28,25,31,30,29,26,
		35,37,33,39,32,34,36,38,11,12, 9,14,10,13, 8,15,21,23,18,16,17,19,22,20,46,41,43,40,45,47,42,44,61,56,59,60,57,58,63,62,50,54,51,55,53,49,48,52,24,28,27,30,29,31,25,26, 2,6, 1, 5, 0, 4, 3, 7,
		35,37,34,32,36,38,39,33,18,21,23,16,20,19,17,22, 1, 2,6, 5, 7, 4, 0, 3,43,46,41,47,40,45,42,44,27,24,28,31,30,29,26,25,59,61,56,58,60,57,63,62,51,50,54,49,55,53,52,48,9,11,12,14,15,13,10, 8,
		35,37,38,36,39,33,32,34, 6, 1, 2,5, 3, 4, 7, 0,12, 9,11,14, 8,13,15,10,41,43,46,45,47,40,42,44,54,51,50,53,49,55,48,52,28,27,24,29,31,30,26,25,56,59,61,57,58,60,62,63,23,18,21,16,22,19,20,17,
		35,38,32,39,34,33,37,36,59,61,56,62,57,63,58,60,43,46,41,42,40,47,44,45,27,24,28,26,29,30,25,31,51,50,54,53,48,49,55,52, 9,11,12,13,15, 8,10,14,18,21,23,19,20,22,16,17, 1, 2,6, 5, 0, 3, 7, 4,
		35,38,33,34,37,36,39,32,41,43,46,42,45,47,40,44, 6, 1, 2,5, 4, 3, 0, 7,28,27,24,30,26,29,25,31,23,18,21,22,19,20,17,16,54,51,50,49,53,48,55,52,12, 9,11,8,13,15,14,10,56,59,61,62,60,63,57,58,
		35,38,36,37,39,32,34,33, 2,6, 1, 5, 7, 3, 4, 0,61,56,59,62,58,63,60,57,24,28,27,29,30,26,25,31,11,12, 9,15, 8,13,10,14,21,23,18,20,22,19,17,16,50,54,51,48,49,53,52,55,46,41,43,42,44,47,45,40,
		35,39,36,34,33,37,32,38,51,50,54,52,55,49,48,53, 9,11,12,14,10, 8,15,13, 1, 2,6, 4, 7, 3, 5, 0,43,46,41,44,45,40,47,42,18,21,23,17,22,20,19,16,27,24,28,29,26,31,25,30,59,61,56,62,60,58,63,57,
		35,39,37,33,32,38,34,36,12, 9,11,14,13, 8,10,15,56,59,61,62,57,58,60,63, 6, 1, 2,3, 4, 7, 5, 0,28,27,24,31,29,26,30,25,41,43,46,40,44,45,47,42,23,18,21,20,17,22,16,19,54,51,50,52,53,49,55,48,
		35,39,38,32,34,36,33,37,61,56,59,62,63,58,57,60,50,54,51,52,48,49,53,55, 2,6, 1, 7, 3, 4, 5, 0,21,23,18,22,20,17,19,16,24,28,27,26,31,29,30,25,46,41,43,45,40,44,42,47,11,12, 9,14,15, 8,13,10,
		36,32,34,38,33,35,37,39,10, 8,15,13,11,9,14,12,25,28,24,27,30,31,26,29,54,52,51,48,53,55,50,49,20,19,17,16,23,21,18,22, 2,1, 5, 0, 3, 7, 4, 6,60,58,63,56,62,59,57,61,45,44,47,40,46,43,42,41,
		36,32,35,33,37,39,38,34,24,25,28,27,29,31,30,26,47,45,44,40,41,43,46,42,51,54,52,55,48,53,50,49,63,60,58,59,56,62,61,57,17,20,19,21,16,23,18,22, 5, 2,1, 7, 0, 3, 6, 4,15,10, 8,13,12, 9,11,14,
		36,32,39,37,38,34,33,35,44,47,45,40,42,43,41,46, 8,15,10,13,14, 9,12,11,52,51,54,53,55,48,50,49, 1, 5, 2,3, 7, 0, 4, 6,58,63,60,62,59,56,61,57,19,17,20,23,21,16,22,18,28,24,25,27,26,31,29,30,
		36,33,32,35,37,34,39,38,25,28,24,27,31,30,29,26,60,58,63,57,56,62,59,61,10, 8,15,11,9,14,13,12, 2,1, 5, 3, 0, 4, 7, 6,54,52,51,55,49,48,53,50,45,44,47,42,41,46,40,43,20,19,17,22,16,23,21,18,
		36,33,34,37,39,38,35,32,63,60,58,57,61,62,56,59,17,20,19,22,18,23,16,21,15,10, 8,14,11,9,13,12,47,45,44,46,42,41,43,40, 5, 2,1, 4, 3, 0, 7, 6,51,54,52,48,55,49,50,53,24,25,28,27,26,30,31,29,
		36,33,38,39,35,32,37,34,19,17,20,22,21,23,18,16,28,24,25,27,29,30,26,31,8,15,10, 9,14,11,13,12,52,51,54,49,48,55,53,50,44,47,45,41,46,42,43,40, 1, 5, 2,0, 4, 3, 6, 7,58,63,60,57,59,62,61,56,
		36,34,35,39,32,38,33,37,54,52,51,50,48,53,55,49,10, 8,15,13, 9,11,12,14, 2,1, 5, 7, 4, 0, 6, 3,45,44,47,46,43,42,41,40,20,19,17,23,16,18,21,22,25,28,24,30,31,26,27,29,60,58,63,57,59,61,56,62,
		36,34,37,33,39,35,32,38,58,63,60,57,56,61,62,59,52,51,54,50,55,53,49,48,1, 5, 2,4, 0, 7, 6, 3,19,17,20,16,18,23,21,22,28,24,25,31,26,30,29,27,44,47,45,43,42,46,40,41,8,15,10,13,12,11,14, 9,
		36,34,38,32,33,37,39,35,15,10, 8,13,14,11,9,12,63,60,58,57,62,61,59,56, 5, 2,1, 0, 7, 4, 6, 3,24,25,28,26,30,31,29,27,47,45,44,42,46,43,41,40,17,20,19,18,23,16,22,21,51,54,52,50,49,53,48,55,
		36,35,33,32,37,38,34,39,28,24,25,27,30,29,31,26, 1, 5, 2,6, 0, 4, 3, 7,19,17,20,21,23,18,22,16,44,47,45,46,41,43,42,40, 8,15,10,11,12, 9,14,13,58,63,60,61,56,59,57,62,52,51,54,50,49,48,55,53,
		36,35,38,37,34,39,32,33, 2,1, 5, 6, 7, 4, 0, 3,54,52,51,50,53,48,49,55,20,19,17,18,21,23,22,16,60,58,63,59,61,56,62,57,45,44,47,43,46,41,42,40,10, 8,15, 9,11,12,13,14,25,28,24,27,26,29,30,31,
		36,35,39,34,32,33,37,38,51,54,52,50,55,48,53,49,24,25,28,27,31,29,26,30,17,20,19,23,18,21,22,16,15,10, 8,12, 9,11,14,13,63,60,58,56,59,61,62,57,47,45,44,41,43,46,40,42, 5, 2,1, 6, 3, 4, 7, 0,
		36,37,32,39,38,35,34,33,47,45,44,40,43,41,42,46, 5, 2,1, 6, 7, 0, 3, 4,24,25,28,29,31,30,27,26,17,20,19,16,21,18,23,22,51,54,52,53,49,55,48,50,15,10, 8,11,14,12,13, 9,63,60,58,57,59,56,62,61,
		36,37,33,34,39,32,38,35,60,58,63,57,62,56,61,59,45,44,47,40,42,41,46,43,25,28,24,31,30,29,27,26,54,52,51,49,55,53,48,50,10, 8,15,14,12,11,9,13,20,19,17,21,18,16,22,23, 2,1, 5, 6, 3, 0, 4, 7,
		36,37,35,38,34,33,39,32, 1, 5, 2,6, 4, 0, 7, 3,58,63,60,57,61,56,59,62,28,24,25,30,29,31,27,26, 8,15,10,12,11,14, 9,13,19,17,20,18,16,21,23,22,52,51,54,55,53,49,50,48,44,47,45,40,46,41,43,42,
		36,38,32,34,33,39,35,37, 8,15,10,13, 9,14,11,12,19,17,20,22,23,21,16,18,44,47,45,42,43,41,40,46,58,63,60,59,62,61,56,57,52,51,54,48,49,53,55,50,28,24,25,29,30,26,27,31,1, 5, 2,6, 3, 7, 0, 4,
		36,38,37,35,34,32,33,39, 5, 2,1, 6, 0, 7, 4, 3,15,10, 8,13,11,14,12, 9,47,45,44,43,41,42,40,46,51,54,52,49,53,48,55,50,24,25,28,30,26,29,31,27,63,60,58,62,61,59,57,56,17,20,19,22,16,21,18,23,
		36,38,39,33,35,37,34,32,20,19,17,22,18,21,23,16, 2,1, 5, 6, 4, 7, 3, 0,45,44,47,41,42,43,40,46,25,28,24,26,29,30,31,27,60,58,63,61,59,62,56,57,54,52,51,53,48,49,50,55,10, 8,15,13,12,14, 9,11,
		36,39,33,38,35,34,32,37,17,20,19,22,23,18,21,16,51,54,52,50,48,55,49,53,63,60,58,61,62,56,57,59, 5, 2,1, 3, 4, 7, 0, 6,15,10, 8,9,12,14,11,13,24,25,28,31,29,26,27,30,47,45,44,40,46,42,41,43,
		36,39,34,35,32,37,38,33,52,51,54,50,53,55,48,49,44,47,45,40,43,42,46,41,58,63,60,56,61,62,57,59,28,24,25,26,31,29,30,27, 1, 5, 2,7, 3, 4, 0, 6, 8,15,10,14, 9,12,13,11,19,17,20,22,16,18,23,21,
		36,39,37,32,38,33,35,34,45,44,47,40,41,42,43,46,20,19,17,22,21,18,16,23,60,58,63,62,56,61,57,59,10, 8,15,12,14, 9,11,13,25,28,24,29,26,31,30,27, 2,1, 5, 4, 7, 3, 6, 0,54,52,51,50,49,55,53,48,
		37,32,33,38,39,36,35,34,55,53,49,48,51,52,50,54,40,45,47,44,41,46,42,43,60,57,58,62,59,56,63,61,30,29,31,27,25,24,28,26, 5, 1, 6, 3, 7, 0, 4, 2,12,11,14,10,13, 8,9,15,18,16,21,23,17,19,22,20,
		37,32,34,35,38,33,39,36,16,21,18,23,22,19,20,17,53,49,55,48,50,52,54,51,57,58,60,59,56,62,63,61,1, 6, 5, 7, 0, 3, 4, 2,11,14,12,13, 8,10,15, 9,29,31,30,25,24,27,26,28,45,47,40,44,42,46,43,41,
		37,32,36,39,35,34,38,33,47,40,45,44,43,46,41,42,21,18,16,23,20,19,17,22,58,60,57,56,62,59,63,61,14,12,11,8,10,13,15, 9,31,30,29,24,27,25,28,26, 6, 5, 1, 0, 3, 7, 2,4,49,55,53,48,54,52,51,50,
		37,33,35,39,34,36,32,38,11,14,12, 9,10,15,13, 8,57,58,60,63,56,59,61,62, 1, 6, 5, 4, 3, 0, 2,7,29,31,30,27,28,25,24,26,45,47,40,46,42,41,43,44,16,21,18,19,22,17,23,20,53,49,55,48,54,51,50,52,
		37,33,36,34,32,38,39,35,60,57,58,63,62,59,56,61,55,53,49,48,52,51,54,50, 5, 1, 6, 0, 4, 3, 2,7,18,16,21,17,19,22,20,23,30,29,31,25,27,28,24,26,40,45,47,41,46,42,44,43,12,11,14, 9, 8,15,10,13,
		37,33,38,32,39,35,34,36,49,55,53,48,50,51,52,54,14,12,11,9,13,15, 8,10, 6, 5, 1, 3, 0, 4, 2,7,47,40,45,42,41,46,43,44,21,18,16,22,17,19,20,23,31,30,29,28,25,27,26,24,58,60,57,63,61,59,62,56,
		37,34,33,36,32,35,38,39,57,58,60,63,59,56,62,61,16,21,18,23,19,22,17,20,11,14,12,10,15,13, 9, 8,45,47,40,42,46,43,41,44, 1, 6, 5, 0, 7, 4, 3, 2,53,49,55,50,52,54,48,51,29,31,30,26,27,28,25,24,
		37,34,35,32,38,39,36,33,18,16,21,23,20,22,19,17,30,29,31,26,24,28,27,25,12,11,14,13,10,15, 9, 8,55,53,49,54,50,52,51,48,40,45,47,43,42,46,41,44, 5, 1, 6, 4, 0, 7, 2,3,60,57,58,63,61,56,59,62,
		37,34,39,38,36,33,32,35,31,30,29,26,25,28,24,27,58,60,57,63,62,56,61,59,14,12,11,15,13,10, 9, 8,6, 5, 1, 7, 4, 0, 3, 2,49,55,53,52,54,50,51,48,47,40,45,46,43,42,44,41,21,18,16,23,17,22,20,19,
		37,35,32,34,38,36,33,39,21,18,16,23,19,20,22,17, 6, 5, 1, 2,0, 3, 7, 4,47,40,45,43,46,41,44,42,31,30,29,27,24,28,25,26,58,60,57,59,61,56,62,63,49,55,53,51,50,54,48,52,14,12,11,9, 8,10,13,15,
		37,35,36,38,33,39,34,32, 1, 6, 5, 2,4, 3, 0, 7,11,14,12, 9,15,10, 8,13,45,47,40,41,43,46,44,42,53,49,55,54,51,50,52,48,29,31,30,28,27,24,25,26,57,58,60,56,59,61,63,62,16,21,18,23,17,20,19,22,
		37,35,39,33,34,32,38,36,12,11,14, 9,13,10,15, 8,18,16,21,23,22,20,17,19,40,45,47,46,41,43,44,42,60,57,58,61,56,59,62,63,55,53,49,50,54,51,52,48,30,29,31,24,28,27,26,25, 5, 1, 6, 2,7, 3, 4, 0,
		37,36,34,33,32,39,35,38,58,60,57,63,56,62,59,61,47,40,45,44,46,43,42,41,31,30,29,25,28,24,26,27,49,55,53,54,52,51,50,48,14,12,11,10, 8,15,13, 9,21,18,16,20,19,17,23,22, 6, 5, 1, 2,7, 4, 0, 3,
		37,36,38,35,33,34,32,39, 5, 1, 6, 2,0, 4, 3, 7,60,57,58,63,59,62,61,56,30,29,31,28,24,25,26,27,12,11,14, 8,15,10,13, 9,18,16,21,19,17,20,22,23,55,53,49,52,51,54,48,50,40,45,47,44,42,43,41,46,
		37,36,39,32,35,38,33,34,45,47,40,44,41,43,46,42, 1, 6, 5, 2,3, 4, 7, 0,29,31,30,24,25,28,26,27,16,21,18,17,20,19,22,23,53,49,55,51,54,52,50,48,11,14,12,15,10, 8,9,13,57,58,60,63,61,62,56,59,
		37,38,32,33,39,34,36,35,53,49,55,48,52,50,51,54,29,31,30,26,25,24,27,28,16,21,18,22,19,20,23,17,11,14,12, 8,13,15,10, 9,57,58,60,62,61,59,56,63,45,47,40,43,41,42,44,46, 1, 6, 5, 2,7, 0, 3, 4,
		37,38,34,39,36,35,33,32,30,29,31,26,28,24,25,27, 5, 1, 6, 2,4, 0, 7, 3,18,16,21,20,22,19,23,17,40,45,47,42,43,41,46,44,12,11,14,15, 8,13,10, 9,60,57,58,59,62,61,63,56,55,53,49,48,54,50,52,51,
		37,38,35,36,33,32,39,34, 6, 5, 1, 2,3, 0, 4, 7,49,55,53,48,51,50,54,52,21,18,16,19,20,22,23,17,58,60,57,61,59,62,56,63,47,40,45,41,42,43,46,44,14,12,11,13,15, 8,9,10,31,30,29,26,27,24,28,25,
		37,39,32,36,35,33,34,38,40,45,47,44,46,41,43,42,12,11,14, 9,10,13, 8,15,55,53,49,51,52,50,48,54, 5, 1, 6, 7, 3, 4, 0, 2,60,57,58,56,61,62,59,63,18,16,21,22,20,17,23,19,30,29,31,26,27,25,24,28,
		37,39,33,35,34,38,36,32,14,12,11,9,15,13,10, 8,31,30,29,26,28,25,27,24,49,55,53,50,51,52,48,54,21,18,16,17,22,20,19,23, 6, 5, 1, 4, 7, 3, 0, 2,58,60,57,62,56,61,63,59,47,40,45,44,42,41,46,43,
		37,39,38,34,36,32,35,33,29,31,30,26,24,25,28,27,45,47,40,44,43,41,42,46,53,49,55,52,50,51,48,54,57,58,60,61,62,56,59,63,16,21,18,20,17,22,19,23, 1, 6, 5, 3, 4, 7, 2,0,11,14,12, 9, 8,13,15,10,
		38,32,35,39,33,37,34,36,59,62,61,56,57,60,63,58,48,49,53,55,50,54,51,52, 6, 2,5, 3, 7, 0, 1, 4,20,22,19,23,21,16,18,17,29,30,26,27,25,24,28,31,42,43,41,47,44,40,46,45,15, 8,13,10,11,12, 9,14,
		38,32,36,34,39,35,33,37, 8,13,15,10, 9,12,14,11,62,61,59,56,63,60,58,57, 2,5, 6, 7, 0, 3, 1, 4,30,26,29,25,24,27,28,31,43,41,42,44,40,47,45,46,22,19,20,21,16,23,17,18,49,53,48,55,51,54,52,50,
		38,32,37,33,34,36,39,35,53,48,49,55,52,54,50,51,13,15, 8,10,14,12,11,9, 5, 6, 2,0, 3, 7, 1, 4,41,42,43,40,47,44,45,46,19,20,22,16,23,21,18,17,26,29,30,24,27,25,31,28,61,59,62,56,58,60,57,63,
		38,33,32,37,34,35,36,39,48,49,53,55,54,50,52,51,42,43,41,46,47,44,40,45,59,62,61,57,60,63,56,58,29,30,26,25,27,28,24,31,6, 2,5, 0, 4, 3, 7, 1,15, 8,13, 9,14,11,10,12,20,22,19,17,23,21,16,18,
		38,33,35,34,36,39,37,32,41,42,43,46,45,44,47,40,19,20,22,17,18,21,23,16,61,59,62,63,57,60,56,58,13,15, 8,11,9,14,12,10,26,29,30,28,25,27,24,31,5, 6, 2,3, 0, 4, 1, 7,53,48,49,55,51,50,54,52,
		38,33,39,36,37,32,34,35,22,19,20,17,16,21,18,23,49,53,48,55,52,50,51,54,62,61,59,60,63,57,56,58,2,5, 6, 4, 3, 0, 7, 1, 8,13,15,14,11,9,12,10,30,26,29,27,28,25,31,24,43,41,42,46,40,44,45,47,
		38,34,32,36,39,37,35,33,13,15, 8,10,12,14, 9,11,26,29,30,31,24,27,25,28,53,48,49,52,54,50,55,51,19,20,22,23,16,18,21,17, 5, 6, 2,7, 4, 0, 3, 1,61,59,62,57,63,58,56,60,41,42,43,46,40,47,44,45,
		38,34,33,35,36,32,39,37,42,43,41,46,44,47,45,40,15, 8,13,10, 9,14,11,12,48,49,53,54,50,52,55,51,6, 2,5, 4, 0, 7, 3, 1,59,62,61,63,58,57,60,56,20,22,19,16,18,23,17,21,29,30,26,31,25,27,28,24,
		38,34,37,39,35,33,36,32,30,26,29,31,28,27,24,25,43,41,42,46,45,47,40,44,49,53,48,50,52,54,55,51,62,61,59,58,57,63,60,56,22,19,20,18,23,16,21,17, 2,5, 6, 0, 7, 4, 1, 3, 8,13,15,10,11,14,12, 9,
		38,35,34,33,36,37,32,39,43,41,42,46,47,45,44,40, 2,5, 6, 1, 0, 7, 4, 3,30,26,29,28,27,24,31,25,22,19,20,23,18,21,16,17,49,53,48,54,51,50,52,55, 8,13,15,12, 9,11,10,14,62,61,59,56,58,57,63,60,
		38,35,37,36,32,39,33,34, 6, 2,5, 1, 3, 7, 0, 4,59,62,61,56,60,57,58,63,29,30,26,24,28,27,31,25,15, 8,13,11,12, 9,14,10,20,22,19,21,23,18,16,17,48,49,53,50,54,51,55,52,42,43,41,46,40,45,47,44,
		38,35,39,32,33,34,36,37,61,59,62,56,63,57,60,58,41,42,43,46,44,45,40,47,26,29,30,27,24,28,31,25,53,48,49,51,50,54,52,55,13,15, 8,9,11,12,14,10,19,20,22,18,21,23,17,16, 5, 6, 2,1, 4, 7, 3, 0,
		38,36,33,39,37,35,32,34,19,20,22,17,21,18,16,23, 5, 6, 2,1, 3, 0, 4, 7,41,42,43,45,44,47,46,40,26,29,30,25,28,24,27,31,61,59,62,60,58,63,57,56,53,48,49,54,52,51,55,50,13,15, 8,10,11,9,14,12,
		38,36,34,32,39,33,37,35,15, 8,13,10,14, 9,12,11,20,22,19,17,16,18,23,21,42,43,41,44,47,45,46,40,59,62,61,58,63,60,57,56,48,49,53,52,51,54,50,55,29,30,26,28,24,25,31,27, 6, 2,5, 1, 4, 0, 7, 3,
		38,36,35,37,32,34,39,33, 2,5, 6, 1, 7, 0, 3, 4, 8,13,15,10,12, 9,11,14,43,41,42,47,45,44,46,40,49,53,48,51,54,52,50,55,30,26,29,24,25,28,27,31,62,61,59,63,60,58,56,57,22,19,20,17,23,18,21,16,
		38,37,33,32,34,39,35,36,49,53,48,55,50,52,54,51,30,26,29,31,27,28,25,24,22,19,20,16,21,18,17,23, 8,13,15,11,14,12, 9,10,62,61,59,57,58,60,63,56,43,41,42,45,47,40,46,44, 2,5, 6, 1, 4, 3, 0, 7,
		38,37,36,35,32,33,34,39, 5, 6, 2,1, 0, 3, 7, 4,53,48,49,55,54,52,51,50,19,20,22,21,18,16,17,23,61,59,62,58,60,57,63,56,41,42,43,47,40,45,44,46,13,15, 8,14,12,11,10, 9,26,29,30,31,25,28,24,27,
		38,37,39,34,35,36,32,33,29,30,26,31,24,28,27,25, 6, 2,5, 1, 7, 3, 4, 0,20,22,19,18,16,21,17,23,42,43,41,40,45,47,44,46,15, 8,13,12,11,14, 9,10,59,62,61,60,57,58,56,63,48,49,53,55,51,52,50,54,
		38,39,32,35,33,36,37,34,62,61,59,56,60,63,57,58,22,19,20,17,21,16,23,18,8,13,15, 9,12,14,10,11,43,41,42,40,44,45,47,46, 2,5, 6, 3, 4, 7, 0, 1,49,53,48,52,50,51,55,54,30,26,29,31,25,24,27,28,
		38,39,34,37,35,32,33,36,26,29,30,31,27,24,28,25,61,59,62,56,57,63,58,60,13,15, 8,12,14, 9,10,11,5, 6, 2,4, 7, 3, 0, 1,53,48,49,50,51,52,54,55,41,42,43,44,45,40,46,47,19,20,22,17,23,16,18,21,
		38,39,36,33,37,34,35,32,20,22,19,17,18,16,21,23,29,30,26,31,28,24,25,27,15, 8,13,14, 9,12,10,11,48,49,53,51,52,50,54,55,42,43,41,45,40,44,47,46, 6, 2,5, 7, 3, 4, 1, 0,59,62,61,56,58,63,60,57,
		39,32,34,33,35,38,36,37, 7, 3, 4, 0, 2,6, 1, 5,56,61,62,59,63,58,57,60,26,31,29,27,25,24,30,28,14, 9,12,10,13, 8,15,11,22,20,17,23,21,16,18,19,51,52,50,53,55,48,54,49,45,40,44,47,43,42,46,41,
		39,32,37,36,33,34,35,38,40,44,45,47,46,42,41,43, 3, 4, 7, 0, 1, 6, 5, 2,31,29,26,25,24,27,30,28,20,17,22,21,16,23,18,19,52,50,51,55,48,53,49,54, 9,12,14,13, 8,10,11,15,61,62,56,59,57,58,60,63,
		39,32,38,35,36,37,33,34,62,56,61,59,60,58,63,57,44,45,40,47,41,42,43,46,29,26,31,24,27,25,30,28,50,51,52,48,53,55,49,54,12,14, 9, 8,10,13,15,11,17,22,20,16,23,21,19,18,4, 7, 3, 0, 5, 6, 2,1,
		39,33,32,34,35,37,38,36, 3, 4, 7, 0, 6, 1, 2,5, 9,12,14,11,13, 8,10,15,40,44,45,46,42,41,47,43,52,50,51,48,55,49,53,54,31,29,26,27,28,25,24,30,61,62,56,60,63,57,59,58,20,17,22,19,21,16,23,18,
		39,33,36,38,34,32,35,37,17,22,20,19,23,16,18,21,4, 7, 3, 0, 2,1, 5, 6,44,45,40,42,41,46,47,43,29,26,31,28,25,27,24,30,62,56,61,63,57,60,58,59,50,51,52,55,49,48,54,53,12,14, 9,11,10, 8,15,13,
		39,33,37,35,38,36,34,32,14, 9,12,11,15, 8,13,10,22,20,17,19,18,16,21,23,45,40,44,41,46,42,47,43,56,61,62,57,60,63,58,59,51,52,50,49,48,55,53,54,26,31,29,25,27,28,30,24, 7, 3, 4, 0, 5, 1, 6, 2,
		39,34,33,32,35,36,37,38,4, 7, 3, 0, 1, 2,6, 5,50,51,52,54,55,49,48,53,17,22,20,23,16,18,19,21,62,56,61,57,63,58,60,59,44,45,40,46,43,42,41,47,12,14, 9,15,13,10,11,8,29,26,31,30,28,25,27,24,
		39,34,36,35,37,38,32,33,52,50,51,54,53,49,55,48,31,29,26,30,24,25,28,27,20,17,22,18,23,16,19,21,9,12,14,10,15,13, 8,11,61,62,56,58,57,63,60,59,40,44,45,42,46,43,47,41,3, 4, 7, 0, 5, 2,1, 6,
		39,34,38,37,32,33,35,36,26,31,29,30,27,25,24,28,7, 3, 4, 0, 6, 2,5, 1,22,20,17,16,18,23,19,21,45,40,44,43,42,46,41,47,14, 9,12,13,10,15, 8,11,56,61,62,63,58,57,59,60,51,52,50,54,48,49,53,55,
		39,35,32,38,36,34,37,33,56,61,62,59,58,63,60,57,51,52,50,54,53,55,48,49, 7, 3, 4, 2,6, 1, 0, 5,22,20,17,21,23,18,16,19,26,31,29,24,28,27,25,30,45,40,44,46,41,43,47,42,14, 9,12,11,10,13, 8,15,
		39,35,33,37,38,32,36,34, 9,12,14,11,8,13,15,10,61,62,56,59,60,63,57,58,3, 4, 7, 6, 1, 2,0, 5,31,29,26,28,27,24,25,30,40,44,45,41,43,46,42,47,20,17,22,23,18,21,19,16,52,50,51,54,48,55,49,53,
		39,35,34,36,37,33,38,32,50,51,52,54,49,55,53,48,12,14, 9,11,15,13,10, 8,4, 7, 3, 1, 2,6, 0, 5,44,45,40,43,46,41,42,47,17,22,20,18,21,23,16,19,29,26,31,27,24,28,30,25,62,56,61,59,57,63,58,60,
		39,36,32,37,33,38,34,35,44,45,40,47,42,41,46,43,17,22,20,19,16,23,21,18,62,56,61,60,58,63,59,57,12,14, 9,10, 8,15,13,11,29,26,31,25,28,24,27,30, 4, 7, 3, 2,1, 5, 0, 6,50,51,52,54,48,53,55,49,
		39,36,35,34,37,32,33,38,51,52,50,54,55,53,49,48,45,40,44,47,46,41,43,42,56,61,62,58,63,60,59,57,26,31,29,28,24,25,27,30, 7, 3, 4, 1, 5, 2,6, 0,14, 9,12, 8,15,10,11,13,22,20,17,19,21,23,18,16,
		39,36,38,33,34,35,37,32,20,17,22,19,18,23,16,21,52,50,51,54,49,53,48,55,61,62,56,63,60,58,59,57, 3, 4, 7, 5, 2,1, 6, 0, 9,12,14,15,10, 8,13,11,31,29,26,24,25,28,30,27,40,44,45,47,43,41,42,46,
		39,37,34,38,32,36,33,35,31,29,26,30,25,24,27,28,40,44,45,47,42,46,43,41,52,50,51,53,49,55,54,48,61,62,56,57,58,60,63,59,20,17,22,16,21,18,23,19, 3, 4, 7, 1, 6, 5, 0, 2,9,12,14,11,10,15,13, 8,
		39,37,35,33,38,34,32,36,12,14, 9,11,13,15, 8,10,29,26,31,30,27,24,28,25,50,51,52,49,55,53,54,48,17,22,20,21,18,16,23,19, 4, 7, 3, 6, 5, 1, 2,0,62,56,61,58,60,57,59,63,44,45,40,47,43,46,41,42,
		39,37,36,32,33,35,38,34,45,40,44,47,41,46,42,43,14, 9,12,11,8,15,10,13,51,52,50,55,53,49,54,48,7, 3, 4, 5, 1, 6, 2,0,56,61,62,60,57,58,63,59,22,20,17,18,16,21,19,23,26,31,29,30,28,24,25,27,
		39,38,33,36,34,37,32,35,22,20,17,19,16,18,23,21,26,31,29,30,25,27,28,24,14, 9,12,15, 8,13,11,10,51,52,50,48,49,53,55,54,45,40,44,42,43,41,46,47, 7, 3, 4, 6, 2,5, 0, 1,56,61,62,59,57,60,63,58,
		39,38,35,32,36,33,34,37,61,62,56,59,63,60,58,57,20,17,22,19,23,18,21,16, 9,12,14, 8,13,15,11,10,40,44,45,43,41,42,46,47, 3, 4, 7, 2,5, 6, 1, 0,52,50,51,49,53,48,54,55,31,29,26,30,28,27,24,25,
		39,38,37,34,32,35,36,33,29,26,31,30,24,27,25,28,62,56,61,59,58,60,57,63,12,14, 9,13,15, 8,11,10, 4, 7, 3, 5, 6, 2,1, 0,50,51,52,53,48,49,55,54,44,45,40,41,42,43,47,46,17,22,20,19,21,18,16,23,
		40,41,42,47,45,46,43,44, 0, 5, 3, 6, 2,7, 1, 4,56,60,61,57,63,59,58,62,48,55,52,53,49,54,50,51,32,37,39,38,33,35,34,36,24,30,31,26,25,28,29,27,16,18,20,19,17,23,22,21,8,12, 9,13,11,15,10,14,
		40,41,44,43,47,42,45,46,12, 9, 8,13,10,15,14,11,5, 3, 0, 6, 1, 7, 4, 2,55,52,48,49,54,53,50,51,30,31,24,25,28,26,29,27,18,20,16,17,23,19,21,22,37,39,32,33,35,38,36,34,60,61,56,57,58,59,62,63,
		40,41,46,45,43,44,47,42,61,56,60,57,62,59,63,58,9, 8,12,13,14,15,11,10,52,48,55,54,53,49,50,51,20,16,18,23,19,17,21,22,39,32,37,35,38,33,34,36,31,24,30,28,26,25,27,29, 3, 0, 5, 6, 4, 7, 2,1,
		40,42,43,45,44,46,41,47,18,20,16,22,19,21,17,23,55,52,48,50,54,49,51,53,30,31,24,29,26,28,27,25,37,39,32,38,34,33,35,36,60,61,56,59,58,63,62,57,12, 9, 8,15,10,11,13,14, 5, 3, 0, 6, 4, 2,1, 7,
		40,42,46,44,41,47,45,43,48,55,52,50,53,49,54,51,0, 5, 3, 6, 7, 2,4, 1,24,30,31,28,29,26,27,25, 8,12, 9,11,15,10,14,13,32,37,39,33,38,34,35,36,56,60,61,63,59,58,57,62,16,18,20,22,23,21,19,17,
		40,42,47,41,45,43,44,46, 3, 0, 5, 6, 1, 2,7, 4,20,16,18,22,17,21,23,19,31,24,30,26,28,29,27,25,61,56,60,58,63,59,62,57, 9, 8,12,10,11,15,14,13,39,32,37,34,33,38,36,35,52,48,55,50,51,49,53,54,
		40,43,41,44,47,46,42,45, 9, 8,12,13,15,14,10,11,31,24,30,27,28,26,25,29,61,56,60,62,59,63,57,58,39,32,37,38,35,34,33,36,52,48,55,49,51,54,53,50, 3, 0, 5, 2,1, 4, 6, 7,20,16,18,22,23,19,17,21,
		40,43,45,42,44,41,47,46,16,18,20,22,17,19,21,23, 8,12, 9,13,10,14,11,15,56,60,61,59,63,62,57,58,48,55,52,51,54,49,53,50, 0, 5, 3, 1, 4, 2,7, 6,32,37,39,35,34,38,36,33,24,30,31,27,25,26,29,28,
		40,43,46,47,42,45,44,41,30,31,24,27,29,26,28,25,18,20,16,22,21,19,23,17,60,61,56,63,62,59,57,58,5, 3, 0, 4, 2,1, 7, 6,37,39,32,34,38,35,33,36,55,52,48,54,49,51,50,53,12, 9, 8,13,11,14,15,10,
		40,44,42,46,41,43,47,45,55,52,48,50,49,54,53,51,12, 9, 8,13,15,10,11,14,18,20,16,19,21,17,22,23,60,61,56,58,59,62,63,57,30,31,24,28,25,29,26,27, 5, 3, 0, 1, 7, 4, 6, 2,37,39,32,36,38,34,33,35,
		40,44,43,41,47,45,46,42, 8,12, 9,13,14,10,15,11,32,37,39,36,35,34,38,33,16,18,20,17,19,21,22,23, 0, 5, 3, 4, 1, 7, 2,6,56,60,61,62,58,59,63,57,24,30,31,29,28,25,27,26,48,55,52,50,51,54,49,53,
		40,44,45,47,46,42,41,43,39,32,37,36,33,34,35,38,52,48,55,50,53,54,51,49,20,16,18,21,17,19,22,23,31,24,30,25,29,28,26,27, 3, 0, 5, 7, 4, 1, 2,6,61,56,60,59,62,58,57,63, 9, 8,12,13,11,10,14,15,
		40,45,41,46,43,42,44,47,56,60,61,57,59,63,62,58,16,18,20,22,19,17,23,21,0, 5, 3, 2,7, 1, 6, 4,24,30,31,25,26,29,28,27,48,55,52,54,51,53,49,50, 8,12, 9,10,14,11,13,15,32,37,39,36,38,33,35,34,
		40,45,42,43,44,47,46,41,20,16,18,22,21,17,19,23,39,32,37,36,34,33,38,35, 3, 0, 5, 1, 2,7, 6, 4, 9, 8,12,11,10,14,15,13,31,24,30,29,25,26,28,27,52,48,55,53,54,51,50,49,61,56,60,57,58,63,59,62,
		40,45,47,44,46,41,43,42,37,39,32,36,35,33,34,38,60,61,56,57,62,63,58,59, 5, 3, 0, 7, 1, 2,6, 4,55,52,48,51,53,54,49,50,12, 9, 8,14,11,10,15,13,30,31,24,26,29,25,27,28,18,20,16,22,23,17,21,19,
		40,46,44,42,41,45,43,47,52,48,55,50,54,53,49,51,61,56,60,57,59,62,58,63,39,32,37,33,34,35,36,38,3, 0, 5, 4, 7, 2,1, 6,20,16,18,19,23,21,17,22, 9, 8,12,14,15,11,13,10,31,24,30,27,25,29,28,26,
		40,46,45,41,43,47,42,44,60,61,56,57,63,62,59,58,30,31,24,27,26,29,25,28,37,39,32,35,33,34,36,38,12, 9, 8,11,14,15,10,13, 5, 3, 0, 2,4, 7, 1, 6,18,20,16,21,19,23,22,17,55,52,48,50,51,53,54,49,
		40,46,47,43,42,44,41,45,24,30,31,27,28,29,26,25,48,55,52,50,49,53,51,54,32,37,39,34,35,33,36,38,16,18,20,23,21,19,17,22, 8,12, 9,15,11,14,10,13, 0, 5, 3, 7, 2,4, 6, 1,56,60,61,57,58,62,63,59,
		40,47,41,42,45,44,46,43, 5, 3, 0, 6, 7, 1, 2,4,37,39,32,36,33,35,38,34,12, 9, 8,10,15,14,13,11,18,20,16,23,17,21,19,22,55,52,48,53,51,49,54,50,60,61,56,62,63,58,57,59,30,31,24,27,25,28,26,29,
		40,47,43,46,42,41,45,44,31,24,30,27,26,28,29,25, 3, 0, 5, 6, 2,1, 4, 7, 9, 8,12,15,14,10,13,11,52,48,55,51,49,53,54,50,61,56,60,63,58,62,59,57,20,16,18,17,21,23,22,19,39,32,37,36,38,35,34,33,
		40,47,44,45,46,43,42,41,32,37,39,36,34,35,33,38,24,30,31,27,29,28,25,26, 8,12, 9,14,10,15,13,11,56,60,61,58,62,63,59,57,16,18,20,21,23,17,19,22,48,55,52,49,53,51,50,54, 0, 5, 3, 6, 4, 1, 7, 2,
		41,40,43,44,42,47,46,45, 9,12,13, 8,15,10,11,14, 0, 6, 5, 3, 4, 2,1, 7,49,54,53,55,52,48,51,50,25,28,26,30,31,24,27,29,17,23,19,18,20,16,22,21,33,35,38,37,39,32,34,36,57,56,61,60,63,62,59,58,
		41,40,45,46,44,43,42,47,56,61,57,60,59,62,58,63,12,13, 9, 8,11,10,14,15,54,53,49,52,48,55,51,50,23,19,17,20,16,18,22,21,35,38,33,39,32,37,36,34,28,26,25,31,24,30,29,27, 6, 5, 0, 3, 1, 2,7, 4,
		41,40,47,42,46,45,44,43, 5, 0, 6, 3, 7, 2,4, 1,61,57,56,60,58,62,63,59,53,49,54,48,55,52,51,50,38,33,35,32,37,39,36,34,26,25,28,24,30,31,27,29,19,17,23,16,18,20,21,22,13, 9,12, 8,14,10,15,11,
		41,42,40,47,46,43,45,44, 0, 6, 5, 3, 2,4, 7, 1,33,35,38,34,37,39,32,36, 9,12,13,15,10,11,8,14,17,23,19,20,18,22,16,21,49,54,53,48,50,55,52,51,57,56,61,59,58,63,60,62,25,28,26,29,30,31,24,27,
		41,42,43,46,45,44,47,40,38,33,35,34,36,39,37,32,26,25,28,29,27,31,30,24,13, 9,12,11,15,10, 8,14,61,57,56,63,59,58,62,60,19,17,23,22,20,18,16,21,53,49,54,55,48,50,51,52, 5, 0, 6, 3, 1, 4, 2,7,
		41,42,44,45,47,40,46,43,28,26,25,29,24,31,27,30, 6, 5, 0, 3, 7, 4, 1, 2,12,13, 9,10,11,15, 8,14,54,53,49,50,55,48,52,51,56,61,57,58,63,59,62,60,23,19,17,18,22,20,21,16,35,38,33,34,32,39,36,37,
		41,43,44,40,42,46,45,47,13, 9,12, 8,11,15,10,14,38,33,35,34,39,36,32,37,19,17,23,18,16,22,21,20, 5, 0, 6, 1, 4, 2,7, 3,61,57,56,59,63,62,58,60,26,25,28,27,31,30,29,24,53,49,54,51,50,52,55,48,
		41,43,46,42,45,47,40,44,35,38,33,34,37,36,39,32,54,53,49,51,48,52,50,55,23,19,17,22,18,16,21,20,28,26,25,30,27,31,24,29, 6, 5, 0, 2,1, 4, 7, 3,56,61,57,62,59,63,60,58,12,13, 9, 8,14,15,11,10,
		41,43,47,45,40,44,42,46,49,54,53,51,55,52,48,50, 9,12,13, 8,10,15,14,11,17,23,19,16,22,18,21,20,57,56,61,63,62,59,58,60,25,28,26,31,30,27,24,29, 0, 6, 5, 4, 2,1, 3, 7,33,35,38,34,32,36,37,39,
		41,44,40,43,42,45,47,46,12,13, 9, 8,10,11,15,14,28,26,25,29,31,24,30,27,56,61,57,59,62,58,60,63,35,38,33,32,39,36,37,34,54,53,49,55,50,52,48,51,6, 5, 0, 7, 4, 1, 3, 2,23,19,17,21,20,16,18,22,
		41,44,45,42,47,46,43,40,25,28,26,29,27,24,31,30,17,23,19,21,22,16,20,18,57,56,61,58,59,62,60,63, 0, 6, 5, 1, 7, 4, 2,3,33,35,38,36,32,39,37,34,49,54,53,52,55,50,51,48,9,12,13, 8,14,11,10,15,
		41,44,46,47,43,40,42,45,19,17,23,21,18,16,22,20,13, 9,12, 8,15,11,14,10,61,57,56,62,58,59,60,63,53,49,54,50,52,55,48,51,5, 0, 6, 4, 1, 7, 2,3,38,33,35,39,36,32,34,37,26,25,28,29,30,24,27,31,
		41,45,42,44,47,43,40,46,26,25,28,29,31,27,24,30,53,49,54,51,55,48,50,52,38,33,35,36,39,37,34,32,19,17,23,20,22,16,18,21,13, 9,12,10,14,11,15, 8,5, 0, 6, 2,7, 1, 3, 4,61,57,56,60,63,59,58,62,
		41,45,43,47,40,46,44,42,54,53,49,51,52,48,55,50,56,61,57,60,62,59,63,58,35,38,33,37,36,39,34,32, 6, 5, 0, 1, 2,7, 4, 3,23,19,17,16,20,22,18,21,12,13, 9,11,10,14, 8,15,28,26,25,29,30,27,31,24,
		41,45,46,40,44,42,47,43,57,56,61,60,58,59,62,63,25,28,26,29,24,27,30,31,33,35,38,39,37,36,34,32, 9,12,13,14,11,10,15, 8,0, 6, 5, 7, 1, 2,4, 3,17,23,19,22,16,20,21,18,49,54,53,51,50,48,52,55,
		41,46,40,45,44,47,43,42,61,57,56,60,62,58,59,63,19,17,23,21,16,18,20,22, 5, 0, 6, 7, 2,4, 3, 1,26,25,28,30,24,27,31,29,53,49,54,52,50,48,55,51,13, 9,12,15,11,14, 8,10,38,33,35,34,32,37,39,36,
		41,46,42,43,45,40,44,47,33,35,38,34,39,37,36,32,57,56,61,60,59,58,63,62, 0, 6, 5, 2,4, 7, 3, 1,49,54,53,50,48,52,55,51,9,12,13,11,14,15,10, 8,25,28,26,24,27,30,29,31,17,23,19,21,20,18,22,16,
		41,46,47,44,43,42,45,40,23,19,17,21,22,18,16,20,35,38,33,34,36,37,32,39, 6, 5, 0, 4, 7, 2,3, 1,12,13, 9,14,15,11,10, 8,28,26,25,27,30,24,31,29,54,53,49,48,52,50,51,55,56,61,57,60,63,58,62,59,
		41,47,42,40,46,44,43,45, 6, 5, 0, 3, 4, 7, 2,1,23,19,17,21,18,22,20,16,28,26,25,24,31,27,29,30,56,61,57,63,58,62,59,60,12,13, 9,15,14,10,11,8,35,38,33,36,37,32,34,39,54,53,49,51,50,55,48,52,
		41,47,44,46,43,45,40,42,17,23,19,21,16,22,18,20,49,54,53,51,52,55,50,48,25,28,26,27,24,31,29,30,33,35,38,32,36,37,39,34,57,56,61,62,63,58,59,60, 9,12,13,10,15,14, 8,11,0, 6, 5, 3, 1, 7, 4, 2,
		41,47,45,43,40,42,46,44,53,49,54,51,48,55,52,50, 5, 0, 6, 3, 2,7, 1, 4,26,25,28,31,27,24,29,30,13, 9,12,14,10,15,11,8,38,33,35,37,32,36,39,34,61,57,56,58,62,63,60,59,19,17,23,21,20,22,16,18,
		42,40,41,47,43,45,46,44, 0, 3, 6, 5, 2,1, 4, 7,18,22,20,16,23,19,17,21,26,28,29,31,24,30,25,27,58,63,59,61,56,60,57,62,10,11,15, 9, 8,12,13,14,34,33,38,39,32,37,35,36,50,55,48,52,54,53,49,51,
		42,40,44,46,47,41,43,45,55,48,50,52,49,53,51,54, 3, 6, 0, 5, 4, 1, 7, 2,28,29,26,24,30,31,25,27,11,15,10, 8,12, 9,13,14,33,38,34,32,37,39,36,35,63,59,58,56,60,61,62,57,22,20,18,16,17,19,21,23,
		42,40,45,43,46,44,47,41,20,18,22,16,21,19,23,17,48,50,55,52,51,53,54,49,29,26,28,30,31,24,25,27,38,34,33,37,39,32,36,35,59,58,63,60,61,56,57,62,15,10,11,12, 9, 8,14,13, 6, 0, 3, 5, 7, 1, 2,4,
		42,41,45,44,40,47,43,46,26,28,29,25,31,24,30,27, 0, 3, 6, 5, 1, 2,7, 4,10,11,15,12,13, 9,14, 8,50,55,48,54,53,49,51,52,58,63,59,56,61,57,60,62,18,22,20,23,19,17,16,21,34,33,38,35,37,36,39,32,
		42,41,46,43,44,45,40,47,33,38,34,35,39,36,32,37,28,29,26,25,30,24,27,31,11,15,10,13, 9,12,14, 8,63,59,58,61,57,56,60,62,22,20,18,19,17,23,21,16,55,48,50,53,49,54,52,51,3, 6, 0, 5, 7, 2,4, 1,
		42,41,47,40,43,46,44,45, 6, 0, 3, 5, 4, 2,1, 7,38,34,33,35,32,36,37,39,15,10,11,9,12,13,14, 8,20,18,22,17,23,19,21,16,48,50,55,49,54,53,51,52,59,58,63,57,56,61,62,60,29,26,28,25,27,24,31,30,
		42,43,40,45,46,41,44,47,18,22,20,16,19,23,21,17,34,33,38,35,39,32,37,36, 0, 3, 6, 2,1, 4, 5, 7,10,11,15, 8,9,13,12,14,26,28,29,30,27,31,24,25,50,55,48,49,51,54,52,53,58,63,59,62,61,56,60,57,
		42,43,41,46,44,47,45,40,38,34,33,35,36,32,39,37,59,58,63,62,57,56,61,60, 6, 0, 3, 4, 2,1, 5, 7,48,50,55,54,49,51,53,52,15,10,11,13, 8,9,12,14,29,26,28,31,30,27,25,24,20,18,22,16,17,23,19,21,
		42,43,47,44,45,40,46,41,63,59,58,62,60,56,57,61,22,20,18,16,21,23,17,19, 3, 6, 0, 1, 4, 2,5, 7,28,29,26,27,31,30,24,25,55,48,50,51,54,49,53,52,11,15,10, 9,13, 8,14,12,33,38,34,35,37,32,36,39,
		42,44,41,45,40,46,47,43,28,29,26,25,24,30,31,27,55,48,50,52,53,49,54,51,33,38,34,39,36,32,35,37,22,20,18,17,19,21,23,16,11,15,10,12, 8,13, 9,14, 3, 6, 0, 4, 1, 7, 5, 2,63,59,58,62,61,57,56,60,
		42,44,43,47,45,41,40,46,59,58,63,62,56,57,60,61,29,26,28,25,31,30,27,24,38,34,33,36,32,39,35,37,15,10,11,8,13,12, 9,14, 6, 0, 3, 1, 7, 4, 2,5,20,18,22,19,21,17,16,23,48,50,55,52,54,49,51,53,
		42,44,46,40,47,43,45,41,50,55,48,52,51,49,53,54,58,63,59,62,60,57,61,56,34,33,38,32,39,36,35,37, 0, 3, 6, 7, 4, 1, 2,5,18,22,20,21,17,19,23,16,10,11,15,13,12, 8,14, 9,26,28,29,25,27,30,24,31,
		42,45,43,40,46,47,41,44,22,20,18,16,23,21,19,17,11,15,10,14, 9,13, 8,12,63,59,58,60,56,57,62,61,55,48,50,54,51,53,49,52, 3, 6, 0, 2,7, 1, 4, 5,33,38,34,36,39,37,35,32,28,29,26,25,27,31,30,24,
		42,45,44,41,40,43,46,47,29,26,28,25,30,31,24,27,20,18,22,16,19,21,17,23,59,58,63,56,57,60,62,61,6, 0, 3, 7, 1, 2,4, 5,38,34,33,39,37,36,32,35,48,50,55,51,53,54,52,49,15,10,11,14, 8,13,12, 9,
		42,45,47,46,41,44,40,43,10,11,15,14,12,13, 9, 8,26,28,29,25,24,31,27,30,58,63,59,57,60,56,62,61,34,33,38,37,36,39,32,35,50,55,48,53,54,51,49,52, 0, 3, 6, 1, 2,7, 5, 4,18,22,20,16,17,21,23,19,
		42,46,40,44,47,45,41,43,48,50,55,52,53,51,49,54,15,10,11,14,12, 9, 8,13,20,18,22,21,19,23,16,17,59,58,63,61,60,57,56,62,29,26,28,24,27,30,31,25, 6, 0, 3, 2,4, 7, 5, 1,38,34,33,35,37,39,32,36,
		42,46,43,41,44,40,47,45,34,33,38,35,32,39,36,37,50,55,48,52,49,51,54,53,18,22,20,19,23,21,16,17,26,28,29,27,30,24,31,25, 0, 3, 6, 4, 7, 2,1, 5,58,63,59,60,57,61,62,56,10,11,15,14, 8,9,13,12,
		42,46,45,47,41,43,44,40,11,15,10,14,13, 9,12, 8,33,38,34,35,36,39,37,32,22,20,18,23,21,19,16,17, 3, 6, 0, 7, 2,4, 1, 5,63,59,58,57,61,60,56,62,28,29,26,30,24,27,25,31,55,48,50,52,54,51,53,49,
		42,47,40,41,43,44,45,46, 3, 6, 0, 5, 1, 4, 2,7,63,59,58,62,56,60,61,57,55,48,50,49,53,51,52,54,33,38,34,37,32,36,39,35,28,29,26,31,27,24,30,25,22,20,18,21,23,17,16,19,11,15,10,14, 8,12, 9,13,
		42,47,44,43,45,46,41,40,58,63,59,62,57,60,56,61,10,11,15,14,13,12, 8,9,50,55,48,51,49,53,52,54,18,22,20,17,21,23,19,16,34,33,38,36,37,32,39,35,26,28,29,24,31,27,25,30, 0, 3, 6, 5, 7, 4, 1, 2,
		42,47,46,45,41,40,43,44,15,10,11,14, 9,12,13, 8,6, 0, 3, 5, 2,4, 7, 1,48,50,55,53,51,49,52,54,29,26,28,27,24,31,30,25,20,18,22,23,17,21,19,16,38,34,33,32,36,37,35,39,59,58,63,62,61,60,57,56,
		43,40,42,45,41,44,46,47,18,16,22,20,19,17,23,21,9,13, 8,12,11,15,10,14,59,63,62,56,60,61,58,57,51,54,49,48,55,52,50,53, 1, 4, 2,0, 5, 3, 6, 7,35,34,38,32,37,39,33,36,27,31,30,24,28,29,26,25,
		43,40,44,41,46,47,45,42, 8,9,13,12,14,15,11,10,30,27,31,24,25,29,28,26,62,59,63,61,56,60,58,57,38,35,34,39,32,37,36,33,49,51,54,52,48,55,50,53, 2,1, 4, 3, 0, 5, 7, 6,22,18,16,20,21,17,19,23,
		43,40,47,46,45,42,41,44,31,30,27,24,26,29,25,28,16,22,18,20,23,17,21,19,63,62,59,60,61,56,58,57, 4, 2,1, 5, 3, 0, 6, 7,34,38,35,37,39,32,36,33,54,49,51,55,52,48,53,50,13, 8,9,12,10,15,14,11,
		43,41,40,44,46,42,47,45, 9,13, 8,12,15,11,14,10,35,34,38,33,32,37,39,36,18,16,22,19,17,23,20,21,1, 4, 2,5, 0, 6, 3, 7,59,63,62,61,57,56,60,58,27,31,30,26,25,28,24,29,51,54,49,53,48,55,52,50,
		43,41,42,46,47,45,44,40,38,35,34,33,36,37,32,39,49,51,54,53,50,55,48,52,22,18,16,23,19,17,20,21,30,27,31,28,26,25,29,24, 2,1, 4, 6, 5, 0, 3, 7,62,59,63,56,61,57,58,60, 8,9,13,12,10,11,15,14,
		43,41,45,47,44,40,46,42,54,49,51,53,52,55,50,48,13, 8,9,12,14,11,10,15,16,22,18,17,23,19,20,21,63,62,59,57,56,61,60,58,31,30,27,25,28,26,29,24, 4, 2,1, 0, 6, 5, 7, 3,34,38,35,33,39,37,36,32,
		43,42,44,47,40,45,41,46,59,63,62,58,56,60,61,57,18,16,22,20,17,19,21,23, 1, 4, 2,3, 6, 0, 7, 5,27,31,30,28,29,26,25,24,51,54,49,55,48,50,52,53, 9,13, 8,11,15,10,12,14,35,34,38,33,39,36,32,37,
		43,42,45,40,41,46,47,44,22,18,16,20,23,19,17,21,38,35,34,33,37,36,39,32, 2,1, 4, 0, 3, 6, 7, 5, 8,9,13,10,11,15,14,12,30,27,31,26,28,29,25,24,49,51,54,50,55,48,53,52,62,59,63,58,57,60,56,61,
		43,42,46,41,47,44,40,45,34,38,35,33,32,36,37,39,63,62,59,58,61,60,57,56, 4, 2,1, 6, 0, 3, 7, 5,54,49,51,48,50,55,52,53,13, 8,9,15,10,11,14,12,31,30,27,29,26,28,24,25,16,22,18,20,21,19,23,17,
		43,44,41,40,46,45,42,47,13, 8,9,12,11,14,15,10, 4, 2,1, 7, 0, 6, 5, 3,54,49,51,52,55,50,53,48,31,30,27,28,25,29,26,24,16,22,18,19,21,17,23,20,34,38,35,36,32,39,33,37,63,62,59,58,57,56,61,60,
		43,44,45,46,42,47,40,41,1, 4, 2,7, 3, 6, 0, 5,59,63,62,58,60,56,57,61,51,54,49,50,52,55,53,48,35,34,38,39,36,32,37,33,27,31,30,29,28,25,26,24,18,16,22,17,19,21,20,23, 9,13, 8,12,10,14,11,15,
		43,44,47,42,40,41,46,45,62,59,63,58,61,56,60,57, 8,9,13,12,15,14,10,11,49,51,54,55,50,52,53,48,22,18,16,21,17,19,23,20,38,35,34,32,39,36,37,33,30,27,31,25,29,28,24,26, 2,1, 4, 7, 5, 6, 3, 0,
		43,45,40,42,41,47,44,46,16,22,18,20,17,23,19,21,54,49,51,53,55,52,48,50,31,30,27,26,29,25,24,28,34,38,35,39,37,36,32,33,63,62,59,56,57,60,61,58,13, 8,9,14,11,10,12,15, 4, 2,1, 7, 5, 3, 0, 6,
		43,45,46,44,42,40,41,47, 2,1, 4, 7, 0, 3, 6, 5,22,18,16,20,19,23,21,17,30,27,31,29,25,26,24,28,62,59,63,57,60,56,61,58,8,9,13,11,10,14,15,12,38,35,34,37,36,39,33,32,49,51,54,53,48,52,50,55,
		43,45,47,41,44,46,42,40,51,54,49,53,50,52,55,48,1, 4, 2,7, 6, 3, 5, 0,27,31,30,25,26,29,24,28,9,13, 8,10,14,11,15,12,35,34,38,36,39,37,32,33,59,63,62,60,56,57,58,61,18,16,22,20,21,23,17,19,
		43,46,40,47,45,44,42,41,30,27,31,24,29,25,26,28,2,1, 4, 7, 3, 0, 5, 6, 8,9,13,14,15,11,12,10,49,51,54,48,52,50,55,53,62,59,63,60,57,61,56,58,22,18,16,19,23,21,20,17,38,35,34,33,39,32,37,36,
		43,46,41,42,47,40,45,44,35,34,38,33,37,32,36,39,27,31,30,24,26,25,28,29, 9,13, 8,15,11,14,12,10,59,63,62,57,61,60,56,58,18,16,22,23,21,19,17,20,51,54,49,52,50,48,53,55, 1, 4, 2,7, 5, 0, 6, 3,
		43,46,44,45,42,41,47,40, 4, 2,1, 7, 6, 0, 3, 5,34,38,35,33,36,32,39,37,13, 8,9,11,14,15,12,10,16,22,18,21,19,23,17,20,54,49,51,50,48,52,55,53,63,62,59,61,60,57,58,56,31,30,27,24,28,25,29,26,
		43,47,41,45,44,42,40,46,49,51,54,53,55,50,52,48,62,59,63,58,56,61,57,60,38,35,34,36,37,32,33,39, 2,1, 4, 5, 6, 3, 0, 7,22,18,16,17,21,23,19,20, 8,9,13,15,14,10,12,11,30,27,31,24,28,26,25,29,
		43,47,42,44,40,46,45,41,63,62,59,58,60,61,56,57,31,30,27,24,29,26,28,25,34,38,35,32,36,37,33,39,13, 8,9,10,15,14,11,12, 4, 2,1, 3, 5, 6, 0, 7,16,22,18,23,17,21,20,19,54,49,51,53,48,50,55,52,
		43,47,46,40,45,41,44,42,27,31,30,24,25,26,29,28,51,54,49,53,52,50,48,55,35,34,38,37,32,36,33,39,18,16,22,21,23,17,19,20, 9,13, 8,14,10,15,11,12, 1, 4, 2,6, 3, 5, 7, 0,59,63,62,58,57,61,60,56,
		44,40,41,43,45,47,42,46,12, 8,13, 9,10,14,11,15,39,36,32,37,38,33,35,34,17,19,21,16,18,20,23,22, 4, 1, 7, 0, 5, 3, 6, 2,62,58,59,56,60,61,57,63,29,28,25,24,30,31,26,27,50,52,55,48,53,49,54,51,
		44,40,46,42,43,41,45,47,52,55,50,48,54,49,51,53, 8,13,12, 9,11,14,15,10,19,21,17,18,20,16,23,22,58,59,62,60,61,56,57,63,28,25,29,30,31,24,27,26, 1, 7, 4, 5, 3, 0, 2,6,36,32,39,37,35,33,34,38,
		44,40,47,45,42,46,43,41,32,39,36,37,34,33,38,35,55,50,52,48,51,49,53,54,21,17,19,20,16,18,23,22,25,29,28,31,24,30,27,26, 7, 4, 1, 3, 0, 5, 6, 2,59,62,58,61,56,60,63,57,13,12, 8,9,15,14,10,11,
		44,41,42,45,46,47,40,43,28,25,29,26,24,27,30,31,19,21,17,23,20,18,22,16,58,59,62,57,56,61,63,60, 1, 7, 4, 0, 6, 5, 3, 2,36,32,39,33,35,38,34,37,52,55,50,49,54,53,48,51,8,13,12, 9,15,10,11,14,
		44,41,43,40,45,42,46,47,13,12, 8,9,11,10,14,15,25,29,28,26,30,27,31,24,59,62,58,56,61,57,63,60,32,39,36,35,38,33,34,37,55,50,52,54,53,49,51,48,7, 4, 1, 6, 5, 0, 2,3,21,17,19,23,22,18,16,20,
		44,41,47,46,40,43,45,42,17,19,21,23,16,18,20,22,12, 8,13, 9,14,10,15,11,62,58,59,61,57,56,63,60,50,52,55,53,49,54,51,48,4, 1, 7, 5, 0, 6, 3, 2,39,36,32,38,33,35,37,34,29,28,25,26,31,27,24,30,
		44,42,40,46,43,47,41,45,55,50,52,48,49,51,54,53,59,62,58,63,61,56,60,57,32,39,36,34,33,38,37,35, 7, 4, 1, 0, 3, 6, 5, 2,21,17,19,18,22,20,16,23,13,12, 8,10,11,15, 9,14,25,29,28,26,31,24,30,27,
		44,42,45,41,46,40,43,47,29,28,25,26,30,24,27,31,50,52,55,48,54,51,53,49,39,36,32,33,38,34,37,35,17,19,21,22,20,18,16,23,12, 8,13,11,15,10,14, 9, 4, 1, 7, 3, 6, 0, 2,5,62,58,59,63,60,56,57,61,
		44,42,47,43,41,45,46,40,58,59,62,63,57,56,61,60,28,25,29,26,27,24,31,30,36,32,39,38,34,33,37,35, 8,13,12,15,10,11,14, 9, 1, 7, 4, 6, 0, 3, 5, 2,19,21,17,20,18,22,23,16,52,55,50,48,53,51,49,54,
		44,43,40,41,45,46,47,42, 8,13,12, 9,14,11,10,15, 1, 7, 4, 2,5, 3, 0, 6,52,55,50,54,49,51,48,53,28,25,29,31,30,27,24,26,19,21,17,16,22,18,20,23,36,32,39,34,38,35,37,33,58,59,62,63,60,61,56,57,
		44,43,42,47,41,40,45,46,59,62,58,63,56,61,57,60,13,12, 8,9,10,11,15,14,55,50,52,49,51,54,48,53,21,17,19,22,18,16,20,23,32,39,36,38,35,34,33,37,25,29,28,30,27,31,26,24, 7, 4, 1, 2,0, 3, 6, 5,
		44,43,46,45,47,42,41,40, 4, 1, 7, 2,6, 3, 5, 0,62,58,59,63,57,61,60,56,50,52,55,51,54,49,48,53,39,36,32,35,34,38,33,37,29,28,25,27,31,30,24,26,17,19,21,18,16,22,23,20,12, 8,13, 9,15,11,14,10,
		44,45,40,47,42,41,46,43,39,36,32,37,33,38,34,35,29,28,25,26,24,30,31,27,12, 8,13,10,14,11,9,15,62,58,59,60,56,57,61,63,17,19,21,20,22,16,18,23,50,52,55,54,51,53,48,49, 4, 1, 7, 2,0, 5, 3, 6,
		44,45,41,42,46,43,47,40,25,29,28,26,27,30,24,31,7, 4, 1, 2,6, 5, 0, 3,13,12, 8,11,10,14, 9,15,55,50,52,53,54,51,49,48,59,62,58,57,60,56,61,63,21,17,19,16,20,22,23,18,32,39,36,37,35,38,33,34,
		44,45,43,46,47,40,42,41,1, 7, 4, 2,3, 5, 6, 0,36,32,39,37,34,38,35,33, 8,13,12,14,11,10, 9,15,19,21,17,22,16,20,18,23,52,55,50,51,53,54,49,48,58,59,62,56,57,60,63,61,28,25,29,26,31,30,27,24,
		44,46,41,47,40,42,43,45,19,21,17,23,18,20,16,22,52,55,50,48,49,54,53,51,28,25,29,24,27,30,26,31,36,32,39,35,33,34,38,37,58,59,62,61,60,57,56,63, 8,13,12,11,14,15, 9,10, 1, 7, 4, 2,0, 6, 5, 3,
		44,46,42,40,43,45,47,41,50,52,55,48,51,54,49,53, 4, 1, 7, 2,3, 6, 0, 5,29,28,25,30,24,27,26,31,12, 8,13,15,11,14,10, 9,39,36,32,34,35,33,38,37,62,58,59,57,61,60,63,56,17,19,21,23,22,20,18,16,
		44,46,45,43,47,41,40,42, 7, 4, 1, 2,5, 6, 3, 0,21,17,19,23,16,20,22,18,25,29,28,27,30,24,26,31,59,62,58,60,57,61,56,63,13,12, 8,14,15,11,10, 9,32,39,36,33,34,35,37,38,55,50,52,48,53,54,51,49,
		44,47,43,42,41,46,40,45,62,58,59,63,61,57,56,60,17,19,21,23,18,16,22,20, 4, 1, 7, 6, 3, 5, 2,0,29,28,25,31,27,24,30,26,50,52,55,49,53,51,54,48,12, 8,13,14,10,15, 9,11,39,36,32,37,35,34,38,33,
		44,47,45,40,42,43,41,46,36,32,39,37,38,34,33,35,58,59,62,63,56,57,60,61,1, 7, 4, 3, 5, 6, 2,0,52,55,50,53,51,49,54,48,8,13,12,10,15,14,11,9,28,25,29,27,24,31,26,30,19,21,17,23,22,16,20,18,
		44,47,46,41,40,45,42,43,21,17,19,23,20,16,18,22,32,39,36,37,33,34,35,38,7, 4, 1, 5, 6, 3, 2,0,13,12, 8,15,14,10,11,9,25,29,28,24,31,27,30,26,55,50,52,51,49,53,48,54,59,62,58,63,60,57,61,56,
		45,40,43,42,47,44,41,46,16,20,22,18,17,21,23,19,37,36,39,32,38,35,34,33, 1, 2,7, 3, 0, 5, 4, 6,11,10,14, 9, 8,12,13,15,29,25,26,31,24,30,27,28,53,54,51,52,48,55,49,50,57,60,56,61,62,59,63,58,
		45,40,44,47,41,46,42,43,39,37,36,32,33,35,38,34,56,57,60,61,58,59,62,63, 7, 1, 2,5, 3, 0, 4, 6,51,53,54,55,52,48,50,49,14,11,10,12, 9, 8,13,15,26,29,25,30,31,24,28,27,22,16,20,18,19,21,17,23,
		45,40,46,41,42,43,47,44,60,56,57,61,63,59,58,62,20,22,16,18,23,21,19,17, 2,7, 1, 0, 5, 3, 4, 6,25,26,29,24,30,31,27,28,54,51,53,48,55,52,50,49,10,14,11,8,12, 9,15,13,36,39,37,32,34,35,33,38,
		45,41,40,46,42,44,43,47,56,57,60,61,59,58,63,62,26,29,25,28,30,31,24,27,39,37,36,33,35,38,32,34,14,11,10, 9,12,13, 8,15, 7, 1, 2,0, 6, 5, 3, 4,22,16,20,17,23,19,18,21,51,53,54,49,55,52,48,50,
		45,41,44,42,43,47,46,40,25,26,29,28,27,31,30,24,54,51,53,49,50,52,55,48,36,39,37,38,33,35,32,34,20,22,16,19,17,23,21,18,10,14,11,13, 9,12, 8,15, 2,7, 1, 5, 0, 6, 4, 3,60,56,57,61,62,58,59,63,
		45,41,47,43,46,40,42,44,53,54,51,49,48,52,50,55,57,60,56,61,63,58,62,59,37,36,39,35,38,33,32,34, 1, 2,7, 6, 5, 0, 3, 4,16,20,22,23,19,17,21,18,11,10,14,12,13, 9,15, 8,29,25,26,28,24,31,27,30,
		45,42,40,43,47,46,44,41,20,22,16,18,21,23,17,19,10,14,11,15, 8,12, 9,13,60,56,57,63,59,58,61,62,54,51,53,55,48,50,52,49, 2,7, 1, 3, 6, 0, 5, 4,36,39,37,33,38,34,32,35,25,26,29,28,24,30,31,27,
		45,42,41,44,43,40,47,46,26,29,25,28,31,30,27,24,22,16,20,18,17,23,19,21,56,57,60,59,58,63,61,62, 7, 1, 2,6, 0, 3, 5, 4,39,37,36,38,34,33,35,32,51,53,54,48,50,55,49,52,14,11,10,15, 9,12,13, 8,
		45,42,46,47,44,41,43,40,11,10,14,15,13,12, 8,9,29,25,26,28,27,30,24,31,57,60,56,58,63,59,61,62,37,36,39,34,33,38,35,32,53,54,51,50,55,48,52,49, 1, 2,7, 0, 3, 6, 4, 5,16,20,22,18,19,23,21,17,
		45,43,41,47,46,44,40,42,54,51,53,49,52,50,48,55, 2,7, 1, 4, 5, 0, 6, 3,25,26,29,27,31,30,28,24,10,14,11,9,13, 8,12,15,36,39,37,35,34,38,33,32,60,56,57,59,63,62,61,58,20,22,16,18,19,17,23,21,
		45,43,42,40,47,41,46,44,22,16,20,18,23,17,21,19,51,53,54,49,48,50,55,52,26,29,25,31,30,27,28,24,39,37,36,34,38,35,33,32,56,57,60,63,62,59,58,61,14,11,10,13, 8,9,15,12, 7, 1, 2,4, 6, 0, 3, 5,
		45,43,44,46,40,42,47,41,1, 2,7, 4, 3, 0, 5, 6,16,20,22,18,21,17,19,23,29,25,26,30,27,31,28,24,57,60,56,62,59,63,58,61,11,10,14, 8,9,13,12,15,37,36,39,38,35,34,32,33,53,54,51,49,55,50,52,48,
		45,44,42,41,43,46,40,47,29,25,26,28,30,27,31,24, 1, 2,7, 4, 0, 3, 6, 5,11,10,14,13,12, 8,15, 9,53,54,51,55,50,52,48,49,57,60,56,59,62,58,63,61,16,20,22,21,17,19,18,23,37,36,39,32,34,33,38,35,
		45,44,46,43,40,47,41,42, 7, 1, 2,4, 5, 3, 0, 6,39,37,36,32,35,33,34,38,14,11,10, 8,13,12,15, 9,22,16,20,19,21,17,23,18,51,53,54,52,55,50,48,49,56,57,60,58,59,62,61,63,26,29,25,28,24,27,30,31,
		45,44,47,40,41,42,43,46,36,39,37,32,38,33,35,34,25,26,29,28,31,27,24,30,10,14,11,12, 8,13,15, 9,60,56,57,62,58,59,63,61,20,22,16,17,19,21,23,18,54,51,53,50,52,55,49,48,2,7, 1, 4, 6, 3, 5, 0,
		45,46,41,40,42,47,44,43,57,60,56,61,58,63,59,62,11,10,14,15,12,13, 9, 8,53,54,51,48,52,50,49,55,16,20,22,19,23,21,17,18,37,36,39,33,34,35,38,32,29,25,26,27,30,24,28,31,1, 2,7, 4, 6, 5, 0, 3,
		45,46,43,44,40,41,42,47, 2,7, 1, 4, 0, 5, 3, 6,60,56,57,61,59,63,62,58,54,51,53,52,50,48,49,55,36,39,37,34,35,33,38,32,25,26,29,30,24,27,31,28,20,22,16,23,21,19,18,17,10,14,11,15, 9,13, 8,12,
		45,46,47,42,44,43,40,41,14,11,10,15, 8,13,12, 9, 7, 1, 2,4, 3, 5, 6, 0,51,53,54,50,48,52,49,55,26,29,25,24,27,30,31,28,22,16,20,21,19,23,17,18,39,37,36,35,33,34,32,38,56,57,60,61,62,63,58,59,
		45,47,40,44,41,43,46,42,37,36,39,32,35,38,33,34,53,54,51,49,52,48,55,50,16,20,22,17,21,23,18,19,29,25,26,24,31,27,30,28,1, 2,7, 5, 6, 3, 0, 4,57,60,56,63,58,62,61,59,11,10,14,15, 9, 8,12,13,
		45,47,42,46,44,40,41,43,10,14,11,15,12, 8,13, 9,36,39,37,32,33,38,34,35,20,22,16,21,23,17,18,19, 2,7, 1, 6, 3, 5, 0, 4,60,56,57,58,62,63,59,61,25,26,29,31,27,24,28,30,54,51,53,49,55,48,50,52,
		45,47,43,41,46,42,44,40,51,53,54,49,50,48,52,55,14,11,10,15,13, 8,9,12,22,16,20,23,17,21,18,19,56,57,60,62,63,58,59,61,26,29,25,27,24,31,30,28,7, 1, 2,3, 5, 6, 4, 0,39,37,36,32,34,38,35,33,
		46,40,41,45,47,43,44,42,61,60,57,56,62,63,58,59,24,27,30,31,25,28,26,29,35,33,34,37,39,32,38,36,11,14,15,12, 9, 8,13,10, 2,4, 7, 5, 3, 0, 6, 1,21,19,23,18,20,16,17,22,50,48,52,55,49,54,53,51,
		46,40,42,44,45,41,47,43,48,52,50,55,53,54,51,49,60,57,61,56,58,63,59,62,33,34,35,39,32,37,38,36, 4, 7, 2,3, 0, 5, 6, 1,19,23,21,20,16,18,22,17,14,15,11,9, 8,12,10,13,27,30,24,31,26,28,29,25,
		46,40,43,47,44,42,45,41,30,24,27,31,29,28,25,26,52,50,48,55,51,54,49,53,34,35,33,32,37,39,38,36,23,21,19,16,18,20,22,17,15,11,14, 8,12, 9,13,10, 7, 2,4, 0, 5, 3, 1, 6,57,61,60,56,59,63,62,58,
		46,41,43,42,40,45,47,44,35,33,34,38,37,39,32,36,61,60,57,56,63,62,59,58,2,4, 7, 0, 6, 5, 1, 3,50,48,52,49,54,53,51,55,11,14,15, 9,12,13, 8,10,24,27,30,25,28,26,31,29,21,19,23,17,16,22,18,20,
		46,41,44,47,42,43,40,45,19,23,21,17,18,22,20,16,33,34,35,38,32,39,36,37, 4, 7, 2,6, 5, 0, 1, 3,14,15,11,12,13, 9, 8,10,27,30,24,28,26,25,29,31,48,52,50,54,53,49,55,51,60,57,61,56,59,62,58,63,
		46,41,45,40,47,44,42,43,57,61,60,56,58,62,63,59,23,21,19,17,20,22,16,18,7, 2,4, 5, 0, 6, 1, 3,30,24,27,26,25,28,29,31,52,50,48,53,49,54,51,55,15,11,14,13, 9,12,10, 8,34,35,33,38,36,39,37,32,
		46,42,41,43,40,44,45,47,33,34,35,38,39,32,37,36,48,52,50,55,54,53,49,51,19,23,21,18,22,20,17,16,27,30,24,26,28,29,25,31,4, 7, 2,0, 3, 6, 5, 1,60,57,61,58,63,59,56,62,14,15,11,10,12,13, 9, 8,
		46,42,44,40,45,47,43,41,50,48,52,55,51,53,54,49,11,14,15,10, 8,13,12, 9,21,19,23,20,18,22,17,16,61,60,57,59,58,63,62,56,24,27,30,29,26,28,25,31,2,4, 7, 6, 0, 3, 1, 5,35,33,34,38,36,32,39,37,
		46,42,47,45,43,41,40,44,15,11,14,10, 9,13, 8,12,34,35,33,38,37,32,36,39,23,21,19,22,20,18,17,16, 7, 2,4, 3, 6, 0, 5, 1,57,61,60,63,59,58,62,56,30,24,27,28,29,26,31,25,52,50,48,55,49,53,51,54,
		46,43,42,41,40,47,44,45,34,35,33,38,32,37,39,36,30,24,27,31,28,29,26,25,15,11,14, 9,13, 8,10,12,57,61,60,59,63,62,58,56,23,21,19,18,16,22,20,17,52,50,48,51,54,49,55,53, 7, 2,4, 1, 3, 6, 0, 5,
		46,43,45,44,41,42,40,47, 2,4, 7, 1, 0, 6, 5, 3,35,33,34,38,39,37,36,32,11,14,15,13, 8,9,10,12,21,19,23,16,22,18,20,17,50,48,52,54,49,51,53,55,61,60,57,63,62,59,56,58,24,27,30,31,26,29,25,28,
		46,43,47,40,44,45,41,42,27,30,24,31,25,29,28,26, 4, 7, 2,1, 5, 6, 3, 0,14,15,11,8,9,13,10,12,48,52,50,49,51,54,53,55,60,57,61,62,59,63,58,56,19,23,21,22,18,16,17,20,33,34,35,38,36,37,32,39,
		46,44,40,42,45,43,41,47,52,50,48,55,54,51,53,49, 7, 2,4, 1, 0, 5, 3, 6,30,24,27,29,28,25,31,26,15,11,14,12, 8,13, 9,10,34,35,33,39,36,32,37,38,57,61,60,62,58,59,56,63,23,21,19,17,16,18,20,22,
		46,44,43,45,41,47,42,40, 4, 7, 2,1, 6, 5, 0, 3,19,23,21,17,22,18,16,20,27,30,24,25,29,28,31,26,60,57,61,59,62,58,63,56,14,15,11,13,12, 8,9,10,33,34,35,32,39,36,38,37,48,52,50,55,49,51,54,53,
		46,44,47,41,42,40,45,43,21,19,23,17,20,18,22,16,50,48,52,55,53,51,49,54,24,27,30,28,25,29,31,26,35,33,34,36,32,39,37,38,61,60,57,58,59,62,63,56,11,14,15, 8,13,12,10, 9, 2,4, 7, 1, 3, 5, 6, 0,
		46,45,40,41,47,42,43,44,60,57,61,56,63,58,62,59,14,15,11,10, 9, 8,12,13,48,52,50,53,54,51,55,49,19,23,21,16,20,22,18,17,33,34,35,37,36,39,32,38,27,30,24,29,25,26,31,28,4, 7, 2,1, 3, 0, 5, 6,
		46,45,42,47,43,44,41,40,11,14,15,10,13, 8,9,12, 2,4, 7, 1, 6, 0, 3, 5,50,48,52,51,53,54,55,49,24,27,30,26,29,25,28,31,21,19,23,22,16,20,18,17,35,33,34,39,37,36,38,32,61,60,57,56,59,58,63,62,
		46,45,44,43,41,40,47,42, 7, 2,4, 1, 5, 0, 6, 3,57,61,60,56,62,58,59,63,52,50,48,54,51,53,55,49,34,35,33,36,39,37,32,38,30,24,27,25,26,29,28,31,23,21,19,20,22,16,17,18,15,11,14,10,12, 8,13, 9,
		46,47,40,43,44,41,42,45,24,27,30,31,28,25,29,26,21,19,23,17,18,20,16,22,61,60,57,62,63,58,56,59, 2,4, 7, 3, 5, 6, 0, 1,35,33,34,32,36,37,39,38,50,48,52,53,51,49,55,54,11,14,15,10,12, 9, 8,13,
		46,47,41,44,42,45,43,40,23,21,19,17,22,20,18,16,15,11,14,10,13, 9,12, 8,57,61,60,58,62,63,56,59,52,50,48,49,53,51,54,55, 7, 2,4, 6, 3, 5, 0, 1,34,35,33,37,32,36,38,39,30,24,27,31,26,25,28,29,
		46,47,45,42,43,40,44,41,14,15,11,10, 8,9,13,12,27,30,24,31,29,25,26,28,60,57,61,63,58,62,56,59,33,34,35,36,37,32,39,38,48,52,50,51,49,53,54,55, 4, 7, 2,5, 6, 3, 1, 0,19,23,21,17,16,20,22,18,
		47,40,42,41,44,45,43,46, 3, 5, 6, 0, 1, 7, 4, 2,32,36,37,39,38,34,33,35,10,15,14,12, 9, 8,11,13,23,17,21,18,20,16,22,19,53,51,49,55,52,48,50,54,62,63,58,60,61,56,59,57,27,24,31,30,29,26,28,25,
		47,40,45,44,43,46,41,42,37,32,36,39,35,34,38,33,31,27,24,30,25,26,29,28,14,10,15, 8,12, 9,11,13,58,62,63,56,60,61,57,59,21,23,17,16,18,20,22,19,49,53,51,48,55,52,54,50, 6, 3, 5, 0, 2,7, 1, 4,
		47,40,46,43,41,42,44,45,24,31,27,30,28,26,25,29, 5, 6, 3, 0, 4, 7, 2,1,15,14,10, 9, 8,12,11,13,51,49,53,52,48,55,50,54,63,58,62,61,56,60,57,59,17,21,23,20,16,18,19,22,36,37,32,39,33,34,35,38,
		47,41,40,42,44,46,45,43, 5, 6, 3, 0, 7, 4, 1, 2,17,21,23,19,20,16,18,22,24,31,27,28,26,25,30,29,63,58,62,56,61,57,60,59,15,14,10,12,13, 9, 8,11,36,37,32,35,38,33,39,34,51,49,53,54,52,48,55,50,
		47,41,43,45,42,40,44,46,49,53,51,54,55,48,50,52, 6, 3, 5, 0, 1, 4, 2,7,31,27,24,26,25,28,30,29,14,10,15,13, 9,12, 8,11,37,32,36,38,33,35,34,39,58,62,63,61,57,56,59,60,21,23,17,19,18,16,22,20,
		47,41,46,44,45,43,42,40,23,17,21,19,22,16,20,18,53,51,49,54,50,48,52,55,27,24,31,25,28,26,30,29,32,36,37,33,35,38,34,39,62,63,58,57,56,61,60,59,10,15,14, 9,12,13,11,8,3, 5, 6, 0, 2,4, 7, 1,
		47,42,41,40,44,43,46,45, 6, 3, 5, 0, 4, 1, 7, 2,58,62,63,59,61,57,56,60,49,53,51,55,48,50,54,52,37,32,36,33,38,34,35,39,31,27,24,28,29,26,25,30,21,23,17,22,20,18,19,16,14,10,15,11,13, 9,12, 8,
		47,42,43,44,46,45,40,41,63,58,62,59,60,57,61,56,15,14,10,11,8,9,13,12,51,49,53,50,55,48,54,52,17,21,23,18,22,20,16,19,36,37,32,34,33,38,35,39,24,31,27,26,28,29,30,25, 5, 6, 3, 0, 2,1, 4, 7,
		47,42,45,46,40,41,44,43,10,15,14,11,12, 9, 8,13, 3, 5, 6, 0, 7, 1, 2,4,53,51,49,48,50,55,54,52,27,24,31,29,26,28,25,30,23,17,21,20,18,22,16,19,32,36,37,38,34,33,39,35,62,63,58,59,56,57,60,61,
		47,43,40,46,41,45,42,44,31,27,24,30,26,25,28,29,49,53,51,54,48,55,52,50,37,32,36,35,34,38,39,33,21,23,17,18,16,22,20,19,14,10,15, 9,13, 8,12,11,6, 3, 5, 1, 4, 2,0, 7,58,62,63,59,56,60,61,57,
		47,43,44,42,46,40,41,45,62,63,58,59,61,60,57,56,27,24,31,30,28,25,29,26,32,36,37,34,38,35,39,33,10,15,14,13, 8,9,12,11,3, 5, 6, 4, 2,1, 7, 0,23,17,21,16,22,18,19,20,53,51,49,54,52,55,50,48,
		47,43,45,41,42,44,46,40,51,49,53,54,50,55,48,52,63,58,62,59,57,60,56,61,36,37,32,38,35,34,39,33, 5, 6, 3, 2,1, 4, 7, 0,17,21,23,22,18,16,20,19,15,14,10, 8,9,13,11,12,24,31,27,30,29,25,26,28,
		47,44,40,45,43,42,46,41,32,36,37,39,34,38,35,33,62,63,58,59,60,61,56,57, 3, 5, 6, 1, 7, 4, 0, 2,53,51,49,52,55,50,48,54,10,15,14, 8,13,12, 9,11,27,24,31,28,25,29,30,26,23,17,21,19,18,20,16,22,
		47,44,41,46,45,40,43,42,17,21,23,19,16,20,22,18,36,37,32,39,35,38,33,34, 5, 6, 3, 7, 4, 1, 0, 2,15,14,10,13,12, 8,9,11,24,31,27,25,29,28,26,30,51,49,53,55,50,52,54,48,63,58,62,59,56,61,57,60,
		47,44,42,43,46,41,45,40,58,62,63,59,57,61,60,56,21,23,17,19,22,20,18,16, 6, 3, 5, 4, 1, 7, 0, 2,31,27,24,29,28,25,26,30,49,53,51,50,52,55,48,54,14,10,15,12, 8,13,11,9,37,32,36,39,33,38,34,35,
		47,45,41,43,42,46,40,44,53,51,49,54,48,50,55,52,10,15,14,11,9,12,13, 8,23,17,21,22,16,20,19,18,62,63,58,56,57,60,61,59,27,24,31,26,29,25,28,30, 3, 5, 6, 7, 1, 2,0, 4,32,36,37,39,33,35,38,34,
		47,45,44,40,43,41,42,46,36,37,32,39,38,35,34,33,51,49,53,54,55,50,52,48,17,21,23,16,20,22,19,18,24,31,27,29,25,26,28,30, 5, 6, 3, 1, 2,7, 4, 0,63,58,62,57,60,56,59,61,15,14,10,11,13,12, 8,9,
		47,45,46,42,40,44,43,41,14,10,15,11,8,12, 9,13,37,32,36,39,34,35,33,38,21,23,17,20,22,16,19,18,6, 3, 5, 2,7, 1, 4, 0,58,62,63,60,56,57,61,59,31,27,24,25,26,29,30,28,49,53,51,54,52,50,48,55,
		47,46,42,45,40,43,41,44,15,14,10,11,9, 8,12,13,24,31,27,30,26,28,29,25,63,58,62,60,57,61,59,56,36,37,32,33,34,35,38,39,51,49,53,48,52,50,55,54, 5, 6, 3, 4, 7, 2,0, 1,17,21,23,19,18,22,20,16,
		47,46,43,40,41,44,45,42,27,24,31,30,25,28,26,29,23,17,21,19,16,22,18,20,62,63,58,61,60,57,59,56, 3, 5, 6, 2,4, 7, 1, 0,32,36,37,35,33,34,38,39,53,51,49,50,48,52,54,55,10,15,14,11,13, 8,9,12,
		47,46,44,41,45,42,40,43,21,23,17,19,20,22,16,18,14,10,15,11,12, 8,13, 9,58,62,63,57,61,60,59,56,49,53,51,52,50,48,55,54, 6, 3, 5, 7, 2,4, 1, 0,37,32,36,34,35,33,39,38,31,27,24,30,29,28,25,26,
		48,49,50,54,52,51,55,53, 0, 6, 4, 2,7, 1, 3, 5, 8,15,14, 9,12,11,10,13,24,29,27,31,25,30,26,28,16,20,19,21,17,22,18,23,56,59,60,58,57,62,61,63,40,42,46,47,41,43,44,45,32,38,33,37,39,35,34,36,
		48,49,51,52,55,53,54,50,14, 8,15, 9,13,11,12,10,33,32,38,37,36,35,39,34,27,24,29,30,31,25,26,28,46,40,42,43,47,41,45,44,19,16,20,22,21,17,18,23,60,56,59,62,58,57,63,61,4, 0, 6, 2,5, 1, 7, 3,
		48,49,53,55,54,50,52,51,38,33,32,37,34,35,36,39, 6, 4, 0, 2,3, 1, 5, 7,29,27,24,25,30,31,26,28,59,60,56,57,62,58,61,63,42,46,40,41,43,47,45,44,20,19,16,17,22,21,23,18,15,14, 8,9,10,11,13,12,
		48,50,51,53,49,54,52,55,24,29,27,26,31,25,30,28,0, 6, 4, 2,1, 7, 5, 3,56,59,60,62,61,58,63,57,32,38,33,39,35,34,36,37,16,20,19,17,21,18,22,23, 8,15,14,12,11,10, 9,13,40,42,46,44,43,45,47,41,
		48,50,54,49,52,55,53,51,4, 0, 6, 2,3, 7, 1, 5,46,40,42,44,41,45,43,47,60,56,59,58,62,61,63,57,14, 8,15,10,12,11,13, 9,33,32,38,34,39,35,36,37,19,16,20,18,17,21,23,22,27,24,29,26,28,25,31,30,
		48,50,55,52,53,51,49,54,42,46,40,44,47,45,41,43,29,27,24,26,30,25,28,31,59,60,56,61,58,62,63,57,20,19,16,21,18,17,22,23,15,14, 8,11,10,12,13, 9,38,33,32,35,34,39,37,36, 6, 4, 0, 2,5, 7, 3, 1,
		48,51,52,49,55,54,50,53,15,14, 8,9,12,13,11,10,59,60,56,63,58,61,57,62,20,19,16,22,17,18,23,21,38,33,32,39,36,35,34,37, 6, 4, 0, 7, 5, 1, 3, 2,42,46,40,45,47,43,44,41,29,27,24,26,28,31,30,25,
		48,51,53,50,49,52,55,54,27,24,29,26,30,31,25,28,14, 8,15, 9,11,13,10,12,19,16,20,17,18,22,23,21,4, 0, 6, 5, 1, 7, 3, 2,46,40,42,47,43,45,41,44,33,32,38,36,35,39,37,34,60,56,59,63,57,61,62,58,
		48,51,54,55,50,53,49,52,56,59,60,63,62,61,58,57,24,29,27,26,25,31,28,30,16,20,19,18,22,17,23,21,40,42,46,43,45,47,41,44,32,38,33,35,39,36,34,37, 0, 6, 4, 1, 7, 5, 2,3, 8,15,14, 9,10,13,12,11,
		48,52,49,51,55,50,53,54, 8,15,14, 9,11,12,13,10,40,42,46,44,47,41,43,45, 0, 6, 4, 7, 1, 3, 2,5,56,59,60,57,58,61,62,63,24,29,27,30,28,31,25,26,32,38,33,34,36,39,37,35,16,20,19,23,21,17,22,18,
		48,52,50,55,53,54,51,49,46,40,42,44,45,41,47,43,19,16,20,23,18,17,21,22, 4, 0, 6, 3, 7, 1, 2,5,33,32,38,39,34,36,35,37,60,56,59,61,57,58,62,63,27,24,29,31,30,28,26,25,14, 8,15, 9,10,12,11,13,
		48,52,54,53,51,49,55,50,20,19,16,23,22,17,18,21,15,14, 8,9,13,12,10,11,6, 4, 0, 1, 3, 7, 2,5,29,27,24,28,31,30,25,26,38,33,32,36,39,34,35,37,59,60,56,58,61,57,63,62,42,46,40,44,43,41,45,47,
		48,53,50,51,49,55,54,52,29,27,24,26,25,30,31,28,38,33,32,37,35,34,39,36,42,46,40,47,45,41,44,43,15,14, 8,10,11,13,12, 9,59,60,56,62,57,61,58,63, 6, 4, 0, 3, 1, 5, 2,7,20,19,16,23,21,18,17,22,
		48,53,52,54,51,50,49,55,19,16,20,23,17,18,22,21,27,24,29,26,31,30,28,25,46,40,42,45,41,47,44,43,60,56,59,57,61,62,58,63, 4, 0, 6, 1, 5, 3, 7, 2,14, 8,15,11,13,10, 9,12,33,32,38,37,39,34,36,35,
		48,53,55,49,54,52,51,50,32,38,33,37,36,34,35,39,16,20,19,23,22,18,21,17,40,42,46,41,47,45,44,43, 0, 6, 4, 5, 3, 1, 7, 2,8,15,14,13,10,11,12, 9,56,59,60,61,62,57,63,58,24,29,27,26,28,30,25,31,
		48,54,49,50,52,53,51,55, 6, 4, 0, 2,1, 3, 7, 5,20,19,16,23,17,22,21,18,38,33,32,34,35,36,37,39,42,46,40,43,41,45,47,44,29,27,24,31,28,25,30,26,15,14, 8,13,12,10, 9,11,59,60,56,63,57,62,58,61,
		48,54,53,52,51,55,50,49,16,20,19,23,18,22,17,21,56,59,60,63,61,62,57,58,32,38,33,36,34,35,37,39, 8,15,14,10,13,12,11,9,40,42,46,45,43,41,47,44,24,29,27,25,31,28,26,30, 0, 6, 4, 2,5, 3, 1, 7,
		48,54,55,51,50,49,52,53,60,56,59,63,58,62,61,57, 4, 0, 6, 2,7, 3, 5, 1,33,32,38,35,36,34,37,39,27,24,29,28,25,31,30,26,14, 8,15,12,10,13,11,9,46,40,42,41,45,43,44,47,19,16,20,23,21,22,18,17,
		48,55,49,53,54,51,50,52,33,32,38,37,35,36,34,39,60,56,59,63,62,58,57,61,14, 8,15,13,11,12, 9,10,19,16,20,21,22,18,17,23,27,24,29,25,28,30,31,26, 4, 0, 6, 7, 3, 5, 2,1,46,40,42,44,43,47,41,45,
		48,55,51,54,50,52,53,49,59,60,56,63,61,58,62,57,42,46,40,44,45,47,43,41,15,14, 8,12,13,11,9,10, 6, 4, 0, 5, 7, 3, 1, 2,20,19,16,18,21,22,17,23,29,27,24,30,25,28,26,31,38,33,32,37,39,36,35,34,
		48,55,52,50,53,49,54,51,40,42,46,44,41,47,45,43,32,38,33,37,34,36,39,35, 8,15,14,11,12,13, 9,10,24,29,27,28,30,25,31,26, 0, 6, 4, 3, 5, 7, 1, 2,16,20,19,22,18,21,23,17,56,59,60,63,57,58,61,62,
		49,48,52,51,53,55,50,54, 8,14, 9,15,11,13,10,12,38,37,33,32,39,34,36,35,30,31,25,27,24,29,28,26,43,47,41,46,40,42,44,45,22,21,17,19,16,20,23,18,62,58,57,60,56,59,61,63, 2,6, 0, 4, 3, 7, 1, 5,
		49,48,54,50,51,52,53,55, 6, 0, 2,4, 1, 7, 5, 3,14, 9, 8,15,10,13,12,11,31,25,30,24,29,27,28,26,21,17,22,16,20,19,23,18,58,57,62,56,59,60,63,61,47,41,43,40,42,46,45,44,37,33,38,32,36,34,35,39,
		49,48,55,53,50,54,51,52,33,38,37,32,35,34,39,36, 0, 2,6, 4, 5, 7, 3, 1,25,30,31,29,27,24,28,26,57,62,58,59,60,56,63,61,41,43,47,42,46,40,44,45,17,22,21,20,19,16,18,23, 9, 8,14,15,12,13,11,10,
		49,50,48,54,51,55,52,53, 0, 2,6, 4, 7, 5, 1, 3,17,22,21,18,20,19,16,23,33,38,37,35,34,39,32,36,41,43,47,46,42,44,40,45,25,30,31,24,26,29,27,28,9, 8,14,11,10,12,15,13,57,62,58,61,59,60,56,63,
		49,50,53,52,54,48,51,55,62,58,57,61,56,60,63,59, 2,6, 0, 4, 1, 5, 3, 7,38,37,33,34,39,35,32,36,30,31,25,26,29,24,27,28,8,14, 9,10,12,11,13,15,43,47,41,42,44,46,45,40,22,21,17,18,16,19,23,20,
		49,50,55,51,52,53,54,48,21,17,22,18,23,19,20,16,58,57,62,61,63,60,59,56,37,33,38,39,35,34,32,36,14, 9, 8,12,11,10,13,15,47,41,43,44,46,42,40,45,31,25,30,29,24,26,28,27, 6, 0, 2,4, 3, 5, 7, 1,
		49,51,48,52,53,54,55,50,14, 9, 8,15,13,10,11,12,47,41,43,45,40,42,46,44, 6, 0, 2,1, 7, 5, 4, 3,58,57,62,59,56,63,60,61,31,25,30,27,26,24,29,28,37,33,38,35,39,36,32,34,21,17,22,18,16,20,19,23,
		49,51,50,55,52,48,53,54,17,22,21,18,19,20,23,16, 9, 8,14,15,11,10,12,13, 0, 2,6, 7, 5, 1, 4, 3,25,30,31,26,24,27,29,28,33,38,37,39,36,35,34,32,57,62,58,56,63,59,61,60,41,43,47,45,46,42,44,40,
		49,51,54,53,55,50,52,48,43,47,41,45,44,42,40,46,22,21,17,18,23,20,16,19, 2,6, 0, 5, 1, 7, 4, 3,38,37,33,36,35,39,34,32,62,58,57,63,59,56,60,61,30,31,25,24,27,26,28,29, 8,14, 9,15,12,10,13,11,
		49,52,50,53,54,55,48,51,58,57,62,61,60,63,56,59,31,25,30,28,29,24,26,27,21,17,22,23,19,20,18,16,47,41,43,46,44,40,42,45,37,33,38,34,36,39,35,32, 6, 0, 2,7, 1, 3, 4, 5,14, 9, 8,15,12,11,10,13,
		49,52,51,48,53,50,54,55, 9, 8,14,15,10,11,13,12,57,62,58,61,56,63,59,60,17,22,21,19,20,23,18,16,33,38,37,36,39,34,35,32, 0, 2,6, 1, 3, 7, 5, 4,41,43,47,44,40,46,45,42,25,30,31,28,26,24,27,29,
		49,52,55,54,48,51,53,50,30,31,25,28,27,24,29,26, 8,14, 9,15,13,11,12,10,22,21,17,20,23,19,18,16, 2,6, 0, 3, 7, 1, 5, 4,43,47,41,40,46,44,42,45,38,37,33,39,34,36,32,35,62,58,57,61,59,63,60,56,
		49,53,48,55,50,52,54,51,38,37,33,32,34,39,35,36,62,58,57,61,60,56,59,63, 8,14, 9,11,13,10,15,12,22,21,17,16,19,23,20,18,30,31,25,29,26,27,24,28,2,6, 0, 1, 5, 3, 4, 7,43,47,41,45,46,40,42,44,
		49,53,51,54,55,48,50,52,47,41,43,45,42,40,44,46,37,33,38,32,35,39,36,34,14, 9, 8,13,10,11,15,12,31,25,30,26,27,29,24,28,6, 0, 2,5, 3, 1, 7, 4,21,17,22,19,23,16,18,20,58,57,62,61,59,56,63,60,
		49,53,52,50,54,51,55,48,57,62,58,61,63,56,60,59,41,43,47,45,44,40,46,42, 9, 8,14,10,11,13,15,12, 0, 2,6, 3, 1, 5, 7, 4,17,22,21,23,16,19,20,18,25,30,31,27,29,26,28,24,33,38,37,32,36,39,34,35,
		49,54,50,48,51,53,55,52, 2,6, 0, 4, 5, 1, 7, 3,43,47,41,45,42,44,46,40,62,58,57,56,60,63,61,59, 8,14, 9,12,10,13,11,15,38,37,33,35,36,34,39,32,22,21,17,23,20,16,18,19,30,31,25,28,26,29,24,27,
		49,54,52,55,48,50,51,53,31,25,30,28,24,29,27,26, 6, 0, 2,4, 7, 1, 3, 5,58,57,62,60,63,56,61,59,37,33,38,36,34,35,39,32,21,17,22,20,16,23,19,18,14, 9, 8,10,13,12,15,11,47,41,43,45,46,44,40,42,
		49,54,53,51,55,52,48,50,41,43,47,45,40,44,42,46,25,30,31,28,27,29,26,24,57,62,58,63,56,60,61,59,17,22,21,16,23,20,19,18,9, 8,14,13,12,10,11,15,33,38,37,34,35,36,32,39, 0, 2,6, 4, 3, 1, 5, 7,
		49,55,51,50,52,54,48,53,22,21,17,18,20,23,19,16,30,31,25,28,24,27,26,29,43,47,41,44,42,40,45,46,62,58,57,59,63,60,56,61,2,6, 0, 7, 3, 5, 1, 4, 8,14, 9,13,11,12,15,10,38,37,33,32,36,35,39,34,
		49,55,53,48,50,51,52,54,37,33,38,32,39,35,34,36,21,17,22,18,19,23,16,20,47,41,43,42,40,44,45,46, 6, 0, 2,3, 5, 7, 1, 4,14, 9, 8,11,12,13,10,15,58,57,62,63,60,59,61,56,31,25,30,28,26,27,29,24,
		49,55,54,52,48,53,50,51,25,30,31,28,29,27,24,26,33,38,37,32,34,35,36,39,41,43,47,40,44,42,45,46, 9, 8,14,12,13,11,10,15,57,62,58,60,59,63,56,61,0, 2,6, 5, 7, 3, 4, 1,17,22,21,18,16,23,20,19,
		50,48,49,54,55,52,51,53, 0, 4, 2,6, 7, 3, 5, 1,42,44,46,40,43,47,41,45,58,62,61,60,56,59,57,63,10,12,11,14, 8,15, 9,13,34,39,35,33,32,38,37,36,18,17,21,19,16,20,22,23,26,29,24,27,30,31,25,28,
		50,48,52,55,51,53,54,49,46,42,44,40,45,47,43,41,24,26,29,27,28,31,30,25,61,58,62,59,60,56,57,63,21,18,17,20,19,16,23,22,11,10,12,15,14, 8,9,13,35,34,39,38,33,32,36,37, 2,0, 4, 6, 1, 3, 7, 5,
		50,48,53,51,54,49,55,52,29,24,26,27,25,31,28,30, 4, 2,0, 6, 5, 3, 1, 7,62,61,58,56,59,60,57,63,39,35,34,32,38,33,37,36,17,21,18,16,20,19,23,22,12,11,10, 8,15,14,13, 9,44,46,42,40,41,47,45,43,
		50,49,51,55,53,52,48,54,17,21,18,22,19,23,16,20,62,61,58,57,59,56,63,60,39,35,34,37,33,38,36,32,12,11,10,14, 9, 8,15,13,44,46,42,47,41,43,45,40,29,24,26,31,25,30,27,28,4, 2,0, 6, 1, 7, 5, 3,
		50,49,52,53,48,54,55,51,58,62,61,57,60,56,59,63, 0, 4, 2,6, 3, 7, 1, 5,34,39,35,38,37,33,36,32,26,29,24,30,31,25,28,27,10,12,11,8,14, 9,15,13,42,44,46,43,47,41,40,45,18,17,21,22,20,23,19,16,
		50,49,54,48,55,51,53,52, 2,0, 4, 6, 5, 7, 3, 1,21,18,17,22,16,23,20,19,35,34,39,33,38,37,36,32,46,42,44,41,43,47,45,40,24,26,29,25,30,31,28,27,11,10,12, 9, 8,14,13,15,61,58,62,57,63,56,60,59,
		50,51,48,53,54,52,49,55,24,26,29,27,31,28,25,30,35,34,39,36,38,33,32,37,46,42,44,45,47,43,40,41,11,10,12,14,15, 9, 8,13,61,58,62,56,63,59,60,57, 2,0, 4, 7, 5, 1, 6, 3,21,18,17,22,20,19,16,23,
		50,51,52,54,49,55,53,48,39,35,34,36,37,33,38,32,17,21,18,22,23,19,20,16,44,46,42,43,45,47,40,41,4, 2,0, 1, 7, 5, 3, 6,12,11,10, 9,14,15, 8,13,62,61,58,59,56,63,57,60,29,24,26,27,30,28,31,25,
		50,51,55,49,53,48,54,52,18,17,21,22,16,19,23,20,26,29,24,27,25,28,30,31,42,44,46,47,43,45,40,41,58,62,61,63,59,56,60,57, 0, 4, 2,5, 1, 7, 3, 6,10,12,11,15, 9,14,13, 8,34,39,35,36,32,33,37,38,
		50,52,53,49,48,55,51,54,61,58,62,57,59,60,56,63,46,42,44,40,47,45,41,43,11,10,12, 8,9,15,13,14, 2,0, 4, 1, 3, 7, 5, 6,21,18,17,19,20,23,16,22,24,26,29,28,31,30,27,25,35,34,39,36,32,37,38,33,
		50,52,54,51,49,53,48,55,34,39,35,36,38,37,33,32,58,62,61,57,56,60,63,59,10,12,11,9,15, 8,13,14,18,17,21,20,23,19,16,22,26,29,24,31,30,28,25,27, 0, 4, 2,3, 7, 1, 6, 5,42,44,46,40,41,45,43,47,
		50,52,55,48,51,54,49,53,44,46,42,40,43,45,47,41,39,35,34,36,33,37,32,38,12,11,10,15, 8,9,13,14,29,24,26,30,28,31,25,27, 4, 2,0, 7, 1, 3, 5, 6,17,21,18,23,19,20,22,16,62,61,58,57,63,60,59,56,
		50,53,49,52,48,51,54,55,62,61,58,57,56,59,60,63,29,24,26,27,31,25,30,28,17,21,18,19,23,16,22,20,44,46,42,41,47,45,43,40,39,35,34,38,32,37,33,36, 4, 2,0, 5, 3, 1, 6, 7,12,11,10,13,14, 9, 8,15,
		50,53,51,48,54,55,52,49,26,29,24,27,28,25,31,30,10,12,11,13,15, 9,14, 8,18,17,21,16,19,23,22,20, 0, 4, 2,1, 5, 3, 7, 6,42,44,46,45,41,47,43,40,34,39,35,37,38,32,36,33,58,62,61,57,63,59,56,60,
		50,53,55,54,52,49,48,51,11,10,12,13, 8,9,15,14,61,58,62,57,60,59,63,56,21,18,17,23,16,19,22,20,35,34,39,32,37,38,33,36, 2,0, 4, 3, 1, 5, 7, 6,46,42,44,47,45,41,40,43,24,26,29,27,30,25,28,31,
		50,54,48,49,55,53,52,51,4, 2,0, 6, 3, 5, 7, 1,12,11,10,13, 8,15,14, 9,29,24,26,25,31,28,27,30,17,21,18,20,16,23,19,22,62,61,58,60,63,56,59,57,44,46,42,45,43,41,40,47,39,35,34,36,32,38,33,37,
		50,54,51,52,49,48,55,53,35,34,39,36,33,38,37,32, 2,0, 4, 6, 7, 5, 1, 3,24,26,29,31,28,25,27,30,61,58,62,63,56,60,59,57,46,42,44,43,41,45,47,40,21,18,17,16,23,20,22,19,11,10,12,13,14,15, 9, 8,
		50,54,53,55,52,51,49,48,10,12,11,13, 9,15, 8,14,34,39,35,36,37,38,32,33,26,29,24,28,25,31,27,30,42,44,46,41,45,43,47,40,18,17,21,23,20,16,19,22,58,62,61,56,60,63,57,59, 0, 4, 2,6, 1, 5, 3, 7,
		50,55,48,52,51,49,53,54,42,44,46,40,47,43,45,41,18,17,21,22,19,16,20,23, 0, 4, 2,7, 3, 5, 6, 1,34,39,35,32,33,37,38,36,58,62,61,59,63,60,56,57,26,29,24,25,28,30,27,31,10,12,11,13,14, 8,15, 9,
		50,55,49,51,53,54,52,48,21,18,17,22,23,16,19,20,11,10,12,13, 9, 8,14,15, 2,0, 4, 5, 7, 3, 6, 1,24,26,29,30,25,28,31,27,35,34,39,37,32,33,38,36,61,58,62,60,59,63,57,56,46,42,44,40,41,43,47,45,
		50,55,54,53,52,48,51,49,12,11,10,13,15, 8,9,14,44,46,42,40,45,43,41,47, 4, 2,0, 3, 5, 7, 6, 1,62,61,58,63,60,59,56,57,29,24,26,28,30,25,31,27,39,35,34,33,37,32,36,38,17,21,18,22,20,16,23,19,
		51,48,49,52,54,55,53,50,14,15, 9, 8,13,12,10,11,56,63,59,60,57,62,58,61,22,17,18,20,19,16,21,23,39,36,35,38,33,32,37,34, 7, 5, 1, 6, 4, 0, 2,3,45,47,43,42,46,40,41,44,26,24,27,29,25,30,31,28,
		51,48,50,53,52,49,54,55,24,27,26,29,31,30,28,25,15, 9,14, 8,10,12,11,13,17,18,22,19,16,20,21,23, 5, 1, 7, 4, 0, 6, 2,3,47,43,45,46,40,42,44,41,36,35,39,33,32,38,34,37,63,59,56,60,58,62,61,57,
		51,48,55,54,53,50,52,49,59,56,63,60,61,62,57,58,27,26,24,29,28,30,25,31,18,22,17,16,20,19,21,23,43,45,47,40,42,46,44,41,35,39,36,32,38,33,37,34, 1, 7, 5, 0, 6, 4, 3, 2,9,14,15, 8,11,12,13,10,
		51,49,52,48,54,53,50,55, 9,14,15, 8,10,13,12,11,43,45,47,41,46,44,40,42, 1, 7, 5, 6, 0, 2,3, 4,59,56,63,58,57,62,61,60,27,26,24,31,25,30,28,29,35,39,36,37,33,38,34,32,18,22,17,21,23,19,20,16,
		51,49,53,54,50,55,48,52,47,43,45,41,42,44,46,40,17,18,22,21,16,19,23,20, 5, 1, 7, 2,6, 0, 3, 4,36,35,39,38,37,33,32,34,63,59,56,62,58,57,61,60,24,27,26,30,31,25,29,28,15, 9,14, 8,11,13,10,12,
		51,49,55,50,48,52,54,53,22,17,18,21,20,19,16,23,14,15, 9, 8,12,13,11,10, 7, 5, 1, 0, 2,6, 3, 4,26,24,27,25,30,31,28,29,39,36,35,33,38,37,32,34,56,63,59,57,62,58,60,61,45,47,43,41,40,44,42,46,
		51,50,49,55,48,53,52,54,17,18,22,21,19,16,20,23,24,27,26,29,30,31,25,28,47,43,45,42,44,46,41,40,63,59,56,58,62,61,57,60, 5, 1, 7, 0, 4, 2,6, 3,15, 9,14,10,12,11,8,13,36,35,39,34,38,37,33,32,
		51,50,53,48,52,54,55,49,26,24,27,29,28,31,30,25,39,36,35,34,32,37,38,33,45,47,43,46,42,44,41,40,14,15, 9,11,10,12,13, 8,56,63,59,61,58,62,57,60, 7, 5, 1, 2,0, 4, 3, 6,22,17,18,21,23,16,19,20,
		51,50,54,52,55,49,48,53,35,39,36,34,33,37,32,38,18,22,17,21,20,16,23,19,43,45,47,44,46,42,41,40, 1, 7, 5, 4, 2,0, 6, 3, 9,14,15,12,11,10,13, 8,59,56,63,62,61,58,60,57,27,26,24,29,25,31,28,30,
		51,52,48,49,54,50,55,53,15, 9,14, 8,12,10,13,11,36,35,39,34,33,32,38,37,24,27,26,31,30,28,29,25,47,43,45,40,46,44,42,41,17,18,22,20,23,19,16,21,63,59,56,61,57,58,60,62, 5, 1, 7, 3, 4, 0, 6, 2,
		51,52,50,54,55,53,49,48,39,36,35,34,37,32,33,38,7, 5, 1, 3, 2,0, 4, 6,26,24,27,28,31,30,29,25,56,63,59,58,61,57,62,60,45,47,43,44,40,46,42,41,22,17,18,19,20,23,21,16,14,15, 9, 8,11,10,12,13,
		51,52,53,55,49,48,54,50, 1, 7, 5, 3, 6, 0, 2,4, 9,14,15, 8,13,10,11,12,27,26,24,30,28,31,29,25,18,22,17,23,19,20,16,21,59,56,63,57,58,61,62,60,43,45,47,46,44,40,41,42,35,39,36,34,38,32,37,33,
		51,53,48,50,52,55,49,54,27,26,24,29,30,28,31,25, 1, 7, 5, 3, 0, 6, 4, 2,59,56,63,61,62,57,60,58,35,39,36,38,32,37,33,34,18,22,17,19,23,16,20,21,9,14,15,13,10,11,8,12,43,45,47,41,40,42,46,44,
		51,53,54,49,50,48,52,55,45,47,43,41,46,42,44,40,26,24,27,29,31,28,25,30,56,63,59,62,57,61,60,58,22,17,18,23,16,19,20,21,14,15, 9,10,11,13,12, 8,39,36,35,32,37,38,34,33, 7, 5, 1, 3, 4, 6, 2,0,
		51,53,55,52,49,54,50,48,5, 1, 7, 3, 2,6, 0, 4,47,43,45,41,44,42,40,46,63,59,56,57,61,62,60,58,15, 9,14,11,13,10,12, 8,36,35,39,37,38,32,33,34,17,18,22,16,19,23,21,20,24,27,26,29,25,28,30,31,
		51,54,48,55,53,49,50,52,56,63,59,60,62,57,61,58,45,47,43,41,42,46,40,44,14,15, 9,13,12,10, 8,11,7, 5, 1, 4, 6, 2,0, 3,22,17,18,16,23,20,19,21,26,24,27,31,28,25,29,30,39,36,35,34,38,33,32,37,
		51,54,49,53,50,52,55,48,43,45,47,41,44,46,42,40,35,39,36,34,37,33,38,32, 9,14,15,10,13,12, 8,11,27,26,24,25,31,28,30,29, 1, 7, 5, 2,4, 6, 0, 3,18,22,17,20,16,23,21,19,59,56,63,60,58,57,62,61,
		51,54,52,50,55,48,53,49,36,35,39,34,32,33,37,38,63,59,56,60,61,57,58,62,15, 9,14,12,10,13, 8,11,17,18,22,23,20,16,19,21,24,27,26,28,25,31,30,29, 5, 1, 7, 6, 2,4, 3, 0,47,43,45,41,40,46,44,42,
		51,55,50,49,48,54,53,52,18,22,17,21,16,20,19,23,59,56,63,60,62,61,58,57,35,39,36,33,37,32,34,38,9,14,15,11,12,13,10, 8,43,45,47,42,40,44,46,41,27,26,24,28,30,25,29,31,1, 7, 5, 3, 4, 2,0, 6,
		51,55,52,53,49,50,48,54, 7, 5, 1, 3, 0, 2,6, 4,22,17,18,21,19,20,23,16,39,36,35,37,32,33,34,38,45,47,43,40,44,42,46,41,26,24,27,30,25,28,31,29,14,15, 9,12,13,11,8,10,56,63,59,60,58,61,57,62,
		51,55,54,48,53,52,49,50,63,59,56,60,57,61,62,58,5, 1, 7, 3, 6, 2,4, 0,36,35,39,32,33,37,34,38,24,27,26,25,28,30,31,29,15, 9,14,13,11,12,10, 8,47,43,45,44,42,40,41,46,17,18,22,21,23,20,16,19,
		52,48,51,49,50,55,54,53,15, 8,9,14,12,11,10,13,46,44,40,42,43,45,47,41,7, 1, 3, 0, 6, 4, 5, 2,57,58,61,56,59,60,63,62,30,28,31,24,29,27,26,25,34,36,39,32,38,33,35,37,23,19,20,16,18,22,17,21,
		52,48,53,54,49,51,50,55,19,20,23,16,17,22,21,18,8,9,15,14,10,11,13,12, 1, 3, 7, 6, 4, 0, 5, 2,28,31,30,29,27,24,26,25,36,39,34,38,33,32,37,35,58,61,57,59,60,56,62,63,44,40,46,42,47,45,41,43,
		52,48,55,50,54,53,49,51,40,46,44,42,41,45,43,47,20,23,19,16,21,22,18,17, 3, 7, 1, 4, 0, 6, 5, 2,39,34,36,33,32,38,37,35,61,57,58,60,56,59,63,62,31,30,28,27,24,29,25,26, 9,15, 8,14,13,11,12,10,
		52,49,48,51,50,53,55,54, 8,9,15,14,11,10,12,13,58,61,57,62,59,60,56,63,19,20,23,17,22,21,16,18,36,39,34,33,38,37,32,35, 1, 3, 7, 0, 2,6, 4, 5,44,40,46,41,43,47,42,45,28,31,30,25,29,27,24,26,
		52,49,53,50,55,54,51,48,57,58,61,62,63,60,59,56,30,28,31,25,26,27,29,24,23,19,20,21,17,22,16,18,46,44,40,47,41,43,45,42,34,36,39,37,33,38,32,35, 7, 1, 3, 6, 0, 2,5, 4,15, 8,9,14,13,10,11,12,
		52,49,54,55,51,48,50,53,31,30,28,25,24,27,26,29, 9,15, 8,14,12,10,13,11,20,23,19,22,21,17,16,18,3, 7, 1, 2,6, 0, 4, 5,40,46,44,43,47,41,45,42,39,34,36,38,37,33,35,32,61,57,58,62,56,60,63,59,
		52,50,48,55,54,51,53,49,46,44,40,42,45,43,41,47,34,36,39,35,32,38,33,37,15, 8,9,12,11,10,14,13,30,28,31,29,24,26,27,25, 7, 1, 3, 4, 2,0, 6, 5,23,19,20,17,21,18,16,22,57,58,61,62,56,59,60,63,
		52,50,49,53,55,48,54,51,58,61,57,62,60,59,63,56,44,40,46,42,41,43,47,45, 8,9,15,11,10,12,14,13, 1, 3, 7, 2,0, 4, 6, 5,19,20,23,21,18,17,22,16,28,31,30,24,26,29,25,27,36,39,34,35,33,38,37,32,
		52,50,51,54,53,49,55,48,39,34,36,35,37,38,32,33,61,57,58,62,63,59,56,60, 9,15, 8,10,12,11,14,13,20,23,19,18,17,21,22,16,31,30,28,26,29,24,27,25, 3, 7, 1, 0, 4, 2,5, 6,40,46,44,42,47,43,45,41,
		52,51,49,48,50,54,53,55, 9,15, 8,14,10,12,11,13,39,34,36,35,38,37,33,32,31,30,28,24,27,26,25,29,40,46,44,47,43,45,41,42,20,23,19,17,18,22,21,16,61,57,58,63,59,56,62,60, 3, 7, 1, 5, 2,6, 0, 4,
		52,51,54,50,53,55,48,49,36,39,34,35,32,37,38,33, 1, 3, 7, 5, 4, 6, 2,0,28,31,30,26,24,27,25,29,58,61,57,56,63,59,60,62,44,40,46,45,47,43,41,42,19,20,23,22,17,18,16,21,8,9,15,14,13,12,10,11,
		52,51,55,53,48,49,50,54, 7, 1, 3, 5, 0, 6, 4, 2,15, 8,9,14,11,12,13,10,30,28,31,27,26,24,25,29,23,19,20,18,22,17,21,16,57,58,61,59,56,63,60,62,46,44,40,43,45,47,42,41,34,36,39,35,33,37,32,38,
		52,53,50,49,55,51,48,54,61,57,58,62,59,63,60,56, 3, 7, 1, 5, 0, 4, 2,6,39,34,36,37,38,32,35,33,31,30,28,29,26,27,24,25, 9,15, 8,11,13,10,12,14,40,46,44,45,41,47,42,43,20,23,19,16,18,17,21,22,
		52,53,51,55,48,54,49,50, 1, 3, 7, 5, 6, 4, 0, 2,19,20,23,16,22,17,18,21,36,39,34,32,37,38,35,33,44,40,46,47,45,41,43,42,28,31,30,27,29,26,24,25, 8,9,15,10,11,13,14,12,58,61,57,62,56,63,59,60,
		52,53,54,48,49,50,55,51,23,19,20,16,21,17,22,18,57,58,61,62,60,63,56,59,34,36,39,38,32,37,35,33,15, 8,9,13,10,11,12,14,46,44,40,41,47,45,43,42,30,28,31,26,27,29,25,24, 7, 1, 3, 5, 2,4, 6, 0,
		52,54,48,53,49,55,51,50,20,23,19,16,22,21,17,18,31,30,28,25,27,24,29,26,40,46,44,41,45,43,42,47,61,57,58,56,60,63,59,62, 3, 7, 1, 6, 2,4, 0, 5, 9,15, 8,12,10,13,14,11,39,34,36,35,33,32,38,37,
		52,54,50,51,53,48,49,55,34,36,39,35,38,32,37,33,23,19,20,16,17,21,18,22,46,44,40,45,43,41,42,47, 7, 1, 3, 2,4, 6, 0, 5,15, 8,9,10,13,12,11,14,57,58,61,60,63,56,62,59,30,28,31,25,29,24,26,27,
		52,54,55,49,51,50,53,48,28,31,30,25,26,24,27,29,36,39,34,35,37,32,33,38,44,40,46,43,41,45,42,47, 8,9,15,13,12,10,11,14,58,61,57,63,56,60,59,62, 1, 3, 7, 4, 6, 2,5, 0,19,20,23,16,18,21,22,17,
		52,55,49,54,51,53,48,50,30,28,31,25,27,26,24,29, 7, 1, 3, 5, 6, 0, 2,4,57,58,61,63,60,59,62,56,34,36,39,33,37,32,38,35,23,19,20,22,18,21,17,16,15, 8,9,11,12,13,14,10,46,44,40,42,47,41,43,45,
		52,55,50,48,54,49,51,53,44,40,46,42,43,41,45,47,28,31,30,25,24,26,29,27,58,61,57,60,59,63,62,56,19,20,23,18,21,22,17,16, 8,9,15,12,13,11,10,14,36,39,34,37,32,33,35,38,1, 3, 7, 5, 2,0, 4, 6,
		52,55,53,51,48,50,54,49, 3, 7, 1, 5, 4, 0, 6, 2,40,46,44,42,45,41,47,43,61,57,58,59,63,60,62,56, 9,15, 8,13,11,12,10,14,39,34,36,32,33,37,38,35,20,23,19,21,22,18,16,17,31,30,28,25,29,26,27,24,
		53,48,49,55,52,54,50,51,38,32,37,33,34,36,39,35,19,23,16,20,21,17,22,18,41,47,45,40,42,46,43,44, 5, 3, 1, 0, 6, 4, 2,7,13,10,11,8,15,14, 9,12,61,62,57,56,59,60,58,63,26,27,29,24,31,25,30,28,
		53,48,51,50,55,49,52,54,27,29,26,24,30,25,28,31,32,37,38,33,39,36,35,34,47,45,41,42,46,40,43,44,10,11,13,15,14, 8,9,12,62,57,61,59,60,56,63,58,3, 1, 5, 6, 4, 0, 7, 2,23,16,19,20,22,17,18,21,
		53,48,54,52,50,51,55,49,16,19,23,20,18,17,21,22,29,26,27,24,28,25,31,30,45,41,47,46,40,42,43,44,57,61,62,60,56,59,63,58,1, 5, 3, 4, 0, 6, 2,7,11,13,10,14, 8,15,12, 9,37,38,32,33,35,36,34,39,
		53,49,50,52,51,54,48,55,62,57,61,58,56,63,59,60,47,45,41,43,46,42,44,40,10,11,13, 9, 8,14,12,15, 3, 1, 5, 0, 2,6, 4, 7,23,16,19,17,22,21,18,20,27,29,26,25,30,31,24,28,32,37,38,33,35,34,39,36,
		53,49,54,51,48,55,52,50,41,47,45,43,40,42,46,44,38,32,37,33,36,34,35,39,13,10,11,14, 9, 8,12,15,26,27,29,31,25,30,28,24, 5, 3, 1, 6, 0, 2,4, 7,19,23,16,21,17,22,20,18,61,62,57,58,60,63,56,59,
		53,49,55,48,52,50,51,54,37,38,32,33,39,34,36,35,57,61,62,58,59,63,60,56,11,13,10, 8,14, 9,12,15,16,19,23,22,21,17,18,20,29,26,27,30,31,25,28,24, 1, 5, 3, 2,6, 0, 7, 4,45,41,47,43,44,42,40,46,
		53,50,48,51,55,54,49,52,29,26,27,24,25,28,30,31,11,13,10,12,14, 8,15, 9,16,19,23,18,17,21,20,22, 1, 5, 3, 0, 4, 2,6, 7,45,41,47,42,44,46,40,43,37,38,32,34,39,35,33,36,57,61,62,58,60,56,59,63,
		53,50,52,49,51,48,55,54,61,62,57,58,59,56,63,60,26,27,29,24,30,28,31,25,19,23,16,17,21,18,20,22,41,47,45,44,46,42,40,43,38,32,37,39,35,34,36,33, 5, 3, 1, 4, 2,0, 7, 6,13,10,11,12,15, 8,9,14,
		53,50,54,55,49,52,51,48,10,11,13,12, 9, 8,14,15,62,57,61,58,63,56,60,59,23,16,19,21,18,17,20,22,32,37,38,35,34,39,36,33, 3, 1, 5, 2,0, 4, 6, 7,47,45,41,46,42,44,43,40,27,29,26,24,31,28,25,30,
		53,51,49,54,48,50,55,52,47,45,41,43,42,46,40,44,27,29,26,24,25,30,31,28,62,57,61,56,63,59,58,60,23,16,19,22,17,18,21,20,10,11,13,14,15, 9, 8,12,32,37,38,39,36,35,33,34, 3, 1, 5, 7, 0, 2,6, 4,
		53,51,50,48,55,52,54,49,26,27,29,24,28,30,25,31,5, 3, 1, 7, 4, 2,0, 6,61,62,57,59,56,63,58,60,38,32,37,35,39,36,34,33,19,23,16,18,22,17,21,20,13,10,11,9,14,15,12, 8,41,47,45,43,44,46,42,40,
		53,51,52,55,54,49,48,50, 1, 5, 3, 7, 6, 2,4, 0,45,41,47,43,40,46,44,42,57,61,62,63,59,56,58,60,11,13,10,15, 9,14, 8,12,37,38,32,36,35,39,34,33,16,19,23,17,18,22,20,21,29,26,27,24,31,30,28,25,
		53,52,48,54,50,49,51,55,19,23,16,20,17,21,18,22,61,62,57,58,56,59,60,63,38,32,37,34,36,39,33,35,13,10,11,15, 8,9,14,12,41,47,45,46,44,40,42,43,26,27,29,30,28,31,24,25, 5, 3, 1, 7, 0, 6, 4, 2,
		53,52,49,50,51,55,54,48,57,61,62,58,63,59,56,60, 1, 5, 3, 7, 2,6, 0, 4,37,38,32,39,34,36,33,35,29,26,27,31,30,28,25,24,11,13,10, 9,15, 8,14,12,45,41,47,40,46,44,43,42,16,19,23,20,22,21,17,18,
		53,52,55,51,54,48,50,49, 3, 1, 5, 7, 4, 6, 2,0,23,16,19,20,18,21,22,17,32,37,38,36,39,34,33,35,47,45,41,44,40,46,42,43,27,29,26,28,31,30,25,24,10,11,13, 8,9,15,12,14,62,57,61,58,60,59,63,56,
		53,54,51,49,48,52,50,55,45,41,47,43,46,40,42,44,16,19,23,20,17,18,22,21,1, 5, 3, 6, 2,4, 7, 0,37,38,32,35,36,34,39,33,57,61,62,56,60,63,59,58,29,26,27,28,25,31,24,30,11,13,10,12,15, 9,14, 8,
		53,54,52,48,50,55,49,51,23,16,19,20,21,18,17,22,10,11,13,12, 8,9,15,14, 3, 1, 5, 4, 6, 2,7, 0,27,29,26,31,28,25,30,24,32,37,38,34,35,36,39,33,62,57,61,63,56,60,58,59,47,45,41,43,44,40,46,42,
		53,54,55,50,49,51,48,52,13,10,11,12,14, 9, 8,15,41,47,45,43,42,40,44,46, 5, 3, 1, 2,4, 6, 7, 0,61,62,57,60,63,56,59,58,26,27,29,25,31,28,30,24,38,32,37,36,34,35,33,39,19,23,16,20,22,18,21,17,
		53,55,48,49,52,51,54,50,32,37,38,33,36,39,34,35, 3, 1, 5, 7, 6, 4, 0, 2,27,29,26,30,25,28,24,31,62,57,61,60,59,63,56,58,47,45,41,40,44,42,46,43,23,16,19,18,21,22,20,17,10,11,13,12,15,14, 8,9,
		53,55,50,54,49,48,52,51,11,13,10,12, 8,14, 9,15,37,38,32,33,34,39,35,36,29,26,27,25,28,30,24,31,45,41,47,44,42,40,46,43,16,19,23,21,22,18,17,20,57,61,62,59,63,60,58,56, 1, 5, 3, 7, 0, 4, 2,6,
		53,55,51,52,54,50,49,48,5, 3, 1, 7, 2,4, 6, 0,13,10,11,12, 9,14,15, 8,26,27,29,28,30,25,24,31,19,23,16,22,18,21,17,20,61,62,57,63,60,59,56,58,41,47,45,42,40,44,43,46,38,32,37,33,35,39,36,34,
		54,48,50,49,53,52,55,51,4, 6, 2,0, 3, 1, 5, 7,16,23,20,19,21,18,17,22,34,35,36,38,33,32,39,37,43,41,45,42,46,40,44,47,31,28,25,29,27,24,26,30,13,12,10,15,14, 8,11,9,63,56,60,59,61,58,62,57,
		54,48,51,55,49,50,53,52,56,60,63,59,62,58,57,61,6, 2,4, 0, 5, 1, 7, 3,35,36,34,33,32,38,39,37,28,25,31,27,24,29,26,30,12,10,13,14, 8,15, 9,11,41,45,43,46,40,42,47,44,23,20,16,19,17,18,22,21,
		54,48,52,53,55,51,49,50,20,16,23,19,22,18,21,17,60,63,56,59,57,58,61,62,36,34,35,32,38,33,39,37,10,13,12, 8,15,14, 9,11,45,43,41,40,42,46,44,47,25,31,28,24,29,27,30,26, 2,4, 6, 0, 7, 1, 3, 5,
		54,49,48,50,53,51,52,55, 6, 2,4, 0, 1, 5, 3, 7,41,45,43,47,46,40,42,44,56,60,63,62,58,57,59,61,12,10,13, 8,14, 9,15,11,35,36,34,38,37,33,32,39,23,20,16,22,21,17,19,18,28,25,31,30,27,24,29,26,
		54,49,51,53,52,55,50,48,43,41,45,47,44,40,46,42,31,28,25,30,26,24,27,29,63,56,60,57,62,58,59,61,16,23,20,17,22,21,18,19,13,12,10, 9, 8,14,15,11,34,35,36,33,38,37,39,32, 4, 6, 2,0, 7, 5, 1, 3,
		54,49,55,52,50,48,53,51,25,31,28,30,29,24,26,27, 2,4, 6, 0, 3, 5, 7, 1,60,63,56,58,57,62,59,61,36,34,35,37,33,38,32,39,20,16,23,21,17,22,18,19,10,13,12,14, 9, 8,11,15,45,43,41,47,42,40,44,46,
		54,50,49,48,53,55,51,52, 2,4, 6, 0, 5, 3, 1, 7,10,13,12,11,14, 9, 8,15,25,31,28,29,24,26,30,27,20,16,23,17,21,18,22,19,60,63,56,62,61,58,57,59,45,43,41,44,46,42,47,40,36,34,35,39,37,33,38,32,
		54,50,52,51,48,49,53,55,34,35,36,39,38,33,32,37, 4, 6, 2,0, 1, 3, 7, 5,31,28,25,24,26,29,30,27,63,56,60,61,58,62,57,59,43,41,45,46,42,44,40,47,16,23,20,21,18,17,19,22,13,12,10,11,8,9,15,14,
		54,50,55,53,51,52,48,49,12,10,13,11,15, 9,14, 8,35,36,34,39,32,33,37,38,28,25,31,26,29,24,30,27,41,45,43,42,44,46,40,47,23,20,16,18,17,21,22,19,56,60,63,58,62,61,59,57, 6, 2,4, 0, 7, 3, 5, 1,
		54,51,50,52,48,55,49,53,35,36,34,39,33,32,38,37,56,60,63,59,58,62,61,57,12,10,13,15, 9,14,11,8,23,20,16,17,18,22,21,19,28,25,31,24,27,26,29,30, 6, 2,4, 5, 1, 7, 0, 3,41,45,43,47,42,44,46,40,
		54,51,53,49,52,50,48,55,45,43,41,47,46,44,40,42,36,34,35,39,38,32,37,33,10,13,12, 9,14,15,11,8,25,31,28,27,26,24,29,30, 2,4, 6, 1, 7, 5, 3, 0,20,16,23,18,22,17,19,21,60,63,56,59,61,62,57,58,
		54,51,55,48,49,53,52,50,63,56,60,59,57,62,58,61,43,41,45,47,40,44,42,46,13,12,10,14,15, 9,11,8,4, 6, 2,7, 5, 1, 3, 0,16,23,20,22,17,18,21,19,31,28,25,26,24,27,30,29,34,35,36,39,37,32,33,38,
		54,52,49,55,50,51,48,53,31,28,25,30,24,26,29,27,34,35,36,39,33,38,37,32,43,41,45,44,40,46,47,42,13,12,10, 8,9,15,14,11,63,56,60,58,61,57,62,59, 4, 6, 2,1, 3, 7, 0, 5,16,23,20,19,17,22,21,18,
		54,52,51,50,48,53,55,49,36,34,35,39,32,38,33,37,20,16,23,19,18,22,17,21,45,43,41,46,44,40,47,42, 2,4, 6, 7, 1, 3, 5, 0,10,13,12,15, 8,9,14,11,60,63,56,57,58,61,59,62,25,31,28,30,27,26,24,29,
		54,52,53,48,55,49,50,51,23,20,16,19,21,22,18,17,28,25,31,30,29,26,27,24,41,45,43,40,46,44,47,42,56,60,63,61,57,58,62,59, 6, 2,4, 3, 7, 1, 5, 0,12,10,13, 9,15, 8,11,14,35,36,34,39,37,38,32,33,
		54,53,48,52,55,50,51,49,16,23,20,19,18,21,22,17,13,12,10,11,15,14, 8,9, 4, 6, 2,3, 1, 5, 0, 7,31,28,25,27,29,26,24,30,34,35,36,32,37,38,33,39,63,56,60,62,57,61,59,58,43,41,45,47,42,46,40,44,
		54,53,49,51,52,48,55,50,41,45,43,47,40,46,44,42,23,20,16,19,22,21,17,18,6, 2,4, 1, 5, 3, 0, 7,35,36,34,37,38,32,33,39,56,60,63,57,61,62,58,59,28,25,31,29,26,27,30,24,12,10,13,11,8,14, 9,15,
		54,53,50,55,51,49,52,48,10,13,12,11,9,14,15, 8,45,43,41,47,44,46,42,40, 2,4, 6, 5, 3, 1, 0, 7,60,63,56,61,62,57,58,59,25,31,28,26,27,29,24,30,36,34,35,38,32,37,39,33,20,16,23,19,17,21,18,22,
		54,55,48,51,49,52,50,53,60,63,56,59,58,57,62,61,25,31,28,30,24,29,27,26,20,16,23,22,18,21,19,17,45,43,41,42,40,44,46,47,36,34,35,33,37,32,38,39, 2,4, 6, 3, 5, 7, 0, 1,10,13,12,11,8,15,14, 9,
		54,55,52,49,50,53,51,48,28,25,31,30,26,29,24,27,12,10,13,11,9,15, 8,14,23,20,16,21,22,18,19,17, 6, 2,4, 7, 3, 5, 1, 0,41,45,43,44,42,40,46,47,35,36,34,32,33,37,39,38,56,60,63,59,61,57,58,62,
		54,55,53,50,51,48,49,52,13,12,10,11,14,15, 9, 8,63,56,60,59,62,57,61,58,16,23,20,18,21,22,19,17,34,35,36,37,32,33,38,39, 4, 6, 2,5, 7, 3, 1, 0,43,41,45,40,44,42,47,46,31,28,25,30,27,29,26,24,
		55,48,50,52,49,53,51,54,42,40,44,46,47,41,43,45,33,37,32,38,39,35,34,36,11,12,13, 8,15,14,10, 9,28,30,25,24,29,27,26,31,3, 5, 7, 0, 6, 4, 2,1,22,18,21,16,20,19,17,23,63,60,59,56,62,61,58,57,
		55,48,53,49,51,54,52,50,32,33,37,38,36,35,39,34,59,63,60,56,57,61,62,58,13,11,12,14, 8,15,10, 9,21,22,18,19,16,20,23,17,25,28,30,27,24,29,26,31,7, 3, 5, 4, 0, 6, 1, 2,44,42,40,46,45,41,47,43,
		55,48,54,51,52,50,49,53,60,59,63,56,58,61,57,62,40,44,42,46,43,41,45,47,12,13,11,15,14, 8,10, 9, 5, 7, 3, 6, 4, 0, 2,1,18,21,22,20,19,16,23,17,30,25,28,29,27,24,31,26,37,32,33,38,34,35,36,39,
		55,49,48,53,51,50,54,52,33,37,32,38,35,39,36,34,22,18,21,17,16,20,19,23,42,40,44,47,41,43,46,45, 3, 5, 7, 6, 0, 2,4, 1,11,12,13,14, 9, 8,15,10,63,60,59,58,57,62,56,61,28,30,25,31,24,29,27,26,
		55,49,50,51,54,52,53,48,21,22,18,17,23,20,16,19,25,28,30,31,26,29,24,27,44,42,40,43,47,41,46,45,59,63,60,62,58,57,61,56, 7, 3, 5, 2,6, 0, 4, 1,13,11,12, 8,14, 9,10,15,32,33,37,38,34,39,35,36,
		55,49,52,54,53,48,51,50,30,25,28,31,27,29,26,24,37,32,33,38,36,39,34,35,40,44,42,41,43,47,46,45,12,13,11,9, 8,14,15,10,60,59,63,57,62,58,61,56, 5, 7, 3, 0, 2,6, 1, 4,18,21,22,17,19,20,23,16,
		55,50,51,49,54,53,48,52,18,21,22,17,16,23,20,19,12,13,11,10,14,15, 9, 8,5, 7, 3, 2,0, 4, 1, 6,30,25,28,24,26,29,27,31,37,32,33,35,34,39,36,38,60,59,63,61,58,62,56,57,40,44,42,46,45,47,43,41,
		55,50,52,48,49,51,54,53,44,42,40,46,43,47,41,45,21,22,18,17,20,23,19,16, 7, 3, 5, 0, 4, 2,1, 6,32,33,37,34,39,35,36,38,59,63,60,58,62,61,57,56,25,28,30,26,29,24,31,27,13,11,12,10, 9,15, 8,14,
		55,50,53,54,48,52,49,51,11,12,13,10, 8,15,14, 9,42,40,44,46,41,47,45,43, 3, 5, 7, 4, 2,0, 1, 6,63,60,59,62,61,58,57,56,28,30,25,29,24,26,27,31,33,37,32,39,35,34,38,36,22,18,21,17,19,23,16,20,
		55,51,48,54,52,53,50,49,59,63,60,56,61,57,58,62, 7, 3, 5, 1, 4, 0, 6, 2,32,33,37,36,35,39,38,34,25,28,30,24,27,26,29,31,13,11,12,15, 9,14, 8,10,44,42,40,47,43,45,46,41,21,22,18,17,19,16,20,23,
		55,51,49,50,54,48,52,53,22,18,21,17,20,16,23,19,63,60,59,56,58,57,62,61,33,37,32,35,39,36,38,34,11,12,13, 9,14,15, 8,10,42,40,44,43,45,47,41,46,28,30,25,27,26,24,31,29, 3, 5, 7, 1, 6, 0, 2,4,
		55,51,53,52,50,49,54,48,5, 7, 3, 1, 2,0, 4, 6,18,21,22,17,23,16,19,20,37,32,33,39,36,35,38,34,40,44,42,45,47,43,41,46,30,25,28,26,24,27,29,31,12,13,11,14,15, 9,10, 8,60,59,63,56,62,57,61,58,
		55,52,48,50,49,54,53,51,40,44,42,46,41,43,47,45,30,25,28,31,29,27,24,26,60,59,63,58,61,57,56,62,18,21,22,19,20,23,16,17,12,13,11,8,9,15,14,10,37,32,33,36,39,34,38,35, 5, 7, 3, 1, 6, 4, 0, 2,
		55,52,51,53,50,48,49,54, 7, 3, 5, 1, 0, 4, 2,6,44,42,40,46,47,43,45,41,59,63,60,61,57,58,56,62,13,11,12, 9,15, 8,14,10,32,33,37,39,34,36,35,38,21,22,18,20,23,19,17,16,25,28,30,31,24,27,26,29,
		55,52,54,49,53,51,50,48,28,30,25,31,26,27,29,24, 3, 5, 7, 1, 2,4, 6, 0,63,60,59,57,58,61,56,62,33,37,32,34,36,39,35,38,22,18,21,23,19,20,16,17,11,12,13,15, 8,9,10,14,42,40,44,46,45,43,41,47,
		55,53,49,48,51,52,50,54,37,32,33,38,39,36,35,34, 5, 7, 3, 1, 0, 2,6, 4,30,25,28,27,29,26,31,24,60,59,63,62,57,61,58,56,40,44,42,47,45,41,43,46,18,21,22,23,16,19,17,20,12,13,11,10, 9, 8,14,15,
		55,53,52,51,50,54,48,49, 3, 5, 7, 1, 4, 2,0, 6,11,12,13,10,15, 8,9,14,28,30,25,26,27,29,31,24,22,18,21,19,23,16,20,17,63,60,59,61,62,57,58,56,42,40,44,41,47,45,46,43,33,37,32,38,34,36,39,35,
		55,53,54,50,48,49,51,52,13,11,12,10,14, 8,15, 9,32,33,37,38,35,36,34,39,25,28,30,29,26,27,31,24,44,42,40,45,41,47,43,46,21,22,18,16,19,23,20,17,59,63,60,57,61,62,56,58,7, 3, 5, 1, 6, 2,4, 0,
		55,54,49,52,53,50,48,51,25,28,30,31,29,26,27,24,13,11,12,10, 8,14, 9,15,21,22,18,23,20,16,17,19, 7, 3, 5, 6, 2,4, 0, 1,44,42,40,41,45,43,47,46,32,33,37,35,36,34,38,39,59,63,60,56,62,58,57,61,
		55,54,50,53,48,51,52,49,12,13,11,10,15,14, 8,9,60,59,63,56,61,58,62,57,18,21,22,16,23,20,17,19,37,32,33,34,35,36,39,38,5, 7, 3, 4, 6, 2,0, 1,40,44,42,43,41,45,46,47,30,25,28,31,24,26,29,27,
		55,54,51,48,52,49,53,50,63,60,59,56,57,58,61,62,28,30,25,31,27,26,24,29,22,18,21,20,16,23,17,19,42,40,44,45,43,41,47,46,33,37,32,36,34,35,39,38,3, 5, 7, 2,4, 6, 1, 0,11,12,13,10, 9,14,15, 8,
		56,57,58,59,63,62,61,60, 0, 7, 6, 1, 3, 4, 5, 2,16,22,23,17,21,19,18,20, 8,14,12,15, 9,11,10,13,48,51,54,52,49,53,50,55,32,39,35,34,33,37,36,38,24,26,28,30,25,29,31,27,40,45,41,46,47,44,42,43,
		56,57,60,61,59,58,63,62,45,41,40,46,42,44,43,47, 7, 6, 0, 1, 5, 4, 2,3,14,12, 8,9,11,15,10,13,39,35,32,33,37,34,36,38,26,28,24,25,29,30,27,31,51,54,48,49,53,52,55,50,22,23,16,17,18,19,20,21,
		56,57,62,63,61,60,59,58,23,16,22,17,20,19,21,18,41,40,45,46,43,44,47,42,12, 8,14,11,15, 9,10,13,28,24,26,29,30,25,27,31,54,48,51,53,52,49,50,55,35,32,39,37,34,33,38,36, 6, 0, 7, 1, 2,4, 3, 5,
		56,58,59,57,63,61,60,62, 6, 0, 7, 1, 5, 3, 4, 2,28,24,26,31,25,27,29,30,35,32,39,34,37,36,38,33,23,16,22,18,21,19,20,17,41,40,45,42,47,44,43,46,54,48,51,50,49,52,55,53,12, 8,14,10,13, 9,15,11,
		56,58,61,63,60,62,57,59,26,28,24,31,30,27,25,29,14,12, 8,10,11,9,13,15,39,35,32,36,34,37,38,33,51,54,48,52,50,49,53,55,22,23,16,19,18,21,20,17,45,41,40,44,42,47,46,43, 7, 6, 0, 1, 2,3, 5, 4,
		56,58,62,60,57,59,63,61,8,14,12,10,15, 9,11,13, 0, 7, 6, 1, 4, 3, 2,5,32,39,35,37,36,34,38,33,40,45,41,47,44,42,43,46,48,51,54,49,52,50,53,55,16,22,23,21,19,18,17,20,24,26,28,31,29,27,30,25,
		56,59,57,58,63,60,62,61,7, 6, 0, 1, 4, 5, 3, 2,51,54,48,55,49,53,52,50,45,41,40,42,44,43,46,47,26,28,24,29,25,27,30,31,14,12, 8,15,13, 9,11,10,22,23,16,20,21,18,17,19,39,35,32,38,33,37,34,36,
		56,59,60,63,62,61,58,57,48,51,54,55,50,53,49,52,32,39,35,38,36,37,33,34,40,45,41,43,42,44,46,47,16,22,23,18,20,21,19,17,24,26,28,27,29,25,30,31,8,14,12, 9,15,13,10,11,0, 7, 6, 1, 2,5, 4, 3,
		56,59,61,62,58,57,63,60,35,32,39,38,34,37,36,33, 6, 0, 7, 1, 3, 5, 2,4,41,40,45,44,43,42,46,47,12, 8,14,13, 9,15,11,10,23,16,22,21,18,20,19,17,28,24,26,25,27,29,31,30,54,48,51,55,52,53,50,49,
		56,60,58,62,57,61,59,63,14,12, 8,10, 9,11,15,13,45,41,40,46,44,42,47,43,26,28,24,30,27,25,31,29,22,23,16,18,19,20,21,17,39,35,32,37,33,36,34,38,7, 6, 0, 5, 4, 2,1, 3,51,54,48,55,52,50,49,53,
		56,60,61,57,59,63,62,58,40,45,41,46,43,42,44,47,48,51,54,55,53,50,52,49,24,26,28,25,30,27,31,29, 0, 7, 6, 2,5, 4, 3, 1,16,22,23,20,18,19,21,17,32,39,35,36,37,33,38,34, 8,14,12,10,13,11,9,15,
		56,60,63,59,62,58,57,61,54,48,51,55,49,50,53,52,12, 8,14,10,15,11,13, 9,28,24,26,27,25,30,31,29,35,32,39,33,36,37,34,38,6, 0, 7, 4, 2,5, 3, 1,23,16,22,19,20,18,17,21,41,40,45,46,47,42,43,44,
		56,61,57,60,59,62,58,63,41,40,45,46,44,43,42,47,35,32,39,38,37,34,33,36,23,16,22,20,19,21,17,18,54,48,51,52,53,50,49,55,12, 8,14, 9,13,11,15,10, 6, 0, 7, 3, 5, 2,1, 4,28,24,26,31,29,30,25,27,
		56,61,62,59,58,63,60,57,39,35,32,38,36,34,37,33,26,28,24,31,27,30,29,25,22,23,16,21,20,19,17,18,7, 6, 0, 2,3, 5, 4, 1,51,54,48,50,52,53,49,55,14,12, 8,11,9,13,10,15,45,41,40,46,47,43,44,42,
		56,61,63,58,60,57,59,62,24,26,28,31,25,30,27,29,40,45,41,46,42,43,47,44,16,22,23,19,21,20,17,18,8,14,12,13,11,9,15,10, 0, 7, 6, 5, 2,3, 4, 1,48,51,54,53,50,52,55,49,32,39,35,38,33,34,36,37,
		56,62,59,61,58,60,57,63,32,39,35,38,37,36,34,33, 8,14,12,10, 9,15,13,11,48,51,54,50,53,49,55,52,24,26,28,29,27,30,25,31,40,45,41,44,47,43,42,46, 0, 7, 6, 4, 3, 2,1, 5,16,22,23,17,18,20,21,19,
		56,62,60,58,57,63,61,59,12, 8,14,10,11,15, 9,13,23,16,22,17,19,20,18,21,54,48,51,49,50,53,55,52, 6, 0, 7, 2,4, 3, 5, 1,28,24,26,30,29,27,25,31,41,40,45,43,44,47,46,42,35,32,39,38,33,36,37,34,
		56,62,63,57,61,59,58,60,22,23,16,17,21,20,19,18,39,35,32,38,34,36,33,37,51,54,48,53,49,50,55,52,45,41,40,47,43,44,42,46, 7, 6, 0, 3, 2,4, 5, 1,26,28,24,27,30,29,31,25,14,12, 8,10,13,15,11,9,
		56,63,57,62,61,58,60,59,16,22,23,17,19,21,20,18,24,26,28,31,30,25,29,27, 0, 7, 6, 3, 4, 5, 1, 2,32,39,35,33,34,36,37,38,8,14,12,11,13,15, 9,10,40,45,41,42,43,47,46,44,48,51,54,55,52,49,53,50,
		56,63,58,61,60,59,62,57,28,24,26,31,27,25,30,29,54,48,51,55,50,49,52,53, 6, 0, 7, 5, 3, 4, 1, 2,41,40,45,47,42,43,44,46,35,32,39,36,33,34,37,38,12, 8,14,15,11,13,10, 9,23,16,22,17,18,21,19,20,
		56,63,59,60,62,57,61,58,51,54,48,55,53,49,50,52,22,23,16,17,20,21,18,19, 7, 6, 0, 4, 5, 3, 1, 2,14,12, 8,13,15,11,9,10,45,41,40,43,47,42,44,46,39,35,32,34,36,33,38,37,26,28,24,31,29,25,27,30,
		57,56,59,58,62,63,60,61,7, 0, 1, 6, 4, 3, 2,5,23,17,16,22,18,20,21,19,15, 9,11,8,14,12,13,10,52,49,53,48,51,54,55,50,34,33,37,32,39,35,38,36,30,25,29,24,26,28,27,31,46,41,45,40,43,42,44,47,
		57,56,61,60,58,59,62,63,41,45,46,40,44,42,47,43, 0, 1, 7, 6, 2,3, 5, 4, 9,11,15,14,12, 8,13,10,33,37,34,39,35,32,38,36,25,29,30,26,28,24,31,27,49,53,52,51,54,48,50,55,17,16,23,22,21,20,19,18,
		57,56,63,62,60,61,58,59,16,23,17,22,19,20,18,21,45,46,41,40,47,42,43,44,11,15, 9,12, 8,14,13,10,29,30,25,28,24,26,31,27,53,52,49,54,48,51,55,50,37,34,33,35,32,39,36,38,1, 7, 0, 6, 5, 3, 4, 2,
		57,58,56,59,62,61,63,60, 0, 1, 7, 6, 3, 2,4, 5,49,53,52,50,51,54,48,55,41,45,46,44,42,47,40,43,25,29,30,28,26,31,24,27, 9,11,15, 8,10,14,12,13,17,16,23,19,18,21,22,20,33,37,34,36,39,35,32,38,
		57,58,60,63,59,56,62,61,37,34,33,36,32,35,38,39, 1, 7, 0, 6, 4, 2,5, 3,45,46,41,42,47,44,40,43,11,15, 9,10,14, 8,12,13,16,23,17,18,21,19,20,22,29,30,25,26,31,28,27,24,53,52,49,50,48,54,55,51,
		57,58,61,62,63,60,59,56,52,49,53,50,55,54,51,48,34,33,37,36,38,35,39,32,46,41,45,47,44,42,40,43,23,17,16,21,19,18,20,22,30,25,29,31,28,26,24,27,15, 9,11,14, 8,10,13,12, 7, 0, 1, 6, 5, 2,3, 4,
		57,59,58,56,62,60,61,63, 1, 7, 0, 6, 2,4, 3, 5,29,30,25,27,26,31,28,24,37,34,33,32,35,38,36,39,16,23,17,21,18,20,19,22,45,46,41,44,43,42,47,40,53,52,49,55,51,48,50,54,11,15, 9,13,10,14, 8,12,
		57,59,60,62,61,63,56,58,25,29,30,27,24,31,26,28,9,11,15,13,12,14,10, 8,33,37,34,38,32,35,36,39,49,53,52,48,55,51,54,50,17,16,23,20,21,18,19,22,41,45,46,42,44,43,40,47, 0, 1, 7, 6, 5, 4, 2,3,
		57,59,63,61,56,58,62,60,15, 9,11,13, 8,14,12,10, 7, 0, 1, 6, 3, 4, 5, 2,34,33,37,35,38,32,36,39,46,41,45,43,42,44,47,40,52,49,53,51,48,55,54,50,23,17,16,18,20,21,22,19,30,25,29,27,28,31,24,26,
		57,60,56,61,58,63,59,62,45,46,41,40,42,47,44,43,37,34,33,36,35,32,39,38,16,23,17,19,20,18,22,21,53,52,49,48,54,55,51,50,11,15, 9,14,10,12, 8,13, 1, 7, 0, 4, 2,5, 6, 3,29,30,25,27,28,24,26,31,
		57,60,62,59,61,56,58,63,30,25,29,27,26,24,31,28,46,41,45,40,44,47,43,42,23,17,16,20,18,19,22,21,15, 9,11,10,12,14, 8,13, 7, 0, 1, 2,5, 4, 3, 6,52,49,53,54,55,48,50,51,34,33,37,36,39,32,38,35,
		57,60,63,58,59,62,61,56,33,37,34,36,38,32,35,39,25,29,30,27,31,24,28,26,17,16,23,18,19,20,22,21,0, 1, 7, 5, 4, 2,3, 6,49,53,52,55,48,54,51,50, 9,11,15,12,14,10,13, 8,41,45,46,40,43,47,42,44,
		57,61,59,63,56,60,58,62, 9,11,15,13,14,12, 8,10,41,45,46,40,42,44,43,47,25,29,30,24,31,26,27,28,17,16,23,21,20,19,18,22,33,37,34,35,39,38,32,36, 0, 1, 7, 2,3, 5, 6, 4,49,53,52,50,48,55,51,54,
		57,61,60,56,58,62,63,59,46,41,45,40,47,44,42,43,52,49,53,50,54,55,48,51,30,25,29,26,24,31,27,28,7, 0, 1, 5, 2,3, 4, 6,23,17,16,19,21,20,18,22,34,33,37,38,35,39,36,32,15, 9,11,13,10,12,14, 8,
		57,61,62,58,63,59,56,60,53,52,49,50,51,55,54,48,11,15, 9,13, 8,12,10,14,29,30,25,31,26,24,27,28,37,34,33,39,38,35,32,36, 1, 7, 0, 3, 5, 2,4, 6,16,23,17,20,19,21,22,18,45,46,41,40,43,44,47,42,
		57,62,56,63,60,59,61,58,23,17,16,22,20,18,19,21,30,25,29,27,24,26,28,31,7, 0, 1, 4, 3, 2,6, 5,34,33,37,39,32,38,35,36,15, 9,11,12,10, 8,14,13,46,41,45,44,47,43,40,42,52,49,53,50,48,51,54,55,
		57,62,58,61,63,56,60,59,49,53,52,50,54,51,55,48,17,16,23,22,19,18,21,20, 0, 1, 7, 3, 2,4, 6, 5, 9,11,15,10, 8,12,14,13,41,45,46,47,43,44,42,40,33,37,34,32,38,39,36,35,25,29,30,27,28,26,31,24,
		57,62,59,60,61,58,63,56,29,30,25,27,31,26,24,28,53,52,49,50,55,51,48,54, 1, 7, 0, 2,4, 3, 6, 5,45,46,41,43,44,47,42,40,37,34,33,38,39,32,35,36,11,15, 9, 8,12,10,13,14,16,23,17,22,21,18,20,19,
		57,63,58,60,59,61,56,62,34,33,37,36,35,38,32,39,15, 9,11,13,14, 8,10,12,52,49,53,55,54,51,50,48,30,25,29,28,31,24,26,27,46,41,45,42,43,47,44,40, 7, 0, 1, 3, 4, 5, 6, 2,23,17,16,22,21,19,18,20,
		57,63,61,59,56,62,60,58,11,15, 9,13,12, 8,14,10,16,23,17,22,20,19,21,18,53,52,49,51,55,54,50,48,1, 7, 0, 5, 3, 4, 2,6,29,30,25,24,28,31,26,27,45,46,41,47,42,43,40,44,37,34,33,36,39,38,35,32,
		57,63,62,56,60,58,59,61,17,16,23,22,18,19,20,21,33,37,34,36,32,38,39,35,49,53,52,54,51,55,50,48,41,45,46,43,47,42,44,40, 0, 1, 7, 4, 5, 3, 2,6,25,29,30,31,24,28,27,26, 9,11,15,13,10, 8,12,14,
		58,56,57,59,61,63,62,60, 0, 6, 1, 7, 3, 5, 2,4,26,31,28,24,29,30,25,27,34,37,36,35,32,39,33,38,18,21,19,23,16,22,17,20,42,47,44,41,40,45,46,43,50,49,52,54,48,51,53,55,10,14, 8,12,11,15, 9,13,
		58,56,60,62,59,57,61,63,14, 8,10,12, 9,15,13,11,6, 1, 0, 7, 2,5, 4, 3,37,36,34,32,39,35,33,38,47,44,42,40,45,41,46,43,49,52,50,48,51,54,55,53,21,19,18,16,22,23,20,17,31,28,26,24,25,30,27,29,
		58,56,63,61,62,60,59,57,28,26,31,24,27,30,29,25, 8,10,14,12,13,15,11,9,36,34,37,39,35,32,33,38,52,50,49,51,54,48,55,53,19,18,21,22,23,16,17,20,44,42,47,45,41,40,43,46, 1, 0, 6, 7, 4, 5, 3, 2,
		58,57,59,56,61,62,60,63, 1, 0, 6, 7, 2,3, 5, 4,52,50,49,53,48,55,51,54,44,42,47,41,45,46,43,40,28,26,31,25,29,30,27,24, 8,10,14, 9,11,15,13,12,19,18,21,17,16,23,20,22,36,34,37,33,38,32,35,39,
		58,57,62,61,60,63,56,59,49,52,50,53,54,55,48,51,37,36,34,33,39,32,38,35,47,44,42,46,41,45,43,40,21,19,18,23,17,16,22,20,31,28,26,30,25,29,27,24,14, 8,10,15, 9,11,12,13, 6, 1, 0, 7, 4, 3, 2,5,
		58,57,63,60,56,59,61,62,34,37,36,33,35,32,39,38,0, 6, 1, 7, 5, 3, 4, 2,42,47,44,45,46,41,43,40,10,14, 8,11,15, 9,13,12,18,21,19,16,23,17,22,20,26,31,28,29,30,25,24,27,50,49,52,53,51,55,54,48,
		58,59,56,57,61,60,63,62, 6, 1, 0, 7, 5, 2,3, 4,21,19,18,20,16,22,23,17,14, 8,10, 9,15,13,12,11,49,52,50,51,48,55,54,53,37,36,34,35,38,32,39,33,31,28,26,27,29,25,24,30,47,44,42,43,40,45,41,46,
		58,59,60,61,63,62,57,56,18,21,19,20,17,22,16,23,42,47,44,43,46,45,40,41,10,14, 8,13, 9,15,12,11,26,31,28,25,27,29,30,24,50,49,52,55,51,48,54,53,34,37,36,32,35,38,33,39, 0, 6, 1, 7, 4, 2,5, 3,
		58,59,62,63,57,56,61,60,44,42,47,43,41,45,46,40, 1, 0, 6, 7, 3, 2,4, 5, 8,10,14,15,13, 9,12,11,36,34,37,38,32,35,39,33,28,26,31,29,25,27,30,24,52,50,49,48,55,51,53,54,19,18,21,20,23,22,17,16,
		58,60,57,63,56,62,59,61,37,36,34,33,32,39,35,38,14, 8,10,12,15, 9,11,13,49,52,50,54,55,48,53,51,31,28,26,25,30,27,29,24,47,44,42,45,40,46,41,43, 6, 1, 0, 2,5, 4, 7, 3,21,19,18,20,23,17,16,22,
		58,60,61,59,63,57,56,62,19,18,21,20,16,17,22,23,36,34,37,33,35,39,38,32,52,50,49,55,48,54,53,51,44,42,47,40,46,45,41,43, 1, 0, 6, 5, 4, 2,3, 7,28,26,31,30,27,25,24,29, 8,10,14,12,11,9,13,15,
		58,60,62,56,59,61,63,57,10,14, 8,12,13, 9,15,11,18,21,19,20,22,17,23,16,50,49,52,48,54,55,53,51,0, 6, 1, 4, 2,5, 3, 7,26,31,28,27,25,30,29,24,42,47,44,46,45,40,43,41,34,37,36,33,38,39,32,35,
		58,61,56,63,62,57,60,59,26,31,28,24,30,29,27,25,50,49,52,53,54,48,51,55, 0, 6, 1, 3, 5, 2,7, 4,42,47,44,40,41,46,45,43,34,37,36,39,38,35,32,33,10,14, 8,9,13,11,12,15,18,21,19,20,23,16,22,17,
		58,61,57,62,60,59,63,56,52,50,49,53,55,48,54,51,19,18,21,20,17,16,23,22, 1, 0, 6, 2,3, 5, 7, 4, 8,10,14,11,9,13,15,12,44,42,47,46,40,41,45,43,36,34,37,35,39,38,33,32,28,26,31,24,25,29,30,27,
		58,61,59,60,63,56,62,57,21,19,18,20,22,16,17,23,31,28,26,24,27,29,25,30, 6, 1, 0, 5, 2,3, 7, 4,37,36,34,38,35,39,32,33,14, 8,10,13,11,9,15,12,47,44,42,41,46,40,43,45,49,52,50,53,51,48,55,54,
		58,62,56,60,59,63,57,61,8,10,14,12,15,13, 9,11,44,42,47,43,45,41,40,46,28,26,31,27,30,29,24,25,19,18,21,23,22,17,16,20,36,34,37,32,38,39,35,33, 1, 0, 6, 3, 2,4, 7, 5,52,50,49,53,51,54,48,55,
		58,62,61,57,60,56,59,63,50,49,52,53,48,54,55,51,10,14, 8,12, 9,13,11,15,26,31,28,30,29,27,24,25,34,37,36,38,39,32,35,33, 0, 6, 1, 2,4, 3, 5, 7,18,21,19,22,17,23,20,16,42,47,44,43,40,41,46,45,
		58,62,63,59,57,61,60,56,47,44,42,43,46,41,45,40,49,52,50,53,55,54,51,48,31,28,26,29,27,30,24,25, 6, 1, 0, 4, 3, 2,5, 7,21,19,18,17,23,22,16,20,37,36,34,39,32,38,33,35,14, 8,10,12,11,13,15, 9,
		58,63,59,62,57,60,56,61,42,47,44,43,45,46,41,40,34,37,36,33,32,35,38,39,18,21,19,17,22,16,20,23,50,49,52,51,55,54,48,53,10,14, 8,15,11,13, 9,12, 0, 6, 1, 5, 3, 4, 7, 2,26,31,28,24,25,27,29,30,
		58,63,60,57,56,61,62,59,36,34,37,33,39,35,32,38,28,26,31,24,30,27,25,29,19,18,21,16,17,22,20,23, 1, 0, 6, 4, 5, 3, 2,7,52,50,49,54,51,55,48,53, 8,10,14,13,15,11,12, 9,44,42,47,43,40,46,45,41,
		58,63,61,56,62,59,57,60,31,28,26,24,29,27,30,25,47,44,42,43,41,46,40,45,21,19,18,22,16,17,20,23,14, 8,10,11,13,15, 9,12, 6, 1, 0, 3, 4, 5, 2,7,49,52,50,55,54,51,53,48,37,36,34,33,38,35,39,32,
		59,56,58,57,60,63,61,62, 6, 7, 1, 0, 5, 4, 2,3,48,55,51,54,52,50,49,53,42,44,43,45,41,40,47,46,29,25,27,26,28,24,31,30,15,13, 9,14,12, 8,10,11,20,21,18,22,23,16,19,17,38,32,35,39,36,34,37,33,
		59,56,62,61,57,58,60,63,32,35,38,39,37,34,33,36, 7, 1, 6, 0, 2,4, 3, 5,44,43,42,41,40,45,47,46,13, 9,15,12, 8,14,10,11,21,18,20,23,16,22,17,19,25,27,29,28,24,26,30,31,55,51,48,54,49,50,53,52,
		59,56,63,60,61,62,57,58,51,48,55,54,53,50,52,49,35,38,32,39,33,34,36,37,43,42,44,40,45,41,47,46,18,20,21,16,22,23,17,19,27,29,25,24,26,28,31,30, 9,15,13, 8,14,12,11,10, 1, 6, 7, 0, 3, 4, 5, 2,
		59,57,56,58,60,62,63,61,7, 1, 6, 0, 4, 2,5, 3,25,27,29,30,28,24,26,31,32,35,38,37,34,33,39,36,21,18,20,16,23,17,22,19,44,43,42,45,46,41,40,47,55,51,48,53,52,49,54,50,13, 9,15,11,12, 8,14,10,
		59,57,61,63,58,56,60,62, 9,15,13,11,14, 8,10,12, 1, 6, 7, 0, 5, 2,3, 4,35,38,32,34,33,37,39,36,43,42,44,46,41,45,40,47,51,48,55,52,49,53,50,54,18,20,21,23,17,16,19,22,27,29,25,30,26,24,31,28,
		59,57,62,60,63,61,58,56,29,25,27,30,31,24,28,26,15,13, 9,11,10, 8,12,14,38,32,35,33,37,34,39,36,48,55,51,49,53,52,50,54,20,21,18,17,16,23,22,19,42,44,43,41,45,46,47,40, 6, 7, 1, 0, 3, 2,4, 5,
		59,58,57,56,60,61,62,63, 1, 6, 7, 0, 2,5, 4, 3,18,20,21,19,23,17,16,22, 9,15,13,14, 8,10,11,12,51,48,55,49,52,50,53,54,35,38,32,37,36,34,33,39,27,29,25,31,28,26,30,24,43,42,44,47,46,41,45,40,
		59,58,61,60,62,63,56,57,21,18,20,19,22,17,23,16,44,43,42,47,40,41,46,45,13, 9,15,10,14, 8,11,12,25,27,29,26,31,28,24,30,55,51,48,50,49,52,53,54,32,35,38,34,37,36,39,33, 7, 1, 6, 0, 3, 5, 2,4,
		59,58,63,62,56,57,60,61,42,44,43,47,45,41,40,46, 6, 7, 1, 0, 4, 5, 3, 2,15,13, 9, 8,10,14,11,12,38,32,35,36,34,37,33,39,29,25,27,28,26,31,24,30,48,55,51,52,50,49,54,53,20,21,18,19,16,17,22,23,
		59,60,56,63,61,58,62,57,48,55,51,54,50,52,53,49,20,21,18,19,22,23,16,17, 6, 7, 1, 5, 4, 2,0, 3,15,13, 9,12,14,10, 8,11,42,44,43,40,46,45,41,47,38,32,35,37,33,36,39,34,29,25,27,30,26,28,24,31,
		59,60,57,62,63,56,61,58,25,27,29,30,24,28,31,26,55,51,48,54,53,52,49,50, 7, 1, 6, 4, 2,5, 0, 3,44,43,42,46,45,40,41,47,32,35,38,33,36,37,34,39,13, 9,15,14,10,12,11,8,21,18,20,19,16,23,17,22,
		59,60,58,61,62,57,63,56,18,20,21,19,17,23,22,16,27,29,25,30,31,28,26,24, 1, 6, 7, 2,5, 4, 0, 3,35,38,32,36,37,33,34,39, 9,15,13,10,12,14, 8,11,43,42,44,45,40,46,47,41,51,48,55,54,49,52,50,53,
		59,61,56,62,57,63,58,60,35,38,32,39,34,33,37,36, 9,15,13,11,8,14,12,10,51,48,55,53,50,52,54,49,27,29,25,26,24,31,28,30,43,42,44,41,46,40,45,47, 1, 6, 7, 5, 2,3, 0, 4,18,20,21,19,16,22,23,17,
		59,61,60,58,62,56,57,63,20,21,18,19,23,22,17,16,38,32,35,39,37,33,36,34,48,55,51,50,52,53,54,49,42,44,43,46,40,41,45,47, 6, 7, 1, 2,3, 5, 4, 0,29,25,27,24,31,26,30,28,15,13, 9,11,12,14,10, 8,
		59,61,63,57,58,60,62,56,13, 9,15,11,10,14, 8,12,21,18,20,19,17,22,16,23,55,51,48,52,53,50,54,49, 7, 1, 6, 3, 5, 2,4, 0,25,27,29,31,26,24,28,30,44,43,42,40,41,46,47,45,32,35,38,39,36,33,34,37,
		59,62,58,63,56,61,57,60,44,43,42,47,41,40,45,46,32,35,38,39,34,37,36,33,21,18,20,22,17,23,19,16,55,51,48,49,50,53,52,54,13, 9,15, 8,12,10,14,11,7, 1, 6, 2,4, 3, 0, 5,25,27,29,30,26,31,28,24,
		59,62,60,57,63,58,56,61,27,29,25,30,28,31,24,26,43,42,44,47,45,40,46,41,18,20,21,17,23,22,19,16, 9,15,13,12,10, 8,14,11,1, 6, 7, 4, 3, 2,5, 0,51,48,55,50,53,49,54,52,35,38,32,39,36,37,33,34,
		59,62,61,56,57,60,63,58,38,32,35,39,33,37,34,36,29,25,27,30,24,31,26,28,20,21,18,23,22,17,19,16, 6, 7, 1, 3, 2,4, 5, 0,48,55,51,53,49,50,52,54,15,13, 9,10, 8,12,11,14,42,44,43,47,46,40,41,45,
		59,63,57,61,58,62,56,60,15,13, 9,11,8,10,14,12,42,44,43,47,41,45,46,40,29,25,27,31,24,28,30,26,20,21,18,16,17,22,23,19,38,32,35,34,36,33,37,39, 6, 7, 1, 4, 5, 3, 0, 2,48,55,51,54,49,53,52,50,
		59,63,60,56,61,57,58,62,55,51,48,54,52,53,50,49,13, 9,15,11,14,10,12, 8,25,27,29,24,28,31,30,26,32,35,38,36,33,34,37,39, 7, 1, 6, 5, 3, 4, 2,0,21,18,20,17,22,16,19,23,44,43,42,47,46,45,40,41,
		59,63,62,58,56,60,61,57,43,42,44,47,40,45,41,46,51,48,55,54,50,53,49,52,27,29,25,28,31,24,30,26, 1, 6, 7, 3, 4, 5, 2,0,18,20,21,22,16,17,23,19,35,38,32,33,34,36,39,37, 9,15,13,11,12,10, 8,14,
		60,56,57,61,63,59,58,62,45,40,46,41,42,43,47,44,54,55,48,51,52,49,53,50,25,30,27,24,26,28,29,31,2,5, 4, 0, 7, 6, 1, 3,20,18,19,16,22,23,17,21,36,37,33,32,39,35,34,38,10,12,14, 8,15, 9,11,13,
		60,56,59,63,58,62,61,57,48,54,55,51,50,49,52,53,14,10,12, 8,13, 9,15,11,27,25,30,28,24,26,29,31,33,36,37,35,32,39,38,34, 4, 2,5, 6, 0, 7, 1, 3,19,20,18,23,16,22,21,17,46,45,40,41,44,43,42,47,
		60,56,62,58,61,57,63,59,12,14,10, 8,11,9,13,15,40,46,45,41,47,43,44,42,30,27,25,26,28,24,29,31,18,19,20,22,23,16,17,21,37,33,36,39,35,32,38,34, 5, 4, 2,7, 6, 0, 3, 1,55,48,54,51,53,49,50,52,
		60,57,58,63,62,59,56,61,37,33,36,34,32,38,39,35,30,27,25,29,28,26,31,24,18,19,20,17,16,23,21,22, 5, 4, 2,0, 1, 7, 6, 3,55,48,54,49,53,52,50,51,12,14,10, 9,11,15, 8,13,40,46,45,41,44,42,47,43,
		60,57,59,62,56,61,63,58,25,30,27,29,24,26,28,31,45,40,46,41,43,42,44,47,20,18,19,23,17,16,21,22,10,12,14,15, 9,11,13, 8,2,5, 4, 7, 0, 1, 6, 3,54,55,48,52,49,53,51,50,36,37,33,34,35,38,32,39,
		60,57,61,56,63,58,62,59,46,45,40,41,47,42,43,44,33,36,37,34,39,38,35,32,19,20,18,16,23,17,21,22,48,54,55,53,52,49,50,51,14,10,12,11,15, 9,13, 8,4, 2,5, 1, 7, 0, 3, 6,27,25,30,29,31,26,24,28,
		60,58,56,62,61,59,57,63,14,10,12, 8,9,13,11,15,19,20,18,21,23,16,22,17,48,54,55,50,49,52,51,53, 4, 2,5, 0, 6, 1, 7, 3,27,25,30,26,31,28,24,29,46,45,40,42,47,44,41,43,33,36,37,34,35,32,39,38,
		60,58,59,61,57,63,62,56,18,19,20,21,17,16,23,22,37,33,36,34,38,32,35,39,55,48,54,52,50,49,51,53,40,46,45,44,42,47,43,41,5, 4, 2,1, 0, 6, 7, 3,30,27,25,28,26,31,29,24,12,14,10, 8,15,13, 9,11,
		60,58,63,57,62,56,61,59,36,37,33,34,39,32,38,35,10,12,14, 8,11,13,15, 9,54,55,48,49,52,50,51,53,25,30,27,31,28,26,24,29,45,40,46,47,44,42,43,41,2,5, 4, 6, 1, 0, 3, 7,20,18,19,21,22,16,17,23,
		60,59,61,58,57,62,56,63,20,18,19,21,23,17,16,22,25,30,27,29,26,24,31,28,2,5, 4, 1, 6, 7, 3, 0,36,37,33,35,38,32,39,34,10,12,14, 9,15,13,11,8,45,40,46,43,42,44,41,47,54,55,48,51,53,50,52,49,
		60,59,62,57,56,63,58,61,27,25,30,29,28,24,26,31,48,54,55,51,49,50,53,52, 4, 2,5, 7, 1, 6, 3, 0,46,45,40,44,43,42,47,41,33,36,37,32,35,38,39,34,14,10,12,13, 9,15, 8,11,19,20,18,21,22,17,23,16,
		60,59,63,56,58,61,57,62,55,48,54,51,52,50,49,53,18,19,20,21,16,17,22,23, 5, 4, 2,6, 7, 1, 3, 0,12,14,10,15,13, 9,11,8,40,46,45,42,44,43,47,41,37,33,36,38,32,35,34,39,30,27,25,29,31,24,28,26,
		60,61,56,57,63,62,59,58,40,46,45,41,43,47,42,44, 5, 4, 2,3, 7, 6, 0, 1,12,14,10,11,9,13, 8,15,37,33,36,35,39,38,32,34,30,27,25,24,31,26,28,29,55,48,54,50,52,53,51,49,18,19,20,21,22,23,16,17,
		60,61,58,59,57,56,63,62,19,20,18,21,16,23,17,22,46,45,40,41,42,47,44,43,14,10,12, 9,13,11,8,15,27,25,30,31,26,24,28,29,48,54,55,52,53,50,49,51,33,36,37,39,38,35,34,32, 4, 2,5, 3, 0, 6, 1, 7,
		60,61,62,63,59,58,57,56, 2,5, 4, 3, 1, 6, 7, 0,20,18,19,21,17,23,22,16,10,12,14,13,11,9, 8,15,54,55,48,53,50,52,49,51,36,37,33,38,35,39,32,34,25,30,27,26,24,31,29,28,45,40,46,41,44,47,43,42,
		60,62,57,59,56,58,61,63,30,27,25,29,26,28,24,31,12,14,10, 8,9,11,15,13,37,33,36,32,38,39,34,35,55,48,54,53,49,50,52,51,18,19,20,23,22,17,16,21,40,46,45,47,43,44,41,42, 5, 4, 2,3, 0, 1, 7, 6,
		60,62,58,56,61,63,59,57,10,12,14, 8,13,11,9,15, 2,5, 4, 3, 6, 1, 0, 7,36,37,33,39,32,38,34,35,45,40,46,44,47,43,42,41,54,55,48,50,53,49,52,51,20,18,19,17,23,22,21,16,25,30,27,29,31,28,26,24,
		60,62,63,61,59,57,56,58,4, 2,5, 3, 7, 1, 6, 0,27,25,30,29,24,28,31,26,33,36,37,38,39,32,34,35,19,20,18,22,17,23,16,21,46,45,40,43,44,47,42,41,48,54,55,49,50,53,51,52,14,10,12, 8,15,11,13, 9,
		60,63,56,59,58,57,62,61,54,55,48,51,49,52,50,53,36,37,33,34,32,39,35,38,45,40,46,42,43,47,41,44,20,18,19,22,16,17,23,21,25,30,27,28,31,24,26,29,10,12,14,11,13,15, 8,9, 2,5, 4, 3, 0, 7, 6, 1,
		60,63,57,58,62,61,59,56,33,36,37,34,38,39,32,35, 4, 2,5, 3, 1, 7, 0, 6,46,45,40,47,42,43,41,44,14,10,12,15,11,13, 9, 8,19,20,18,17,22,16,23,21,27,25,30,24,28,31,29,26,48,54,55,51,53,52,49,50,
		60,63,61,62,59,56,58,57, 5, 4, 2,3, 6, 7, 1, 0,55,48,54,51,50,52,53,49,40,46,45,43,47,42,41,44,30,27,25,31,24,28,26,29,12,14,10,13,15,11,9, 8,18,19,20,16,17,22,21,23,37,33,36,34,35,39,38,32,
		61,56,58,63,57,60,62,59,26,24,31,28,30,25,29,27,41,46,40,45,47,44,42,43,19,21,20,16,22,23,18,17,13,11,9, 8,14,12,10,15, 5, 2,3, 0, 7, 6, 1, 4,53,50,52,48,51,54,49,55,38,35,39,32,37,36,34,33,
		61,56,59,62,63,58,57,60,35,39,38,32,34,36,33,37,24,31,26,28,29,25,27,30,21,20,19,22,23,16,18,17, 2,3, 5, 7, 6, 0, 1, 4,50,52,53,51,54,48,55,49,11,9,13,14,12, 8,15,10,46,40,41,45,42,44,43,47,
		61,56,60,57,62,59,63,58,40,41,46,45,43,44,47,42,39,38,35,32,33,36,37,34,20,19,21,23,16,22,18,17,52,53,50,54,48,51,55,49, 9,13,11,12, 8,14,10,15, 3, 5, 2,6, 0, 7, 4, 1,31,26,24,28,27,25,30,29,
		61,57,56,60,62,58,59,63,41,46,40,45,44,47,43,42,53,50,52,49,48,51,54,55,26,24,31,30,25,29,28,27, 5, 2,3, 7, 0, 1, 6, 4,19,21,20,23,17,16,22,18,38,35,39,34,33,37,32,36,13,11,9,15, 8,14,12,10,
		61,57,58,62,59,63,60,56,52,53,50,49,55,51,48,54, 9,13,11,15,10,14, 8,12,31,26,24,29,30,25,28,27,39,38,35,37,34,33,36,32, 3, 5, 2,1, 7, 0, 6, 4,20,19,21,16,23,17,18,22,40,41,46,45,42,47,44,43,
		61,57,63,59,60,56,62,58,11,9,13,15,12,14,10, 8,46,40,41,45,43,47,42,44,24,31,26,25,29,30,28,27,21,20,19,17,16,23,22,18,35,39,38,33,37,34,36,32, 2,3, 5, 0, 1, 7, 4, 6,50,52,53,49,54,51,55,48,
		61,58,60,59,56,63,57,62,19,21,20,18,16,22,23,17,26,24,31,28,25,30,27,29, 5, 2,3, 6, 1, 0, 4, 7,38,35,39,37,36,34,33,32,13,11,9,14, 8,10,12,15,41,46,40,47,44,42,45,43,53,50,52,49,54,55,48,51,
		61,58,62,57,59,60,56,63,50,52,53,49,48,55,51,54,21,20,19,18,23,22,17,16, 2,3, 5, 1, 0, 6, 4, 7,11,9,13, 8,10,14,12,15,46,40,41,44,42,47,43,45,35,39,38,36,34,37,32,33,24,31,26,28,27,30,29,25,
		61,58,63,56,57,62,59,60,31,26,24,28,29,30,25,27,52,53,50,49,51,55,54,48,3, 5, 2,0, 6, 1, 4, 7,40,41,46,42,47,44,43,45,39,38,35,34,37,36,33,32, 9,13,11,10,14, 8,15,12,20,19,21,18,17,22,16,23,
		61,59,57,63,60,58,56,62, 9,13,11,15,14,10,12, 8,20,19,21,18,16,23,17,22,52,53,50,55,51,48,49,54, 3, 5, 2,7, 1, 6, 0, 4,31,26,24,25,27,29,30,28,40,41,46,44,43,42,45,47,39,38,35,32,37,34,33,36,
		61,59,58,60,56,62,63,57,21,20,19,18,22,23,16,17,35,39,38,32,36,34,37,33,50,52,53,48,55,51,49,54,46,40,41,42,44,43,47,45, 2,3, 5, 6, 7, 1, 0, 4,24,31,26,29,25,27,28,30,11,9,13,15, 8,10,14,12,
		61,59,62,56,63,57,60,58,38,35,39,32,33,34,36,37,13,11,9,15,12,10, 8,14,53,50,52,51,48,55,49,54,26,24,31,27,29,25,30,28,41,46,40,43,42,44,47,45, 5, 2,3, 1, 6, 7, 4, 0,19,21,20,18,17,23,22,16,
		61,60,57,56,62,63,58,59,46,40,41,45,47,43,44,42, 2,3, 5, 4, 0, 1, 7, 6,11,9,13,12,14,10,15, 8,35,39,38,37,33,36,34,32,24,31,26,30,27,25,29,28,50,52,53,55,48,54,49,51,21,20,19,18,17,16,23,22,
		61,60,59,58,56,57,62,63,20,19,21,18,23,16,22,17,40,41,46,45,44,43,42,47, 9,13,11,14,10,12,15, 8,31,26,24,27,25,30,29,28,52,53,50,48,54,55,51,49,39,38,35,33,36,37,32,34, 3, 5, 2,4, 7, 1, 6, 0,
		61,60,63,62,58,59,56,57, 5, 2,3, 4, 6, 1, 0, 7,19,21,20,18,22,16,17,23,13,11,9,10,12,14,15, 8,53,50,52,54,55,48,51,49,38,35,39,36,37,33,34,32,26,24,31,25,30,27,28,29,41,46,40,45,42,43,47,44,
		61,62,56,59,63,60,58,57,39,38,35,32,36,33,34,37, 3, 5, 2,4, 6, 0, 7, 1,40,41,46,43,44,47,45,42, 9,13,11,8,12,10,14,15,20,19,21,22,17,23,16,18,31,26,24,30,29,27,28,25,52,53,50,49,54,48,51,55,
		61,62,57,58,59,56,63,60,53,50,52,49,51,48,55,54,38,35,39,32,34,33,37,36,41,46,40,44,47,43,45,42,19,21,20,17,23,22,16,18,26,24,31,29,27,30,25,28,13,11,9,12,10, 8,15,14, 5, 2,3, 4, 7, 0, 1, 6,
		61,62,60,63,58,57,59,56, 2,3, 5, 4, 1, 0, 6, 7,50,52,53,49,55,48,54,51,46,40,41,47,43,44,45,42,24,31,26,27,30,29,25,28,11,9,13,10, 8,12,14,15,21,20,19,23,22,17,18,16,35,39,38,32,37,33,36,34,
		61,63,56,58,57,59,60,62,24,31,26,28,25,29,30,27,11,9,13,15,14,12, 8,10,35,39,38,34,36,33,32,37,50,52,53,54,51,55,48,49,21,20,19,16,17,22,23,18,46,40,41,43,47,42,45,44, 2,3, 5, 4, 7, 6, 0, 1,
		61,63,59,57,60,62,58,56,13,11,9,15,10,12,14, 8,5, 2,3, 4, 1, 6, 7, 0,38,35,39,33,34,36,32,37,41,46,40,42,43,47,44,45,53,50,52,55,54,51,48,49,19,21,20,22,16,17,18,23,26,24,31,28,27,29,25,30,
		61,63,62,60,58,56,57,59, 3, 5, 2,4, 0, 6, 1, 7,31,26,24,28,30,29,27,25,39,38,35,36,33,34,32,37,20,19,21,17,22,16,23,18,40,41,46,47,42,43,44,45,52,53,50,51,55,54,49,48,9,13,11,15, 8,12,10,14,
		62,56,57,63,59,61,60,58,23,22,17,16,20,21,18,19,32,38,39,35,33,37,34,36,53,49,50,51,54,48,52,55,47,43,44,45,41,40,46,42, 3, 2,4, 7, 6, 0, 1, 5,27,30,29,26,28,24,25,31,10, 8,12,14, 9,11,15,13,
		62,56,58,60,63,57,59,61,8,12,10,14,15,11,13, 9,22,17,23,16,18,21,19,20,49,50,53,54,48,51,52,55, 2,4, 3, 6, 0, 7, 1, 5,30,29,27,28,24,26,31,25,43,44,47,41,40,45,42,46,38,39,32,35,34,37,36,33,
		62,56,61,59,60,58,63,57,39,32,38,35,36,37,33,34,12,10, 8,14,13,11,9,15,50,53,49,48,51,54,52,55,29,27,30,24,26,28,31,25,44,47,43,40,45,41,46,42, 4, 3, 2,0, 7, 6, 5, 1,17,23,22,16,19,21,20,18,
		62,57,60,59,58,61,56,63,30,29,27,25,26,31,28,24,49,50,53,52,48,54,55,51,2,4, 3, 1, 7, 0, 5, 6,43,44,47,45,46,41,40,42,38,39,32,37,34,33,36,35, 8,12,10,11,15, 9,14,13,22,17,23,16,19,20,18,21,
		62,57,61,58,56,63,59,60,53,49,50,52,51,54,48,55,23,22,17,16,21,20,19,18,3, 2,4, 0, 1, 7, 5, 6,10, 8,12, 9,11,15,13,14,47,43,44,41,45,46,40,42,32,38,39,33,37,34,35,36,27,30,29,25,24,31,26,28,
		62,57,63,56,59,60,58,61,17,23,22,16,18,20,21,19,29,27,30,25,28,31,24,26, 4, 3, 2,7, 0, 1, 5, 6,39,32,38,34,33,37,36,35,12,10, 8,15, 9,11,13,14,44,47,43,46,41,45,42,40,50,53,49,52,55,54,51,48,
		62,58,57,61,56,60,63,59,49,50,53,52,54,48,51,55, 8,12,10,14,11,15, 9,13,30,29,27,26,31,28,25,24,38,39,32,34,37,36,33,35, 2,4, 3, 0, 6, 1, 7, 5,22,17,23,18,21,19,16,20,43,44,47,42,45,46,41,40,
		62,58,59,63,61,57,56,60,44,47,43,42,41,46,40,45,50,53,49,52,51,48,55,54,29,27,30,31,28,26,25,24, 4, 3, 2,6, 1, 0, 7, 5,17,23,22,21,19,18,20,16,39,32,38,37,36,34,35,33,12,10, 8,14, 9,15,13,11,
		62,58,60,56,63,59,61,57,10, 8,12,14,13,15,11,9,47,43,44,42,40,46,45,41,27,30,29,28,26,31,25,24,23,22,17,19,18,21,20,16,32,38,39,36,34,37,33,35, 3, 2,4, 1, 0, 6, 5, 7,53,49,50,52,55,48,54,51,
		62,59,56,61,60,57,58,63,32,38,39,35,37,33,36,34,27,30,29,25,26,28,24,31,23,22,17,20,21,18,16,19, 3, 2,4, 6, 7, 1, 0, 5,53,49,50,48,55,51,54,52,10, 8,12,15,13, 9,14,11,47,43,44,42,45,41,40,46,
		62,59,57,60,58,63,61,56,29,27,30,25,31,28,26,24,44,47,43,42,46,41,45,40,17,23,22,18,20,21,16,19,12,10, 8,9,15,13,11,14, 4, 3, 2,1, 6, 7, 0, 5,50,53,49,51,48,55,52,54,39,32,38,35,34,33,37,36,
		62,59,63,58,61,56,60,57,43,44,47,42,40,41,46,45,38,39,32,35,36,33,34,37,22,17,23,21,18,20,16,19,49,50,53,55,51,48,54,52, 8,12,10,13, 9,15,11,14, 2,4, 3, 7, 1, 6, 5, 0,30,29,27,25,24,28,31,26,
		62,60,56,58,63,61,57,59,12,10, 8,14,11,13,15, 9, 4, 3, 2,5, 0, 7, 6, 1,39,32,38,36,37,33,35,34,44,47,43,45,40,46,41,42,50,53,49,54,55,48,51,52,17,23,22,20,18,19,16,21,29,27,30,25,24,26,28,31,
		62,60,59,57,58,56,63,61,27,30,29,25,28,26,31,24,10, 8,12,14,15,13, 9,11,32,38,39,37,33,36,35,34,53,49,50,55,48,54,51,52,23,22,17,18,19,20,21,16,47,43,44,40,46,45,42,41,3, 2,4, 5, 6, 7, 1, 0,
		62,60,61,63,57,59,58,56, 2,4, 3, 5, 1, 7, 0, 6,30,29,27,25,31,26,24,28,38,39,32,33,36,37,35,34,22,17,23,19,20,18,21,16,43,44,47,46,45,40,41,42,49,50,53,48,54,55,52,51,8,12,10,14, 9,13,11,15,
		62,61,58,57,56,59,60,63,50,53,49,52,48,51,54,55,39,32,38,35,37,36,34,33,44,47,43,41,46,40,42,45,17,23,22,19,21,20,18,16,29,27,30,26,24,31,28,25,12,10, 8,13,11,9,14,15, 4, 3, 2,5, 6, 1, 0, 7,
		62,61,59,56,60,63,57,58,38,39,32,35,33,36,37,34, 2,4, 3, 5, 7, 1, 6, 0,43,44,47,40,41,46,42,45, 8,12,10, 9,13,11,15,14,22,17,23,20,19,21,18,16,30,29,27,31,26,24,25,28,49,50,53,52,55,51,48,54,
		62,61,63,60,57,58,56,59, 3, 2,4, 5, 0, 1, 7, 6,53,49,50,52,54,51,55,48,47,43,44,46,40,41,42,45,27,30,29,24,31,26,28,25,10, 8,12,11,9,13,15,14,23,22,17,21,20,19,16,18,32,38,39,35,34,36,33,37,
		62,63,56,57,59,58,61,60,22,17,23,16,21,18,20,19,43,44,47,42,41,40,45,46, 8,12,10,15,11,13,14, 9,30,29,27,24,28,31,26,25,49,50,53,51,55,54,48,52,38,39,32,36,33,34,35,37, 2,4, 3, 5, 6, 0, 7, 1,
		62,63,58,59,61,60,57,56,47,43,44,42,46,40,41,45, 3, 2,4, 5, 1, 0, 6, 7,10, 8,12,13,15,11,14, 9,32,38,39,34,36,33,37,35,27,30,29,31,24,28,26,25,53,49,50,54,51,55,52,48,23,22,17,16,19,18,21,20,
		62,63,60,61,57,56,59,58,4, 3, 2,5, 7, 0, 1, 6,17,23,22,16,20,18,19,21,12,10, 8,11,13,15,14, 9,50,53,49,55,54,51,48,52,39,32,38,33,34,36,37,35,29,27,30,28,31,24,25,26,44,47,43,42,45,40,46,41,
		63,56,60,59,57,62,58,61,54,51,55,48,49,53,52,50,16,17,22,23,18,19,20,21,4, 5, 3, 7, 6, 0, 2,1,13,15,11,14,12, 8,10, 9,43,47,42,45,41,40,46,44,34,36,33,39,35,32,37,38,31,24,28,26,30,27,25,29,
		63,56,61,58,59,60,57,62,24,28,31,26,25,27,29,30,51,55,54,48,52,53,50,49, 5, 3, 4, 6, 0, 7, 2,1,47,42,43,41,40,45,46,44,36,33,34,35,32,39,38,37,15,11,13,12, 8,14, 9,10,17,22,16,23,20,19,21,18,
		63,56,62,57,58,61,59,60,22,16,17,23,21,19,18,20,28,31,24,26,29,27,30,25, 3, 4, 5, 0, 7, 6, 2,1,33,34,36,32,39,35,38,37,11,13,15, 8,14,12,10, 9,42,43,47,40,45,41,44,46,55,54,51,48,50,53,49,52,
		63,57,56,62,58,60,61,59,16,17,22,23,19,18,21,20,34,36,33,37,39,35,32,38,54,51,55,49,53,52,48,50,43,47,42,41,45,46,40,44, 4, 5, 3, 0, 1, 7, 6, 2,31,24,28,25,29,30,26,27,13,15,11,9,14,12, 8,10,
		63,57,59,61,62,56,58,60,15,11,13, 9, 8,12,10,14,17,22,16,23,21,18,20,19,51,55,54,53,52,49,48,50, 5, 3, 4, 1, 7, 0, 6, 2,24,28,31,29,30,25,27,26,47,42,43,45,46,41,44,40,36,33,34,37,32,35,38,39,
		63,57,60,58,61,59,62,56,33,34,36,37,38,35,39,32,11,13,15, 9,10,12,14, 8,55,54,51,52,49,53,48,50,28,31,24,30,25,29,27,26,42,43,47,46,41,45,40,44, 3, 4, 5, 7, 0, 1, 2,6,22,16,17,23,20,18,19,21,
		63,58,56,61,59,62,60,57,28,31,24,26,27,29,25,30,42,43,47,44,40,45,41,46,22,16,17,21,19,18,23,20,11,13,15,14, 8,10,12, 9, 3, 4, 5, 6, 1, 0, 7, 2,55,54,51,49,52,50,48,53,33,34,36,37,32,39,35,38,
		63,58,57,60,61,56,59,62,34,36,33,37,35,39,38,32,31,24,28,26,25,29,30,27,16,17,22,19,18,21,23,20, 4, 5, 3, 1, 0, 6, 7, 2,54,51,55,52,50,49,53,48,13,15,11,8,10,14, 9,12,43,47,42,44,41,45,46,40,
		63,58,62,59,60,57,61,56,47,42,43,44,46,45,40,41,36,33,34,37,38,39,32,35,17,22,16,18,21,19,23,20,51,55,54,50,49,52,53,48,15,11,13,10,14, 8,12, 9, 5, 3, 4, 0, 6, 1, 2,7,24,28,31,26,30,29,27,25,
		63,59,56,60,57,61,62,58,51,55,54,48,53,52,49,50,15,11,13, 9,12, 8,14,10,24,28,31,25,27,29,26,30,36,33,34,32,35,38,39,37, 5, 3, 4, 7, 1, 6, 0, 2,17,22,16,21,18,20,23,19,47,42,43,44,41,40,45,46,
		63,59,58,62,60,56,57,61,42,43,47,44,45,40,46,41,55,54,51,48,49,52,50,53,28,31,24,27,29,25,26,30, 3, 4, 5, 1, 6, 7, 0, 2,22,16,17,18,20,21,19,23,33,34,36,35,38,32,37,39,11,13,15, 9,14, 8,10,12,
		63,59,61,57,62,58,60,56,13,15,11,9,10, 8,12,14,43,47,42,44,46,40,41,45,31,24,28,29,25,27,26,30,16,17,22,20,21,18,19,23,34,36,33,38,32,35,39,37, 4, 5, 3, 6, 7, 1, 2,0,54,51,55,48,50,52,53,49,
		63,60,58,57,61,62,56,59,36,33,34,37,39,38,35,32, 5, 3, 4, 2,0, 6, 1, 7,47,42,43,46,45,40,44,41,15,11,13,14,10,12, 8,9,17,22,16,19,20,18,21,23,24,28,31,27,25,30,26,29,51,55,54,48,50,49,52,53,
		63,60,59,56,57,58,61,62,55,54,51,48,52,49,53,50,33,34,36,37,35,38,32,39,42,43,47,45,40,46,44,41,22,16,17,20,18,19,21,23,28,31,24,25,30,27,29,26,11,13,15,10,12,14, 9, 8,3, 4, 5, 2,1, 6, 7, 0,
		63,60,62,61,56,59,57,58,4, 5, 3, 2,7, 6, 0, 1,54,51,55,48,53,49,50,52,43,47,42,40,46,45,44,41,31,24,28,30,27,25,29,26,13,15,11,12,14,10, 8,9,16,17,22,18,19,20,23,21,34,36,33,37,32,38,39,35,
		63,61,57,59,62,60,56,58,11,13,15, 9,12,10, 8,14, 3, 4, 5, 2,7, 0, 1, 6,33,34,36,38,35,39,37,32,42,43,47,41,46,40,45,44,55,54,51,53,50,52,49,48,22,16,17,19,21,20,23,18,28,31,24,26,30,25,29,27,
		63,61,58,56,59,57,62,60,31,24,28,26,29,25,27,30,13,15,11,9, 8,10,14,12,34,36,33,35,39,38,37,32,54,51,55,50,52,53,49,48,16,17,22,21,20,19,18,23,43,47,42,46,40,41,44,45, 4, 5, 3, 2,1, 0, 6, 7,
		63,61,60,62,56,58,59,57, 5, 3, 4, 2,6, 0, 7, 1,24,28,31,26,27,25,30,29,36,33,34,39,38,35,37,32,17,22,16,20,19,21,18,23,47,42,43,40,41,46,45,44,51,55,54,52,53,50,48,49,15,11,13, 9,14,10,12, 8,
		63,62,57,56,58,59,60,61,17,22,16,23,18,21,19,20,47,42,43,44,45,46,41,40,15,11,13, 8,12,10, 9,14,24,28,31,30,29,27,25,26,51,55,54,49,50,53,52,48,36,33,34,38,39,32,37,35, 5, 3, 4, 2,1, 7, 0, 6,
		63,62,59,58,60,61,56,57,43,47,42,44,40,46,45,41,4, 5, 3, 2,6, 7, 1, 0,13,15,11,10, 8,12, 9,14,34,36,33,32,38,39,35,37,31,24,28,27,30,29,25,26,54,51,55,53,49,50,48,52,16,17,22,23,20,21,18,19,
		63,62,61,60,56,57,58,59, 3, 4, 5, 2,0, 7, 6, 1,22,16,17,23,19,21,20,18,11,13,15,12,10, 8,9,14,55,54,51,50,53,49,52,48,33,34,36,39,32,38,35,37,28,31,24,29,27,30,26,25,42,43,47,44,41,46,40,45,
};