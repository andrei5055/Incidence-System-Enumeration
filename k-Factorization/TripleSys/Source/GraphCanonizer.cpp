#include "GraphCanonizer.h"

#if DEBUGING_CANON
int canonMatrCntr = 0;

#define SAME_FILE_OUTPUT 1  // Output into the same (1) or different (0) file(s)
#define LOG_FILE_PATH	""
#if SAME_FILE_OUTPUT   
#define CANON_MATR_LOG	LOG_FILE_PATH "canonMatrLog.txt"
#define CANON_MATR_FRMT "%s"
#else   
#define CANON_MATR_LOG	"canonMatrLog"
#define CANON_MATR_FRMT LOG_FILE_PATH "%s_%02d.txt"
#endif
#endif

#if PRINT_MATRICES || PRINT_NUM_CUR_GRAPH
extern int numCurrGraph;
#endif

#if PRINT_MATRICES
tchar* ptrGraph[2] = { NULL };
static int nIter = 0;
static bool printFlag = false;
static int hhh;
#endif

CGraphCanonizer::CGraphCanonizer(int nVert) :
	Generators<ushort>(0, "\nVertex orbits and group generators of graph", nVert) {
	Init(nVert);
}

int CGraphCanonizer::Init(int nVert) {
	if (nVert <= 0)
		return -1;

	releaseCanonizerMemory();

	m_v = nVert;
	m_lenGraphMatr = m_v * m_v;

	const int coeff = PRINT_MATRICES ? 3 : 0;
	m_pGraph[0] = new tchar[(2 + coeff) * m_lenGraphMatr];
	m_pGraph[1] = m_pGraph[0] + m_lenGraphMatr;

	m_pNumOrbits = new ushort[(3 + coeff) * m_v];
	m_pGroupOrbits = m_pNumOrbits + m_v;
	m_pOrbits = (m_pNumOrbits + m_v) + m_v;
	m_bUsedFlags = new tchar[m_v * m_v];

	const auto len = nVert * (nVert - 1) / 2;
	m_pLenOrbits = new ushort[3 * len];
	m_pSavedOrbits = m_pLenOrbits + len;
	m_pSavedOrbIdx = m_pSavedOrbits + len;
	return 0;
}

void CGraphCanonizer::initCanonizer() {
	// Initialize orbits
	m_pNumOrbits[0] = 1;
	m_pLenOrbits[0] = m_v;

	// Initiate flags
	memset(m_bUsedFlags, 0, m_v * m_v);
}

void CGraphCanonizer::releaseCanonizerMemory() {
	delete[] m_pGraph[0];
	delete[] m_pNumOrbits;
	delete[] m_bUsedFlags;
	delete[] m_pLenOrbits;
}

ctchar* CGraphCanonizer::canonize_graph(ctchar* pGraph, int* pCanonIndex) {
	initCanonizer();
	if (pGraph)
		memcpy(m_pGraph[0], pGraph, m_lenGraphMatr);

	int i, firstVert = 0;
	i = 0;

	PRINT_ADJ_MATRIX(m_pGraph[0], ++canonMatrCntr, m_v, NULL, CANON_MATR_LOG, (SAME_FILE_OUTPUT && canonMatrCntr? "a" : "w"), CANON_MATR_FRMT);

#if PRINT_NUM_CUR_GRAPH
	printfYellow("**** numCurrGraph = %d ****\n", ++numCurrGraph);
#endif
#if PRINT_MATRICES
#if !PRINT_NUM_CUR_GRAPH
	++numCurrGraph;
#endif
	auto* pInitOrbits = m_pOrbits + m_v;
	auto* pResOrbits = pInitOrbits + m_v;
	if (numCurrGraph == N_MATR) {
		printFlag = true;
		for (int j = m_v; j--;)
			m_pOrbits[j] = j;

		memcpy(ptrGraph[0] = m_pGraph[i] + 2 * m_lenGraphMatr, m_pGraph[0], m_lenGraphMatr * sizeof(m_pGraph[0][0]));
		ptrGraph[1] = ptrGraph[0] + m_lenGraphMatr;
		memcpy(pInitOrbits, m_pOrbits, m_v * sizeof(m_pOrbits[0]));
		PRINT_ADJ_MATRIX(m_pGraph[i], -1, m_v);
	}
	else
		printFlag = false;
#endif
	while (firstVert = canonizeGraph(m_pGraph[i], m_pGraph[1 - i], firstVert)) {
		createGraphOut(m_pGraph[i], m_pGraph[1 - i]);

#if PRINT_MATRICES
		if (printFlag) {
			for (int j = m_v; j--;)
				pResOrbits[j] = pInitOrbits[m_pOrbits[j]];

			memcpy(pInitOrbits, pResOrbits, m_v * sizeof(m_pOrbits[0]));
			// Construct matrix from the original one by using pInitOrbits permutation
			createGraphOut(ptrGraph[0], ptrGraph[1], 0, 0, pInitOrbits);
			PRINT_ADJ_MATRIX(ptrGraph[1], nIter, m_v, pInitOrbits, "bbb");
			// Matrix after nIter iterations of canonization
			// PRINT_ADJ_MATRIX(m_pGraph[1 - i], nIter, m_v);
			// These two matrices should be identical
			assert(!memcmp(ptrGraph[1], m_pGraph[1 - i], m_lenGraphMatr * sizeof(*ptrGraph[1])));
			nIter++;
		}
#endif
		i = 1 - i;
	}

	if (pCanonIndex)
		*pCanonIndex = i;

	PRINT_ADJ_MATRIX(m_pGraph[i], canonMatrCntr, m_v, NULL, CANON_MATR_LOG, "a", CANON_MATR_FRMT);
	return m_pGraph[i];
}

void CGraphCanonizer::initVertexGroupOrbits() {
	setGroupOrder(1);
	setStabilizerLengthAut(m_v);
	for (int i = m_v; i--;)
		m_pGroupOrbits[i] = i;
	memcpy(this->getObject(1), m_pGroupOrbits, m_v * sizeof(m_pGroupOrbits[0]));
}

tchar* CGraphCanonizer::createGraphOut(ctchar* pGraph, tchar* pGraphOut, int startVertex, int endVertex, const ushort* pOrb) const {
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

ushort* CGraphCanonizer::restoreParam(int& i, int iStart, ushort* pLenOrbits) {
	do {
		// Restore previous parameters;
		while (i > iStart && !*(pLenOrbits -= m_v - i--));

		if (i == iStart)
			return NULL;

	} while (m_pSavedOrbIdx[i]++ >= *pLenOrbits);

	return pLenOrbits;
}

int CGraphCanonizer::canonizeGraph(ctchar* pGraph, tchar* pGraphOut, int firstVert) {
#if PRINT_MATRICES
	if (firstVert) {
		// Copy the top part of the adjacency matrix
		memcpy(pGraphOut, pGraph, m_v * (m_v - firstVert));
	}
#endif
#if 0
	startCanonize:
#endif
	const auto defineAut = firstVert > 0;
	int lastUnfixedVertexIndex = 0;
	// Vertex index which will be swapped with the vertex firstVert
	int indVertMax = 1;
	int indVert = 0;
	resetGroup();
	releaseAllObjects();
	// Allocate space for three vectors
	// (a) vertex orbits
	// (b) vertex orbits under the stabilizer of vertex 0
	// (c) trivial permutation
	reserveObjects(3);

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
		if (DO_PRINT(nIter)) {
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
			if (DO_PRINT(nIter)) {
				FOPEN_F(f, "bbb.txt", "a");
				fprintf(f, "ff = %2d  i = %2d  indVert = %2d\n", ff, i, indVert);
				FCLOSE_F(f);
			}
#endif
			if (indVert <= indVertMax) {
				// swapping of i-th and (i+indVert)-th vertices
				m_pOrbits[m_pOrbits[i] = i + indVert] = i;
				++*pLenOrbits;
				break;
			}
		}

		int flag = 0;
		int idxRight = 1;  // Continue, until at least one orbit length > 1
		if (defineAut) {
			const auto iStart = i;
			while (idxRight && i < m_v) {
				if (*pLenOrbits > 1) {
					m_pSavedOrbIdx[i] = 0;
					memset(m_bUsedFlags + i * m_v, 0, m_v);
				}
			canonRow:
#if FFF
				if (DO_PRINT(nIter)) {
					const bool flg = i == FF_ && idxRight == 1;
					sprintf_s(buffer, "ccc_%04d.txt", hhh += flg ? 1 : 0);
					if (flg) {
						FOPEN_F(f, "bbb.txt", "a");
						fprintf(f, "ff = %2d  canonizeMatrixRow: hhh = %3d\n", ff, hhh);
						FCLOSE_F(f);
						printAdjMatrix(pGraphOut, buffer, i + 1);
					}

					FOPEN_F(f, buffer, "a");
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
					auto pOrbit = m_pOrbits + i;
					auto j = m_pSavedOrbIdx[i];
					auto pUsedFlag = m_bUsedFlags + i * m_v;
					if (j > 1) {
						while (pUsedFlag[pOrbit[j]])
							j--;
					}

					// Replace the current orbit’s leading vertex with 
					// the next candidate and mark it as used
					const auto tmp = *pOrbit;
					*pOrbit = *(pOrbit + j);
					pUsedFlag[*(pOrbit + j) = tmp] = 1;
					++*pLenOrbits;
#if 0
					if (nIter == 40 && m_pOrbits[0] == 29 && i == 1)
						i += 0;
#endif
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

				FOPEN_F(f1, "ccc.txt", fff > 1 ? "a" : "w");
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
				if (DO_PRINT(nIter)) {
					printAdjMatrix(pGraphOut, buffer, m_v);
					printAdjMatrix(pGraph, buffer, m_v);
				}
#endif
				flag = memcmp(pGraphLast, pGraph + i * m_v, m_v * (m_v - i));
				if (!flag) {
					addAutomorphism(m_v, m_pOrbits, m_pGroupOrbits, true);
					if (!m_pOrbits[0]) {
						// Set vertex orbits under the stabilizer of vertex 0
						memcpy(this->getObject(1), m_pGroupOrbits, sizeof(m_pGroupOrbits[0]) * m_v);
					}
				}
			}
		}
		else {
#if 0
			while (idxRight && i < m_v) {
				flag = canonizeMatrixRow(pGraph, pGraphOut, i++, &pLenOrbits, idxRight, flag, lastUnfixedVertexIndex);
			}
#else
			while (i < m_v) {
				if (idxRight)
					flag = canonizeMatrixRow(pGraph, pGraphOut, i++, &pLenOrbits, idxRight, flag, lastUnfixedVertexIndex);
				else {
					const auto vertIdx = i++;
					const auto pVertexIn = pGraph + m_pOrbits[vertIdx] * m_v;
					auto pVertOut = pGraphOut + vertIdx * m_v;
					for (int j = 0; j < m_v; j++)
						pVertOut[j] = pVertexIn[m_pOrbits[j]];

					flag = memcmp(pVertOut, pGraph + vertIdx * m_v, m_v);
					if (flag)
						return lastUnfixedVertexIndex;
				}
			}

#if 0
			// This is an attempt to create group during re-check of cannonicity
			// Commented out because it needs some debuging
			// Run CI, to see the problem. 
			firstVert = lastUnfixedVertexIndex;
			goto startCanonize;
#else
			return 0;
#endif
#endif
		}

		if (flag > 0 || !defineAut)
			return lastUnfixedVertexIndex;
	}

	return 0;
}

int CGraphCanonizer::canonizeMatrixRow(ctchar* pGraph, tchar* pGraphOut, int vertIdx,
	ushort** ppLenOrbits, int& idxRight, int flag, int& lastUnfixedVertexIndex) {

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
		memcpy(pLenOrbitsNext, ++pLenOrbits, len);
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
		ASSERT_IF(!lenOrb);
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
			ASSERT_IF(!splitPos || lenOrb == splitPos);
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

void CGraphCanonizer::outAdjMatrix(ctchar* pGraphOut, FILE* f, int endVertex) const {
	const auto v = groupDegree();
	if (!endVertex)
		endVertex = v;

	char buf[1024], * pBuf;
	for (int i = 0; i < endVertex; i++, pGraphOut += v) {
		pBuf = buf;
		for (int j = 0; j < v; j++)
			SPRINTFD(pBuf, buf, "%1d", pGraphOut[j]);

		*(pBuf - v + i) = '*';		// *'s on diagonal
		SPRINTFD(pBuf, buf, " %2d", i); // row #'s from left
		fprintf(f, "%s\n", buf);
	}
}

#if PRINT_MATRICES || DEBUGING_CANON
void CGraphCanonizer::printAdjMatrix(ctchar* pGraphOut, int idx, int endVertex, ushort* pOrb, const char* fName, cchar* mode, cchar *frmt) const {
	char buf[1024];
	snprintf(buf, sizeof(buf), frmt, fName, idx);
	printAdjMatrix(pGraphOut, buf, endVertex, idx, pOrb, mode);
}

void CGraphCanonizer::printAdjMatrix(ctchar* pGraphOut, const char* fileName, int endVertex, int idx, ushort* pOrb, cchar* mode) const {
	FOPEN_F(f, fileName, mode);
	fprintf(f, "\nAdjacency matrix for graph %d\n", idx);

#if !DEBUGING_CANON
	char buf[1024], * pBuf = buf;
	const auto v = groupDegree();
	int iMax = v <= 30 ? v : 20;
	int k = 0;
	if (!pOrb)
		pOrb = vertOrbits();

	const auto jMax = (v + iMax - 1) / iMax;
	for (int j = 0; j++ < jMax;) {
		if (j == jMax && jMax > 1)
			iMax = v % iMax;

		for (int i = 0; i < iMax; i++)
			SPRINTFD(pBuf, buf, "%3d", pOrb[k++]);

		SPRINTFD(pBuf, buf, "\n");
	}
	fprintf(f, "%s\n", buf);

	pBuf = buf;
	for (int i = 0; i < v; i += 10) {
		SPRINTFD(pBuf, buf, "%1d", i / 10);
		for (int j = 0; j < 9; j++)
			SPRINTFD(pBuf, buf, " ");
	}
	fprintf(f, "%s\n", buf);

	pBuf = buf;
	for (int i = 0; i < v; i++)
		SPRINTFD(pBuf, buf, "%1d", i % 10);

	fprintf(f, "%s\n", buf);
#endif
	outAdjMatrix(pGraphOut, f, endVertex);
	FCLOSE_F(f);
}
#endif