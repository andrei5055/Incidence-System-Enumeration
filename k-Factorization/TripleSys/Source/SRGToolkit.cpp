#include "SRGToolkit.h"
#include "SRGSupport.h"
#include <cstring>

#pragma execution_character_set("utf-8")

#define PRINT_NUM_CUR_GRAPH TRACE_GROUP_ORDER
#define PRINT_MATRICES		0

extern short* pGenerator = NULL;
#if PRINT_NUM_CUR_GRAPH
int numCurrGraph;
#endif

#if PRINT_MATRICES
#define PRINT_ADJ_MATRIX(...) printAdjMatrix(__VA_ARGS__)
#define DO_PRINT(nIter)   (printFlag && FFF != -1 && nIter >= FFF)
// Parameters specific to the bug we are trying to fix.
#define N_MATR 3     // Number of matrix to activate the output
#define NUM_GENERATOR	1  // Generator number that is presumably missing from the matrix with the smaller group
#define FFF 22  // should be set equal to index NM of last bbb_NM.txt file
#define FF_ 1 // 2   // 4 - for 26. 2 for 22
tchar* pGraph[2] = { NULL };
static int nIter = 0;
static bool printFlag = false;
static int hhh;
#else
#define PRINT_ADJ_MATRIX(...)
#endif

#define PERMUT			1
#define USE_COMPLEMENT -1   // -1 use compliment graph only when 2*k > v
                            //  0 don't use compliment graph
                            //  1 always use complement graph
                             

// Use PERMUT equal 
//    0 to find the automorphism of original matrix (it will be in vvv.txt)
//    1 to find the automorphism of canonized complement matrix
#if PRINT_NUM_CUR_GRAPH && PRINT_MATRICES && PERMUT
ushort autIni[] = {  // Could be found in vvv.txt
#if N_MATR == 3
#if NUM_GENERATOR == 1
	 0, 34, 18, 17, 32, 47, 13,  7, 23, 43, 20, 41, 31,  6, 35, 37, 29,  3,  2, 25,
	10, 21, 44,  8, 46, 19, 26, 40, 28, 16, 36, 12,  4, 33,  1, 14, 30, 15, 48, 45,
	27, 11, 42,  9, 22, 39, 24,  5, 38
#elif NUM_GENERATOR == 2
	19, 46, 21, 26, 44, 40,  8,  7,  6, 41, 23, 31, 43, 20, 34, 32, 47, 17, 18,  0,
	13,  2, 37, 10, 35, 25,  3, 29, 48, 27, 30, 11, 15, 45, 14, 24, 42, 22, 38, 39,
	 5,  9, 36, 12,  4, 33,  1, 16, 28
#endif
#elif N_MATR == 5
	41, 40, 37, 39, 35, 36, 38, 23, 26, 24, 25, 27, 22, 21, 48, 45, 43, 47, 42, 44,
	46, 13, 12,  7,  9, 10,  8, 11, 29, 28, 33, 32, 31, 30, 34,  4,  5,  2,  6,  3,
	 1,  0, 18, 16, 19, 15, 20, 17, 14
#endif
};

ushort autLost[] = {
#if N_MATR == 3
#if NUM_GENERATOR == 1
	 0,  7,  9, 13, 11, 17, 15,  1,  8,  2, 10,  4, 12,  3, 14,  6, 18,  5, 16, 23,
	29, 27, 43, 19, 30, 25, 45, 21, 47, 20, 24, 33, 41, 31, 48, 35, 39, 46, 40, 36,
	38, 32, 44, 22, 42, 26, 37, 28, 34
#endif
#elif N_MATR == 5
	29, 30,  8,  7, 43, 45, 47,  3,  2, 12, 20, 21,  9, 38, 33, 42, 31, 39, 35, 22,
	10, 11, 19, 40, 32, 37, 36, 41, 34,  0,  1, 16, 24, 14, 28, 18, 26, 25, 13, 17,
	23, 27, 15,  4, 48,  5, 46,  6, 44
#endif
};
#endif

SRGToolkit::SRGToolkit(const kSysParam* p, int nCols, int nRows, int groupSize, const std::string& resFileName, int exploreMatrices) :
	m_pParam(p), m_nCols(nCols), m_nRows(nRows), m_groupSize(groupSize), m_v(groupDegree()), m_nExploreMatrices(exploreMatrices),
	m_len(groupDegree() * (groupDegree() - 1) / 2), m_lenGraphMatr(groupDegree() * groupDegree()),
	m_resFileName(resFileName), Generators<ushort>(0, "\nVertex orbits and group generators of graph", nRows * nCols / groupSize) {
	setOutFileName(NULL);
	const int coeff = PRINT_MATRICES? 3 : 0;
	m_pGraph[0] = new tchar[(2 + coeff) * m_lenGraphMatr];
	m_pGraph[1] = m_pGraph[0] + m_lenGraphMatr;
	m_subgraphVertex = new ushort[m_v];
	m_pOrbits = new ushort[(1 + coeff) * m_v];
	m_pNumOrbits = new ushort[2 * m_v];
	m_pGroupOrbits = m_pNumOrbits + m_v;
	m_pLenOrbits = new ushort [3 * m_len];
	m_pSavedOrbits = m_pLenOrbits + m_len;
	m_pSavedOrbIdx = m_pSavedOrbits + m_len;
	for (int i = 0; i < 2; i++) {
		m_bChekMatr[i] = true;
		m_pGraphParam[i] = new SRGParam();
		m_pMarixStorage[i] = new CBinaryMatrixStorage((int)m_len, 50);
	}

	m_bUsedFlags = new tchar[m_v * m_v];
}

SRGToolkit::~SRGToolkit() { 
	delete[] m_pGraph[0];
	delete[] m_subgraphVertex;
	delete[] m_pNumOrbits;
	delete[] m_pOrbits;
	delete[] outFileName();
	delete[] m_bUsedFlags;
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

bool SRGToolkit::exploreMatrix(ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID, uint srcGroupOrder) {
	int counter = 0;
	for (int i = 0; i < 2; i++) {
		if (!m_bChekMatr[i])
			continue;

		if (exploreMatrixOfType(i, pMatr, pGraphDB+i, sourceMatrID, srcGroupOrder))
			counter++;			
		else 
			if (!(m_nExploreMatrices & 2))
				m_bChekMatr[i] = false;
	}

	return counter > 0;
}

bool SRGToolkit::exploreMatrixOfType(int typeIdx, ctchar* pMatr, GraphDB* pGraphDB, uint sourceMatrID, uint srcGroupOrder) {
	const auto numGroups = m_nCols / m_groupSize;
	const auto pVertexLast = pMatr + m_nRows * m_nCols;

	auto graphParam = m_pGraphParam[typeIdx];
	const auto cond = typeIdx != 0;
	const auto mult = cond ? m_groupSize : numGroups - m_groupSize;
	// Create graph
	auto pAdjacencyMatrix = m_pGraph[0];
	memset(pAdjacencyMatrix, 0, m_lenGraphMatr * sizeof(pAdjacencyMatrix[0]));
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
			int nConnected = (m_nRows - i - 1) * mult;
			auto* pNextVertex = pMatr + (i + 1) * m_nCols;
			while (pNextVertex < pVertexLast) {
				numNextVertex++;
				// We do have two pointers (pVertex and pNextVertex).
				// If one out of groupSize elements is the same, then we have an edge between two vertices.
				if (one_common_element(pVertex, pNextVertex, m_groupSize) == cond) {
					// Add edges to the graph for the vertices intersecting by one element
					pAdjacencyMatrix[numVertex * m_v + numNextVertex] =
						pAdjacencyMatrix[numNextVertex * m_v + numVertex] = 1;

					if (!(--nConnected))
						break;
				}

				// Move to the next vertex
				pNextVertex += m_groupSize;
			}
		}
	}

	const auto graphType = checkSRG(pAdjacencyMatrix, graphParam);
	pGraphDB->setGraphType(graphType);
	int flg = 1;
	switch (graphType) {
	case t_nonregular:
	case t_complete:
		delete graphParam;
		m_pGraphParam[typeIdx] = NULL;
		return false;
	case t_regular:
		flg = 2;
	}

	const bool canonize = m_nExploreMatrices > 0;
	if (!((canonize? m_nExploreMatrices : -m_nExploreMatrices) & flg)) {
		delete graphParam;
		m_pGraphParam[typeIdx] = NULL;
		return false;
	}

	bool rank3 = false;
	tchar* pGraph[2] = { NULL, NULL };
	tchar* pUpperDiag = NULL;
	ctchar* pResGraph = NULL;
	if (canonize) {
		initCanonizer();
		int i, firstVert = 0;
		i = 0;
#if PRINT_NUM_CUR_GRAPH
		printfYellow("**** numCurrGraph = %d ****\n", ++numCurrGraph);
#endif
#if PRINT_MATRICES
		auto* pInitOrbits = m_pOrbits + m_v;
		auto* pResOrbits = pInitOrbits + m_v;
		if (numCurrGraph == N_MATR) {
			printFlag = true;
			for (int j = m_v; j--;)
				m_pOrbits[j] = j;

			memcpy(pGraph[0] = m_pGraph[i] + 2 * m_lenGraphMatr, m_pGraph[0], m_lenGraphMatr * sizeof(m_pGraph[0][0]));
			pGraph[1] = pGraph[0] + m_lenGraphMatr;
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
				createGraphOut(pGraph[0], pGraph[1], 0, 0, pInitOrbits);
				PRINT_ADJ_MATRIX(pGraph[1], nIter, m_v, pInitOrbits, "bbb");
				// Matrix after nIter iterations of canonization
				// PRINT_ADJ_MATRIX(m_pGraph[1 - i], nIter, m_v);
				// These two matrices should be identical
				assert(!memcmp(pGraph[1], m_pGraph[1 - i], m_lenGraphMatr * sizeof(*pGraph[1])));
				nIter++;
			}
#endif
			i = 1 - i;
		}

#if 0
		// Check if canonical graph is really canonical
		// Warning: This call removes all previously determined generators of Aut(G).
		//          For debugging use only.
		ASSERT(canonizeGraph(m_pGraph[i], m_pGraph[1 - i], 0));
#endif
#if PRINT_NUM_CUR_GRAPH && PRINT_MATRICES
		if (printFlag) {
			// Create the reverse permutation for pInitOrbits
			for (int j = m_v; j--;)
				pResOrbits[pInitOrbits[j]] = j;

			// Construct matrix from the resulting one by using pResOrbits permutation
			createGraphOut(pGraph[1], m_pGraph[1 - i], 0, 0, pResOrbits);
			PRINT_ADJ_MATRIX(m_pGraph[1 - i], 0, m_v, pResOrbits, "ddd");
			// The constructed graph must be identical to the original
			assert(!memcmp(pGraph[0], m_pGraph[1 - i], m_lenGraphMatr * sizeof(*pGraph[0])));

			// Save resulting permutation which convert matrix of the resulting graph into it's original state.
			memcpy(pResOrbits + m_v, pResOrbits, m_v * sizeof(*pResOrbits));
		}
#endif
		// Copy elements above the main diagonal into the array.
		pResGraph = m_pGraph[i];
		auto pFrom = m_pGraph[i];
		auto pTo = pUpperDiag = m_pGraph[1 - i];
		for (int j = m_v, i = 0; --j; pTo += j, pFrom += m_v)
			memcpy(pTo, pFrom + ++i, j);

		// Copying vertex orbits and trivial permutation
		auto pntr = this->getObject(0);
		memcpy(pntr, m_pGroupOrbits, m_v * sizeof(*pntr));
		pntr += 2 * (i = m_v);
		rank3 = true;
		while (i--) {
			pntr[i] = i;
			if (rank3)
				rank3 = !m_pGroupOrbits[i];
		}

		// Analyze the stabilizer of first vertex
		pntr -= m_v;
		const auto j = graphParam->k + 1;
		for (i = 1; i < j; i++)
			rank3 &= pntr[i] == 1;

		while (i < m_v && rank3)
			rank3 = pntr[i++] == j;

#if PRINT_NUM_CUR_GRAPH && PRINT_MATRICES
		if (printFlag) {
			// pInitOrbits - permutation, that canonize original matrix
			const auto ppp = pResOrbits + m_v;
			ushort s[49];
#if PERMUT
			// Convert it to the automorphism of the original matrix 
			// Using ChatGPT suggestion https://chatgpt.com/share/68ab681a-aa0c-8010-a03d-7c9afb22af14
			// s = q∘p∘q^−1.
			// Note: q^−1 is in ppp = pResOrbits + m_v

			for (int j = m_v; j--;)
				pResOrbits[j] = ppp[autIni[j]];

			for (int j = m_v; j--;)
				s[j] = pResOrbits[pInitOrbits[j]];

			createGraphOut(pResGraph, pGraph[1], 0, 0, s);
			assert(!memcmp(pResGraph, pGraph[1], m_lenGraphMatr * sizeof(*pGraph[0])));
			PRINT_ADJ_MATRIX(pGraph[1], 0, m_v, s, "vvv");

			ASSERT(canonizeGraph(m_pGraph[i], m_pGraph[1 - i], 0));
#else
			pntr += (NUM_GENERATOR + 1) * m_v;  // pointer to the first non-trivial generator
			// Check if it's an automorphism
			createGraphOut(pResGraph, pGraph[1], 0, 0, pntr);
			assert(!memcmp(pResGraph, pGraph[1], m_lenGraphMatr * sizeof(*pGraph[0])));

			// Multiply this automorphism by pInitOrbits
			for (int j = m_v; j--;)
				pResOrbits[j] = pInitOrbits[pntr[j]];

			for (int j = m_v; j--;)
				pInitOrbits[pResOrbits[j]] = j;

			createGraphOut(pResGraph, pGraph[1], 0, 0, pInitOrbits);
			PRINT_ADJ_MATRIX(pGraph[1], 0, m_v, pInitOrbits, "zzz");
			assert(!memcmp(pGraph[0], pGraph[1], m_lenGraphMatr * sizeof(*pGraph[0])));
			
			// Using ChatGPT suggestion https://chatgpt.com/share/68ab681a-aa0c-8010-a03d-7c9afb22af14
			// multiply q (saved in pResOrbits + m_v) by p^(-1) (which is pInitOrbits)
			// s=q∘p^−1
#define	A	1
#if A
			const auto q = ppp;
			const auto p = pInitOrbits;
#else
			const auto q = pInitOrbits;
			const auto p = ppp;
#endif
			for (int j = m_v; j--;)
				pResOrbits[p[j]] = j;   // pResOrbits is p^(-1)

			for (int j = m_v; j--;)
				s[j] = pResOrbits[q[j]];

			// Obtained permutation (pResOrbits) should be an automorphism of pGraph[0]
			// Let's check it
			createGraphOut(pGraph[0], pGraph[1], 0, 0, s);
			assert(!memcmp(pGraph[0], pGraph[1], m_lenGraphMatr * sizeof(*pGraph[0])));
			PRINT_ADJ_MATRIX(pGraph[1], 0, m_v, s, "vvv");
#endif
		}
#endif
	}

	const char* pGraphDescr = "";
	if (rank3)
		pGraphDescr = "rank 3 graph";
	else
	if (graphType == t_4_vert)
		pGraphDescr = "4-vertex condition";

	char buf[512], *pBuf = buf;
	if (graphType != t_regular)
		SPRINTFD(pBuf, buf, "Strongly regular graphs with parameters: (v,k,λ,μ) = (%d,%2d,%d,%d)",
			m_v, graphParam->k, graphParam->λ, graphParam->μ);
	else
		SPRINTFD(pBuf, buf, "Regular graphs with parameters: (v,k) = (%d,%2d)", m_v, graphParam->k);

	pGraphDB->setTableTitle(buf);
	if (rank3)
		graphParam->m_cntr[4]++;

	bool newGraph = true;
	int prevMatrNumb = m_nPrevMatrNumb;
	if (pUpperDiag) {
		prevMatrNumb = m_pMarixStorage[typeIdx]->numObjects();
		m_pMarixStorage[typeIdx]->updateRepo(pUpperDiag);
		newGraph = prevMatrNumb < m_pMarixStorage[typeIdx]->numObjects() ? 1 : 0;
	}
	else
		m_nPrevMatrNumb++;

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
		SrgSummary srgSummary(m_pParam->strVal[t_ResultFolder], param(t_out_CSV_file), m_pParam->strVal[t_CSV_FileName]);
		pBuf = buf;
#if OUT_SRG_TO_SEPARATE_FILE
		if (!prevMatrNumb) {
			if (graphType != t_regular)
				fprintf(f, "List of SRGs of type %d with parameters (v,k,λ,μ) = (%d,%2d,%d,%d):\n", typeIdx + 1, m_v, graphParam->k, graphParam->λ, graphParam->μ);
			else
				fprintf(f, "List of regular graphs of type %d with parameters (v,k) = (%d,%2d):\n", typeIdx + 1, m_v, graphParam->k);
		}

		SPRINTFD(pBuf, buf, "\nGraph #%d:  |Aut(G)| = %zd", prevMatrNumb + 1, groupOrder());
#else
#if PRINT_NUM_CUR_GRAPH
		SPRINTFD(pBuf, buf, "\n  *** numCurrGraph = %d ***", numCurrGraph);
#endif
		if (graphType != t_regular)
			SPRINTFD(pBuf, buf, "\nSRG #%d of type %d with parameters (v,k,λ,μ) = (%d,%2d,%d,%d): |Aut(G)| = %zd", 
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

		if (graphType != t_regular) {
			const auto rank3graph = rank3 ? 1 : (canonize ? -1 : 0);
			const auto grOrder = canonize ? groupOrder() : 0;
			srgSummary.outSRG_info(m_v, graphParam, graphType, rank3graph, grOrder, m_groupSize, m_nCols / m_groupSize, srcGroupOrder);
		}

		fprintf(f, "%s\n", buf);
		if (pResGraph)
			outAdjMatrix(pResGraph, f);

		FCLOSE_F(f);

		if (pResGraph && groupOrder() > 1)
			makeGroupOutput(NULL, false, false);
	}
	else {
		// Graph is isomorphic to a previously constructed one — adjust statistics to avoid duplicate counting.
		--graphParam->m_cntr[0];
		--graphParam->m_cntr[1];
		if (graphType != t_regular) 
			--graphParam->m_cntr[2];

		if (graphType == t_4_vert)
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

	// Initiate flags
	memset(m_bUsedFlags, 0, m_v * m_v);
}

void SRGToolkit::initVertexGroupOrbits() {
	setGroupOrder(1);
	setStabilizerLengthAut(m_v);
	for (int i = m_v; i--;)
		m_pGroupOrbits[i] = i;
	memcpy(this->getObject(1), m_pGroupOrbits, m_v * sizeof(m_pGroupOrbits[0]));
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
					sprintf_s(buffer, "ccc_%04d.txt", hhh += flg? 1 : 0);
					if (flg) {
						FOPEN_F(f, "bbb.txt", "a");
						fprintf(f, "ff = %2d  canonizeMatrixRow: hhh = %3d\n", ff, hhh);
						FCLOSE_F(f);
						printAdjMatrix(pGraphOut, buffer, i+1);
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

	int nCommon[10];
	bool flag;
	// Check if constructed graph is strongly regular
	const auto graphType = checkSRG(pGraph, graphDegree, nCommon, sizeof(nCommon), flag);
	switch (graphType) {
	case t_regular: 
		if (pGraphParam && !pGraphParam->m_cntr[1]++)
			pGraphParam->k = graphDegree;
	case t_complete:return graphType;

	default: // Check if complementary graph is a complete graph or a set of complete graphs 
		if (graphDegree == m_v - 1 || graphDegree == nCommon[1])
			return t_complete;
	}

#if USE_COMPLEMENT != 0
	if (USE_COMPLEMENT == 1 || 2 * graphDegree > m_v && graphDegree < m_v - 1) {
		graphDegree = m_v - 1 - graphDegree;
		// Compute the complement graph by inverting the adjacency relations.
		auto pVertex = pGraph;
		for (int i = 0; i < m_v; i++, pVertex += m_v) {
			for (int j = 0; j < m_v; j++)
				pVertex[j] = 1 - pVertex[j];

			pVertex[i] = 0;
		}
		checkSRG(pGraph, graphDegree, nCommon, sizeof(nCommon), flag);
		printfYellow_TGO("\ncomplement\n");
	}
	else
		printfYellow_TGO("\nNot complement (because of type)\n");
#else
	printfYellow_TGO("\nNot complement (n/a)\n");
#endif

	if (pGraphParam && !pGraphParam->m_cntr[1]++)
		pGraphParam->k = graphDegree;

	return pGraphParam->updateParam(nCommon, flag);
}

t_graphType SRGToolkit::checkSRG(const tchar *pGraph, int graphDegree, int * nCommon, size_t lenCommon, bool &flag) const {
	memset(nCommon, 0, lenCommon);
	flag = true;
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
						i, j, idx ? "mu" : "lambda", nCommonCurr, nCommon[idx], nCommon[4 * idx + 4], nCommon[4 * idx + 5]);
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
						i, j, idx ? "β" : "α", alpha, nCommon[idx + 2], nCommon[4 * idx + 4], nCommon[4 * idx + 5]);
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

	return t_srg;
}

t_graphType SRGParam::updateParam(int* pCommon, bool flag_4_ver) {
	if (this) {
		if (flag_4_ver) {
			if (α && α != pCommon[2] || β && β != pCommon[3])
				printfRed("Found graph with 4-vertex condition for different (α, β) = (%d, %d) != (%d, %d)\n",
					pCommon[2], pCommon[3], α, β);
			α = pCommon[2];
			β = pCommon[3];
			m_cntr[3]++; // # of graphs satisfying 4-vertex condition
		}

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

		unsigned int n4VertCond;
		int a = -1, b = -1;
		if (graphParam.m_cntr[4]) {
			printfYellow(" • %d - rank 3 graph%s\n", graphParam.m_cntr[4], (graphParam.m_cntr[4] > 1 ? "s" : ""));
			n4VertCond = graphParam.m_cntr[4] - graphParam.m_cntr[3];
			if (n4VertCond) {
				printfRed(" • %d graph%s satisf%s 4-vertex conditions: (α = %d, β = %d)\n",
					n4VertCond, (n4VertCond > 1 ? "s" : ""), (n4VertCond > 1 ? "y" : "ies"), graphParam.α, graphParam.β);
				a = graphParam.α;
				b = graphParam.β;
			}
		}
		else {
			n4VertCond = graphParam.m_cntr[3];
			if (n4VertCond) {
				printfYellow(" • %d graph%s satisf%s 4-vertex conditions: (α = %d, β = %d)\n",
					n4VertCond, (n4VertCond > 1 ? "s" : ""), (n4VertCond > 1 ? "y" : "ies"), graphParam.α, graphParam.β);
				a = graphParam.α;
				b = graphParam.β;
			}
		}
		printfYellow(" • %d %s strongly regular with parameters: (v, k, λ, μ) = (%d,%2d,%d,%d)\n",
			graphParam.m_cntr[2], pntr1, m_v, graphParam.k, graphParam.λ, graphParam.μ);
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
void Generators<ushort>::makeGroupOutput(const CRepository<ushort>* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	this->printTable(this->getObject(0), false, false, this->numObjects(), "");
}

template<>
void Generators<tchar>::makeGroupOutput(const CRepository<tchar>* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	createOrbits(pElemGroup);
	this->printTable(this->getObject(0), false, outToScreen, this->numObjects(), "");

	if (checkNestedGroups) {
		const auto retVal = testNestedGroups((const CGroupInfo*)pElemGroup);
		reportNestedGroupCheckResult(retVal, outToScreen);
	}

	m_bOrbitsCreated = false;
}