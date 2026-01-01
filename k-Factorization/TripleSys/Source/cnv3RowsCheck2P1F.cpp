#include "TripleSys.h"
#include "Table.h"

#define FIND_ERROR 0
#if FIND_ERROR && _DEBUG
// Function to find error in cnv3RowsCheck2P1F defining non-compatible solutions
// between pre-calculated ones. to be 
bool findError(ctchar* sol1, ctchar* sol2, int lenSol, bool assert = false)
{
	// Put here the last rows of non-constructed matrix
	tchar sol[] = {
		0,  4,   1, 12,   2, 13,   3,  9,   5, 10,   6,  8,   7, 11,
		0,  5,   1,  6,   2,  9,   3, 12,   4, 11,   7, 13,   8, 10,
		0,  6,   1, 10,   2, 12,   3, 11,   4, 13,   5,  9,   7,  8,
		0,  7,   1, 11,   2,  4,   3, 13,   5,  6,   8, 12,   9, 10,
		0,  8,   1,  3,   2,  6,   4,  7,   5, 11,   9, 13,  10, 12,
		0,  9,   1,  2,   3,  5,   4,  6,   7, 12,   8, 11,  10, 13,
		0, 10,   1,  9,   2, 11,   3,  7,   4,  8,   5, 13,   6, 12,
		0, 11,   1, 13,   2,  5,   3,  8,   4, 12,   6, 10,   7,  9,
		0, 12,   1,  8,   2, 10,   3,  4,   5,  7,   6, 13,   9, 11,
		0, 13,   1,  7,   2,  8,   3, 10,   4,  9,   5, 12,   6, 11,
	};

	for (int i = 0; i < sizeof(sol) / lenSol; i++) {
		const auto ret = MEMCMP(sol + i * lenSol, sol1, lenSol);
		if (ret > 0)
			return true;

		if (ret == 0) {
			// Found first solution, now search for the second one
			for (int j = i + 1; j < sizeof(sol) / lenSol; j++) {
				const auto ret2 = MEMCMP(sol + j * lenSol, sol2, lenSol);
				if (ret2 > 0)
					return true;

				if (ret2 == 0) {
					printf("Error found at rows %d and %d\n", i, j);
					if (assert) {
						ASSERT_IF(1);
					}
					else
						return false;
				}
			}
		}
	}
	return true;
}
#else
#define findError(sol1, sol2, lenSol, assert)
#endif

CC int alldata::cnv3RowsCheck2P1F(ctchar* p1, ctchar* p1Neighbors, ctchar* p2, ctchar* p2Neighbors, int& mode) const
{
	//if (mode == 0 && (p1[1] > 7 || p2[1] > 7))
	//	return 0; //leo

	// p1, p2 - rows to calculate compatibilities with first 3 rows
	// mode: = 0 - p1 is new,
	//       > 0 - p1 is the same as in prev call and mode is a length of solution to check.

	m_playerIndex = m_numPlayers;
	if (mode == -1) {
		for (int i = 0; i < 3; i++)
			u1fSetTableRow(neighborsPC(i), result(i), true);
		return 0;
	}

	tchar tr[MAX_PLAYER_NUMBER];
	int iRet = 0;
	auto* pTestedTRs = testedTrs();
	const auto* neighbor0 = neighborsPC(0);
	const auto* neighbor1 = neighborsPC(1);
	const auto* neighbor2 = neighborsPC(2);

	if (!mode) {
		pTestedTRs->resetGroupOrder();

		// create all tr to convert (p1, p2) to result(0,1) and to result(1,0))
		for (int j = 0; j < 2; j++) {
			auto* pf0 = j == 0 ? neighbor0 : neighbor1;
			auto* pf1 = j == 0 ? neighbor1 : neighbor0;
			for (int k = 0; k < m_numPlayers; k++) {
				for (int i = 0; i < 3; i++) {
					if (create2P1FTr(tr, k, pf0, pf1, neighborsPC(i), p1Neighbors)) {
						pTestedTRs->isProcessed(tr);
					} else {
						ASSERT_IF(1);
						exit(106);
					}
				}
			}
		}
	} else {
		if (MEMCMP(prevP2(), p2, mode) == 0)
			return mode;
	}

	const int nTrs = pTestedTRs->numObjects();
	for (int itr = 0; itr < nTrs; itr++) {
		tchar* trt = pTestedTRs->getObject(itr);
		if (cnvCheckOneRow(trt, p2, true)) { // if true : p2 with applied trt is a third row and less than result(2)
			iRet = -1;
			mode = m_playerIndex;
			if (m_doNotExitEarlyIfNotCanonical)
				continue;
			goto Ret;
		}
	}

	mode = m_numPlayers;
	for (int j = 0; j < 2; j++) {
		auto* pf0 = j == 0 ? neighbor0 : neighbor1;
		auto* pf1 = j == 0 ? neighbor1 : neighbor0;
		for (int k = 0; k < m_numPlayers; k++) {
			if (create2P1FTr(tr, k, pf0, pf1, p1Neighbors, p2Neighbors)) {
				for (int i = 0; i < 3; i++) {
					if (cnvCheckOneRow(tr, result(i), false)) { // if true : result(i) with applied tr is a third row and less than result(2)
						findError(p1, p2, m_numPlayers, false);
						iRet = -1;
						if (m_doNotExitEarlyIfNotCanonical)
							continue;
						goto Ret;
					}
				}
			}
			else {
				ASSERT_IF(1);
				exit(107);
			}
		}
	}
Ret:
#if 0	
	for (int i = 0; i < 3; i++)
		u1fSetTableRow(neighborsPC(i), result(i));
#endif
	if (iRet < 0) {
		memcpy(prevP2(), p2, m_numPlayers);
		return mode;
	}
	return 0;
}
bool alldata::cnvCheckOneRow(ctchar* tr, ctchar* pRow, bool bCalcLength) const {
	tchar tRow = param(t_CBMP_Graph) == 2 ? 5 : 3;
#if 1
	if (!kmTranslate2AndCheck(m_Km, pRow, tr, m_numPlayers, tRow))
		return false;
	kmSortGroupsByFirstValue(m_Km, m_Ktmp);
#else
	kmTranslate(m_Km, pRow, tr, m_numPlayers);
	(this->*m_pSortGroups)(m_Km, m_numPlayers);
	kmSortGroupsByFirstValue(m_Km, m_Ktmp);
	if (m_Ktmp[1] != tRow)
		return false;
#endif
	const auto* r = result(2);
	if (MEMCMP(m_Ktmp, r, m_numPlayers) != -1)
		return false;

	if (bCalcLength) {
		int minLength = 0;
		for (minLength = 0; minLength < m_numPlayers; minLength++) {
			if (m_Ktmp[minLength] < r[minLength])
				break;
		}
		setPlayerIndexByPos(tr, m_Ktmp, pRow, 0, minLength);
		/**
		if (minLength > m_numPlayers - 3) { // leo
			return -1;
		}
		tchar ttr1[MAX_PLAYER_NUMBER];
		tchar ttr2[MAX_PLAYER_NUMBER];
		for (int i = 0; i < m_numPlayers; i++) {
			ttr1[m_Ktmp[i]] = i;
		}
		int mxLength = 0;
		for (int i = 0; i < m_numPlayers; i++) {
			if (ttr[pRow[i]])
		}
		*/
	}
	//printTable("r3", result(2), 1, 18, 2);
	//printTable("n3", m_Ktmp, 1, 18, 2);
	//setPlayerIndex(tr, 0, 0, m_Ktmp, m_Km, pRow);

	return true;
}
