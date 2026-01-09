#include "TripleSys.h"
#include "Table.h"

#define USE_INTRINSIC		1
#if USE_INTRINSIC
#include <immintrin.h>
#include <iostream>
#include <algorithm>

bool transform_and_find_pair(tchar* pOut, ctchar* pIn, ctchar* pTr); //, tchar v1, tchar v2);
__m128i swap_mask = _mm_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14);
__m128i merge_mask = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
__m128i search_vec;
#endif

#define FIND_ERROR 0
#if FIND_ERROR && _DEBUG
// Function to find error in cnvPrecalcRowsCompCheck defining non-compatible solutions
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

CC int alldata::cnvPrecalcRowsCompCheck(int& mode, ctchar* p1, ctchar* p1Neighbors, ctchar* p2, ctchar* p2Neighbors) const
{
	//if (mode == 0 && (p1[1] > 7 || p2[1] > 7))
	//	return 0; //leo

	// p1, p2 - rows to calculate compatibilities with first precalculate rows rows
	// mode: = 0 - p1 is new,
	//       > 0 - p1 is the same as in prev call and mode is a length of solution to check.
	const auto nRowsPrecalc = param(t_useRowsPrecalculation);
	const tchar tRow = param(t_CBMP_Graph) == 2 ? 5 : 3; // check 3rd row even if param(t_useRowsPrecalculation) == 4 (not 3)
	/**
	tchar out[16];
	tchar inp[16] = { 0,1,2,3,4,5,6,7,8,9,10,15,14,13,12,11 };
	tchar ttr[16] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0 };
	printTable("out0", out, 1, m_numPlayers, 2);
	for (int i = 0; i < 1000000000; i++)
		kmTranslate2AndCheck(out, inp, ttr, m_numPlayers, 11);
	printTable("out1", out, 1, m_numPlayers, 2);
	bool ret = false;
	for (int i = 0; i < 1000000000; i++)
		ret = transform_and_find_pair(out, inp, ttr, 0, 11);
	
	printf("ret=%s\n", ret ? "pair exists" : "not exists");
	printTable("out2", out, 1, m_numPlayers, 2);
	exit(1);
	**/
	m_playerIndex = m_numPlayers;
	if (mode == -1) {
#if USE_INTRINSIC
		ctchar v1 = 0;
		uint16_t target = (static_cast<uint16_t>(tRow) << 8) | v1;
		search_vec = _mm_set1_epi16(target);
#endif
		for (int i = 0; i < nRowsPrecalc; i++)
			u1fSetTableRow(neighborsPC(i), result(i), true);
		return 0;
	}

	tchar tr[MAX_PLAYER_NUMBER];
	int iRet = 0;
	auto* pTestedTRs = testedTrs();
	const auto* neighbors0 = neighborsPC(0);
	const auto* neighbors1 = neighborsPC(1);

	if (!mode) {
		pTestedTRs->resetGroupOrder();

		// create all tr to convert (p1, p2) to result(0,1) and to result(1,0))
		for (int j = 0; j < 2; j++) {
			auto* pf0 = j == 0 ? neighbors0 : neighbors1;
			auto* pf1 = j == 0 ? neighbors1 : neighbors0;
			for (int k = 0; k < m_numPlayers; k++) {
				for (int i = 0; i < nRowsPrecalc; i++) {
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
		if (cnvCheckOneRow(trt, p2, tRow, true)) { // if true : p2 with applied trt is a third row and less than result(2)
			iRet = -1;
			mode = m_playerIndex;
			if (m_doNotExitEarlyIfNotCanonical)
				continue;
			goto Ret;
		}
	}
	/**/
	mode = m_numPlayers;
	for (int j = 0; j < 2; j++) {
		auto* pf0 = j == 0 ? neighbors0 : neighbors1;
		auto* pf1 = j == 0 ? neighbors1 : neighbors0;
		for (int k = 0; k < m_numPlayers; k++) {
			if (create2P1FTr(tr, k, pf0, pf1, p1Neighbors, p2Neighbors)) {
				for (int i = 0; i < nRowsPrecalc; i++) {
					if (cnvCheckOneRow(tr, result(i), tRow, false)) { // if true : result(i) with applied tr is a third row and less than result(2)
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
	/**/
Ret:
	if (iRet < 0) {
		memcpy(prevP2(), p2, m_numPlayers);
		return mode;
	}
	return 0;
}
bool alldata::cnvCheckOneRow(ctchar* tr, ctchar* pRow, ctchar tRow, bool bCalcLength) const {
#if 1
	bool ret;
#if USE_INTRINSIC
	if (m_numPlayers == 16)
		ret = transform_and_find_pair(m_Km, pRow, tr);
	else
#endif
		ret = kmTranslate2AndCheck(m_Km, pRow, tr, m_numPlayers, tRow);
	if (!ret)
		return false;
	kmSortGroupsByFirstValue(m_Km, m_Ktmp);
#else
	kmTranslate(m_Km, pRow, tr, m_numPlayers);
	(this->*m_pSortGroups)(m_Km, m_numPlayers);
	kmSortGroupsByFirstValue(m_Km, m_Ktmp);
	if (m_Ktmp[1] != tRow)
		return false;
#endif
	const auto* r = result(2); // check 3rd row even if param(t_useRowsPrecalculation) == 4
	if (MEMCMP(m_Ktmp, r, m_numPlayers) != -1)
		return false;

	if (bCalcLength) {
		int minLength = 0;
		for (minLength = 0; minLength < m_numPlayers; minLength++) {
			if (m_Ktmp[minLength] < r[minLength])
				break;
		}
		setPlayerIndexByPos(tr, m_Ktmp, pRow, 0, minLength);
	}

	return true;
}
#if USE_INTRINSIC
// Function to translate 16 bytes and sort each internal pair
// Returns true if specific pair (v1, v2) is found among the 8 pairs (we use it to verify if this is expected row)
bool transform_and_find_pair(tchar* pOut, ctchar* pIn, ctchar* pTr) { //, tchar v1, tchar v2) {
	// 1. Load Data
	__m128i inp = _mm_loadu_si128((__m128i*)pIn);
	__m128i trTable = _mm_loadu_si128((__m128i*)pTr);

	// 2. Translate: pIn values 0-15 pick from trTable
	__m128i src = _mm_shuffle_epi8(trTable, inp);

	// 3. Internal Pair Sort: [min, max] for each consecutive pair
	//__m128i swap_mask = _mm_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14);
	__m128i swapped = _mm_shuffle_epi8(src, swap_mask);
	__m128i min_v = _mm_min_epu8(src, swapped);
	__m128i max_v = _mm_max_epu8(src, swapped);

	// Merge: Byte 0=min, Byte 1=max...
	__m128i result = _mm_blendv_epi8(min_v, max_v, merge_mask);
	//	_mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1));

	// 4. Store current result
	_mm_storeu_si128((__m128i*)pOut, result);
	/**
	// 5. Search for specific pair (v1, v2)
	// We treat the target as a 16-bit word: (v2 << 8) | v1
	uint16_t target = (static_cast<uint16_t>(v2) << 8) | v1;
	__m128i search_vec = _mm_set1_epi16(target);
	**/
	// Compare each of the 8 pairs simultaneously
	__m128i cmp = _mm_cmpeq_epi16(result, search_vec);

	// Movemask returns bits; if any pair matched, mask will be non-zero
	return (_mm_movemask_epi8(cmp) != 0);
}

#endif
