#include "TripleSys.h"
#include "Table.h"

#ifndef USE_CUDA
#if USE_INTRINSIC
#include <immintrin.h>
#include <iostream>
#include <algorithm>

void transform16_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut);
void transform32_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut, int count);
bool find_2byte_sequence(short int* data, int nPairs, unsigned short* trSearch);
void sortGroupsI(tchar* pIn, tchar* pOut, int nPairs);
int sortGroupsAndCompare(ctchar* input, ctchar* target, tchar* sorted, int nPairs);

// Byte-level masks for _mm256_blendv_epi8 and final comparison
// 0xFF means "active pair", 0x00 means "padding/ignored"
alignas(32) uint8_t BYTE_MASKS[32];

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
CC void alldata::setReversSearchValues(tchar* tr, int nc, tchar tRow) const {
	tchar rs[2];
	int nrs = 0;
	for (int i = 0; i < nc; i++) {
		if (tr[i] == 0 || tr[i] == tRow) {
			rs[nrs++] = i;
			if (nrs >= 2)
				break;
		}
	}
	ASSERT_IF(nrs != 2);
	tr[nc] = tr[nc + 3] = rs[0];
	tr[nc + 1] = tr[nc + 2] = rs[1];
}

CC bool alldata::notCBMPtr(short int* tr, int nc) const {
	// return true if tr not CBMP tr.
	// Note: looks like this function never return true for p1f. Do we need it for p1f?
	const short int trMask = tr[0] & 0x0101;
	if (trMask != 0x0001 && trMask != 0x0100)
		return true;
	for (int i = 1; i < nc / 2; i++) {
		if ((tr[i] & 0x0101) != trMask)
			return true;
	}
	return false;
}
CC int alldata::cnvPrecalcRowsCompCheck(int& mode, ctchar* p1, ctchar* p1Neighbors, ctchar* p2, ctchar* p2Neighbors) const
{
	// p1, p2 - rows to calculate compatibilities with first 3 rows 
	// mode: = 0: p1 is new,
	//       > 0 && < m_numPlayers: p1 is the same as in prev call and mode is a length of solution to check.
	//       = m_numPlayers: p1 is the same as in prev call, but we can use only prev calculated tr's, not solution length
	alignas(32) tchar tr[(32 + 4 + 31) / 32 * 32]; 
	auto* pTestedTRs = testedTrs();
	const auto nRowsPrecalc = param(t_useRowsPrecalculation);
	const bool bCBMP = param(t_CBMP_Graph) == 2;
	const auto u1fPntr = sysParam()->u1fCycles[0];
	bool bP1f = ((!u1fPntr || u1fPntr[1] == m_numPlayers) && !m_allowUndefinedCycles);
	const tchar tRow = bCBMP ? 5 : 3; // check 3rd row even if param(t_useRowsPrecalculation) == 4 (not 3)

	ASSERT_IF(m_numPlayers > 32);
	
	m_playerIndex = m_numPlayers;
	if (mode == -1) {
		memset(tr, 0, sizeof(tr));
		pTestedTRs->resetGroupOrder();
		for (int i = 0; i < nRowsPrecalc; i++)
			u1fSetTableRow(neighborsPC(i), result(i), true);
		return 0;
	}
	int iRet = 0;
	const auto* neighbors0 = neighborsPC(0);
	const auto* neighbors1 = neighborsPC(1);

	if (mode > 0 && mode < m_numPlayers) {
		if (MEMCMP(prevP2(), p2, mode) == 0)
			return mode;
	}
	else if (mode == 0) {
		pTestedTRs->resetGroupOrder();
		// create 2 rows target: result(0,1) or result(1,0))
		for (int j = 0; j < 2; j++) {
			auto* pf0 = j == 0 ? neighbors0 : neighbors1;
			auto* pf1 = j == 0 ? neighbors1 : neighbors0;
			if (!bP1f) {
				if (create2U1FTr(pTestedTRs, tr, pf0, pf1, neighborsPC(), p1Neighbors, nRowsPrecalc, tRow))
					continue;
				ASSERT_IF(1);
				exit(116);
			}
			for (int k = 0; k < m_numPlayers; k++) {
				for (int i = 0; i < nRowsPrecalc; i++) {
					// create all tr's to convert 2 rows (result(i), p1) to target: (result(0,1) or to result(1,0))
					if (create2P1FTr(tr, k, pf0, pf1, neighborsPC(i), p1Neighbors)) {
						if (bCBMP && notCBMPtr((short int*)tr, m_numPlayers)) {
							exit(222);
							continue;
						}
						setReversSearchValues(tr, m_numPlayers, tRow);
						pTestedTRs->isProcessed(tr);
					}
					else {
						ASSERT_IF(1);
						exit(106);
					}
				}
			}
		}
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
	mode = m_numPlayers;

	// Calculate tr to convert p1 and p2 to the first two rows and check if this tr creates 3 rows matrix that less than original
	// The following code double the time needed to calculate compatibility but add less than 0.01% of non compatible rows
	if (param(t_useCompatibilityCheck) > 1 && bP1f) {
		for (int j = 0; j < 2; j++) {
			// create 2 rows target: result(0,1) or result(1,0))
			auto* pf0 = j == 0 ? neighbors0 : neighbors1;
			auto* pf1 = j == 0 ? neighbors1 : neighbors0;
			for (int k = 0; k < m_numPlayers; k++) {
				// create all tr's to convert 2 rows (p1, p2) to target: (result(0,1) or to result(1,0))
				if (create2P1FTr(tr, k, pf0, pf1, p1Neighbors, p2Neighbors)) {
					if (bCBMP && notCBMPtr((short int*)tr, m_numPlayers)) {
						exit(222);
						continue;
					}
					for (int i = 0; i < nRowsPrecalc; i++) {
						setReversSearchValues(tr, m_numPlayers, tRow);
						if (cnvCheckOneRow(tr, result(i), tRow, false)) { // if true : result(i) with applied tr is a third row and less than result(2)
							/**
							printTable("tr", tr, 1, m_numPlayers, 2);
							printTable("result", result(0), 3, m_numPlayers, 2);
							printTable("p1", p1, 1, m_numPlayers, 2);
							printTable("p2", p2, 1, m_numPlayers, 2);
							printTable("p3", m_Ktmp, 1, m_numPlayers, 2);
							**/
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
	}
Ret:
	if (iRet < 0) {
		memcpy(prevP2(), p2, m_numPlayers);
		return mode;
	}
	return 0;
}
bool alldata::cnvCheckOneRow(ctchar* tr, ctchar* pRow, ctchar tRow, bool bCalcLength) const {
#if USE_INTRINSIC
	alignas(32) tchar pOut[32];
	if (!find_2byte_sequence((short int*)pRow, m_numPlayers / 2, (unsigned short*)(tr + m_numPlayers)))
		return false;
	if (m_numPlayers == 16)
		transform16_and_sort_each_pair(pRow, tr, pOut);
	else
		transform32_and_sort_each_pair(pRow, tr, pOut, m_numPlayers);
	const auto iRet = sortGroupsAndCompare(pOut, result(2), m_Ktmp, m_nGroups);
	if (iRet != -1)
		return false;
#else
	const auto bRet = kmTranslate2AndCheck(m_Km, pRow, tr, m_numPlayers, tRow);
	if (!bRet)
		return false;
	kmSortGroupsByFirstValue(m_Km, m_Ktmp);
	if (MEMCMP(m_Ktmp, result(2), m_numPlayers) != -1)
		return false;
#endif
	if (bCalcLength) {
		const auto* r = result(2); // check 3rd row even if param(t_useRowsPrecalculation) == 4
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
void transform16_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut) {
	static const __m128i swap_mask = _mm_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14);
	static const __m128i merge_mask = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
	// 1. Load Data
	__m128i inp = _mm_loadu_si128((__m128i*)pIn);
	__m128i trTable = _mm_loadu_si128((__m128i*)pTr);

	// 2. Translate: pIn values 0-15 pick from trTable
	__m128i src = _mm_shuffle_epi8(trTable, inp);

	// 3. Internal Pair Sort: [min, max] for each consecutive pair
	__m128i swapped = _mm_shuffle_epi8(src, swap_mask);
	__m128i min_v = _mm_min_epu8(src, swapped);
	__m128i max_v = _mm_max_epu8(src, swapped);

	// Merge: Byte 0=min, Byte 1=max...
	__m128i result = _mm_blendv_epi8(min_v, max_v, merge_mask);

	// 4. Store current result
	_mm_storeu_si128((__m128i*)pOut, result);
}
// Function to find a 2-byte sequence in a N*16 - byte array using SIMD
// Returns true if found, or false if not found.
bool find_2byte_sequence(short int* data, int ng, unsigned short* trSearch) {
	if (ng <= 8) {
		static const __m128i masks[9] = {
			_mm_setr_epi16(0,0,0,0,0,0,0,0), // N=0
			_mm_setr_epi16(-1,0,0,0,0,0,0,0), // N=1
			_mm_setr_epi16(-1,-1,0,0,0,0,0,0), // N=2
			_mm_setr_epi16(-1,-1,-1,0,0,0,0,0), // N=3
			_mm_setr_epi16(-1,-1,-1,-1,0,0,0,0), // N=4
			_mm_setr_epi16(-1,-1,-1,-1,-1,0,0,0), // N=5
			_mm_setr_epi16(-1,-1,-1,-1,-1,-1,0,0), // N=6
			_mm_setr_epi16(-1,-1,-1,-1,-1,-1,-1,0), // N=7
			_mm_setr_epi16(-1,-1,-1,-1,-1,-1,-1,-1) // N=8
		};

		__m128i inp = _mm_loadu_si128((__m128i*)data);

		if (ng < 8) {
			inp = _mm_and_si128(inp, masks[ng]);
		}
		__m128i cmp = _mm_cmpeq_epi16(inp, _mm_set1_epi16(trSearch[0]));
		if (!_mm_testz_si128(cmp, cmp))
			return true;

		cmp = _mm_cmpeq_epi16(inp, _mm_set1_epi16(trSearch[1]));
		if (!_mm_testz_si128(cmp, cmp))
			return true;

		return false; // Sequence not found
	}
	ASSERT_IF(ng > 16);
	// 1. Create a 256-bit mask to handle lengths between 18-32 bytes (9-16 elements)
	// For performance, you can pre-align this or use a simple bitmask logic.
	// If ng is always > 8 and <= 16, we can often skip the mask if the buffer is padded.

	// Load the full 32-byte block
	__m256i inp = _mm256_loadu_si256((__m256i*)data);

	// 2. Optional: Masking if data is not padded to 32 bytes
	// If your memory is not 32-byte allocated, use a mask to avoid out-of-bounds reads.
	if (ng < 16) {
		// Create a mask with -1 (0xFFFF) for the first 'ng' shorts and 0 for the rest
		// In practice, for small fixed ranges, a lookup table like your masks[9] is fastest.
		alignas(32) static const short m_data[32] = {
			-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
			 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0
		};
		// Load mask shifted by (16 - ng) to get 'ng' consecutive -1s
		__m256i mask = _mm256_loadu_si256((__m256i*) & m_data[16 - ng]);
		inp = _mm256_and_si256(inp, mask);
	}

	// 3. Compare all 16 elements against trSearch[0] and trSearch[1]
	__m256i cmp0 = _mm256_cmpeq_epi16(inp, _mm256_set1_epi16(trSearch[0]));
	__m256i cmp1 = _mm256_cmpeq_epi16(inp, _mm256_set1_epi16(trSearch[1]));

	// 4. Combine results: Check if either trSearch[0] OR trSearch[1] exists
	__m256i combined = _mm256_or_si256(cmp0, cmp1);

	// 5. Test if any match bit is set
	// _mm256_testz returns 1 if the intersection is zero (no match found)
	return !_mm256_testz_si256(combined, combined);
}

/**
 * Applies a 32-byte transition (lookup table) to 32 bytes of input data.
 * @param pInp: Pointer to 32 bytes of input (indices 0-31).
 * @param pTr: Pointer to 32 bytes of transition values.
 * @param pOut: Pointer where the 32-byte result will be stored.
 * Warning: all arrays must have available 32 bytes
 
* @param pIn, pTr, pOut: Pointers to 32-byte buffers.
* @param count: Number of valid bytes to process (0-32).
*/
void transform32_and_sort_each_pair(ctchar* pIn, ctchar* pTr, tchar* pOut, int count) {
	// 1. Setup Dynamic Mask (0xFF for valid bytes, 0x00 for undefined)
	__m256i indices = _mm256_setr_epi8(
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
	__m256i v_count = _mm256_set1_epi8((char)count);
	__m256i valid_mask = _mm256_cmpgt_epi8(v_count, indices);

	// 2. Load and Clean Input
	// Any byte > 23 in your 24-byte case is forced to 0
	__m256i input = _mm256_and_si256(_mm256_loadu_si256((const __m256i*)pIn), valid_mask);
	__m256i trans = _mm256_loadu_si256((const __m256i*)pTr);

	// 3. Lane-Crossing Transformation (Look up any index 0-31)
	// Low lane (0-15) and High lane (16-31) of 'trans'
	__m256i trans_low = _mm256_permute2x128_si256(trans, trans, 0x00); // Broadcast low 128 to both lanes
	__m256i trans_high = _mm256_permute2x128_si256(trans, trans, 0x11); // Broadcast high 128 to both lanes

	// If input byte < 16, look in trans_low. If >= 16, look in trans_high (offset index by 16)
	__m256i res_low = _mm256_shuffle_epi8(trans_low, input);
	__m256i res_high = _mm256_shuffle_epi8(trans_high, _mm256_sub_epi8(input, _mm256_set1_epi8(16)));

	// Blend: use res_high only where input byte bit 4 (0x10) is set
	__m256i is_high = _mm256_cmpgt_epi8(input, _mm256_set1_epi8(15));
	__m256i transformed = _mm256_blendv_epi8(res_low, res_high, is_high);

	// 4. Byte-Pair Sorting [min, max]
	static const __m256i swap_mask = _mm256_setr_epi8(
		1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14,
		17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30);
	static const __m256i merge_mask = _mm256_setr_epi8(
		0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1,
		0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);

	__m256i swapped = _mm256_shuffle_epi8(transformed, swap_mask);
	__m256i min_v = _mm256_min_epu8(transformed, swapped);
	__m256i max_v = _mm256_max_epu8(transformed, swapped);

	__m256i result = _mm256_blendv_epi8(min_v, max_v, merge_mask);

	// 5. Cleanup result (optional: zero out 24-31 so they don't contain garbage)
	result = _mm256_and_si256(result, valid_mask);

	_mm256_storeu_si256((__m256i*)pOut, result);
}

// Precomputed globally: BYTE_MASKS has 2*n bytes of 0xFF
void  sortGroupsI(tchar* pIn, tchar* pOut, int nPairs) {
	const uint16_t* input = (const uint16_t*)pIn;
	uint16_t* sorted_buf = (uint16_t*)pOut;

	// 1. Single efficient unaligned load
	__m256i v = _mm256_loadu_si256((const __m256i*)input);

	// 2. Extract Keys (Low Byte of each 16-bit lane)
	// We only need to do this once.
	__m256i keys = _mm256_and_si256(v, _mm256_set1_epi16(0x00FF));

	// Get valid mask once
	int valid_lane_mask = _mm256_movemask_epi8(_mm256_loadu_si256((const __m256i*)BYTE_MASKS));

	// 3. Rank-Based Sort
	for (int i = 0; i < nPairs; ++i) {
		uint16_t cur_pair = input[i];
		__m256i cur_key_vec = _mm256_set1_epi16(cur_pair & 0x00FF);

		// REMOVED: No more XOR shifting needed for keys 0-31
		__m256i cmp = _mm256_cmpgt_epi16(cur_key_vec, keys);

		// Each match in 16-bit lane = 2 bits in mask
		int mask = _mm256_movemask_epi8(cmp) & valid_lane_mask;

		// Rank = popcount / 2
		sorted_buf[_mm_popcnt_u32(mask) >> 1] = cur_pair;
	}
}

// Precomputed globally: BYTE_MASKS has 2*n bytes of 0xFF
int sortGroupsAndCompare(ctchar* input_ptr, ctchar* target_ptr, tchar* sorted, int nPairs) {
	const uint16_t* input = (const uint16_t*)input_ptr;
	uint16_t* sorted_buf = (uint16_t*)sorted;

	// 1. Single efficient unaligned load
	__m256i v = _mm256_loadu_si256((const __m256i*)input);

	// 2. Extract Keys (Low Byte of each 16-bit lane)
	// We only need to do this once.
	__m256i keys = _mm256_and_si256(v, _mm256_set1_epi16(0x00FF));

	// Get valid mask once
	int valid_lane_mask = _mm256_movemask_epi8(_mm256_loadu_si256((const __m256i*)BYTE_MASKS));

	// 3. Rank-Based Sort
	for (int i = 0; i < nPairs; ++i) {
		uint16_t cur_pair = input[i];
		__m256i cur_key_vec = _mm256_set1_epi16(cur_pair & 0x00FF);

		// REMOVED: No more XOR shifting needed for keys 0-31
		__m256i cmp = _mm256_cmpgt_epi16(cur_key_vec, keys);

		// Each match in 16-bit lane = 2 bits in mask
		int mask = _mm256_movemask_epi8(cmp) & valid_lane_mask;

		// Rank = popcount / 2
		sorted_buf[_mm_popcnt_u32(mask) >> 1] = cur_pair;
	}

	// 4. Optimized Comparison
	__m256i sorted_v = _mm256_loadu_si256((const __m256i*)sorted_buf);
	__m256i target_v = _mm256_loadu_si256((const __m256i*)target_ptr);
	__m256i eq = _mm256_cmpeq_epi8(sorted_v, target_v);

	int diff_mask = ~_mm256_movemask_epi8(eq) & valid_lane_mask;
	if (diff_mask == 0) return 0;

	int first_idx = _tzcnt_u32(diff_mask);
	uint8_t a = ((uint8_t*)sorted_buf)[first_idx];
	uint8_t b = target_ptr[first_idx];

	return (a < b) ? -1 : 1;
}
#endif
#else
// For the version that uses GPU, the following function needs to be implemented.
CC int alldata::cnvPrecalcRowsCompCheck(int& mode, ctchar* p1, ctchar* p1Neighbors, ctchar* p2, ctchar* p2Neighbors) const {
	return 0;
}
#endif
