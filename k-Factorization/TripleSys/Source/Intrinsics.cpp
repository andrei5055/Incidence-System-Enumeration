#include "Intrinsics.h"

#if USE_INTRINSIC

void multiplyAll_(ll* __restrict a,
	const ll* __restrict b,
	const ll* __restrict c,
	int n)
{
	int i = 0;

	const int vec = 4;          // 4 x 64-bit per AVX2 register
	const int step = 4 * vec;   // unroll by 4 vectors (16 elements)

	for (; i + step <= n; i += step) {
		__m256i b0 = _mm256_loadu_si256((__m256i const*)(b + i + 0));
		__m256i c0 = _mm256_loadu_si256((__m256i const*)(c + i + 0));
		__m256i b1 = _mm256_loadu_si256((__m256i const*)(b + i + 4));
		__m256i c1 = _mm256_loadu_si256((__m256i const*)(c + i + 4));
		__m256i b2 = _mm256_loadu_si256((__m256i const*)(b + i + 8));
		__m256i c2 = _mm256_loadu_si256((__m256i const*)(c + i + 8));
		__m256i b3 = _mm256_loadu_si256((__m256i const*)(b + i + 12));
		__m256i c3 = _mm256_loadu_si256((__m256i const*)(c + i + 12));

		_mm256_storeu_si256((__m256i*)(a + i + 0), _mm256_and_si256(b0, c0));
		_mm256_storeu_si256((__m256i*)(a + i + 4), _mm256_and_si256(b1, c1));
		_mm256_storeu_si256((__m256i*)(a + i + 8), _mm256_and_si256(b2, c2));
		_mm256_storeu_si256((__m256i*)(a + i + 12), _mm256_and_si256(b3, c3));
	}

	for (; i + vec <= n; i += vec) {
		__m256i vb = _mm256_loadu_si256((__m256i const*)(b + i));
		__m256i vc = _mm256_loadu_si256((__m256i const*)(c + i));
		_mm256_storeu_si256((__m256i*)(a + i),
			_mm256_and_si256(vb, vc));
	}

	for (; i < n; ++i)
		a[i] = b[i] & c[i];
}

void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	int i = 0;
	// Unroll to process 8 elements (2 vectors) at a time
	for (; i + 8 <= n; i += 8) {
		__m256i v_b1 = _mm256_loadu_si256((__m256i*) & b[i]);
		__m256i v_c1 = _mm256_loadu_si256((__m256i*) & c[i]);
		_mm256_storeu_si256((__m256i*) & a[i], _mm256_and_si256(v_b1, v_c1));
		__m256i v_b2 = _mm256_loadu_si256((__m256i*) & b[i + 4]);
		__m256i v_c2 = _mm256_loadu_si256((__m256i*) & c[i + 4]);
		_mm256_storeu_si256((__m256i*) & a[i + 4], _mm256_and_si256(v_b2, v_c2));
	}

	// Process remaining 4-element chunks if they exist
	for (; i + 4 <= n; i += 4) {
		__m256i v_b = _mm256_loadu_si256((__m256i*) & b[i]);
		__m256i v_c = _mm256_loadu_si256((__m256i*) & c[i]);
		_mm256_storeu_si256((__m256i*) & a[i], _mm256_and_si256(v_b, v_c));
	}

	// Scalar remainder loop
	for (; i < n; ++i) {
		a[i] = b[i] & c[i];
	}
}


// Byte-level masks for _mm256_blendv_epi8 and final comparison
// 0xFF means "active pair", 0x00 means "padding/ignored"
alignas(32) uint8_t BYTE_MASKS[32];

void initIntrinsics(int numPlayers) {
	char* mask = (char*)BYTE_MASKS;
	if (numPlayers <= 32) {
		memset(mask, 255, numPlayers);
		memset(mask + numPlayers, 0, 32 - numPlayers);
	}
}

// Precomputed globally: BYTE_MASKS has 2*n bytes of 0xFF
void sortGroupsI(tchar* pIn, tchar* pOut, int nPairs) {
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


#else
CC void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	const int n8max = n & (~7);
	int j = 0;
	for (; j < n8max; j += 8) {
		a[j] = b[j] & c[j];
		a[j + 1] = b[j + 1] & c[j + 1];
		a[j + 2] = b[j + 2] & c[j + 2];
		a[j + 3] = b[j + 3] & c[j + 3];
		a[j + 4] = b[j + 4] & c[j + 4];
		a[j + 5] = b[j + 5] & c[j + 5];
		a[j + 6] = b[j + 6] & c[j + 6];
		a[j + 7] = b[j + 7] & c[j + 7];
	}
	for (; j < n; j++)
		a[j] = b[j] & c[j];
}
#endif

