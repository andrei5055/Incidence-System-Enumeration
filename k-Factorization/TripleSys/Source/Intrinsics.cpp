#include "Intrinsics.h"

#if USE_INTRINSIC

// Highly optimized and unrolled AND operation
CC void multiplyAll(ll* __restrict a, const ll* __restrict b, const ll* __restrict c, int n) {
	int i = 0;
	// Process 16 elements (4 vectors) at a time for maximum throughput
	const int batch = 16;
	for (; i + batch <= n; i += batch) {
		const auto vb0 = _mm256_loadu_si256((__m256i const*)(b + i + 0));
		const auto vb1 = _mm256_loadu_si256((__m256i const*)(b + i + 4));
		const auto vb2 = _mm256_loadu_si256((__m256i const*)(b + i + 8));
		const auto vb3 = _mm256_loadu_si256((__m256i const*)(b + i + 12));
		const auto vc0 = _mm256_loadu_si256((__m256i const*)(c + i + 0));
		const auto vc1 = _mm256_loadu_si256((__m256i const*)(c + i + 4));
		const auto vc2 = _mm256_loadu_si256((__m256i const*)(c + i + 8));
		const auto vc3 = _mm256_loadu_si256((__m256i const*)(c + i + 12));

		_mm256_storeu_si256((__m256i*)(a + i + 0), _mm256_and_si256(vb0, vc0));
		_mm256_storeu_si256((__m256i*)(a + i + 4), _mm256_and_si256(vb1, vc1));
		_mm256_storeu_si256((__m256i*)(a + i + 8), _mm256_and_si256(vb2, vc2));
		_mm256_storeu_si256((__m256i*)(a + i + 12), _mm256_and_si256(vb3, vc3));
	}

	// Process remaining 4-element chunks (1 vector)
	for (; i + 4 <= n; i += 4) {
		const auto vb = _mm256_loadu_si256((__m256i const*)(b + i));
		const auto vc = _mm256_loadu_si256((__m256i const*)(c + i));
		_mm256_storeu_si256((__m256i*)(a + i), _mm256_and_si256(vb, vc));
	}

	// Fast scalar remainder
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
	int valid_lane_mask = _mm256_movemask_epi8(_mm256_loadu_si256((__m256i const*)BYTE_MASKS));

	// 3. Rank-Based Sort
	for (int i = 0; i < nPairs; ++i) {
		uint16_t cur_pair = input[i];
		__m256i cur_key_vec = _mm256_set1_epi16(cur_pair & 0x00FF);

		// REMOVED: No more XOR shifting needed for keys 0-31
		__m256i cmp = _mm256_cmpgt_epi16(cur_key_vec, keys);

		// Each match in 16-bit lane = 2 bits in mask
		const int mask = _mm256_movemask_epi8(cmp) & valid_lane_mask;

		// Rank = popcount / 2
		sorted_buf[_mm_popcnt_u32(mask) >> 1] = cur_pair;
	}
}

int sortGroupsAndCompare(ctchar* input_ptr, ctchar* target_ptr, tchar* sorted, int nPairs) {
	const uint16_t* input = (const uint16_t*)input_ptr;
	auto* sorted_buf = (uint16_t*)sorted;

	// 1. Single efficient unaligned load
	const auto v = _mm256_loadu_si256((const __m256i*)input);

	// 2. Extract Keys (Low Byte of each 16-bit lane)
	// We only need to do this once.
	const auto keys = _mm256_and_si256(v, _mm256_set1_epi16(0x00FF));

	// Get valid mask once
	const int valid_lane_mask = _mm256_movemask_epi8(_mm256_loadu_si256((__m256i const*)BYTE_MASKS));

	// 3. Rank-Based Sort
	for (int i = 0; i < nPairs; ++i) {
		uint16_t cur_pair = input[i];
		const auto cur_key_vec = _mm256_set1_epi16(cur_pair & 0x00FF);

		// REMOVED: No more XOR shifting needed for keys 0-31
		const auto cmp = _mm256_cmpgt_epi16(cur_key_vec, keys);

		// Each match in 16-bit lane = 2 bits in mask
		int mask = _mm256_movemask_epi8(cmp) & valid_lane_mask;

		// Rank = popcount / 2
		sorted_buf[_mm_popcnt_u32(mask) >> 1] = cur_pair;
	}

	// 4. Optimized Comparison
	const auto sorted_v = _mm256_loadu_si256((__m256i const*)sorted_buf);
	const auto target_v = _mm256_loadu_si256((__m256i const*)target_ptr);
	const auto eq = _mm256_cmpeq_epi8(sorted_v, target_v);

	int diff_mask = ~_mm256_movemask_epi8(eq) & valid_lane_mask;
	if (diff_mask == 0) return 0;

	int first_idx = _tzcnt_u32(diff_mask);
	uint8_t a = ((uint8_t*)sorted_buf)[first_idx];
	uint8_t b = target_ptr[first_idx];

	return (a < b) ? -1 : 1;
}
#endif

#ifdef USE_CUDA
#define BIT_SCAN_FORWARD64 BitScanForward64_Cuda
// Dummy function written by AI and never tested
CC void BitScanForward64_Cuda(unsigned long* _Index, unsigned long long _Mask) {
	unsigned long index = 0;
	if (_Mask == 0)
		return;
	while (!(_Mask & 1)) {
		_Mask >>= 1;
		index++;
	}
	*_Index = index;
}
#else
#define BIT_SCAN_FORWARD64  _BitScanForward64
#endif
CC bool findAndClearNextBit(uint& first, uint last, ll *pCompSol) {
	const auto lastB = IDX(last);
	auto firstB = first >> SHIFT;
#if USE_INTRINSIC
	// Process 4 words (256 bits) at a time using AVX2
	while (firstB + 4 <= lastB) {
		const auto v = _mm256_loadu_si256((__m256i const*)(pCompSol + firstB));
		if (!_mm256_testz_si256(v, v))
			break;
		firstB += 4;
	}
#endif
	while (firstB < lastB && !pCompSol[firstB])
		firstB++;

	if (firstB >= lastB)
		return false;

	unsigned long iBit;
	BIT_SCAN_FORWARD64(&iBit, (unsigned long long)pCompSol[firstB]);
	if ((first = (firstB << SHIFT) + iBit) >= last)
		return false;

	pCompSol[firstB] ^= (ll)1 << iBit;
	return true;
}
#if !USE_INTRINSIC
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

