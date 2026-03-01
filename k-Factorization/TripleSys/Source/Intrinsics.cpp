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

