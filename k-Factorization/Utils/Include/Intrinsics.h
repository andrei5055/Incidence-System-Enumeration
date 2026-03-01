#pragma once
#include "CudaSupport.h"

#define USE_INTRINSIC		!USE_CUDA

typedef long long ll;

CC void multiplyAll(ll* a, const ll* b, const ll* c, int n);

#if USE_INTRINSIC
#include <immintrin.h> // Header for AVX2 intrinsics
#include <algorithm>

#define NOT_ALL_4_ZEROS(pToA, j)  const __m256i fourValues = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pToA+j)); \
                                  !_mm256_testz_si256(fourValues, fourValues)
#else
#define NOT_ALL_4_ZEROS(pToA, j)  pToA[j] || pToA[j + 1] || pToA[j + 2] || pToA[j + 3]
#endif
