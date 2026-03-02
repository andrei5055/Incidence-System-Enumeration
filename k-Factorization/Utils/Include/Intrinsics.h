#pragma once
#include "CudaSupport.h"
#include "k-SysSupport.h"

#ifdef USE_CUDA
#define USE_INTRINSIC	 0
#else
#define USE_INTRINSIC	 1
#endif

typedef long long ll;

CC void multiplyAll(ll* a, const ll* b, const ll* c, int n);

#if USE_INTRINSIC
#include <immintrin.h> // Header for AVX2 intrinsics
#include <algorithm>

#define NOT_ALL_4_ZEROS(pToA, j)  const __m256i fourValues = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pToA+j)); \
                                  !_mm256_testz_si256(fourValues, fourValues)
void initIntrinsics(int numPlayers);
#else
#define NOT_ALL_4_ZEROS(pToA, j)  pToA[j] || pToA[j + 1] || pToA[j + 2] || pToA[j + 3]
#define initIntrinsics(numPlayers)
#endif

