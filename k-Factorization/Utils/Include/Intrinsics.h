#pragma once
#include "CudaSupport.h"
#include "k-SysSupport.h"

#ifdef USE_CUDA
#define USE_INTRINSIC	 0
#define USE_64_BIT_MASK	 0
#else
#define USE_INTRINSIC	 1
#define USE_64_BIT_MASK	 1
#endif

typedef long long ll;

#if USE_64_BIT_MASK
typedef ll tmask;
#define SHIFT						6
#else
typedef tchar tmask;
#define SHIFT						3
#endif

#define MASK_BIT(idx)				((tmask)1 << ((idx) & ((1<<SHIFT) - 1)))	
#define IDX(n)						((n + (1<<SHIFT) - 1) >> SHIFT)
#define REM(n)						(n % ((tmask)1<<SHIFT))			// remainder from division
#define SET_MASK_BIT(mask, idx)		(mask)[(idx) >> SHIFT] |= MASK_BIT(idx)
#define RESET_MASK_BIT(mask, idx)	(mask)[(idx) >> SHIFT] ^= MASK_BIT(idx)
#define CHECK_MASK_BIT(mask, idx)	((mask)[(idx) >> SHIFT] & MASK_BIT(idx))
#define AVALABLE_PLAYER_MASK_LENGTH  8

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

CC bool findAndClearNextBit(uint& first, uint last, ll* pCompSol);

