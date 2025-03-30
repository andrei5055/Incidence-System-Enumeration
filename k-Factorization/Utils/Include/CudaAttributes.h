#pragma once

#ifdef USE_CUDA
#include "cuda_runtime_api.h"
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>

#define CC __host__ __device__		// CUDA_CALLABLE
#define MIN2(x, y)              ((x) < (y)? (x) : (y))
#define MAX2(x, y)              ((x) > (y)? (x) : (y))
#define MEMCMP(s1, s2, n)       memcmp_gpu(s1, s2, n)
// Copying overlapping array                              
#define MEMMOVE(dest, src, len) { auto* pTmp = new tchar[len];  \
                                  memcpy(pTmp, src, len);       \
                                  memcpy(dest, pTmp, len);      \
                                  delete [] pTmp;               \
                                }
#else
#define CC
#define MIN2(x, y)              min(x, y)
#define MAX2(x, y)              max(x, y)
#define MEMCMP(s1, s2, n)       std::memcmp(s1, s2, n)
#define MEMMOVE(dest, src, len) memmove(dest, src, len)
#endif

#define CUDA_PRINTF(x, ...)
//#define CUDA_PRINTF(x, ...) printf(x, __VA_ARGS__)

#if USE_CUDA
#define ASSERT(c,...)
#define ASSERT_(c,...)
#else
#define ASSERT_FUNC         myAssert
#define ASSERT_(c,...)       if (c) { \
                                 __VA_ARGS__; \
                                 ASSERT_FUNC(c, __FILE__, __LINE__); \
                             }
#if NDEBUG
#define ASSERT(c,...)
#else

#define ASSERT              ASSERT_
#endif
#endif


#ifdef USE_CUDA
CC int memcmp_gpu(const void* s1, const void* s2, size_t n) {
    const unsigned char* p1 = (const unsigned char*)s1;
    const unsigned char* p2 = (const unsigned char*)s2;

    for (size_t i = 0; i < n; ++i) {
        if (p1[i] != p2[i]) {
            return (p1[i] < p2[i]) ? -1 : 1;
        }
    }
    return 0;
}
#endif
