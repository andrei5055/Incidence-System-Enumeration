#pragma once

#ifndef USE_CUDA
#ifdef UTILS_EXPORTS
#    define UTIL_LIBRARY __declspec(dllexport)
#else
#    define UTIL_LIBRARY __declspec(dllimport)
#endif

#ifdef K_SYS_LIBRARY_EXPORTS
#    define K_SYS_LIBRARY_API __declspec(dllexport)
#else
#    define K_SYS_LIBRARY_API __declspec(dllimport)
#endif

#define CC
#define CK
#else
#define CC __host__ __device__		// CUDA_CALLABLE
#if CONSTR_ON_GPU
#define CK	CC __noinline__
#else
#define CK
#endif
#define UTIL_LIBRARY
#define K_SYS_LIBRARY_API
#endif
