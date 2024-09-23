#pragma once
#include <stdio.h>

#ifdef USE_CUDA
#define CONSTR_ON_GPU				0						// 1 - Start using GPU for object construction
#define USE_THREADS					1						// Should be at least 1
#define CANON_ON_GPU                (1 && CONSTR_ON_GPU==0)	// 1 - Start using GPU for canonicity testing
#define NUM_GPU_WORKERS             128
//#include "host_defines.h"
#include "cuda_runtime_api.h"
#define CC __host__ __device__		// CUDA_CALLABLE
#if CONSTR_ON_GPU
#define CK	CC __noinline__
#else
#define CK
#endif
#else
#define USE_THREADS					15 //(TEST? 0 : 15)	// default numer if thread used, If NOT 0, it can be chanaged by THREAD_NUMBER
#define CONSTR_ON_GPU				0
#define CANON_ON_GPU				0
#define NUM_GPU_WORKERS				0
#define CC
#define CK
#endif

typedef unsigned int		uint;
typedef unsigned short		ushort;
typedef unsigned char		uchar;
typedef unsigned long long	ulonglong;

#define countof(x)  (sizeof(x)/sizeof(x[0]))

#ifdef WIN
#define OPEN_FILE(x, y, z)	 fopen_s(&x, y, z)
#define FOPEN(x, y, z)	  	 FILE *x = NULL; if (y && strlen(y)) fopen_s(&x, y, z)
#else
#define sprintf_s(x, y, ...) sprintf(x, __VA_ARGS__)
#define strcpy_s(x, y, z)    strcpy(x, z)
#define memcpy_s(x, y, ...)  memcpy(x, __VA_ARGS__)
#define OPEN_FILE(x, y, z)	 x = fopen(y, z)
#define FOPEN(x, y, z)	  	 FILE *OPEN_FILE(x, y, z)
#endif

#define FCLOSE(file)			if (file) fclose(file)

#define SNPRINTF(x, len, ...)	static_cast<size_t>(snprintf(x, len, __VA_ARGS__))
#define SPRINTF(x, ...)			SNPRINTF(x, sizeof(x), __VA_ARGS__)

