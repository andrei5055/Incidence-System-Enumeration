#pragma once
#include <string>

#ifndef USE_CUDA
#ifdef UTILS_EXPORTS
#    define UTIL_LIBRARY __declspec(dllexport)
#else
#    define UTIL_LIBRARY __declspec(dllimport)
#endif

#ifdef LIBRARY_EXPORTS
#    define LIBRARY_API __declspec(dllexport)
#else
#    define LIBRARY_API __declspec(dllimport)
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
#define LIBRARY_API
#endif


#define RedText "\x1b[1;31m"
#define GreenText "\x1b[1;32m"
#define YellowText "\x1b[1;33m"
#define ResetTextColor "\x1b[0m"
#define printfRed(fmt, ...) printf(RedText fmt ResetTextColor, __VA_ARGS__)
#define printfGreen(fmt, ...) printf(GreenText fmt ResetTextColor, __VA_ARGS__)
#define printfYellow(fmt, ...) printf(YellowText fmt ResetTextColor, __VA_ARGS__)

typedef unsigned char tchar;
typedef const tchar  ctchar;

template<typename T>
CC T* reallocStorageMemory(T** pObjects, size_t lenObj) {
	auto* pNewObjMemory = new T[lenObj];
	memcpy(pNewObjMemory, *pObjects, lenObj >>= 1);
	delete[] * pObjects;
	return (*pObjects = pNewObjMemory) + lenObj;
}

UTIL_LIBRARY int readTable(const std::string& fn, int nRows, int nCols, tchar** pSm, int nmax, int reservedElement = 0, char infoSymb = '\"');
