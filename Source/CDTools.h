#pragma once

#ifdef LIBRARY_EXPORTS
#    define LIBRARY_API __declspec(dllexport)
#else
#    define LIBRARY_API __declspec(dllimport)
#endif

LIBRARY_API void* createCanonizer(int v, int lenGroup);
LIBRARY_API void releaseCanonizer(void *pCanonizer);
LIBRARY_API const unsigned char* runCanonizer(void* pCanonizer, const unsigned char* pMatrix, int k, int dayNumb = 0);