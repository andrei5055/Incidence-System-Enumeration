#pragma once
#include "k-SysSupport.h"

LIBRARY_API void* createCanonizer(int v, int lenGroup);
LIBRARY_API void releaseCanonizer(void *pCanonizer);
LIBRARY_API ctchar* runCanonizer(void* pCanonizer, ctchar *pMatrix, int k);