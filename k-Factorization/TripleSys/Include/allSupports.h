#pragma once
#include "k-SysSupport.h"

typedef signed long long sLongLong;

K_SYS_LIBRARY_API const char* getFileNameAttr(const kSysParam* param, const char** uf = NULL);
K_SYS_LIBRARY_API void initTripleSysData();
void convertLinksToResult(const char* ci, char* res, int np, int gs, bool cmpGraph);
void _printf(FILE* f, bool toScreen, const char* format, const char* pStr = NULL);
void printTableColor(char const* name, ctchar* c, int nl, int nc, int np, int ns = 0, bool makeString = false, ctchar* co=NULL, int* t=NULL);

template<typename T>
CC void printValue(tchar v, double scale) {
	if (v == unset)
		printfGreen(" -1");
	else
		printf("%3d", v);
}

template<typename T>
CC void printValue(int v, double scale) {
	if (scale > 0.0)
		printf(" %4.1f", v * scale);
	else
		printf(" %4d", v);
}

template<typename T>
CC void printValue(double v, double scale) {
	if ((v * scale) != 0)
		printfGreen(" %4.1f", v);
	else
		printf(" %4.1f", v);
}

template<typename T>
CC void printTable(char const* name, const T* c, int nl, int nc, int np, int ns = 0, bool makeString = false, double scale = 1.0) {
	if (name[0])
	{
		printf("%s:", name);
		if (nl > 1)
			printf("\n");
	}
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			const auto v = c[j * nc + i];
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");

			printValue<T>(v, scale);
			if (!makeString && ns == -1)
				printf(",");
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}

int compare_fn(const void* pA, const void* pB);

template<typename T>
CC void elemOrdering(T* pElems, size_t numElem, size_t groupSize) {
	// Ordering elements in the groups od size groupSize
	auto j = numElem + groupSize;
	switch (groupSize) {
	case 2:
		// Ordering groups of pairs.
		for (; j -= 2; pElems += 2) {
			if (pElems[0] > pElems[1]) {
				const auto tmp = pElems[0];
				pElems[0] = pElems[1];
				pElems[1] = tmp;
			}
		}
		return;
	case 3:
		// Ordering groups of triples.
		for (; j -= 3; pElems += 3) {
			const auto tmp0 = pElems[0];
			const auto tmp1 = pElems[1];
			const auto tmp2 = pElems[2];
			if (tmp2 > tmp1) {
				if (tmp0 > tmp1) {
					pElems[0] = tmp1;
					if (tmp2 < tmp0) {
						pElems[1] = tmp2;
						pElems[2] = tmp0;
					}
					else
						pElems[1] = tmp0;
				}
			}
			else {
				if (tmp2 > tmp0) {
					pElems[1] = tmp2;
					pElems[2] = tmp1;
				}
				else {
					pElems[0] = tmp2;
					if (tmp0 < tmp1) {
						pElems[1] = tmp0;
						pElems[2] = tmp1;
					}
					else
						pElems[2] = tmp0;
				}
			}
		}
		return;
	default:
		for (; j -= groupSize; pElems += groupSize)
			qsort(pElems, groupSize, sizeof(T), compare_fn);
	}

	return;
}

bool _CheckMatrix(ctchar* matrix, int nl, int nc, int gs, tchar* links, bool printError, int* errLine, int* errGroup, int* dubLine, int nr);
void printTransformed(int nrows, int ncols, int groupSize, ctchar* tr, ctchar* ttr, ctchar* pImatr, ctchar* pTmatr, int numRow = 0, sLongLong nLoops = 0, int finalKMindex = 0);
void createStartFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs);
CC bool p1fCheck2(ctchar* u1fCycles, ctchar* neigborsi, ctchar* neighborsj, int nc);
CC tchar checkForUnexpectedCycle(ctchar iv, ctchar ic, ctchar nc, ctchar* lnk, ctchar* v);
void myExit(int code);
void reportEOJ(int code);
void setConsoleOutputMode();
CC void kmTranslate(tchar* mo, ctchar* mi, ctchar* tr, int len);
CC bool kmTranslate2AndCheck(tchar* mo, ctchar* mi, ctchar* tr, int len, tchar tRow);
CC void getTT14ForG3(tchar* tt1, tchar* tt2, tchar* tt3, tchar* tt4, ctchar* v, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, int gn);
int getLS(unsigned char* ls, const int n, const int iRandomStart);

