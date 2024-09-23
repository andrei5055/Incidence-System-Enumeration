#pragma once

typedef signed long long sLongLong;

void convertLinksToResult(const char* ci, char* res, int np, int gs);
void _printf(FILE* f, bool toScreen, const char* format, const char* pStr = NULL);
void printTableColor(char const* name, const char* c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const char* c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const int* c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 0.0);
void printTable(char const* name, const double* c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 1.0);
bool _CheckMatrix(const char* matrix, int nl, int nc, int gs, char* links, bool printError, int* errLine, int* errGroup, int* dubLine);
int kmProcessMatrix2(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax);
int kmProcessMatrix3(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax);
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc);
int kmProcessMatrix(char* mo, char* mi, char* tmp, int nr, int nc, int gs, char* tm = NULL, int* pDayMax = NULL);
void kmSortGroupsByFirstValue(char* mo, char* mi, char nr, char nc, char np);
void kmSortRowsBy2ndValue(char* mo, char* mi, char nr, char nc, char* tm);
int factorial(int n);
bool setLinksForOnePlayer(int id, int np, char* lnk, char* p, int ip, char v);
void linksFromMatrix(char* lnk, char* iv, int nr, int np);
void printThreadsStat(sLongLong* cntTotal, sLongLong* cnt, int nMatrices, int nProcessed, int nrowsStart, int nrowsOut, int nThreads, clock_t iTime, bool bPrintSetup);
void printTransformed(int nrows, int ncols, const char* tr, const char* ttr, const char* pImatr, const char* pTmatr, int numRow = 0, sLongLong nLoops = 0, int finalKMindex = 0);
void createFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs, int iset);
void createStartFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs);
void saveStartData(char* fn, char* sm, int nm, int np, int nr, int gs);
int readStartData(char* fn, char* sm, int nm, int np, int nr, int gs);
void Stat(const char* t, int ind, bool bAdd);
void TestStatPrint(const char* hdr, int d);