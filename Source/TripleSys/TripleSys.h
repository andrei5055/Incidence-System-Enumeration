#pragma once
#include <string> 
#include <iostream>

#define nPlayers  10
#define GroupSize 2

#define UseMultiThreading false

#if UseMultiThreading
#define NThreads 12				  // number of threads to calculate matrices		
#define NRowsInStartMatrix 4	  // number of rows in initial set of matrices for threads
#define NRowsInResultMatrix 0	  // number of rows in result, if 0 then calculate full matrix
#define MaxNumberOfStartMatrices 1000000
#else
#define NThreads 1
#define NRowsInStartMatrix 0
#define NRowsInResultMatrix 0
#define MaxNumberOfStartMatrices 1
#endif

#define USE_cnvCheckNew 1		  // 1 - cnvCheckNew (no arrays), 0 - cnvCheck (use two precalculated arrays)

#define LoopsMax 200000000000.
#define ImproveResults 0		  // 0 - with no explanation of the matrix rejection reason;
                                  // 1 - reason of rejection will be explained;
								  // 2 - improve matrix as much as possible.
#define CreateImprovedResult 0    // 2 - cnvCheckNew create improved matrix in m_Km, 0 - not (0works faster) 
#define StartFolder				  "../Logs1/" 	// Name of folder with 'Start Matrices, must present in multithread mode.
#define ResultFolder			  "../Logs2/" 	// Name of folder with the results, "" - no file output.
#define ResultNameFormat		  "R%010d.txt"  // Format for result file name
#define PrintImprovedResults 0	  // Set this value to >= 2 if you want to see improved results on screen.
#define ImprovedResultFolder	  "../ImprovedResults/"	// Name of folder for files with the improved results.
#define ImprovedResultNameFormat  "I%010d.txt"	// Format of the output file name with the improved results.
#define UsePos_1_4_condition 1
#define UseCheckLinksV 1
#define UseCheckLinksH 1
#define ReportInterval 120000
#define StartMatricesFolder			  "../StartMatrices/" 	// Name of folder to save/Load Start Matrices, "" - do not use.


// Conditions to be tested on day 1:
#define USE_STATEMENT_7   1    // players[4] < players[7]  
#define USE_STATEMENT_17  1    // Only the players 4 or 9 could be in position[1, 4].
#define USE_STATEMENT_18  1    // The positions of any two players who were in the same group 
                               // on the first day on the second day must be in ascending order of their numbers.
#define USE_CHANGING_DAY_0_GROUPS 0 // Use permutation of first 3 groups of day 0
#define CHECK_WITH_GROUP  1    // Use new group on completely constructed matrices
#define USE_TRANSLATE_BY_LEO 0
#define USE_ORBTS		  1
#define USE_EQUAL		  1
#define PRINT_MATR_CNTR	  1
#define PRINT_TRANSFORMED 0

#define USE_STATEMENT_19  1
#define CHECK_PERMUTATIONS	1
#define IMPROVE_RESULT_MAX  1000 // Maximal number of attempts to improve matrix 
								 // (works only with ImproveResults=2)
#define OUTPUT_VECTOR_STAT	0    // Output of some statistics for constructed matrices

#define PrintLinksStat 1
#define PrintNVminmax 0
#define UseSS 0
#define nGroups (nPlayers / GroupSize)
#define NDays ((nPlayers - 1) / (GroupSize - 1))
#define unset ((char)(-1))
#define printfRed(fmt, v) printf("\x1b[1;31m" fmt "\x1b[0m", v)
#define printfGreen(fmt, v) printf("\x1b[1;32m" fmt "\x1b[0m", v)
#define printfYellow(fmt, v) printf("\x1b[1;33m" fmt "\x1b[0m", v)


#define FOPEN_W(x, y, z, w)	 FILE *x = w; \
                             if (!x && y && strlen(y)) fopen_s(&x, y, z)
#define FOPEN_F(x, y, z)	 FOPEN_W(x, y, z, NULL)
#define SPRINTFS(x, y, size,...)	 x += sprintf_s(x, size - (x - y), __VA_ARGS__)
#define SPRINTFD(x, y, ...)	 SPRINTFS(x, y, sizeof(y), __VA_ARGS__)
#define FCLOSE_W(f, w)		 if (f != w) fclose(f)
#define FCLOSE_F(f)			 FCLOSE_W(f, NULL)

typedef signed long long sLongLong;
typedef enum {
	eCalcStart,
	eCalcResult,
} eThreadStartMode;

template<typename T>void initArray(T** pPntr, int len, T val = 0) {
	auto *ptr = *pPntr = new T[len];
	while (len--)
		*(ptr+len) = val;
}

class SizeParam {
protected:
	SizeParam(int numDays, int numPlayers, int groupSize) :
		m_numDays(numDays),
		m_numPlayers(numPlayers),
		m_groupSize(groupSize) {
		if (m_numDays < 1 || m_numDays * (m_groupSize - 1) != m_numPlayers - 1 || 
			m_groupSize < 0 || m_numPlayers / m_groupSize * m_groupSize != m_numPlayers)
		{
			printf("*** Incorrect parameters: nPlayers=%d GroupSize=%d, Exit\n", m_numPlayers, m_groupSize);
			exit(1);
		}
	}
	~SizeParam() {
	}

	const int m_numDays;
	const int m_numPlayers;
	const int m_groupSize;
};

class CChecklLink : private SizeParam {
public:
	CChecklLink(int numDays, int numPlayers, int groupSize);
	~CChecklLink();
	bool checkLinks(char *c, int id, bool printLinksStatTime = false);
	bool checkLinks2(char* pLinks, int id, bool printLinksStatTime);
	bool checkLinks27(char *c, int id);
	bool checkLinksH(const char* c, const char* v, int nv, int nvo, int ind1, int ind2, char* vo);
	bool checkLinksH2(const char* c, const char* v, int nv, int nvo, int ind1, int ind2, char* vo);
	void reportCheckLinksData();
private:
	bool checkLinksV(const char* links, const char* v, int nv, int ind, char* vo);

	sLongLong cnt = 0;
	sLongLong tmtotal = 0;
	sLongLong tmtotalFalse = 0;
	sLongLong tmtotalOk = 0;
	sLongLong cntErr = 0;
	sLongLong cntOk = 0;
	sLongLong* counts = NULL;
	sLongLong* tmfalse = NULL;
	sLongLong* tmok = NULL;
	char *faults = NULL;
	char *m_pLinksCopy = NULL;
	char *m_co = NULL;
	char *m_v = NULL;
	char *m_vo = NULL;
#if PrintNVminmax
	char *nvmn = NULL;
	char *nvmx = NULL;
	void setNV_MinMax(int id, int idx, char nv);
#else
#define setNV_MinMax(id, idx, nv)   // empty macros
#endif
};

#ifdef CD_TOOLS
template<typename T, typename S> class CCanonicityChecker;
#define CheckCanon CCanonicityChecker<unsigned __int8, unsigned __int8>
#else
template<typename T> class CCheckerCanon;
#define CheckCanon CCheckerCanon<unsigned __int8>
#endif

class alldata : private SizeParam {
public:
	alldata(int numPlayers, int groupSize=GroupSize, 
		bool useCheckLinksV = UseCheckLinksV, bool useCheckLinksH = UseCheckLinksH, 
		int improveResult = ImproveResults, int createImprovedResult = CreateImprovedResult);

	~alldata();
	bool Run(int threadNumber=0, int iCalcMode=eCalcStart, 
		char* mstart0=NULL, char* mstart=NULL,
		int nrowsStart=0, int nrowsOut=0, sLongLong* pcnt=NULL, bool bPrint=false);
	bool initStartValues(const char* ivc, bool printStartValues=true);
	bool improveMatrix(int improveResult, unsigned char* bResults, const int lenResult, unsigned char** pbRes1 = NULL);
	int cnvCheckKm1(char* tr, int nrows, unsigned char* pOrbits=NULL, int* pDayMax=NULL) const;
	inline void initDayIdx(int nDays) const     {
		memcpy(m_DayIdx, result(), (m_NumDaysToTransform = nDays));
	}
private:
	void Init();
	inline auto numPlayers() const				{ return m_numPlayers; }
	inline auto numDays() const					{ return m_numDays; }
	bool initPrevDay();
	inline char *result(int nDay = 0) const		{ return m_pResults + nDay * m_numPlayers; }
	inline auto *links(int nPlayer = 0) const	{ return m_pLinks + nPlayer * m_numPlayers; }
	void getUnselected(char* v, int nv);
	void getPrevPlayer();
	int getNextPlayer();
	bool initCurrentDay();
	bool unsetLinksForOnePlayer(char* p, int ip) const;
	void setCheckLinks();
	bool processOneDay();
	int checkPlayer1(int iPlayerNumber);
	int checkPlayer3(int iPlayerNumber, int lastPlayer);
	void outputResults(int iDay, const unsigned char* pResult, int cntr = 0) const;
	void outputError() const;
	bool CheckMatrix(const char* matrix, int nl, int nc, int gs, bool printError, int* errDay, int* errGroup, int* dubLine);
	bool cnvCheck();
	bool cnvCheckKm(char* tr, char* tg);
	bool cnvCheckTg(char* tr, char* tg, int ntg, int gsf);
	void cnvInit();
	bool cnvCheckTgNew(char* tr, int gsf, bool bSkipFirst);
	bool cnvCheckNew();
	void testImproveMatrixSpeed();

	inline void addCanonCall(int idx = 0)		{ m_nCanonCalls[idx]++; }
	inline auto canonCalls(int idx) const		{ return m_nCanonCalls[idx]; }
	size_t m_nCanonCalls[2] = 				// Counter of CanonicityChecker::CheckCanonicity calls 
				{ (size_t)-1, (size_t)-1 };	// (total # of calls, # of calls with negative answer) 
	char* maxResult;
	int maxDays;
	int maxDaysPlayers;
	sLongLong nLoops;
	bool noMoreResults;
	char* m_pResults;
	char* m_pLinks;
	char* selPlayers;
	char* tmpPlayers;
	char* indexPlayer;
	char* m_h = NULL;
	char* m_ho = NULL;
	int  iPlayer, iDay;
	int m_groupIndex;     // Index of the group to change if the matrix is not canonical.
	bool bPrevResult;
	size_t m_nLenResults;
	int m_finalKMindex;
	int m_groupSizeFactorial;
	int m_nallTr;
	int m_nallTg;
	int m_bestTr;
	int m_bestTg;
	int m_improveResult;
	int m_createImprovedResult;
	int *m_TrTest;
	int *m_TgTest;
	char* m_Km;
	char* m_Km2;
	char* m_Ktmp;
	char* m_KmSecondRow;
	char* m_Km2ndRowInd;
	char* m_allTr;
	char* m_allTg;
	char* m_trmk;
	char* m_groups;
	mutable int m_NumDaysToTransform;
	mutable unsigned char* m_DayIdx;

	const int m_np2;
	const int m_nGroups;
	const bool m_bCheckLinkV;
	const bool m_bCheckLinkH;
	CChecklLink *m_pCheckLink = NULL;
	CheckCanon *m_pCheckCanon = NULL;
	char ImprovedResultFile[512];
	char ResultFile[512];
	FILE* m_file = NULL;    // File for output of improved matrices.
};

void convertLinksToResult(const char* ci, char* res, int np, int gs);
void _printf(FILE* f, bool toScreen, const char* format, const char* pStr = NULL);
void printTableColor(char const* name, const char *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const char *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const int *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 0.0);
void printTable(char const* name, const double *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 1.0);
bool _CheckMatrix(const char* matrix, int nl, int nc, int gs, char* links, bool printError, int* errLine, int* errGroup, int* dubLine);
int kmProcessMatrix2(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax);
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc);
int kmProcessMatrix(char* mo, char* mi, char* tmp, int nr, int nc, int gs, char* tm = NULL, int* pDayMax = NULL);
void kmSortGroupsByFirstValue(char* mo, char* mi, char nr, char nc, char np);
void kmSortRowsBy2ndValue(char* mo, char* mi, char nr, char nc, char* tm);
int factorial(int n);
void printThreadsStat(sLongLong* cntTotal, sLongLong* cnt, int nMatrices, int nProcessed, int nrowsStart, int nrowsOut, int nThreads, clock_t iTime, bool bPrintSetup);
void printTransformed(int nrows, int ncols, const char* tr, const char* ttr, const char* pImatr, const char* pTmatr, int numRow=0, sLongLong nLoops=0, int finalKMindex=0);
void createFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs, int iset);
void createStartFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs);
void saveStartData(char* fn, char* sm, int nm, int np, int nr, int gs);
int readStartData(char* fn, char* sm, int nm, int np, int nr, int gs);
bool setLinksForOnePlayer(int id, int np, char* lnk, char* p, int ip, char v);
void linksFromMatrix(char* lnk, char* iv, int nr, int np);
