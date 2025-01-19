#pragma once
#include <iostream>
#include <format>
#include <assert.h>
#include "groupInfo.h"

#include "allsupports.h"
#include "MatrixDB.h"
#include "cycles.h"
#include "RowStorage.h"

#define UseRowsPrecalculation  0 // The number of first rows of the matrix from which we start to pre-calculate the remaining rows. 
                                 // ntd: change to parameter, check input params that numbers of rows in input is 3 or less

#define UseMultiThreading   0

#if UseMultiThreading
#define NThreads 10				  // number of threads to calculate matrices		
#define NRowsInStartMatrix 3	  // number of rows in initial set of matrices for threads
#define NRowsInResultMatrix 5	  // number of rows in result, if 0 then calculate full matrix
#define MaxNumberOfStartMatrices 1000000
#define FirstIndexOfStartMatrices 0 // 12220
#else
#define NThreads 1
#define NRowsInStartMatrix 0 // must be 0
#define NRowsInResultMatrix 2
#define MaxNumberOfStartMatrices 1
#define FirstIndexOfStartMatrices 0
#endif
#define NRBase					  2
#define UseP1fCheckGroups         1
#define USE_BINARY_CANONIZER	  0
#define AllowNotP1FRows    0 // Only for 3U1F. if 0 - all rows pairs must be p1f for all common sets combinations.
									// 1 - no check for rows pairs related to number of cycles.
#define Any2RowsConvertToFirst2   1 // Only for 3U1F. (*) If 1, then any two rows must be converted to first 2 rows.
									// If 0, then a. no requirement (*), b. more than (*) first 2 rows pairs used 
									// If 2, then a. no requirement (*), b. same as (*) set of first 2 rows pairs used
#define GenerateSecondRowsFor3U1F 0 // Only for 3U1F. Use 1 in single thread mode only with NRowsInResultMatrix=2. 
									// SW will use on the fly calculated 2nd rows (instead of precalculated '_expected2ndRow3p1f')
									// Calculated rows will be printed at the end
#define UseUniform1Factorization  0
#define U1FName ""

#define Use2RowsCanonization    0 
#define SubmatrixGroupOrderMin  0
#define ResultGroupOrderMin     0
#define LoopsMax 200000000000
#define PrintMatrices			1
#define SavingMatricesToDisk	1
#define ExpectedResult	   -1     // Expected Number of Matrices to be Constructed
#define CheckConstructedMatrices	0
#define UseImproveMatrix	0
#define ImproveResults		0	  // 0 - with no explanation of the matrix rejection reason;
                                  // 1 - reason of rejection will be explained;
								  // 2 - improve matrix as much as possible.
#define CreateImprovedMatrix 0	  // 1 - cnvCheckNew will create improved matrix in m_Km, 0 - not (0 works faster) 
#define StartFolder				  "../Logs1/" 	// Name of folder with 'Start Matrices', must present in multithread mode.
#define ResultFolder			  "../Logs1/" 	// Name of folder with the results, "" - no file output.
#define PrintImprovedResults 0	  // Set this value to >= 2 if you want to see improved results on screen.
#define ImprovedResultFolder	  "../ImprovedResults/"	// Name of folder for files with the improved results.
#define UsePos_1_4_condition 1
#define UseCheckLinksV 1
#define UseCheckLinksT 0
#define StartMatricesFolder		  "../StartMatrices/" 	// Name of folder to save/Load Start Matrices, "" - do not use.


// Conditions to be tested on day 1:
#define USE_STATEMENT_7   1    // players[4] < players[7]  
#define USE_STATEMENT_17  1    // Only the players 4 or 9 could be in position[1, 4].
#define USE_STATEMENT_18  1    // The positions of any two players who were in the same group 
                               // on the first day on the second day must be in ascending order of their numbers.
#define USE_CHANGING_DAY_0_GROUPS 1 // Use permutation of first 3 groups of day 0
#define CHECK_WITH_GROUP  0    // Use new group on completely constructed matrices
#define Use_GroupOrbits	  0
#define USE_TRANSLATE_BY_LEO 0
#define USE_ORBTS		  0
#define USE_EQUAL		  0
#define PRINT_MATR_CNTR	  0
#define PRINT_TRANSFORMED 0

#define USE_STATEMENT_19  1
#define CHECK_PERMUTATIONS	0
#define IMPROVE_RESULT_MAX  1000 // Maximal number of attempts to improve matrix 
								 // (works only with ImproveResults=2)
#define OUTPUT_VECTOR_STAT	0    // Output of some statistics for constructed matrices

#define UseSS 0

#define FOPEN_W(x, y, z, w)	 FILE *x = w; \
                             if (!x && y && strlen(y)) fopen_s(&x, y, z)
#define FOPEN_F(x, y, z)	 FOPEN_W(x, y, z, NULL)
#define SPRINTFS(x, y, size,...)	 x += sprintf_s(x, size - (x - y), __VA_ARGS__)
#define SPRINTFD(x, y, ...)	 SPRINTFS(x, y, sizeof(y), __VA_ARGS__)
#define FCLOSE_W(f, w)		 if (f != w) fclose(f)
#define FCLOSE_F(f)			 FCLOSE_W(f, NULL)

#define GetBit(a, n) (*((char*)a+((n) / 8)) & (1<<(7-((n) & 7))))
#define SetBit(a, n) (*((char*)a+((n) / 8)) |= (1<<(7-((n) & 7))))
#define ResetBit(a, n) (*((char*)a+((n) / 8)) &= ~(1<<(7-((n) & 7))))

#define SWAP(a, b)          a ^= (b ^= (a ^= b))

typedef enum {
	eCalcResult,
	eCalculateRows,
	eCalculateMatrices,
	eDisabled,
} eThreadStartMode;

template<typename T>CC void initArray(T** pPntr, int len, T val = 0) {
	auto *ptr = *pPntr = new T[len];
	while (len--)
		*(ptr+len) = val;
}

#ifdef CD_TOOLS
template<typename T, typename S> class CCanonicityChecker;
#define CheckCanon CCanonicityChecker<unsigned __int8, unsigned __int8>
#else
template<typename T> class CCheckerCanon;
#define CheckCanon CCheckerCanon<unsigned __int8>
#endif

template<typename T> class CGroupOrbits;
#define GroupOrbits CGroupOrbits<unsigned __int8>

typedef struct TrCycles {
	int counter;
	int ncycles;
	tchar start[MAX_CYCLES_PER_SET];
	tchar length[MAX_CYCLES_PER_SET];
	tchar fullPath[MAX_PLAYER_NUMBER * 2];
} TrCycles;


typedef CRepository CBinaryMatrixStorage;

class alldata : private CGroupInfo, CGroupUtilisation, CycleSupport, CChecklLink {
	typedef bool(alldata::*checkU1F)(int);
	typedef void(alldata::*sortGroups)(tchar *, int) const;
	typedef int(alldata::*processMatrix2)(ctchar* mi, ctchar* tr, int nr, tchar ind) const;
	typedef bool(alldata::*checkInvalidCycle)(int ncycles, ctchar* cycles) const;
public:
	CC alldata(const SizeParam& p, const kSysParam* pSysParam, CRowStorage* pRowStorage = NULL, bool useCheckLinksT = UseCheckLinksT,
		int improveResult = ImproveResults, bool createImprovedResult = CreateImprovedMatrix);
	CC ~alldata();
	CC sLongLong Run(int threadNumber=0, int iCalcMode=eCalcResult, tchar* mstart0 = NULL, tchar* mstart = NULL,
		int nrowsOut=0, sLongLong* pcnt=NULL, std::string* pOutResult = NULL, int iThread = 0);
	bool initStartValues(const char* ivc, bool printStartValues=true);
	CC bool improveMatrix(int improveResult, tchar* bResults, const int lenResult, tchar** pbRes1 = NULL);
	CC int cnvCheckKm1(ctchar* tr, int nrows, tchar* pOrbits=NULL);
	inline GroupOrbits* orbits() const			{ return m_pOrbits; }
	CC inline auto numPlayers() const			{ return m_numPlayers; }
	CC inline auto RowStorage() const			{ return m_pRowStorage; }
	CC inline auto groupSize()	const			{ return m_groupSize; }
#if !USE_CUDA
	inline MatrixDB* matrixDB()					{ return &m_matrixDB; }
#endif
	CC inline void initCheckByGroup(tchar iDay, tchar iMode) {
		// Creating the sequences 0,1,2,3,... as the day's indices.
		memcpy(m_DayIdx, result(), (m_NumDaysToTransform = iDay));
		if (!(m_cnvMode = iMode))
			resetGroupOrder();

#if NEED_TO_DEBUG
		orbits()->resetOrbits((ctchar*)result());
#endif
	}
	CC bool p1fCheck3(ctchar* rowi, ctchar* rowj, ctchar* neighborsi, ctchar* neighborsj) const;
	CC inline auto numDaysResult() const		{ return m_numDaysResult; }
	CC inline tchar* result(int nDay = 0) const { return m_pResults + nDay * m_numPlayers; }
private:
	CC void Init();
	inline auto numDays() const					{ return m_numDays; }
	CC bool initPrevDay();
	CC inline tchar* links(int nPlayer = 0) const	{ return m_pLinks + nPlayer * m_numPlayers; }
	CC void getUnselected(tchar* v, int nv) const;
	CC void getPrevPlayer();
	CC int getNextPlayer();
	CC bool initCurrentDay();
	CC bool unsetLinksForOnePlayer(ctchar* p, int ip) const;
	void setCheckLinks();
	CC bool processOneDay();
	CC int checkPlayer1(int iPlayerNumber);
	CC int checkPlayer3(int iPlayerNumber, int lastPlayer);
#if !USE_CUDA && PrintImprovedResults
	void outputResults(int iDay, const unsigned char* pResult, int cntr = 0) const;
#else
#define outputResults(iDay, pResult, ...)
#endif
	void outputError() const;
	bool CheckMatrix(ctchar* matrix, int nl, int nc, int gs, bool printError, int* errDay, int* errGroup, int* dubLine);
	CC bool cnvCheckKm(ctchar* tr, ctchar* tg, int nrows);
	CC void cnvInit();
	CC bool cnvCheckTgNew(ctchar* tr, int nrows, int ngroups);
	CC bool cnvCheckNew(int iMode, int nrows, bool useAutomorphisms = true);
	CC bool canonizator(int iMode, int nrows);
	CC bool cnvCheck45(int nrows);
	CC bool cnvCheck2U1F(int nrows);
	CC bool cnvCheck2P1F(int nrows);
	CC bool cnvCheck3U1F(int nrows);
	void testCanonizatorSpeed();
	void testRightNeighbor(int nr);
	void TestkmProcessMatrix(int nrows, tchar n, ctchar* tr, ctchar* ttr, int icmp) const;
	CC bool matrixStat(ctchar* table, int nr, bool* pNeedOutput = NULL);
	char *matrixStatOutput(char* str, int maxStr) const;
	CC int p3Cycles(TrCycles* trc, int ncr, ctchar* t1, ctchar* t2, ctchar* v, ctchar* res1, ctchar* res2) const;
	CC int u1fGetCycleLength(TrCycles* trc, int ncr, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, int ind = 0) const;
	CC int u1fGetCycleLength3(TrCycles* trc, int ncr, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, int ind = 0) const;
	CC void sortCycles(tchar* cycles, tchar* cyclcesStart, int ncycles) const;
	CC int collectCyclesAndPath(TrCycles* trc);
	CC bool getCyclesAndPath3(TrCycles* trc, ctchar* v, ctchar* t0, ctchar* t1, ctchar* res0, ctchar* res1) const;
	CC int getCyclesAndPath(TrCycles* trc, int ncr, ctchar* tt1, ctchar* tt2, ctchar* tt3 = NULL, ctchar* tt4 = NULL) const;
	CC int checkCurrentResult(bool bPrint, void* pIS_Canonizer = NULL);
	CC int kmProcessMatrix2p1f(tchar* tr, int nr, int ind0, int ind1);
	CC void goBack();
	CC void p1fCheckStartMatrix(int nr);
	CC bool create2P1FTr(tchar* tr, tchar kStart, ctchar* pf0, ctchar* pf1, ctchar* pfi, ctchar* pfj) const;
	CC bool create3U1FTr(tchar* tr, tchar k0Start, tchar k1Start, ctchar* v0, ctchar* v1,
		ctchar* t0, ctchar* t1, ctchar* res1, tchar ir0, tchar ir1, int idir, int iPrint = 0);
	CC bool create3U1FTr1(tchar* tr, tchar k0Start, tchar k1Start, ctchar* v0, ctchar* v1,
		ctchar* t0, ctchar* t1, ctchar* res1, tchar ir0, tchar ir1, int idir, int iPrint = 0) const;
	CC bool createU1FTr(tchar* tr, TrCycles* trCycles01, TrCycles* trCycles,
		ctchar* dir, ctchar* offset, ctchar* start, int iPrint = 0);


	CC int getAllV(tchar* allv, int maxv, tchar ir1, tchar ir2, tchar* pt2 = NULL) const;
	CC int p1fCheck2ndRow() const;
	CC ctchar* expected2ndRow3p1f(int iSet) const;
	CC void updateIndexPlayerMinMax();
	void cnvPrintAuto(ctchar* tr, int nrows);

	inline void addCanonCall(int idx = 0)		{ m_nCanonCalls[idx]++; }
	inline auto canonCalls(int idx) const		{ return m_nCanonCalls[idx]; }
	CC inline bool checkCanonicity() const		{ return iDay == m_matrixCanonInterval; }
	CC void kmSortGroups3(tchar* mi, int nr) const;
	CC void kmSortGroups2(tchar* mi, int nr) const;
	CC void kmSortGroups(tchar* mi, int nr) const;
	CC int kmProcessMatrix(ctchar* mi, ctchar* tm, int nr, tchar ind = 0) const;
	CC int kmProcessMatrix2(ctchar* mi, ctchar* tr, int nr, tchar ind) const;
	CC int kmProcessMatrix3(ctchar* mi, ctchar* tr, int nr, tchar ind) const;
	CC int kmProcessOneNot1stRow3(tchar* mo, ctchar* mi, int mind, tchar* tb, tchar* tc, ctchar* tr, int nr) const;
	CC int kmProcessOneNot1stRow2(ctchar* mi, int mind, tchar* tb, ctchar* tr, int nr, int irow = 2) const;
	CC void setPlayerIndex(ctchar* tr, int iDayMax, int iDayCurrent, ctchar* co, ctchar* ci, ctchar* ciFrom, int nc) const;
	CC void setPlayerIndexByPos(ctchar* tr, ctchar* co, ctchar* ciFrom, int iDayMax, int iDayCurrent, int nc, int ip) const;
	CC void releaseBinaryMatricesStorage();
	bool FindIsomorphicBaseElements(const std::string& fileName);
	void checkCommonValues(ctchar* pBaseValues, int numSets);
    void checkCommonValues();
	void printResultWithHistory(char const* name, int nRows);

	CC void setArraysForLastRow(int nrows);
	CC void adjustPlayerPosition(tchar* path, tchar length, tchar nrows);

	CC bool CycleIsInvalid(int ncycles, ctchar* cycles) const {
		return ncycles != 1 || cycles[0] != m_numPlayers;
	}
	CC bool CycleIsInvalid_27_999(int ncycles, ctchar* cycles) const {
		return ncycles != 3 || cycles[0] != 9 || cycles[1] != 9;
	}

	CC bool CycleIsInvalid_21_669_912(int ncycles, ctchar* cycles) const {
		return (ncycles != 2 || (cycles[0] != 9 && cycles[0] != 12)) &&
			(ncycles != 3 || (cycles[0] != 6 && cycles[0] != 9));
	}
	Cycles* m_cycles;
	size_t m_nCanonCalls[2] = 				// Counter of CanonicityChecker::CheckCanonicity calls 
				{ (size_t)-1, (size_t)-1 };	// (total # of calls, # of calls with negative answer) 
	char* maxResult;
	int maxDays;
	int maxDaysPlayers;
	sLongLong nLoops;
	bool noMoreResults;
	clock_t* m_rowTime;
	tchar* m_pResults;
	tchar* m_pResultsPrev;
	tchar* m_pLinks;
	tchar* selPlayers;
	tchar* tmpPlayers;
	tchar* indexPlayer;
	tchar* m_indexPlayerMin;
	tchar* m_indexPlayerMax;
	tchar* m_h = NULL;
	tchar* m_ho = NULL;
	tchar* m_p3fSecondRows;
	int m_p3fNumSecondRows;
	int m_p3fNumSecondRowsAct;
	int m_p3fV01nv;
	int  iPlayer, iPlayerIni, iDay;
	int m_firstNotSel;
	mutable int m_groupIndex;      // Index of the group to change if the matrix is not canonical.
	mutable int m_playerIndex;     // Index of the player to change if the matrix is not canonical.
	bool bPrevResult;
	int m_improveResult;
	int m_TrInd;
	int m_cnvMode;
	int m_useRowsPrecalculation;
	int m_secondPlayerInRow4;
	int m_secondPlayerInRow4Last;

	tchar* m_Km;
	tchar* m_Km2;
	tchar* m_Ktmp;
	tchar* m_Km2ndRowInd;
	tchar* m_trmk;
	tchar* m_groups;
	tchar* m_tx = NULL;

	mutable int m_numCycles;
	mutable int m_NumDaysToTransform;
	mutable unsigned char* m_DayIdx;

	const int m_nGroups;
	const size_t m_nLenResults;
	const bool m_bCheckLinkV;
	const bool m_bCheckLinkT;
	CheckCanon *m_pCheckCanon = NULL;
	GroupOrbits *m_pOrbits = NULL;
	unsigned int m_p1f_counter = 0;
	bool m_checkForUnexpectedCycle;
	int m_matrixCanonInterval;   // 0 - Canonicity will be verified only for fully constructed matrices.
	                             // != 0 Canonicality will be checked on lines with numbers proportional to m_matrixCanonInterval
#if !USE_CUDA
	int m_finalKMindex;
	std::string ImprovedResultFile;
	std::string ResultFile;
	FILE* m_file = NULL;    // File for output of improved matrices.
	MatrixDB m_matrixDB;
#endif
	void* m_pRes;
	TrCycles m_TrCycles;
	TrCycles m_TrCyclesAll[MAX_CYCLE_SETS];
	checkU1F m_pCheckFunc;
	sortGroups m_pSortGroups;
	processMatrix2 m_pProcessMatrix;
	checkInvalidCycle m_pInvalidCycle;
	CBinaryMatrixStorage** m_ppBinMatrStorage = NULL;
	bool m_bRowStorageOwner;
	CRowStorage* m_pRowStorage = NULL;
	CRowUsage* m_pRowUsage = NULL;
};

inline bool is_number(const std::string& s)
{
	return !s.empty() && std::find_if(s.begin(), s.end(), 
		[](unsigned char c) { return !std::isdigit(c); }) == s.end();
}
