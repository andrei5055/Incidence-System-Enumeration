#pragma once
#include <string> 
using namespace std;
#define nPlayers0 21
#define GroupSize 3
#define LoopsMax 30000
#define ImproveResults 0		// 0 - with no explanation of the matrix rejection reason;
                                // 1 - reason of rejection will be explained;
								// 2 - improve matrix as much as possible.
#define ResultFile		   "../bbb.txt" // Name of output file with the results, "" - no file output.
#define PrintImprovedResults 1	// Set this value to >= 2 if you want to see improved results on screen.
#define ImprovedResultFile "../aaa_.txt"  // Name of output file with the improved results, "" - no file output.
#define UsePos_1_4_condition 1
#define UseCheckLinksV 1
#define UseCheckLinksH 1
#define ReportInterval   120000

// Conditions to be tested on day 1:
#define USE_STATEMENT_7   1    // players[4] < players[7]  
#define USE_STATEMENT_17  1    // Only the players 4 or 9 could be in position[1, 4].
#define USE_STATEMENT_18  1    // The positions of any two players who were in the same group 
                               // on the first day on the second day must be in ascending order of their numbers.
#define USE_STATEMENT_19  1
#define CHECK_PERMUTATIONS	0
#define IMPROVE_RESULT_MAX  1000 // Maximal number of attempts to improve matrix 
								 // (works only with ImproveResults=2)
#define OUTPUT_VECTOR_STAT	0    // Output of some statistics for constructed matrices

#define PrintLinksStat 1
#define PrintLinksStatTime 0 /* 1 - ~50% more cpu required */
#define PrintNVminmax 0
#define UseSS 0
#define nPlayers (nPlayers0/GroupSize*GroupSize)
#define nGroups (nPlayers / GroupSize)
#define unset ((char)(-1))
#define printfRed(fmt, v) printf("\x1b[1;31m" fmt "\x1b[0m", v)
#define printfGreen(fmt, v) printf("\x1b[1;32m" fmt "\x1b[0m", v)
#define printfYellow(fmt, v) printf("\x1b[1;33m" fmt "\x1b[0m", v)


#define FOPEN_W(x, y, z, w)	 FILE *x = w; \
                             if (!x && y && strlen(y)) fopen_s(&x, y, z)
#define FOPEN_F(x, y, z)	 FOPEN_W(x, y, z, NULL)
#define SPRINTFD(x, y, ...)	 x += sprintf_s(x, sizeof(y) - (x - y), __VA_ARGS__)
#define FCLOSE_W(f, w)		 if (f != w) fclose(f)
#define FCLOSE_F(f)			 FCLOSE_W(f, NULL)


template<typename T>void initArray(T** pPntr, int len, T val = 0) {
	auto *ptr = *pPntr = new T[len];
	while (len--)
		*(ptr+len) = val;
}

class SizeParam {
protected:
	SizeParam(int numDays, int numPlayers, int groupSize=0) :
		m_numDays(numDays),
		m_numPlayers(numPlayers),
		m_groupSize(groupSize) {
		m_co = new char[numDays * numPlayers];
		m_lo = new char[numPlayers * numPlayers];
		if (groupSize)
			m_pBuf = new char[groupSize];
	}
	~SizeParam() {
		delete[] m_co;
		delete[] m_lo;
		delete[] m_pBuf;
	}
	void convertLinksToResult(const char *ci);

	const int m_numDays;
	const int m_numPlayers;
	const int m_groupSize;
	char* m_co = NULL;
	char* m_lo = NULL;
	char *m_pBuf = NULL;
};

class CChecklLink : private SizeParam {
public:
	CChecklLink(int numDays, int numPlayers);
	~CChecklLink();
	bool checkLinks(char *c, int id, bool printLinksStatTime = false);
	bool checkLinks27(char *c, int id);
	bool checkLinksH(const char* c, char* st, char *s, int ipos, const char* v, int nv, int nvo, int ind1, int ind2, char* vo, double* counter = NULL);
	void reportCheckLinksData();
private:
	bool checkLinksV(const char* links, const char* v, int nv, int ind, char* vo);

	double cnt = 0;
	double tmtotal = 0.0;
	double tmtotalFalse = 0;
	double tmtotalOk = 0;
	double cntErr = 0;
	double cntOk = 0;
	double *counts = NULL;
	double *tmfalse = NULL;
	double *tmok = NULL;
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
	alldata(int numPlayers, int groupSize=GroupSize, bool useCheckLinksV = UseCheckLinksV, bool useCheckLinksH = UseCheckLinksH);

	~alldata();
	bool Run(int improveResult=0);
	bool initStartValues(const char* ivc, bool printStartValues=true);
	bool improveMatrix(int improveResult, unsigned char* bResults, const int lenResult, unsigned char** pbRes1 = NULL);
private:
	void Init();
	inline auto numPlayers() const				{ return m_numPlayers; }
	inline auto numDays() const					{ return m_numDays; }
	bool initPrevDay();
	inline auto *result(int nDay = 0) const		{ return m_pResults + nDay * m_numPlayers; }
	inline auto *links(int nPlayer = 0) const	{ return m_pLinks + nPlayer * m_numPlayers; }
	void getUnselected(char* v, int nv);
	void getPrevPlayer();
	int getNextPlayer();
	bool initCurrentDay();
	bool setLinksForOnePlayer(char* p, int ip, char v) const;
	bool unsetLinksForOnePlayer(char* p, int ip) const;
	void setCheckLinks();
	bool processOneDay();
	int checkPlayer1(int iPlayerNumber);
	int checkPlayer3(int iPlayerNumber, int lastPlayer);
	void outputResults(int iDay, const unsigned char* pResult, int cntr = 0) const;
	void outputError() const;
    void sortLinks();
	bool CheckMatrix(const char* matrix, int nl, int nc, bool printError, int* errDay, int* errGroup, int* dubLine);

	inline void addCanonCall(int idx = 0)		{ m_nCanonCalls[idx]++; }
	inline auto canonCalls(int idx) const		{ return m_nCanonCalls[idx]; }
	size_t m_nCanonCalls[2] = 				// Counter of CanonicityChecker::CheckCanonicity calls 
				{ (size_t)-1, (size_t)-1 };	// (total # of calls, # of calls with negative answer) 
	char* maxResult;
	int maxDays;
	int maxDaysPlayers;
	int nLoops;
	bool noMoreResults;
	char* m_pResults;
	char* m_pLinks;
	char* selPlayers;
	char* tmpPlayers;
	char* indexPlayer;
	char* m_h = NULL;
	char* m_ho = NULL;
	int  iPlayer, iDay;
	bool bPrevResult;
	size_t m_nLenResults;

	const int m_np2;
	const int m_nGroups;
	const bool m_bCheckLinkV;
	const bool m_bCheckLinkH;
	CChecklLink *m_pCheckLink = NULL;
	CheckCanon *m_pCheckCanon = NULL;
	FILE* m_file = NULL;    // File for output of improved matrices.
};


int getLastSixIndex(alldata* s);
void _printf(FILE* f, bool toScreen, const char* format, const char* pStr = NULL);
void printTableColor(char const* name, const char *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const char *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const int *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 0.0);
void printTable(char const* name, const double *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 1.0);
bool _CheckMatrix(const char* matrix, int nl, int nc, char* links, bool printError, int* errLine, int* errGroup, int* dubLine);
void initPrevDay(alldata* s);
int processLastSix(alldata* s);

int compareMatrix(const char *result, int ncolumns, char* transition);


