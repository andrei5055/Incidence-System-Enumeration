#pragma once
#include <string> 
using namespace std;
#define nPlayers0 15
#define GroupSize 3
#define LoopsMax 1
#define UseCheckLinks 0
#define PrintLinksStat 0
#define PrintLinksStatTime 0 /* requered ~50% more cpu */
#define PrintNVminmax 0
#define UseLastSixAsGroup 1 /* Value 1 makes it 6 times faster */
#define UseLastSix (GroupSize == 3 && nPlayers > 9 && UseLastSixAsGroup != 0)
#define nPlayers (nPlayers0/GroupSize*GroupSize)
#define nGroups (nPlayers / GroupSize)
#define unset ((char)(-1))
#define printfRed(fmt, v) printf("\x1b[1;31m" fmt "\x1b[0m", v)
#define printfGreen(fmt, v) printf("\x1b[1;32m" fmt "\x1b[0m", v)
#define printfYellow(fmt, v) printf("\x1b[1;33m" fmt "\x1b[0m", v)

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
		if (groupSize)
			m_pBuf = new char[groupSize];
	}
	~SizeParam()			{ 
		delete[] m_co; 
		delete[] m_pBuf;
	}
	void convertLinksToResult(const char *ci);

	const int m_numDays;
	const int m_numPlayers;
	const int m_groupSize;
	char *m_co = NULL;
	char *m_pBuf = NULL;
};

class CChecklLink : private SizeParam {
public:
	CChecklLink(int numDays, int numPlayers);
	~CChecklLink();
	bool checkLinks(char *c, int id, bool printLinksStatTime = false);
	bool checkLinks27(char *c, int id);
private:
	bool checkLinksV(const char* c, const char* v, int nv, int ind, char* vo);

	double cnt = 0;
	double tmtotal = 0.0;
	double tmtotalFalse = 0;
	double tmtotalOk = 0;
	double cntErr = 0;
	double cntOk = 0;
	int icnt = 0;
	double *counts = NULL;
	double *tmfalse = NULL;
	double *tmok = NULL;
	char *faults = NULL;
	char *cb = NULL;
	char *m_co = NULL;
#if PrintNVminmax
	char *nvmn = NULL;
	char *nvmx = NULL;
	void setNV_MinMax(int id, int idx, char nv);
#else
#define setNV_MinMax(id, idx, nv)   // empty macros
#endif
};

template<typename T, typename S> class CCanonicityChecker;

class alldata : private SizeParam {
public:
	alldata(int numPlayers, int groupSize=GroupSize, bool useCheckLinks=UseCheckLinks);

	~alldata();
	bool Run();
	bool initStartValues(const char* ivc, bool printStartValues=true);
private:
	void Init();
	inline auto numPlayers() const				{ return m_numPlayers; }
	inline auto numDays() const					{ return m_numDays; }
	void initPrevDayProcess();
	inline auto *result(int nDay = 0) const		{ return m_pResults + nDay * m_numPlayers; }
	inline auto *links(int nPlayer = 0) const	{ return m_pLinks + nPlayer * m_numPlayers; }
	int getLastSixIndex(const char* resDay);
	void getLastSix(char* v);
	int processLastSix();
	void getPrevPlayer();
	int getNextPlayer();
	void initCurrentDay();
	bool setLinksAndDevCounts(char* p, int ip, char iset);
	void setCheckLinks();

	char* maxResult;
	int maxDays;
	int maxDaysPlayers;
	int nLoops;
	bool noResult;
	char* m_pResults;
	char* m_pLinks;
	char* selPlayers;
	char* tmpPlayers;
	char* indexPlayer;
	char* index6;
	int  iPlayer, iDay;
	int  iPrevResult;
	size_t m_nLenResults;

	const int m_np2;
	const int m_np3;
	const int m_nGroups;
	CChecklLink *m_pCheckLink = NULL;
	CCanonicityChecker<unsigned __int8, unsigned __int8> *m_pCheckCanon = NULL;
};

int getLastSixIndex(alldata* s);
void printTableColor(char const* name, const char *c, int nl, int nc, int ns, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const char *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false);
void printTable(char const* name, const int *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 0.0);
void printTable(char const* name, const double *c, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, double scale = 1.0);
void initPrevDayProcess(alldata* s);
int processLastSix(alldata* s);

int compareMatrix(const char *result, int ncolumns, char* transition);


