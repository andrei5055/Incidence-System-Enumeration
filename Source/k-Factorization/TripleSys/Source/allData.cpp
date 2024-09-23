#include "TripleSys.h"
#ifdef CD_TOOLS
#include "CanonicityChecker.h"
#else
#include "CheckCanon.h"
#endif
using namespace std;

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinksV, bool useCheckLinksH, bool useCheckLinksT,
	int improveResult, bool createImprovedMatrix) :
	SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
	m_np2(numPlayers* numPlayers),
	m_nGroups(numPlayers / groupSize),
	m_bCheckLinkV(groupSize == 3 && useCheckLinksV),
	m_bCheckLinkH(groupSize <= 3 && useCheckLinksH),
	m_bCheckLinkT(groupSize == 3 && useCheckLinksT),
	m_improveResult(improveResult),
	m_createImprovedMatrix(createImprovedMatrix) {
	m_nLenResults = m_numDays * numPlayers;
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[5 * m_numPlayers];
	tmpPlayers = selPlayers + m_numPlayers;
	indexPlayer = tmpPlayers + m_numPlayers;
	m_h = indexPlayer + m_numPlayers;
	m_ho = m_h + m_numPlayers;
	ImprovedResultFile[0] = '\0';
	ResultFile[0] = '\0';
	m_groupSizeFactorial = factorial(m_groupSize);
	m_TrInd = 0;
	m_nallTr = 1; 
	m_nallTg = m_groupSizeFactorial;
	for (int j = 2; j <= m_nGroups; j++)
	{
		m_nallTr *= j; m_nallTg *= m_groupSizeFactorial;
	}
	m_nTr = m_nallTr * m_nallTg; // number of all transitions (without iterations for days)
	m_nTrBytes = (m_nTr + 7) / 8; // length of transitions bitmask array in bytes
#if UseTrMask
	m_TrMask = new char[m_nTrBytes * m_numDays];
	if (!m_TrMask)
	{
		printf("*** No memory for m_TrMask, exit (m_nallTr=%d m_nallTg=%d m_numDays=%d)\n", m_nallTr, m_nallTg, m_numDays);
		exit(1);
	}
	memset(m_TrMask, 0, m_nTrBytes * m_numDays);
#endif
#if USE_cnvCheckNew == 0
	m_allTr = new char[m_nallTr * m_nGroups];
	m_allTg = new char[m_nallTg * m_nGroups];
	m_TrTest = new int[m_nallTr];
	memset(m_TrTest, 0, m_nallTr * sizeof(m_TrTest[0]));
	m_TgTest = new int[m_nallTg];
	memset(m_TgTest, 0, m_nallTr * sizeof(m_TgTest[0]));
#endif
	m_finalKMindex = 0;
	m_bestTr = 0;
	m_bestTg = 0;
	m_Km = new char[m_numPlayers * m_numPlayers]; // m_Km can be used for sort and needs m_numPlayers rows
	m_KmSecondRow = m_Km + m_numPlayers;
	m_Ktmp = new char[m_numPlayers * m_numPlayers]; // m_Ktmp can be used for sort and needs m_numPlayers rows
	m_Km2ndRowInd = new char[m_numPlayers];
	memset(m_Km2ndRowInd, 0, m_numPlayers);
	m_trmk = new char[m_numPlayers];
	m_groups = new char[m_groupSizeFactorial * m_groupSize];
	m_pLinks = new char[m_np2];
	m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers, m_groupSize);
	m_DayIdx = new unsigned char[m_numDays];
#if !USE_EQUAL && !CHECK_WITH_GROUP
	auto nDays = m_numDays;
	while (nDays--) m_DayIdx[nDays] = nDays;
#endif

#ifdef CD_TOOLS
	m_pCheckCanon = new CCanonicityChecker<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize, t_kSystems);
#else
	//m_pCheckCanon = new CCheckerCanon<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize);
	m_pCheckCanon = new CCheckerCanon<unsigned char>(m_numDays, numPlayers, groupSize);
#endif

	Init();
	if (!strchr(ImprovedResultFile, '_')) {
		FOPEN_F(f, ImprovedResultFile, "w");
		m_file = f;
	}

	m_pOrbits = new CGroupOrbits<unsigned char>(m_numPlayers);
}

alldata::~alldata() {
	delete[] m_pResults;
	delete[] selPlayers;
	delete[] m_pLinks;
#if UseTrMask
	delete[] m_TrMask;
#endif
#if USE_cnvCheckNew == 0
	delete[] m_allTr;
	delete[] m_allTg;
	delete[] m_TrTest;
	delete[] m_TgTest;
#endif
	delete[] m_Km;
	delete[] m_Km2ndRowInd;
	delete[] m_Ktmp;
	delete[] m_trmk;
	delete[] m_groups;
	delete[] m_DayIdx;
	delete m_pCheckLink;
	delete m_pCheckCanon;
	delete m_pOrbits;
	FCLOSE_F(m_file);
}

void alldata::Init() {
	memset(links(), unset, m_numPlayers * m_numPlayers);
	memset(result(), 0, m_nLenResults);
	maxDays = -1;
	nLoops = 0;
	noMoreResults = false;
	iDay = 0;
	bPrevResult = false; // can be false, or true to go to prev day
	cnvInit();
}

