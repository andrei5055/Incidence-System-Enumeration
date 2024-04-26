#include "TripleSys.h"
#ifdef CD_TOOLS
#include "../CanonicityChecker.h"
#else
#include "CheckCanon.h"
#endif

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinksV, bool useCheckLinksH, int improveResult, int createImprovedResult) :
	SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
	m_np2(numPlayers* numPlayers),
	m_nGroups(numPlayers / groupSize),
	m_bCheckLinkV(groupSize == 3 && useCheckLinksV),
	m_bCheckLinkH(groupSize <= 3 && useCheckLinksH),
	m_improveResult(improveResult),
	m_createImprovedResult(createImprovedResult) {
	m_nLenResults = m_numDays * numPlayers;
	maxResult = new char[m_nLenResults];
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[5 * m_numPlayers];
	tmpPlayers = selPlayers + m_numPlayers;
	indexPlayer = tmpPlayers + m_numPlayers;
	m_h = indexPlayer + m_numPlayers;
	m_ho = m_h + m_numPlayers;
	ImprovedResultFile[0] = '\0';
	ResultFile[0] = '\0';
	m_groupSizeFactorial = factorial(m_groupSize);
	int n = 1, m = m_groupSizeFactorial;
	for (int j = 2; j <= m_nGroups; j++)
	{
		n *= j; m *= m_groupSizeFactorial;
	}
#if USE_cnvCheckNew
	m = n = 2;
#endif
	m_nallTr = n;
	m_nallTg = m;
	m_finalKMindex = 0;
	m_allTr = new char[n * m_nGroups];
	m_allTg = new char[m * m_nGroups];
	m_bestTr = 0;
	m_bestTg = 0;
	m_TrTest = new int[n];
	memset(m_TrTest, 0, n * sizeof(m_TrTest[0]));
	m_TgTest = new int[m];
	memset(m_TgTest, 0, m * sizeof(m_TgTest[0]));
	m_Km = new char[m_numPlayers * m_numDays];
	m_KmSecondRow = m_Km + m_numPlayers;
	m_Ktmp = new char[m_numPlayers * m_numPlayers];
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
}

alldata::~alldata() {
	delete[] maxResult;
	delete[] m_pResults;
	delete[] selPlayers;
	delete[] m_pLinks;
	delete[] m_allTr;
	delete[] m_allTg;
	delete[] m_TrTest;
	delete[] m_TgTest;
	delete[] m_Km;
	delete[] m_Km2ndRowInd;
	delete[] m_Ktmp;
	delete[] m_trmk;
	delete[] m_groups;
	delete[] m_DayIdx;
	delete m_pCheckLink;
	delete m_pCheckCanon;
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

