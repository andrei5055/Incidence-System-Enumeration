#include "TripleSys.h"
#ifdef CD_TOOLS
#include "../CanonicityChecker.h"
#else
#include "CheckCanon.h"
#endif
using namespace std;

CC alldata::alldata(const SizeParam& p, const kSysParam* pSysParam,
	bool useCheckLinksT, int improveResult, bool createImprovedMatrix) :
	CGroupInfo(pSysParam->val[t_numPlayers], 200), CGroupUtilisation(pSysParam),
	CycleSupport(pSysParam->val[t_numPlayers]), CChecklLink(p, pSysParam),
	m_nGroups(numPlayers() / m_groupSize),
	m_bCheckLinkV(m_groupSize == 3 && param(t_useCheckLinksV)),
	m_bCheckLinkT(m_groupSize == 3 && useCheckLinksT),
	m_nLenResults(m_numDays * numPlayers()) {
	m_improveResult = improveResult;
	m_rowTime = new clock_t[m_numDays];
	m_rowTime[0] = 0;
	m_pResults = new tchar[m_nLenResults];
	m_pResultsPrev = new tchar[m_nLenResults];
	int numPlayers64 = (m_numPlayers + 7) / 8 * 8;

	selPlayers = new tchar[7 * numPlayers64];
	tmpPlayers = selPlayers + numPlayers64;
	indexPlayer = tmpPlayers + numPlayers64;
	m_h = indexPlayer + numPlayers64;
	m_ho = m_h + numPlayers64;
	m_indexPlayerMin = m_ho + numPlayers64;
	m_indexPlayerMax = m_indexPlayerMin + numPlayers64;
	for (int i = 0; i < m_numPlayers; i++)
	{
		m_indexPlayerMax[i] = (i % m_groupSize) == 0 ? i : m_numPlayers - 1;
		m_indexPlayerMin[i] = i / m_groupSize + (i % m_groupSize);
	}

#if !USE_CUDA
	ImprovedResultFile[0] = ResultFile[0] = '\0';
	m_finalKMindex = 0;
#endif
	m_TrInd = 0;
	m_tx = new tchar[m_numPlayers * (m_numPlayers - 2) / 2];

	const auto np2 = m_numPlayers * m_numPlayers;
	m_Km = new tchar[np2];		// m_Km can be used for sort and needs m_numPlayers rows
	m_Ktmp = new tchar[np2];	// m_Ktmp can be used for sort and needs m_numPlayers rows
	m_Km2ndRowInd = new tchar[m_numPlayers];
	memset(m_Km2ndRowInd, 0, m_numPlayers);
	m_trmk = new tchar[m_numPlayers];
	m_groups = new tchar[m_groupSizeFactorial * m_groupSize];
	m_pLinks = new tchar[np2];
	m_p3fSecondRows = new tchar[m_numPlayers * MAX_3PF_SECOND_ROWS];
	m_p3fNumSecondRows = 0;
	m_p3fNumSecondRowsAct = 0;
	m_DayIdx = new tchar[m_numDays];
	m_numCycles = 0;
	m_firstNotSel = 0;

	iPlayerIni = m_numPlayers - 1;
	if (m_numPlayers > m_groupSize)
		iPlayerIni -= m_groupSize;

#if !USE_EQUAL && !CHECK_WITH_GROUP
	auto nDays = m_numDays;
	while (nDays--) m_DayIdx[nDays] = nDays;
#endif

	if (param(t_useCheckLinksV) || param(t_useCheckLinksV)) {
#ifdef CD_TOOLS
		m_pCheckCanon = new CCanonicityChecker<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize, t_kSystems);
#else
		m_pCheckCanon = new CCheckerCanon<unsigned char>(m_numDays, numPlayers(), m_groupSize);
#endif
	}

	Init();
#if !USE_CUDA
	if (ImprovedResultFile.find('_') == string::npos) {
		FOPEN_F(f, ImprovedResultFile.c_str(), "w");
		m_file = f;
	}
#endif
	if (m_p1f && m_groupSize == 3)
		InitCycleSupport(m_nGroups);

	m_pCheckP1F = m_p1f? ((m_groupSize == 2) ? &alldata::cnvCheck2P1F : &alldata::cnvCheck3P1F) : NULL;
	m_pSortGroups = m_groupSize == 2 ? &alldata::kmSortGroups2 : (m_groupSize == 3 ? &alldata::kmSortGroups3 : &alldata::kmSortGroups);

	m_pProcessMatrix = createImprovedMatrix || m_groupSize > 3 ? &alldata::kmProcessMatrix : (m_groupSize == 2 ? &alldata::kmProcessMatrix2 : &alldata::kmProcessMatrix3);
	if (m_groupSize <= 3 && numPlayers() >= 10 && param(t_useCheckLinksH)) {
		// NOTE: the use of &alldata::checkLinksH is the cause of the warning:
		// CUDACOMPILE : ptxas warning : Stack size for entry function '_Z10initKernelPP7alldata9SizeParamPK9kSysParami' 
		// cannot be statically determined
		m_pCheckLinksH = m_groupSize == 2 ? &alldata::checkLinksH2 : &alldata::checkLinksH;
	}

	expected2ndRow3p1f(-1);

	typedef struct {
		int numPlayer;
		tchar u1f[20];
		int size;
		alldata::checkInvalidCycle pFunc;
	} SpecialCaseFuncDescr;


	SpecialCaseFuncDescr invCyclesFunc[2] = { 
		{27, {1, 9, 9, 9}, 4, &alldata::CycleIsInvalid_27_999},
		{21, {2, 6, 6, 9, 0, 9, 12}, 7, &alldata::CycleIsInvalid_21_669_912}
	};

	// Assign the special case function for determination of the invalid cycles.
	m_pInvalidCycle = &alldata::CycleIsInvalid;
	const auto u1f_in = pSysParam->u1f[0];
	if (u1f_in) {
		for (int i = countof(invCyclesFunc); i--;) {
			const auto& specCase = invCyclesFunc[i];
			if (numPlayers() != specCase.numPlayer)
				continue;

			const auto& u1f = specCase.u1f;
			int j = specCase.size;
			while (j-- && u1f[j] == u1f_in[j]);
			if (j < 0) {
				m_pInvalidCycle = specCase.pFunc;
				expected2ndRow3p1f(-2);
				break;
			}
		}
	}
			               
#if Use_GroupOrbits
	m_pOrbits = new CGroupOrbits<unsigned char>(m_numPlayers);
#endif
	m_cycles = new Cycles[m_numDays * (m_numDays - 1)];
#if !USE_CUDA
#if 0   // Testing the construction of all permutations
	tchar perm[] = { 0, 1, 2, 3, 4, 5 };
	int i = 0;
	char buffer[256];
	const auto len = countof(buffer);
	do {
		auto* pBuff = buffer;
		SPRINTFS(pBuff, buffer, len, "%4d: ", ++i);
		for (int j = 0; j < 6; j++)
			SPRINTFS(pBuff, buffer, len, " %2d", perm[j]);

		printf("%s\n", buffer);
	} while (nextPerm(perm, countof(perm)));
#endif

#if 0  // Testing the construction of all permutations directions/indices/permut_cycle combinations
	tchar lenCycles[] = { 15}, startCycles[] = { 0 };
	//tchar lenCycles[] = { 6, 6, 9 }, startCycles[] = { 9, 15, 0 } /*{ 0, 6, 12 }*/;
	const auto numCycles = countof(lenCycles);
	
	ctchar *pDir, *pStartOut;
	auto pIdx = InitCycleMapping(lenCycles, startCycles, numCycles, 3, &pDir, &pStartOut);
	int k = 0;
	char buffer[256];
	const auto len = countof(buffer);
	do {
		auto *pBuff = buffer;
		SPRINTFS(pBuff, buffer, len, "%4d: ", ++k);
		for (int j = 0; j < numCycles; j++)
			SPRINTFS(pBuff, buffer, len, " %2d", pDir[j]);

		SPRINTFS(pBuff, buffer, len, "    ");
		for (int j = 0; j < numCycles; j++)
			SPRINTFS(pBuff, buffer, len, " %2d", pIdx[j]);

		SPRINTFS(pBuff, buffer, len, "    ");
		for (int j = 0; j < numCycles; j++)
			SPRINTFS(pBuff, buffer, len, " %2d", pStartOut[j]);

		SPRINTFS(pBuff, buffer, len, "    ");
		for (int j = 0; j < numCycles; j++)
			SPRINTFS(pBuff, buffer, len, " %2d", permCycles()[j]);
		
		printf("%s\n", buffer);
	} while (ProceedToNextMapping());
#endif
#endif
}

CC alldata::~alldata() {
	delete[] m_rowTime;
	delete[] m_pResults;
	delete[] m_pResultsPrev;
	delete[] selPlayers;
	delete[] m_pLinks;
	delete[] m_Km;
	delete[] m_Km2ndRowInd;
	delete[] m_Ktmp;
	delete[] m_trmk;
	delete[] m_groups;
	delete[] m_DayIdx;
	delete[] m_p3fSecondRows;
	delete m_pCheckCanon;
	delete m_pOrbits;
	delete[] m_tx;
	delete[] m_cycles;
#if !USE_CUDA
	FCLOSE_F(m_file);
#endif
}

CC void alldata::Init() {
	memset(links(), unset, m_numPlayers * m_numPlayers);
	memset(result()+ m_numPlayers, 0, m_nLenResults-m_numPlayers);
	for (int i = 0; i < m_numPlayers; i++)
		result()[i] = i;
	memcpy(m_pResultsPrev, result(), m_nLenResults);
	initCheckByGroup(1, 0);
	p1fSetTableRow(p1ftable(), result());

	maxDays = -1;
	nLoops = iDay = 0;
	noMoreResults = bPrevResult = false; // can be false, or true to go to prev day
	cnvInit();
}

#if !USE_CUDA
bool alldata::FindIsomorphicBaseElements(const string& fn) {
	// Testing the equivalence of the sets of base element
	auto const lenSet = m_numPlayers / 3;
	CStorage<tchar> baseElements(2, lenSet);
	const auto nSets = readTable(fn, -1, lenSet, baseElements.getObjectsPntr(), 1, 2);
	if (!nSets) {
		printfRed("Cannot read file \"%s\"\n", fn.c_str());
		return false;
	}

	checkCommonValues(baseElements.getObject(), nSets);

	myExit(1);
	return true;
}
#endif