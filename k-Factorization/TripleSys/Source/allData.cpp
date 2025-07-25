#include "TripleSys.h"
#include "CheckCanon.h"

using namespace std;

CC alldata::alldata(const SizeParam& p, const kSysParam* pSysParam, int createSecondRow, CRowStorage* pRowStorage,
	bool useCheckLinksT, int improveResult, bool createImprovedMatrix) :
	CGroupInfo(pSysParam->val[t_numPlayers], 200), CGroupUtilisation(pSysParam),
	CycleSupport(pSysParam->val[t_numPlayers], !pSysParam->completeGraph()), CChecklLink(p, pSysParam),
	m_nGroups(numPlayers() / m_groupSize),
	m_bCheckLinkV(m_groupSize == 3 && param(t_useCheckLinksV)),
	m_bCheckLinkT(m_groupSize == 3 && useCheckLinksT),
	m_nLenResults(m_numDays * numPlayers()) {
	m_improveResult = improveResult;
	if (!(m_numDaysResult = pSysParam->val[t_nRowsInResultMatrix]))
		m_numDaysResult = m_numDays;
	m_createSecondRow = createSecondRow;
	m_rowTime = new int[m_numDaysResult];
	m_rowTime[0] = 0;
	m_pResults = new tchar[m_nLenResults];
	m_pResultsPrev = new tchar[m_nLenResults * 2];
	m_pResultsPrev2 = m_pResultsPrev + m_nLenResults;
	const auto numPlayers64 = (m_numPlayers + 7) / 8 * 8;

	selPlayers = new tchar[7 * numPlayers64];
	tmpPlayers = selPlayers + numPlayers64;
	indexPlayer = tmpPlayers + numPlayers64;
	m_h = indexPlayer + numPlayers64;
	m_ho = m_h + numPlayers64;
	m_indexPlayerMin = m_ho + numPlayers64;
	m_indexPlayerMax = m_indexPlayerMin + numPlayers64;
	m_groupSizeRemainder = new tchar[m_numPlayers];
	for (int i = 0; i < m_numPlayers; i++)
	{
		m_groupSizeRemainder[i] = i % m_groupSize;
		m_indexPlayerMax[i] = m_groupSizeRemainder[i] == 0 ? i : m_numPlayers - 1;
		m_indexPlayerMin[i] = i / m_groupSize + m_groupSizeRemainder[i];
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
	m_pLinks = new tchar[np2 * 2];
	m_pSecondRowsDB = NULL;
	m_DayIdx = new tchar[m_numDays];
	m_numCycles = 0;
	m_firstNotSel = 0;
	bool bCBMP = !pSysParam->completeGraph();
	m_matrixCanonInterval = param(t_matrixCanonInterval);
	m_checkForUnexpectedCycle = !m_allowUndefinedCycles && m_groupSize == 2 && m_numPlayers > 4 &&
		param(t_u1f) && (!pSysParam->u1fCycles[0] || pSysParam->u1fCycles[0][1] > 4); // assume all cycle sets are sorted
	if (param(t_rejectCycleLength) == 4)
		m_checkForUnexpectedCycle = true; // reject cycle length 4

	iPlayerIni = m_numPlayers - 1;
	if (m_numPlayers > m_groupSize)
		iPlayerIni -= m_groupSize;

#if !USE_EQUAL && !CHECK_WITH_GROUP
	auto nDays = m_numDays;
	while (nDays--) m_DayIdx[nDays] = nDays;
#endif

	if (param(t_useImproveMatrix)) {
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
	m_maxCommonVSets = m_groupSize == 2 ? 1 : (bCBMP ? m_groupSizeFactorial : (int)pow(6.0, m_nGroups / 3.0) + 1);
	InitCycleSupport(m_nGroups, m_maxCommonVSets);

	const auto u1fPntr = sysParam()->u1fCycles[0];
	
	const auto nPreconstructedRows = param(t_useRowsPrecalculation);
	if (nPreconstructedRows) {
		if (m_bRowStorageOwner = (pRowStorage == NULL))
			m_pRowStorage = new CRowStorage(sysParam(), m_numPlayers, 30000, this);
		else
			m_pRowStorage = pRowStorage;

		if (pRowStorage || !(param(t_MultiThreading) == 2 && param(t_useRowsPrecalculation)))
			m_pRowUsage = new CRowUsage(m_pRowStorage);
	}
	m_TrCyclesAll = new TrCycles[MAX_CYCLE_SETS];
	m_TrCyclesPair = new TrCycles[MAX_CYCLE_SETS];
	m_TrCyclesFirst2Rows = new TrCycles[MAX_CYCLE_SETS];
	bool bp1f = ((!u1fPntr || u1fPntr[1] == m_numPlayers) && !m_allowUndefinedCycles);
	m_pCheckFunc = NULL;
	if (m_use2RowsCanonization) {
		if (m_groupSize == 2) {
			if (param(t_useFastCanonizerForG2))
				m_pCheckFunc = &alldata::cnvCheck3U1F;
			else 
				m_pCheckFunc = (bp1f && !bCBMP) ? &alldata::cnvCheck2P1F : &alldata::cnvCheck2U1F;
		}
		else if (m_groupSize == 3 || bCBMP)
			m_pCheckFunc = &alldata::cnvCheck3U1F;
	}
	else if (!bCBMP &&
		((m_numPlayers == 16 && m_groupSize == 4) || (m_numPlayers == 25 && m_groupSize == 5)))
		m_pCheckFunc = &alldata::cnvCheck45;

	if (m_createSecondRow)
		m_pCheckFunc = &alldata::cnvCheck3U1F;
	if (param(t_useFastCanonizerForG2) == 2)
		m_pSortGroups = NULL;

	m_pSortGroups = m_groupSize == 2 ? &alldata::kmSortGroups2 : (m_groupSize == 3 ? &alldata::kmSortGroups3 : &alldata::kmSortGroups);

	if (bCBMP)
		m_pProcessMatrix = &alldata::kmProcessMatrix;
	else
		m_pProcessMatrix = createImprovedMatrix || m_groupSize > 3 ? &alldata::kmProcessMatrix : (m_groupSize == 2 ? &alldata::kmProcessMatrix2 : &alldata::kmProcessMatrix3);

#if !USE_CUDA
	const auto* binaryCanonRows = sysParam()->strVal[t_binaryCanonizer];
	if (binaryCanonRows) {
		m_ppBinMatrStorage = new CBinaryMatrixStorage * [m_numDays + 1];
		memset(m_ppBinMatrStorage, 0, (m_numDays + 1) * sizeof(m_ppBinMatrStorage[0]));
		const auto numGroups = m_numPlayers / m_groupSize;
		string s(*binaryCanonRows);
		size_t pos = 0;
		bool flagOK = true;
		bool lastStorageNeeded = false;
		while (pos <= s.size()) {
			int idx1 = -1;
			pos = s.find(",");
			auto line = s.substr(0, pos);
			trim(line);
			if (!is_number(line)) {
				const auto pos1 = s.find("-");
				if (pos1 != string::npos) {
					auto line1 = line.substr(pos1 + 1);
					if (is_number(line1)) {
						idx1 = atoi(line1.c_str());
						line = line.substr(0, pos1);
					}
				}
				if (idx1 < 0) {
					flagOK = false;
					break;
				}
			} 
			s.erase(0, pos + 1);
			auto idx = atoi(line.c_str());
			lastStorageNeeded |= idx >= 1;
			if (idx1 < 0)
				idx1 = idx;
			else
				if (idx > idx1) {
					flagOK = false;
					break;
				}
					

			if (idx1 < 2 || idx >= m_numDays)
				continue;

			while (idx <= idx1) {
				if (idx > 1 && !m_ppBinMatrStorage[idx]) {
					const auto lenBinaryMatrix = m_numPlayers * idx * numGroups;
					m_ppBinMatrStorage[idx] = new CBinaryMatrixStorage(lenBinaryMatrix, 50 * lenBinaryMatrix);
				}
				idx++;
			}
		}

		if (flagOK && lastStorageNeeded) {
			const auto idx = m_numDays;
			if (!m_ppBinMatrStorage[idx]) {
				const auto lenBinaryMatrix = m_numPlayers * m_numDays * numGroups;
				m_ppBinMatrStorage[idx] = new CBinaryMatrixStorage(lenBinaryMatrix, 50 * lenBinaryMatrix);
			}
		}
		else {
			releaseBinaryMatricesStorage();
			if (!flagOK)
				printfRed("*** Unexpected value for \"UseBinaryCanonizer\": \"%s\"\n", binaryCanonRows->c_str());
		}
	}
#endif

#if Use_GroupOrbits
	m_pOrbits = new CGroupOrbits<unsigned char>(m_numPlayers);
#endif
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
	delete[] m_groupSizeRemainder;
	delete[] m_DayIdx;
	delete m_pCheckCanon;
	delete m_pOrbits;
	delete[] m_tx;
	delete[] m_TrCyclesAll;
	delete[] m_TrCyclesPair;
	delete[] m_TrCyclesFirst2Rows;
	if (m_bRowStorageOwner)
		delete m_pRowStorage;
	delete m_pRowUsage;
	releaseBinaryMatricesStorage();
#if !USE_CUDA
	FCLOSE_F(m_file);
#endif
}

CC void alldata::Init() {
	memset(links(), unset, m_numPlayers * m_numPlayers);
	memset(result() + m_numPlayers, 0, m_nLenResults - m_numPlayers);
	for (int i = 0; i < m_numPlayers; i++)
		result()[i] = i;
	memcpy(m_pResultsPrev, result(), m_nLenResults);
	memcpy(m_pResultsPrev2, result(), m_nLenResults);
	initCheckByGroup(1, 0);
	u1fSetTableRow(neighbors(), result());

	maxDays = -1;
	nLoops = iDay = 0;
	noMoreResults = bPrevResult = false; // can be false, or true to go to prev day
	cnvInit();
}

CC void alldata::releaseBinaryMatricesStorage() {
	if (!m_ppBinMatrStorage)
		return;

	for (int i = 0; i < m_numDays; i++)
		delete m_ppBinMatrStorage[i];

	delete[] m_ppBinMatrStorage;
	m_ppBinMatrStorage = NULL;
}

#if !USE_CUDA
bool alldata::FindIsomorphicBaseElements(const string& fn) {
	// Testing the equivalence of the sets of base element
	auto const lenSet = m_numPlayers / 3;
	CStorage<tchar> baseElements(2, lenSet);
	int reserved = 2;
	const auto nSets = readTable(fn, -1, lenSet, 1, 0, baseElements.getObjectsPntr(), reserved, param(t_nMaxNumberOfStartMatrices));
	if (!nSets) {
		printfRed("Cannot read file \"%s\"\n", fn.c_str());
		return false;
	}

	checkCommonValues(baseElements.getObject(), nSets);

	myExit(1);
	return true;
}
#endif
