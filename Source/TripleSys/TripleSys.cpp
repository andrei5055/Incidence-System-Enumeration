#include <iostream>
#include "TripleSys.h"
#include "Table.h"

#ifdef CD_TOOLS
   #include "../CanonicityChecker.h"
#else
   #include "CheckCanon.h"
#endif

Table<char>* pRes = NULL;

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinksV, bool useCheckLinksH) :
	SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
	m_np2(numPlayers * numPlayers),
	m_nGroups(numPlayers / groupSize),
	m_bCheckLinkV(groupSize == 3 && useCheckLinksV),
	m_bCheckLinkH(groupSize == 3 && useCheckLinksH) {
	m_nLenResults = m_numDays * numPlayers;
	maxResult = new char[m_nLenResults];
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[5 * m_numPlayers];
	tmpPlayers = selPlayers + m_numPlayers;
	indexPlayer = tmpPlayers + m_numPlayers;
	m_h = indexPlayer + m_numPlayers;
	m_ho = m_h + m_numPlayers;
	m_groupSizeFactorial = factorial(m_groupSize);
	int n = 1, m = m_groupSizeFactorial;
	for (int j = 2; j <= m_nGroups;j++)
	{
		n *= j; m *= m_groupSizeFactorial;
	}
	//n = m = 50;
	m_nallTr = n;
	m_nallTg = m;
	m_finalKMindex = 0;
	m_allTr = new char[n * m_nGroups];
	m_allTg = new char[m * m_nGroups];
	m_bestTr = 0;
	m_bestTg = 0;
	m_Km = new char[m_numPlayers * m_numDays * 2];
	m_trmk = new char[m_numPlayers];
	m_groups = new char[m_groupSizeFactorial * m_groupSize];
	m_pLinks = new char[m_np2];
	m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers, m_groupSize);

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
	delete[] m_Km;
	delete[] m_trmk;
	delete[] m_groups;
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

void _printf(FILE* f, bool toScreen, const char* format, const char* pStr) {
	if (f)
		fprintf(f, format, pStr);

	if (toScreen)
		printf(format, pStr);
}

void alldata::outputResults(int iDay, const unsigned char *pResult, int cntr) const
{
	char buffer[256];
	const bool toScreen = PrintImprovedResults > 1;
	FOPEN_W(f, ImprovedResultFile, canonCalls(1) || cntr? "a" : "w", m_file);

	const unsigned char* pDayPerm = NULL;
	auto flag = true;
	if (cntr) {	
		flag = m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainMatr);
		if (flag) {
			pDayPerm = pResult + iDay * numPlayers();
			sprintf_s(buffer, "Improved Result #%d for %d days  m_groupIndex = %d:\n", cntr, iDay, m_groupIndex);
			_printf(f, toScreen, buffer);
		}

		if (m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainTxt))
			_printf(f, toScreen, "%s\n", m_pCheckCanon->comment());
	}
	else {
		sprintf_s(buffer, "Initial Result #%zd (%zd):\n", canonCalls(0), pRes->counter());
		_printf(f, toScreen, buffer);
	}

	if (flag)
		outMatrix(pResult, iDay, numPlayers(), m_groupSize, 0, f, false, toScreen, cntr, pDayPerm);

	FCLOSE_W(f, m_file);
}

void alldata::outputError() const {
	extern char lastError[];
	FOPEN_W(f, ImprovedResultFile, "a", m_file);
	_printf(f, false, lastError);
	FCLOSE_W(f, m_file);
}

bool alldata::Run(int improveResult) {
	// Input parameter:
	//      improveResult: 
	//           0 (default) - do not try to improve given "result"; 
	//		   !=0 - m_pCheckCanon will return the permutations of the sets of players 
	//               and days that improve given "results";
	//          >1 - try to improve the “results” as much as possible.
	clock_t rTime, mTime, iTime = clock();
	unsigned char* bResults = NULL;
	const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
	rTime = mTime = iTime;
	if (improveResult)
		bResults = new unsigned char[(improveResult > 1 ? 2 : 1) * lenResult];
	Table<char> Result("Result table", m_numDays, m_numPlayers, 0, GroupSize, true, true);
	pRes = &Result;
	if (iDay == m_numDays)
	{
		char* bRes1 = NULL;
		const auto bResults_1 = new unsigned char[2 * lenResult];
		m_pCheckCanon->setPreordered(false);
		kmFullSort(result(), iDay, m_numPlayers, m_groupSize);
		bRes1 = m_Km;
#if 1   // Change to 0 to use "improveMatrix"
		while(!cnvCheck())
#else
        if (improveMatrix(2, (unsigned char *)bResults_1, lenResult, (unsigned char **)& bRes1))
#endif
		{
			printTable("Result improved", (const char*)bRes1, iDay, m_numPlayers, 0, m_groupSize, true);
			memcpy(result(0), bRes1, m_nLenResults);
			int errLine = 0, errGroup = 0, dubLine = 0;
			if (!_CheckMatrix(result(0), iDay, m_numPlayers, links(0), true, &errLine, &errGroup, &dubLine))
			{
				printf("Dublicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
				abort();
			}
			printTableColor("Links improved", links(0), numPlayers(), numPlayers());
		}
		delete[] bResults_1;
		m_pCheckCanon->setPreordered(true);
	}
	const auto numDaysAdj = CalcOnlyNFirstLines != 0 ? CalcOnlyNFirstLines : m_numDays;
	while (nLoops < LoopsMax)
	{
		clock_t dTime = clock();
		mTime = clock();
		while (iDay < numDaysAdj || bPrevResult)
		{
			if (iDay < 0)
			{
				noMoreResults = true;
				break;
			}
			if (bPrevResult)
			{
				if (!initPrevDay())
					continue;
			}
			else if (!initCurrentDay())
				continue;
ProcessOneDay:
			if (!processOneDay())
			{
				bPrevResult = true;
				continue;
			}

			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			clock_t cTime = clock();
			if (maxDays < iDay || cTime - rTime > ReportInterval)
			{
				if (1)//CalcOnlyNFirstLines == 0)
				{
					m_pCheckLink->reportCheckLinksData();
					printf("Current result for matrix %.0f: rows=%d, build time=%d, time since start=%d\n",
						nLoops + 1, iDay + 1, cTime - mTime, cTime - iTime);
					printTable("Current result", result(), iDay + 1, m_numPlayers, 0, m_groupSize, true);
				}
				rTime = cTime;
				maxDays = iDay;
				memcpy(maxResult, result(0), m_nLenResults);
			}
			iDay++;
#if UseSS != 2
			if (iDay > 1)
			{
#if 1 // set to 0 to disable improvement 
				bPrevResult = improveMatrix(improveResult, bResults, lenResult);
				if (bPrevResult)
				{
					if (m_groupIndex / m_nGroups < iDay - 1)
						printf("*** More then one day back (from group %d to %d) request\n", iDay * m_nGroups, m_groupIndex);
					else if (m_groupIndex >= iDay * m_nGroups)
					{
						printf("*** After current day 'back' request\n");
						abort();
					}
					else
					{
						if (iDay == m_numDays && m_groupIndex > (iDay - 1) * m_nGroups - 2)
							m_groupIndex = (iDay - 1) * m_nGroups - 2;
						iDay--;
						bPrevResult = false;
						while(iDay * m_numPlayers + iPlayer >= (m_groupIndex + 1) * m_groupSize)
						{
							getPrevPlayer();
						}
						goto ProcessOneDay;
					}
				}
#endif
			}
#endif
		}

		if (noMoreResults)
		{
			memcpy(result(0), maxResult, m_nLenResults);
			printf("no more results\n");
			break;
		}
		if (CalcOnlyNFirstLines > 0)
			memcpy(maxResult, result(0), m_nLenResults);
		nLoops++;
		if (iDay < numDaysAdj)
			abort();
		if (CalcOnlyNFirstLines == 0 || nLoops >= LoopsMax)
		{
			//printf("*** cnvCheck starts\n");
			if (CHECK_WITH_GROUP && USE_TRANSLATE_BY_LEO == 0 || cnvCheck())
			{
				//report result
				clock_t cTime = clock();
				printf("Result %.0f: matrix build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);
				Result.printTable(result(), true, ResultFile);

				//memcpy(m_Km + m_finalKMindex * m_numDays * m_numPlayers, result(), m_numDays * m_numPlayers);
				m_finalKMindex++;
			}
		}
		bPrevResult = true;
	}
	//imStatReport();
	printf("\nA total of %d non-isomorphic matrices (%dx%d, group size=%d) were selected (from %.0f calculated)\n",
		m_finalKMindex, m_numPlayers, m_numDays, m_groupSize, nLoops);
	printf("Total time = %d ms\n", clock() - iTime);
	if (CalcOnlyNFirstLines > 0 && nLoops < LoopsMax)
	{
		clock_t cTime = clock();
		printf("Result %.0f: matrix build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);
		printTable("Last result", result(), m_numDays, m_numPlayers, 0, m_groupSize, true);
	}

	delete[] bResults;
	return true;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
