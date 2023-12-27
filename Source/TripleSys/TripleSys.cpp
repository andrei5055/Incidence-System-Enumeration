#include <iostream>
#include "TripleSys.h"
#include "Table.h"

#ifdef CD_TOOLS
   #include "../CanonicityChecker.h"
#else
   #include "CheckCanon.h"
#endif

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinksV, bool useCheckLinksH) :
	SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
	m_np2(numPlayers * numPlayers),
	m_nGroups(numPlayers / groupSize),
	m_bCheckLinkV(groupSize == 3 && useCheckLinksV),
	m_bCheckLinkH(groupSize == 3 && useCheckLinksH) {
	m_nLenResults = m_numDays * numPlayers;
	maxResult = new char[m_nLenResults];
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[m_numPlayers];
	tmpPlayers = new char[m_numPlayers];
	indexPlayer = new char[m_numPlayers];
	m_pLinks = new char[m_np2];
	index6 = new char[m_numDays];
	m_h = new char[m_numPlayers];
	m_ho = new char[m_numPlayers];
	m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers);

#ifdef CD_TOOLS
	m_pCheckCanon = new CCanonicityChecker<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize, t_kSystems);
#else
	m_pCheckCanon = new CCheckerCanon<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize);
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
	delete[] tmpPlayers;
	delete[] indexPlayer;
	delete[] m_pLinks;
	delete[] index6;
	delete[] m_h;
	delete[] m_ho;
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

	iDay++;
	const unsigned char* pDayPerm = NULL;
	auto flag = true;
	if (cntr) {	
		flag = m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainMatr);
		if (flag) {
			pDayPerm = pResult + iDay * numPlayers();
			sprintf_s(buffer, "Improved Result #%d for %d days:\n", cntr, iDay);
			_printf(f, toScreen, buffer);
		}

		if (m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainTxt))
			_printf(f, toScreen, "%s\n", m_pCheckCanon->comment());
	}
	else {
		sprintf_s(buffer, "Initial Result #%zd:\n", canonCalls(0));
		_printf(f, toScreen, buffer);
	}

	if (flag)
		outMatrix(pResult, iDay, numPlayers(), m_groupSize, 0, f, false, toScreen, cntr, pDayPerm);

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
#if 1
	if (improveResult)
		bResults = new unsigned char[(improveResult > 1? 2 : 1) * lenResult];
#else
	bResults = new unsigned char[2 * lenResult];
#endif
	Table<char> Result("Result table", m_numDays, m_numPlayers, 0, GroupSize, true, true);

	if (iDay == m_numDays)
	{
		unsigned char* bRes1 = NULL;
		const auto bResults_1 = new unsigned char[2 * lenResult];
		iDay = iDay - 1;
		if (improveMatrix(2, bResults_1, lenResult, &bRes1))
		{
			printTable("Result improved", (const char*)bRes1, iDay + 1, m_numPlayers, 0, 3, true);
			memcpy(result(0), bRes1, m_numPlayers * (iDay + 1));
			memset(links(0), unset, m_numPlayers * m_numPlayers);
			for (int j = 0; j <= iDay; j++)
			{
				char* c = links(j);
				for (int i = 0; i < m_numPlayers; i = i + 3)
				{
					if (!setLinksForOnePlayer(result(j), i + 1, 1) ||
						!setLinksForOnePlayer(result(j), i + 2, 1))
						abort();
				}
			}
		}
		delete[] bResults_1;
		iDay = iDay + 1;
	}

	while (nLoops < LoopsMax)
	{
		mTime = clock();
		while (iDay < m_numDays || bPrevResult)
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

			if (!processOneDay())
			{
				bPrevResult = true;
				continue;
			}

			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			clock_t cTime = clock();
			if (maxDays < iDay || cTime - rTime > ReportInterval)
			{
				/**/
				m_pCheckLink->reportCheckLinksData();
				printf("Current result for matrix %d: days=%d, build time=%d, time since start=%d\n", 
					nLoops + 1, iDay + 1, cTime - mTime, cTime - iTime);
				Result.printTable(result());
				//printTable("Links", links[0], m_numPlayers, m_numPlayers);
				/**/
				rTime = cTime;
				maxDays = iDay;
				memcpy(maxResult, result(0), m_nLenResults);
				sortLinks();
			}
#if UseSS == 0
			if (iDay > 0)
			{
				bPrevResult = improveMatrix(improveResult, bResults, lenResult);
			}
#endif
			iDay++;
		}

		if (noMoreResults)
		{
			printf("no more results\n");
			printTableColor("Links", links(), m_numPlayers, m_numPlayers);
			break;
		}

		nLoops++;
		if (iDay < m_numDays)
			abort();
		//report result
		clock_t cTime = clock();
		printf("Result %d: matrix build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);
		//printTable("Links", links(), m_numPlayers, m_numPlayers);
		Result.printTable(result(), true, ResultFile);
		if (nLoops >= LoopsMax)
			break;
		bPrevResult = true;
	}
	printf("\nA total of %d %d-configurations were built.\n", nLoops, GroupSize);
	printf("Total time = %d ms\n", clock() - iTime);
	if (nLoops == 1)
	{
		if (memcmp(maxResult, result(), m_nLenResults) != 0)
			printTable("'Maximum days' Result", maxResult, m_numDays, m_numPlayers);
		printTableColor("Links", links(), m_numPlayers, m_numPlayers);
		convertLinksToResult(links());
		if (memcmp(m_co, result(), m_nLenResults) != 0)
			printTable("Result from link (different than result)", m_co, m_numDays, m_numPlayers);
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
