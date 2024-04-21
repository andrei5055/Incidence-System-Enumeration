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
	ImprovedResultFile[0] = '\0';
	ResultFile[0] = '\0';
	m_groupSizeFactorial = factorial(m_groupSize);
	int n = 1, m = m_groupSizeFactorial;
	for (int j = 2; j <= m_nGroups;j++)
	{
		n *= j; m *= m_groupSizeFactorial;
	}
	n = m = 50;
	m_nallTr = n;
	m_nallTg = m;
	m_finalKMindex = 0;
	m_allTr = new char[n * m_nGroups];
	m_allTg = new char[m * m_nGroups];
	m_bestTr = 0;
	m_bestTg = 0;
	m = n = 2;
	m_TrTest = new int[n];
	memset(m_TrTest, 0, n * sizeof(m_TrTest[0]));
	m_TgTest = new int[m];
	memset(m_TgTest, 0, m * sizeof(m_TgTest[0]));
	m_Km = new char[m_numPlayers * m_numDays * 2];
	m_KmSecondRow = m_Km + m_numPlayers;
	m_Km2 = m_Km + m_numPlayers * m_numDays;
	m_trmk = new char[m_numPlayers];
	m_groups = new char[m_groupSizeFactorial * m_groupSize];
	m_pLinks = new char[m_np2];
	m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers, m_groupSize);
	m_DayIdx = new unsigned char[m_numDays];

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

bool alldata::Run(int threadNumber, int iStartStopMode, int improveResult, 
	char* mStart0, char* mStart, char* mStop, int nrows, int mStep, double* pcnt, bool bPrint) {
	// Input parameter:
	//      improveResult: 
	//           0 (default) - do not try to improve given "result"; 
	//		   !=0 - m_pCheckCanon will return the permutations of the sets of players 
	//               and days that improve given "results";
	//          >1 - try to improve the “results” as much as possible.
	clock_t rTime, mTime, iTime = clock();
	double dNumMatricesMax = 0;
	int nBytesInStartMatrix = nrows * m_numPlayers;
	double startStopMatrixCount = 0;
	unsigned char* bResults = NULL;
	const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
	if (iStartStopMode == eCalcStartStop)
		dNumMatricesMax = *pcnt;
	rTime = mTime = iTime;
	if (improveResult)
		bResults = new unsigned char[(improveResult > 1 ? 2 : 1) * lenResult];
	Table<char> Result("Result table", m_numDays, m_numPlayers, 0, GroupSize, true, true);
	if (strlen(ImprovedResultFilePrefix) > 0)
		sprintf_s(ImprovedResultFile, "%s%04d.txt", ImprovedResultFilePrefix, threadNumber);
	if (strlen(ResultFilePrefix) > 0)
		sprintf_s(ResultFile, "%s%04d.txt", ResultFilePrefix, threadNumber);

	pRes = &Result;
	if (mStart0 != NULL)
	{
		linksFromMatrix(mStart0, nrows);
		iDay = nrows;
	}
	if (iDay == m_numDays)
	{
		char* bRes1 = NULL;
		const auto bResults_1 = new unsigned char[2 * lenResult];
		m_pCheckCanon->setPreordered(false);
		kmFullSort(result(), iDay, m_numPlayers, m_groupSize);
		bRes1 = m_Km;
#if 0   // Change to 0 to use "improveMatrix"
		while(!cnvCheckNew()) // cnvCheck (if 2 playes in group) works only with full matrix
#else
        if (improveMatrix(2, (unsigned char *)bResults_1, lenResult, (unsigned char **)& bRes1))
#endif
		{
			printTable("Result improved", (const char*)bRes1, iDay, m_numPlayers, 0, m_groupSize, true);
			memcpy(result(0), bRes1, m_nLenResults);
			int errLine = 0, errGroup = 0, dubLine = 0;
			if (!CheckMatrix(result(0), iDay, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
			{
				printf("Duplicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
				abort();
			}
			//printTableColor("Links improved", links(0), numPlayers(), numPlayers());
		}
		delete[] bResults_1;
		m_pCheckCanon->setPreordered(true);
	}
	const auto numDaysAdj = iStartStopMode == eCalcResult ? m_numDays : nrows;
	while (nLoops < LoopsMax)
	{
		if (threadNumber == 2)
			nLoops = nLoops;
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
			if (bPrint)
			{
				clock_t cTime = clock();
				if (maxDays < iDay || cTime - rTime > ReportInterval)
				{
					if (1)//nrows == 0)
					{
						//m_pCheckLink->reportCheckLinksData();
						printf("Thread %d: Current result for matrix %.0f: rows=%d, build time=%d, time since start=%d\n",
							threadNumber, nLoops + 1, iDay + 1, cTime - mTime, cTime - iTime);
						printTable("Current result", result(), iDay + 1, m_numPlayers, 0, m_groupSize, true);
					}
					rTime = cTime;
					maxDays = iDay;
					memcpy(maxResult, result(0), m_nLenResults);
				}
			}
			iDay++;
#if UseSS != 2
			if (iDay > 1)
			{
				if (iStartStopMode == eCalcResult && iDay == nrows && memcmp(result(0), mStop, nBytesInStartMatrix) > 0)
				{
					noMoreResults = true;
					break;
				}
#if 1 // set to 0 to disable improvement 
				bPrevResult = improveMatrix(improveResult, bResults, lenResult);
				if (bPrevResult)
				{
					if (m_groupIndex / m_nGroups < iDay - 1)
					{
						printf("*** Thread %d: More than one day back (from group %d to %d) request. Task exit.\n", 
							threadNumber, iDay * m_nGroups, m_groupIndex);
						printTable("Input matrix", result(), iDay, m_numPlayers, 0, m_groupSize, true);
						//m_groupIndex = iDay * m_nGroups - 2;
						//exit(0);
					}
					/* else */ if (m_groupIndex >= iDay * m_nGroups)
					{
						printf("*** Thread %d: After current day 'back' request\n", threadNumber);
						printTable("Input matrix", result(), iDay, m_numPlayers, 0, m_groupSize, true);
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
			if (bPrint)
				printf("no more results\n");
			break;
		}
		if (iDay < numDaysAdj)
			abort();
		if (iStartStopMode == eCalcStartStop)
		{
			if (startStopMatrixCount >= dNumMatricesMax)
				break;
			memcpy(mStop, result(), nBytesInStartMatrix);
			if (startStopMatrixCount == 0)
				memcpy(mStart, result(), nBytesInStartMatrix);
			startStopMatrixCount++;
			if (startStopMatrixCount == mStep)
			{
				dNumMatricesMax -= mStep;
				startStopMatrixCount = 0;
				mStart += nBytesInStartMatrix;
				mStop += nBytesInStartMatrix;
			}
		}
		nLoops++;
		if (iStartStopMode == eCalcResult)
		{
			//printf("*** cnvCheck starts\n");
			if (CHECK_WITH_GROUP || cnvCheckNew())
			{
				if (bPrint)
				{
					//report result
					const clock_t cTime = clock();
#if PRINT_MATR_CNTR
					extern int matr_cntr, perm_cntr;
					printf("Result %.0f (%3d, %3d)\n", nLoops, matr_cntr, perm_cntr);
#else
					printf("Result %.0f: matrix build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);
#endif
					Result.printTable(result(), true, ResultFile);
				}
				m_finalKMindex++;
				if (pcnt != NULL)
					*pcnt = -m_finalKMindex - 1;
			}
			if (pcnt != NULL)
				*(pcnt + 1) = nLoops;
		}
		bPrevResult = true;
	}
	//imStatReport();
	if (iStartStopMode == eCalcResult)
	{
		if (bPrint)
		{
			printf("\nThread %d: %d non-isomorphic matrices (%d,%d,%d) were selected (from %.0f generated)\n",
				threadNumber, m_finalKMindex, m_numPlayers, m_numDays, m_groupSize, nLoops);
			printf("Thread execution time = %d ms\n", clock() - iTime);
		}
		if (pcnt != NULL)
		{
			*pcnt = m_finalKMindex;
			*(pcnt + 1) = nLoops;
		}
	}
	else
	{
		if (bPrint)
		    printf("%.0f matrices, time since start=%d\n", nLoops, clock() - iTime);
		if (pcnt != NULL)
			*pcnt = nLoops;
	}
	delete[] bResults;
	return true;
}

void printTransformed(int nrows, int ncols, const char* tr, const char* ttr, const char *pImatr, const char* pTmatr, int numRow, const double nLoops, int finalKMindex)
{
	printf("Calculated Matrix %.0f can be improved. See below (nKm=%d row=%d)\n", nLoops, finalKMindex, numRow);
	printTable("Tr source", tr, 1, ncols);
	printTable("Tr actual", ttr, 1, ncols);
	printTable("Original", pImatr, nrows, ncols);
	printTable("Translated", pTmatr, nrows, ncols);
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
