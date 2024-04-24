#include <iostream>
#include "TripleSys.h"
#include "Table.h"

#ifdef CD_TOOLS
   #include "../CanonicityChecker.h"
#else
   #include "CheckCanon.h"
#endif

#if 0
#ifdef _MSC_VER
#include <intrin.h> /* for rdtscp and clflush */
#pragma optimize("gt",on)
#else
#include <x86intrin.h> /* for rdtscp and clflush */
#endif
#endif
//#include <windows.h>

Table<char>* pRes = NULL;

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

bool alldata::Run(int threadNumber, int iCalcMode, 
	char* mStart0, char* mStart, int nrowsStart, int nrowsOut, sLongLong* pcnt, bool bPrint) {
	// Input parameter:
	//      improveResult: 
	//           0 (default) - do not try to improve given "result"; 
	//		   !=0 - m_pCheckCanon will return the permutations of the sets of players 
	//               and days that improve given "results";
	//          >1 - try to improve the “results” as much as possible.
	clock_t rTime, mTime, iTime = clock();
	int nMatricesMax = 0;
	int nBytesInStartMatrix = nrowsStart * m_numPlayers;
	int startMatrixCount = 0;
	unsigned char* bResults = NULL;
	const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
	int minRows = nrowsStart;
	int numDaysAdj = nrowsOut;
	rTime = mTime = iTime;

	//HANDLE hwnd;
	//hwnd = GetCurrentProcess(); // current process handle
	//FILETIME time1, time2, dum1, dum2, dum3;
	//LARGE_INTEGER time1, time2, dum1, dum2, dum3;

	if (iCalcMode == eCalcStart)
	{
		minRows = 0;
		numDaysAdj = nrowsOut;
		nMatricesMax = (int)(*pcnt);
	}
	else if (iCalcMode == eCalcResult)
	{
		minRows = nrowsStart;
		numDaysAdj = nrowsOut;
	}
	if (m_improveResult)
		bResults = new unsigned char[(m_improveResult > 1 ? 2 : 1) * lenResult];

	Table<char> Result("Result table", m_numDays, m_numPlayers, 0, GroupSize, true, true);

	createFolderAndFileName(ImprovedResultFile, sizeof(ImprovedResultFile), ImprovedResultFolder, ImprovedResultNameFormat,
		m_numPlayers, numDaysAdj, m_groupSize, threadNumber);
	createFolderAndFileName(ResultFile, sizeof(ResultFile), ResultFolder, ResultNameFormat,
		m_numPlayers, numDaysAdj, m_groupSize, threadNumber);

	pRes = &Result;
	if (mStart0 != NULL)
	{
		iDay = nrowsStart;
		memcpy(result(), mStart0, nrowsStart * m_numPlayers);
		linksFromMatrix(links(), mStart0, nrowsStart, m_numPlayers);
	}
	else if (iDay > 0) //== m_numDays && threadNumber == 0)
	{
		//uint64_t time1, time2;
		unsigned int junk = 0;
		char* bRes1 = NULL;
		const auto bResults_1 = new unsigned char[2 * lenResult];
		kmFullSort(m_Ktmp, result(), iDay, m_numPlayers, m_groupSize);
		bRes1 = m_Km;
		//time1 = __rdtscp(&junk);

		//GetProcessTimes(hwnd, &dum1, &dum2, &dum3, &time1);
		//QueryPerformanceCounter(&time1);
		//GetThreadTimes(hwnd, &dum1, &dum2, &dum3, &time1);

		clock_t tTime = clock();
		int improveResult = m_improveResult;
		int createImprovedResult = m_createImprovedResult;
		m_createImprovedResult = 2;
		m_improveResult = 2;
		int nTests = 100;
		for (int i = 0; i < nTests; i++)
		{
#if 0   // Change to 0 to use "improveMatrix"
			if (!cnvCheckNew())
#else
			bRes1 = 0;
			if (improveMatrix(m_improveResult, (unsigned char*)bResults_1, lenResult, (unsigned char**)&bRes1) && bRes1 != 0)
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
				printTableColor("Links improved", links(0), m_numPlayers, m_numPlayers);
			}
			else
			{
				m_createImprovedResult = createImprovedResult;
				m_improveResult = improveResult;
			}
		}
		m_createImprovedResult = createImprovedResult;
		m_improveResult = improveResult;
		//time2 = __rdtscp(&junk);
		//GetProcessTimes(hwnd, &dum1, &dum2, &dum3, &time2);
		//QueryPerformanceCounter(&time2);
		//GetThreadTimes(hwnd, &dum1, &dum2, &dum3, &time2);
		printf("+++ %.1f ms needed per one improvement check\n", (double(clock() - tTime)) / nTests);
		//printf(" %.1f ms (%.1f cpu) needed per one improvement check\n", (double(clock() - tTime)) / nTests);
		//	(time2.QuadPart - time1.QuadPart) / nTests / 1000.0);
		//	(double(time2.dwLowDateTime - time1.dwLowDateTime)) / nTests / 1000.0);
		delete[] bResults_1;
	}
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
			if (iDay < minRows && iCalcMode == eCalcResult)
			{
				noMoreResults = true;
				break;
			}
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
					if (1)//nrowsStart == 0)
					{
						//m_pCheckLink->reportCheckLinksData();
						printf("Thread %d: Current result for matrix %zd: rows=%d, build time=%d, time since start=%d\n",
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
#if 1 // set to 0 to disable improvement
				bPrevResult = improveMatrix(m_improveResult, bResults, lenResult);
#if 1
				if (!bPrevResult && (
					(!CHECK_WITH_GROUP && iDay == m_numDays)
					|| (iCalcMode == eCalcStart && iDay <= nrowsStart)
					|| (iDay == numDaysAdj && numDaysAdj != m_numDays)
				//	|| (iDay < m_numDays)
					))
				{
#if USE_cnvCheckNew
					bPrevResult = !cnvCheckNew();
#else
					bPrevResult = !cnvCheck();
#endif
					/**
					if (CHECK_WITH_GROUP && bPrevResult && iDay == m_numDays)
					{
						printTable("Input matrix", result(), iDay, m_numPlayers, 0, m_groupSize, true);
						printf("*** improvement!!!***\n");
						exit(0);
					}**/
					if (bPrevResult)
					    m_groupIndex = iDay * m_nGroups - 2;
				}
#endif
				if (bPrevResult)
				{
#if 0
					if (m_groupIndex < (iDay - 1) * m_nGroups - 1)
					{
						printf("*** Thread %d: More than one day back (from group %d to %d) request.\n", 
							threadNumber, iDay * m_nGroups, m_groupIndex);
						printTable("Input matrix", result(), iDay, m_numPlayers, 0, m_groupSize, true);
						//m_groupIndex = iDay * m_nGroups - 2;
						//exit(0);
					}
#endif
					if (m_groupIndex >= iDay * m_nGroups)
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
#if NEW_CODE && !USE_TRANSLATE_BY_LEO
						while (iDay > m_groupIndex / m_nGroups)
						{
							while (iPlayer >= 0)
							{
								getPrevPlayer();
							}
							initPrevDay();
						}
#endif
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
		if (iCalcMode == eCalcStart)
		{
			memcpy(mStart, result(), nBytesInStartMatrix);
			startMatrixCount++;
			mStart += nBytesInStartMatrix;
			if (startMatrixCount >= nMatricesMax)
				break;
		}
		nLoops++;
		if (iCalcMode == eCalcResult)
		{
			//printf("*** cnvCheck starts\n");
			//if (CHECK_WITH_GROUP || cnvCheckNew())
			{
				if (bPrint)
				{
					//report result
					clock_t cTime = clock();
					printf("Result %zd: matrix build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);
				}
				Result.printTable(result(), true, ResultFile, bPrint, numDaysAdj);
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
	if (iCalcMode == eCalcResult)
	{
		if (bPrint)
		{
			printf("\nThread %d: %d non-isomorphic matrices (%d,%d,%d) were selected (from %zd generated)\n",
				threadNumber, m_finalKMindex, m_numPlayers, numDaysAdj, m_groupSize, nLoops);
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
		    printf("%zd matrices, time since start=%d\n", nLoops, clock() - iTime);
		if (pcnt != NULL)
			*pcnt = nLoops;
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
