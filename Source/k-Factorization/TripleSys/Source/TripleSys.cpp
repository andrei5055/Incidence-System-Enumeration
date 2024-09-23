#include <iostream>
#include "TripleSys.h"
#include "Table.h"

#ifdef CD_TOOLS
   #include "CanonicityChecker.h"
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

	TableAut Result("|Aut(M)|", m_numDays, m_numPlayers, 0, GroupSize, true, true);
	Result.allocateBuffer(32);

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
	else if (iDay > 0)
	{
		testImproveMatrixSpeed();
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
							threadNumber, nLoops + 1, iDay + 1, cTime - rTime, cTime - iTime);
						printTable("Current result", result(), iDay + 1, m_numPlayers, 0, m_groupSize, true);
#if 0
						TestStatPrint("iDay=%d:", iDay);
#endif
					}
					rTime = cTime;
					maxDays = iDay;
				}
			}
			iDay++;
#if UseSS != 2
			if (iDay > 1)
			{
				if (m_bCheckLinkT && !m_pCheckLink->checkLinksT(links(), iDay))
				{
					m_groupIndex = iDay * m_nGroups - 2;
					bPrevResult = true;
				}
#if 1 // set to 0 to disable improvement by improveMatrix and cnvCheckNew/cnvCheck
				if (!bPrevResult)
				bPrevResult = improveMatrix(m_improveResult, bResults, lenResult);
#if 1
				if (!bPrevResult && (
					(!CHECK_WITH_GROUP && iDay == m_numDays)
				//	|| (iDay > m_numDays - 2)
					|| (iDay == numDaysAdj && numDaysAdj != m_numDays)
					))
				{
#if USE_cnvCheckNew
					bPrevResult = !cnvCheckNew();
#else
					bPrevResult = !cnvCheck();
#endif
				}
#endif
				if (bPrevResult)
				{
					if (m_groupIndex >= iDay * m_nGroups)
					{
						printf("*** Thread %d: After current day 'back' request (m_groupIndex=%d)\n", threadNumber, m_groupIndex);
						printTable("Input matrix", result(), iDay, m_numPlayers, 0, m_groupSize, true);
						abort();
					}
					else
					{
						if (iDay == m_numDays && m_groupIndex > (iDay - 1) * m_nGroups - 2)
							m_groupIndex = (iDay - 1) * m_nGroups - 2;
						iDay--;
#if UseTrMask == 1
						if (iDay > 0)
							memcpy(m_TrMask + iDay * m_nTrBytes, m_TrMask + (iDay - 1) * m_nTrBytes, m_nTrBytes);
#endif
						while (iDay > m_groupIndex / m_nGroups)
						{
							while (iPlayer >= 0)
							{
								getPrevPlayer();
							}
							initPrevDay();
						}
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
			m_finalKMindex++;
			{
				if (bPrint)
				{
					//report result
					clock_t cTime = clock();
					printf("%5zd: Result Matrix, build time=%d, time since start=%d\n", nLoops, cTime - mTime, cTime - iTime);

				}
				Result.setGroupOrder(reportMatrixStats(bPrint));
				Result.printTable(result(), true, ResultFile, bPrint, numDaysAdj);
#if 0
				TestStatPrint("iDay=%d:", iDay);
#endif
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
			printf("\nThread %d: %d non-isomorphic matrices (%d,%d,%d) created\n",
				threadNumber, m_finalKMindex, m_numPlayers, numDaysAdj, m_groupSize);
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
