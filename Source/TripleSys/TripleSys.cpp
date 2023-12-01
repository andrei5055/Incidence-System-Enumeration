// TripleSys.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <windows.h>
#include "TripleSys.h"
#include "CanonicityChecker.h"

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinks) :
		SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
		m_np2(numPlayers * numPlayers),
		m_np3(m_np2 * numPlayers),
		m_nGroups(numPlayers / groupSize) {
	m_nLenResults = m_numDays * numPlayers;
	maxResult = new char[m_nLenResults];
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[m_numPlayers];
	tmpPlayers = new char[m_numPlayers];
	indexPlayer = new char[m_numPlayers];
	m_pLinks = new char[m_np2];
	index6 = new char[m_numDays];
	if (groupSize == 3 && useCheckLinks)
		m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers);

	m_pCheckCanon = new CCanonicityChecker<unsigned char, unsigned char>(m_numDays, numPlayers, 3);
	Init();
}

alldata::~alldata() {
	delete[] maxResult;
	delete[] m_pResults;
	delete[] selPlayers;
	delete[] tmpPlayers;
	delete[] indexPlayer;
	delete[] m_pLinks;
	delete[] index6;
	delete m_pCheckLink;
	delete m_pCheckCanon;
}

void alldata::Init() {
	int lastRep = GetTickCount();
	memset(links(), unset, m_numPlayers * m_numPlayers);
	memset(result(), 0, m_nLenResults);
	maxDays = -1;
	nLoops = 0;
	noResult = false;
	iDay = 0;
	iPrevResult = unset;
}

bool alldata::Run() {
	int iTime = GetTickCount();
	while (nLoops < LoopsMax)
	{
		while (iDay < m_numDays)
		{
			if (iPrevResult >= 0)
			{
				iDay = iPrevResult;
				initPrevDayProcess();
				iPrevResult = unset;
				if (iDay < 0)
				{
					noResult = true;
					goto finish;
				}
			}
			else
				initCurrentDay();

			while (iPlayer < m_numPlayers)
			{
				int iPlayerNumber = 0;
				if (iPlayer < 0)
				{
					int k = -1;
					if (iDay <= 0)
					{
						noResult = true;
						goto finish;
					}
					else
					{
						memset(result(iDay), 0, m_numPlayers);
						iDay--;
						iPrevResult = iDay;
						goto runday;
					}
					abort();
				}

				if (iPlayer >= m_numPlayers)
					abort();

				if (UseLastSix && iPlayer == m_numPlayers - 6)
				{
					iPlayerNumber = processLastSix();
					if (iPlayer == m_numPlayers)
						goto nextDay;
				}
				else
					iPlayerNumber = getNextPlayer();
				
				if (iPlayerNumber >= m_numPlayers)
				{
					//printf("\n%d  go back", iDay);
					getPrevPlayer();
					continue;
				}
				indexPlayer[iPlayer] = iPlayerNumber;
				tmpPlayers[iPlayer] = iPlayerNumber;
				if ((iPlayer + 1) < m_numPlayers)
					indexPlayer[(iPlayer + 1)] = 0;
				selPlayers[iPlayerNumber] = iPlayer;
				iPlayer++;
			}
nextDay:
			if (noResult)
				break;

			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			if (maxDays < iDay)
			{
				/**/
				printf("day %d  Time = %d\n", iDay, GetTickCount() - iTime);
				printTable("Result", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
				//printTable("Links", links[0], m_numPlayers, m_numPlayers);
				/**/
				maxDays = iDay;
				memcpy(maxResult, result(0), m_nLenResults);
			}
			if (iDay > 0)
			{
				if (!m_pCheckCanon->CheckCanonicity((unsigned char *)result(), iDay + 1))
				{
					// get new matrix
					iPrevResult = iDay - 1;
					continue;
				}
			}
			iDay++;
		runday:;
		}
	finish:

		nLoops++;

		if (noResult)
		{
			printf("no more results\n");
			break;
		}

		//report result
		printf("Result %d, Time = %d\n", nLoops, GetTickCount() - iTime);
		//printTable("Links", links(), m_numPlayers, m_numPlayers);
		printTable("Result table", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
		if (nLoops >= LoopsMax)
			break;
		iDay--;
		iPrevResult = iDay;
	}
	printf("Total time = %d ms\n", GetTickCount() - iTime);
	if (nLoops == 1)
	{
		if (memcmp(maxResult, result(), m_nLenResults) != 0)
			printTable("'Maximum days' Result", maxResult, m_numDays, m_numPlayers);
		printTable("Links", links(), m_numPlayers, m_numPlayers);
		convertLinksToResult(links());
		if (memcmp(m_co, result(), m_nLenResults) != 0)
			printTable("Result from link (different than result)", m_co, m_numDays, m_numPlayers);
	}

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
