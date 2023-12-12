#include <iostream>
#include "TripleSys.h"
#ifdef CD_TOOLS
   #include "../CanonicityChecker.h"
#else
   #include "CanonicityChecker.h"
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
	m_pCheckCanon = new CCanonicityChecker<unsigned char, unsigned char>(m_numDays, numPlayers, groupSize);
#endif

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
	delete[] m_h;
	delete[] m_ho;
	delete m_pCheckLink;
	delete m_pCheckCanon;
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

void _printf(FILE* f, const char* format, const char* pStr = NULL) {
	if (f)
		fprintf(f, format, pStr);

	printf(format, pStr);
}

void alldata::outputResults(int iDay, int cntr, const unsigned char *pPlayersDayPerm) const
{
	static unsigned int cntr_init;
	char buffer[256];
	const auto lenBuf = sizeof(buffer);
	
	FILE* f = NULL;
	if (strlen(ImprovedResultFile))
		fopen_s(&f, ImprovedResultFile, cntr_init? "a" : "w");

	if (!cntr) {
		sprintf_s(buffer, "Initial Result #%d:\n", ++cntr_init);
		_printf(f, buffer);
		for (int j = 0; j <= iDay; j++) {
			const auto* pRes = result() + j * m_numPlayers;
			char* pBuf = buffer;
			pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), " \"");
			for (int i = 0; i < m_numPlayers; i++) {
				if (!(i % m_groupSize))
					pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), "  %3d", *pRes++);
				else
					pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), "%3d", *pRes++);
			}

			_printf(f, "%s \\n\"\n", buffer);
		}
	}
	else {
		sprintf_s(buffer, "Improved Result #%d:\n", cntr);
		_printf(f, buffer);
		const auto* pDayPerm = pPlayersDayPerm + m_numPlayers;
		bool flag = true;
		for (int j = 0; j <= iDay; j++) {
			const int day = pDayPerm[j];
			const auto* pPlayers = pPlayersDayPerm;
			const auto *pRes = result() + day * m_numPlayers;
			char* pBuf = buffer;
			pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), " \"");
			for (int i = 0; i < m_numPlayers; i++) {
				flag &= j != 0 || pRes[*pPlayers] == i;
				if (!(i % m_groupSize))
					pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), "  %3d", pRes[*pPlayers++]);
				else
					pBuf += sprintf_s(pBuf, lenBuf - (pBuf - buffer), "%3d", pRes[*pPlayers++]);
			}

			_printf(f, "%s \\n\"\n", buffer);
		}

		assert(flag);
	}

	if (f)
		fclose(f);
}

bool alldata::Run(int improveResult) {
	// Input parameter:
	//      improveResult: 
	//           0 (default) - do not try to improve given "result"; 
	//		   !=0 - m_pCheckCanon will return the permutations of the sets of players 
	//               and days that improve given "results";
	//          >1 - try to improve the “results” as much as possible.
	clock_t iTime = clock();
	unsigned char* bResults = NULL;
	const auto lenResult = m_numDays * (m_numPlayers + m_numDays);
	if (improveResult)
		bResults = new unsigned char[(improveResult > 1? 2 : 1) * lenResult];

	while (nLoops < LoopsMax)
	{
		while (iDay < m_numDays || bPrevResult)
		{
			if (iDay < 0)
			{
				noMoreResults = true;
				break;
			}
			if (bPrevResult)
			{
				/**/
				if (nLoops == 1)
				{
					//printf("%d ", iDay);

					//printTable("", result(), m_numDays, m_numPlayers);
					//printTable("Links", links(), m_numPlayers, m_numPlayers);
				}/**/

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

			if (maxDays < iDay)
			{
				/**/
				printf("day %d  Time = %d\n", iDay, clock() - iTime);
				printTable("Result", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
				//printTable("Links", links[0], m_numPlayers, m_numPlayers);
				/**/
				maxDays = iDay;
				memcpy(maxResult, result(0), m_nLenResults);
			}
			if (iDay > 0)
			{
#if 0
				for (int j = 2; j <= iDay; j++) {
					const auto* pntr = result() + m_numPlayers * (j+1) - 1;
					if (*pntr == 14 && (*(pntr - 1) == 11 && *(pntr - 2) == 8))
						nLoops += 0;
				}
/*
NOt a canonical one:
Initial Result:

	0  1  2    3  4  5    6  7  8    9 10 11   12 13 14
	0  3  6    1  7 12    2  8  9    4 10 14    5 11 13
	0  4  7    1  3  8    2 10 13    5  9 14    6 11 12
	0  5  8    1  4 11    2  7 14    3 10 12    6  9 13
	0  9 12    1  5 10    2  4  6    3  7 13    8 11 14
*/
				/*
				Initial Result :
				"    0  1  2    3  4  5    6  7  8    9 10 11   12 13 14 \n"
				"    0  3  6    1  4 12    2  7  9    5 10 14    8 11 13 \n"
				Improved Result #1:
				"    0  1  2    3  4  5    6  7  8    9 10 11   12 14 13 \n"
				"    0  3  6    1  4  9    2  7 12    8 10 13    5 11 14 \n"
				*/
#endif
				static int fff = 0; fff++;
				if (!m_pCheckCanon->CheckCanonicity((unsigned char *)result(), iDay+1, bResults))
				{
					if (PrintImprovedResults || improveResult > 1) {
						int cntr = 0;
						auto* bRes1 = bResults;
						auto* bRes2 = bResults + lenResult;
						do {
							if (PrintImprovedResults) {
								if (!cntr) {
									// Output of initial results
									outputResults(iDay);
								}

								outputResults(iDay, ++cntr, bRes1);
							}

							if (true || improveResult == 1) // Not ready yet
								break;

							// Swap the the best results buffers
							auto* bRes = bRes1;
							bRes1 = bRes2;
							bRes2 = bRes;
						} while (!m_pCheckCanon->CheckCanonicity(bRes2, iDay+1, bRes1));
					}

					// get new matrix
					bPrevResult = true;
				}
			}
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
		printf("Result %d, Time = %d\n", nLoops, clock() - iTime);
		//printTable("Links", links(), m_numPlayers, m_numPlayers);
		printTable("Result table", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
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
