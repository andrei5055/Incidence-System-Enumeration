#include <iostream>
#include "TripleSys.h"

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
void alldata::testImproveMatrixSpeed()
{
	// calculate average time for one call to cnvCheckNew or improveMatrix
	// 
	//HANDLE hwnd;
	//hwnd = GetCurrentProcess(); // current process handle
	//FILETIME time1, time2, dum1, dum2, dum3;
	//LARGE_INTEGER time1, time2, dum1, dum2, dum3;
	//uint64_t time1, time2;
    //unsigned int junk = 0;

	unsigned char* bResults = NULL;
	const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
	char* bRes1 = NULL;
	const auto bResults_1 = new unsigned char[2 * lenResult];

	// sort matrix from data.h
	kmProcessMatrix(m_Km, result(), m_Ktmp, iDay, m_numPlayers, m_groupSize);
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
	int nTests = 1;
	for (int i = 0; i < nTests; i++)
	{
#if 1   // Change to 0 to use "improveMatrix"
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
