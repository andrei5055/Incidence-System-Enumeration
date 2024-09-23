#include <iostream>
#include "TripleSys.h"

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
int alldata::reportMatrixStats(bool bPrint) const
{
	int nbrv = 0;
    int nbrt = 0;
#if UseTrMask == 2
	int nbits[32];
	assert (m_numDays <= sizeof(nbits) / sizeof(nbits[0]));

	memset(nbits, 0, m_numDays * sizeof(nbits[0]));
	int nbrMin = m_numDays;
	int nbrMax = 0;

	SetBit(m_TrMask, 0);
	for (int i = 0; i < m_nTr; i++)
	{
		int nbr = 0;
		for (int j = 0; j < m_numDays; j++)
		{
			if (GetBit(m_TrMask + j * m_nTrBytes, i))
			{
				nbits[j]++;
				nbr++;
			}
		}
		if (nbr > 0 && bPrint)
		{
			printf("%10d: ", i);
			for (int j = 0; j < m_numDays; j++)
			{
				if (GetBit(m_TrMask + j * m_nTrBytes, i))
					printf("1");
				else
					printf(".");
			}
			printf("\n");
			nbrMin = std::min(nbr, nbrMin);
			nbrMax = std::max(nbr, nbrMax);
			nbrt++;
		}
	}
	nbrv = nbits[0];
	for (int j = 1; j < m_numDays; j++)
	{
		if (nbits[j] > 0 && nbrv != nbits[j])
		{
			assert(nbrv == nbits[j]);
			//nbrv = std::max(nbrv, nbits[j]);
		}
	}
	if (bPrint)
	{
		printfYellow2("%5d: |Aut(M)| = %d", m_finalKMindex, nbrt);
		printfYellow2(", full matrix size=%dx%d", m_nTr, m_numDays);
		printfYellow(", non empty dimension=%dx", nbrv);
		if (nbrMin == nbrMax)
			printfYellow("%d\n", nbrMin);
		else
			printfYellow2("(%d-%d)\n", nbrMin, nbrMax);
#if 0
		printf("      column/count: ");
		for (int j = 0; j < m_numDays; j++)
		{
			printf("%d/%d ", j, nbits[j]);
		}
		printf("\n");
#endif
	}
#endif
	return nbrt;
}
void alldata::TestkmProcessMatrix(int nrows, unsigned char n, char* tr, char* ttr, int icmp, int* pDayMax) const
{
//  test for kmProcessMatrix
	char* res = result();
	int dayMax = nrows - 1;
	int icmp2 = kmProcessMatrix(m_Km, res, m_Ktmp, nrows, m_numPlayers, m_groupSize, ttr, &dayMax);
	if (icmp != icmp2 || (pDayMax != NULL && *pDayMax != dayMax))
	{
		printTransformed(nrows, m_numPlayers, tr, ttr, res, m_Km, n, nLoops, m_finalKMindex);
		if (icmp != icmp2 || (pDayMax != NULL && *pDayMax != dayMax))
			printf("ic=%d must be %d\n", icmp, icmp2);
		if (pDayMax != NULL && *pDayMax != dayMax)
			printf("dayMax %d must be %d\n", *pDayMax, dayMax);
		//exit(1);
	}
}
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
	const auto createImprovedMatrix = m_createImprovedMatrix;
	m_createImprovedMatrix = true;
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
			m_createImprovedMatrix = createImprovedMatrix;
			m_improveResult = improveResult;
		}
	}
	m_createImprovedMatrix = createImprovedMatrix;
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
sLongLong __a[16] = { 0 };
const char* __t[16] = { 0 };
void Stat(const char* t, int ind, bool bAdd)
{
	if (ind < 0 || ind >= 16)
		abort();
	__t[ind] = t;
	if (bAdd)__a[ind]++;
}
void TestStatPrint(const char* hdr, int d)
{
	bool np = true;
	printf(hdr, d);
	for (int i = 0; i < 16; i++)
	{
		if (__t[i])
		{
			np = false;
			printf(" %s:%zd", __t[i], __a[i]);
		}
	}
	if (np)
		printf(".");
	else
		printf("\n");
	memset(__a, 0, sizeof(__a)); 
	memset(__t, 0, sizeof(__t));
}
