#include "TopGun.h"

void printTableColor(char const* name, const tchar *c, int nl, int nc, int np, int ns, bool makeString)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			char v = c[j * nc + i];
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			if (v >= 0 && v < 67)
			{
				printf("\x1b[38;5;%dm%2d", 28 + v * 3, v);
				printf("\x1b[0m");
			}
			else
				printf("%2d", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}

sLongLong TopGun::printThreadsStat(int nMatrices, int nProccesed, const clock_t& iTime, bool bPrintSetup)
{
	const sLongLong* cntTotal = m_cntTotal;
	const sLongLong* cnt = m_cnt;
	sLongLong sum1 = 0, sum2 = 0;
	for (int i = 0; i < numThreads; i++)
	{
		const auto j = i * 2;
		sLongLong d = cnt[j];
		if (d < 0)
			d = -1 - d;
		sum1 += cntTotal[j] + d;
		sum2 += cntTotal[j + 1] + cnt[j + 1];
	}

	const char* fhdr = getFileNameAttr(paramPtr());
	m_reportInfo = std::format(
		"T = {:.0f}sec: {} (from {}) {}-matrices ({}x{}) processed by {} threads. {} {}-matrices ({}x{}) generated\n",
		(clock() - iTime) / 1000., nProccesed, nMatrices, fhdr, numPlayers(), nRowsStart(), numThreads, sum1, fhdr, numPlayers(), nRowsOut());

	if (bPrintSetup)
	{
		char buffer[256], *pBuf = buffer;
		const auto lenBuf = countof(buffer);
		for (int i = 0; i < numThreads; i++)
		{
			if ((i % 12) == 0)
				SPRINTFS(pBuf, buffer, lenBuf, "\n");

			const auto j = i * 2;
			sLongLong d = cnt[j];
			if (d < 0)
				d = -1 - d;
			SPRINTFS(pBuf, buffer, lenBuf, " %d:%zd", i + 1, cntTotal[j] + d);
		}
		m_reportInfo += std::format("Thread:Matrices generated{}\n", buffer);
	}
	std::cout << m_reportInfo;
	return sum1;
}
void printTransformed(int nrows, int ncols, int groupSize, const tchar* tr, const tchar* ttr, const tchar* pImatr, const tchar* pTmatr, int numRow, sLongLong nLoops, int finalKMindex)
{
	if (nLoops)
		printf("Calculated Matrix %zd can be improved. See below (nKm=%d row=%d)\n", nLoops, finalKMindex, numRow);
	
	printTable("Tr source", tr, 1, ncols, groupSize);
	printTable("Tr actual", ttr, 1, ncols, groupSize);
	printTable("Original", pImatr, nrows, ncols, groupSize);
	printTable("Translated", pTmatr, nrows, ncols, groupSize);
}
