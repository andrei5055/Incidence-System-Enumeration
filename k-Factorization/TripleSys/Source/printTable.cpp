#include "TopGun.h"

void alldata::printPermutationMatrices(const int iMode) const {
	if (m_groupSize != 3 || iMode < 2)
		return;
	int maxStr = 144;
	auto v = getV0();
	tchar cycles[MAX_3PF_SETS * MAX_GROUP_NUMBER];
	tchar pm0[MAX_GROUP_NUMBER * MAX_GROUP_NUMBER * 2];
	tchar pm[MAX_3PF_SETS * MAX_GROUP_NUMBER * MAX_GROUP_NUMBER];
	for (int i = 1; i < iDay; i++) {
		for (int i0 = 0; i0 < i; i0++) {
			if (i0 == i)
				continue;
			if (iMode == 2 && i0 != 0)
				break;
			//if (iMode == 3 && i0 != 1 && i != 4)
			if (iMode == 3 && i0 != i - 1)
				continue;
			auto nv = getAllV(v, MAX_3PF_SETS, i0, i);
			memset(pm0, 0, sizeof(pm0));
			memset(pm, 0, sizeof(pm));
			memset(cycles, unset, sizeof(cycles));
			memset(&m_TrCyclesAll, 0, sizeof(m_TrCyclesAll));
			for (int j = 0; j < nv; j++) {
				for (int k = 0; k < m_nGroups; k++) {
					int m = neighbors(i)[v[j * m_nGroups + k]];
					int n = result(i)[m];
					pm[j * m_nGroups + m / m_groupSize + k * nv * m_nGroups] = 1;
					pm0[m / m_groupSize + k * m_nGroups * 2] = n;
					pm0[m / m_groupSize + m_nGroups + k * m_nGroups * 2] = 1;
				}
				TrCycles trc;
				const auto ncycles = p3Cycles(&trc, 2, neighbors(i0), neighbors(i), v + j * m_nGroups, result(i0), result(i));
				if (ncycles > 0)
					memcpy(cycles + j * m_nGroups, trc.length, ncycles);
			}
			printf("\nResult rows (%d,%d), cycles and corresponding Permutation matrices\n", i0, i);
			printTableColor("", result(i0), 1, m_numPlayers, m_groupSize);
			printTableColor("", result(i), 1, m_numPlayers, m_groupSize);
			printf("\n");
			printTableColor("", pm0, m_nGroups, m_nGroups * 2, m_nGroups);
			int istr = 0;
			int nvMax = 6;
			for (int inv = 0; inv < nv; inv += nvMax) {
				printf("\n");
				printTableColor("c", cycles + inv * m_nGroups, 1, m_nGroups * min((nv - inv), nvMax), m_nGroups);
				for (int k = 0; k < m_nGroups; k++) {
					printTableColor("", pm + (inv + k * nv) * m_nGroups, 1, m_nGroups * min((nv - inv), nvMax), m_nGroups);
				}
			}
		}
	}
}

void printTableColor(char const* name, ctchar* c, int nl, int nc, int np, int ns, bool makeString, ctchar* co, clock_t* t)
{
	int ind = 0;
	if (name && name[0]) {
		printf("%s:", name);
		if (t)
			printfGreen("     Last column: Row change time stamp (sec)\n");
		else if (nl > 1)
			printf("\n");
	}
	for (int j = 0; j < nl; j++, ind += nc)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			ctchar v = c[ind + i];
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			if (co) {
				if (v == co[ind + i])
					nc > 16 ? printf("%2d", v) : printf("%3d", v);
				else
					nc > 16 ? printfYellow("%2d", v) : printfYellow("%3d", v);
			}
			else if (v < 67)
			{
				if (v == 1)
					printf(" 1");
				else {
					printf("\x1b[38;5;%dm%2d", 28 + v * 3, v);
					printf("\x1b[0m");
				}
			}
			else if (v == unset)
				printf("  ");
			else
				printf("%2d", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
		{
			if (makeString)
				printf(" \"");
			if (t)
				printfGreen("//%7d", t[j] / 1000);
			printf("\n");
		}
		else
			printf(" ");
	}
}

void alldata::printResultWithHistory(char const* name, int nRows)
{
	printTableColor(name, result(), nRows, m_numPlayers, m_groupSize, 0, true, m_pResultsPrev, m_rowTime);
	memcpy(m_pResultsPrev, result(), m_nLenResults);
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
	int time[4];
	const int multTime[] = { 60, 60, 24 };
	auto *pUnitTime = "smhd";
	time[0] = (int)((clock() - iTime) / 1000. + 0.5);
	for (int i = 0; i < countof(multTime); i++) {
		time[i + 1] = time[i] / multTime[i];
		time[i] %= multTime[i];
	}
	
	char timeBuf[256], * pTime = timeBuf;
	for (int i = countof(time); i--;)
		if (time[i] || pTime != timeBuf || !i)
			SPRINTFD(pTime, timeBuf, "%2d%c", time[i], pUnitTime[i]);

	m_reportInfo = std::format(
		"T ={}: {} (from {}) {}-matrices ({}x{}) processed by {} threads. {} {}-matrices ({}x{}) generated\n",
		timeBuf, nProccesed, nMatrices, fhdr, numPlayers(), nRowsStart(), numThreads, sum1, fhdr, numPlayers(), nRowsOut());

	if (bPrintSetup)
	{
		char buffer[256], *pBuf = buffer;
		const auto lenBuf = sizeof(buffer);
		for (int i = 0; i < numThreads; i++)
		{
			if ((i % 20) == 0)
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
	outputIntegratedResults(NULL, 3);
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
