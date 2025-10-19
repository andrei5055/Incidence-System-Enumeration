#include "TripleSys.h"  
void aq()
{
#define N 28
	tchar v[N];
	tchar vk[N];
	tchar tr[N];
	int trib[N];
	int cnt[N*2];
	trib[0] = 0;
	trib[1] = 1;
	trib[2] = 1;
	double nsOld = 1;
	for (int i = 3; i < N; i++)
		trib[i] = trib[i - 1] + trib[i - 2] + trib[i - 3];
	for (int i = 0;i < N;i++)
		v[i] = i + 1;
	for (int nv = 3; nv <= 8; nv++) {
		int nvm2 = nv - 1;
		int n2 = 1 << nvm2;
		int ns = 0;
		int nlast = 0;
		memset(cnt, 0, sizeof(cnt));
		for (int j = 0; j < n2; j++) {
			int km = nv;
			int k1 = 0;
			for (int n = 0; n < nv; n++) {
				if (j & (1 << n))
					vk[--km] = n + 1;
				else
					vk[k1++] = n + 1;
			}
			k1 = 0;

			for (int n = 0; n < nv; n++)
				tr[vk[n] - 1] = n + 1;
			
			//printTableColor("tr", tr, 1, nv, 1);
			//printTableColor("vks", vk, 1, nv, 1);
			for (int n = 0; n < nv; n++) {
				if (tr[k1] == 0)
					goto err1;
				int k0 = k1;
				k1 = tr[k1] - 1;
				tr[k0] = 0;
			}
			ns++;
			printf(" %d", j);
			if (vk[1] == nv)
				nlast++;
			//cnt[vk[0] - 1]++;
			for (int ic = 1; ic < nv; ic++)
				cnt[(vk[ic] - vk[ic-1] + nv)]++;
			//printTableColor("vkr", vk, 1, nv, 1);
		err1: continue;
		}
		int is = 0;
		trib[0] = 1;
		for (int ic = nv - 3; ic >= 0; ic--)
			is += trib[ic];

		//printf("N=%2d S=%6d %6d TF=%6d nl=%d", nv, ns, trib[nv - 2], is, nlast);
		//printf("N=%2d S=%6d R=%.3f log2=%.3f", nv, ns, ns / nsOld, log2(ns) - nv + 6);
		printf("\nN=%2d S=%6d R=%.3f log2=%.3f", nv, ns, ns / nsOld, log2(ns) - nv + 6);
		nsOld = ns;
		/**
		for (int ic = 0; ic < sizeof(cnt) / sizeof(cnt[0]); ic++)
			if (cnt[ic])
				printf(" %d:%d", ic - nv, cnt[ic]);**/
		printf("\n");
		//printf("N=%2d S=%6d approximation=%.0f\n", nv, ns, pow(1.74, nv-2));
	}
	exit(0);
}
bool addRow(tchar* allPaths, int maxP, int* pnp, tchar* path, int pathLength)
{
	int np = *pnp;
	if (np < maxP)
	{
		int ip = 0;
		for (ip = 0; ip < np; ip++)
		{
			switch (memcmp(allPaths + ip * pathLength, path, pathLength)) {
			case 0: return false;
			case -1: continue;
			case 1:
				for (int jp = np; jp > ip; jp--)
				{
					memcpy(allPaths + jp * pathLength, allPaths + (jp - 1) * pathLength, pathLength);
				}
				break;
			}
			break;
		}
		np++;
		memcpy(allPaths + ip * pathLength, path, pathLength);
	}
	*pnp = np;
	return true;
}
void formOneSp(tchar* r, tchar* rp, int dir, int nc, int gs)
{
	tchar sp[MAX_PLAYER_NUMBER];
	//printf("   ");
	for (int j = 0; j < nc; j++)
		sp[r[j]] = j;
	for (int j = 0; j < nc; j++)
	{
		tchar pl = sp[j];
		tchar rn;
		if (dir & (1 << (pl / gs)))
			rn = (pl % gs) == 0 ? r[pl + gs - 1] : r[pl - 1];
		else
			rn = (pl % gs) == gs - 1 ? r[pl - gs + 1] : r[pl + 1];
		rp[j] = rn;
	}
}
void alldata::testRightNeighbor(int nr)
{
	int iRow0 = 0, iRow1 = 1;
	tchar mLS[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	memset(mLS, 0, sizeof(mLS));
	tchar mLSRow0[MAX_PLAYER_NUMBER];
	tchar mLSRow1[MAX_PLAYER_NUMBER];
	for (int i = 0; i < nr; i++)
	{
		formOneSp(result(i), mLSRow0, 0, m_numPlayers, m_groupSize);
		memcpy(mLS + mLSRow0[0] * m_numPlayers, mLSRow0, m_numPlayers);
		formOneSp(result(i), mLSRow0, -1, m_numPlayers, m_groupSize);
		memcpy(mLS + mLSRow0[0] * m_numPlayers, mLSRow0, m_numPlayers);
	}
	printTable("\nmln+mrn", mLS + m_numPlayers, m_numPlayers - 1, m_numPlayers, m_groupSize);
#define MaxPaths 300
	tchar path[MAX_PLAYER_NUMBER * 2];
	tchar allPaths[MAX_PLAYER_NUMBER * 2 * MaxPaths];
	memset(allPaths, 0, sizeof(allPaths));
	int pathLength = m_numPlayers * 2;
	int np = 0;
	int nGroupsDirectins = 1 << m_nGroups;
	for (int idir0 = 0; idir0 < nGroupsDirectins; idir0++) // idir0 - bits to define "groups directions", 1 bit for each group
	{
		formOneSp(result(iRow0), mLSRow0, idir0, m_numPlayers, m_groupSize);
		for (int ifirst = 0; ifirst < m_groupSize; ifirst++)
		{
			for (int idir = 0; idir < nGroupsDirectins; idir++) // groups direction bitmask for second row
			{
				if (idir0 == 0x18 && idir == 0x1c && ifirst == 2)
					ifirst = ifirst;
				tchar middlePlayersCounter[MAX_PLAYER_NUMBER];
				formOneSp(result(iRow1), mLSRow1, idir, m_numPlayers, m_groupSize);
				int iPathIndex = 0, ind = result(iRow1)[ifirst];
				tchar* p0 = mLSRow0;
				tchar* p1 = mLSRow1;
				memset(path, unset, sizeof(path));
				tchar tst[MAX_PLAYER_NUMBER];
				memset(tst, unset, sizeof(tst));
				int ierr = 0;
				memset(middlePlayersCounter, 0, sizeof(middlePlayersCounter));
				int middlePlayersPairsCount = 0;
				for (int i = 0; i < m_numPlayers / m_groupSize; i++)
				{
					path[iPathIndex++] = ind;
					if (m_groupSize == 3)
					{
						path[iPathIndex++] = ind = p1[ind];
						middlePlayersCounter[ind]++;
						if (middlePlayersCounter[ind] == 2)
							middlePlayersPairsCount++;
						path[iPathIndex++] = ind = p1[ind];
					}
					path[iPathIndex++] = ind;
					if (m_groupSize == 3)
					{
						path[iPathIndex++] = ind = p0[ind];
						middlePlayersCounter[ind]++;
						if (middlePlayersCounter[ind] == 2)
							middlePlayersPairsCount++;
					}
					path[iPathIndex++] = ind = p0[ind];
					if (ind == path[0])
						break;
				}
				if ((iPathIndex % (m_groupSize * 2)) != 0 || 
					path[iPathIndex - 1] != path[0] || middlePlayersPairsCount != iPathIndex / (m_groupSize * 2))
					continue;
				//printTable("mLSRow0", mLSRow0, 1, m_numPlayers, 0);
				//printTable("mLSRow1", mLSRow1, 1, m_numPlayers, 0);
				if (addRow(allPaths, MaxPaths, &np, path, pathLength))
				{
					printf("%03d dir0=0x%02x dir=0x%02x\npath:", np, idir0, idir);
					printTable("", path, 1, iPathIndex, 0);
				}
				else
				{
					printf("%03d dir0=0x%02x dir=0x%02x\npath:", np, idir0, idir);
					printf("dublicate path ignored\n");
					//printTable("", path, 1, iPathIndex, 0);
				}
			}
		}
	}
	printTable("N&D Paths", allPaths, np, pathLength, m_groupSize);

	bool added = false;
	auto v1 = getV1();
	const int nv1 = getAllV(v1, m_maxCommonVSets, iRow0, iRow1);
	const auto* pV1 = v1;
	int nce = 0;
	for (int iv1 = 0; iv1 < nv1; iv1++, pV1 += m_nGroups)  // Andrei nv1 is equal to 1
	{
		TrCycles trCycles;

		if (!getCyclesAndPath3(&trCycles, pV1, neighbors(iRow1), neighbors(iRow0), result(iRow1), result(iRow0), eNoErrorCheck))
			continue;
		int m = 0;
		for (m = 0; m < pathLength; m += m_groupSize * 2)
		{
			if (trCycles.fullPath[m] == 0 || trCycles.fullPath[m + 1] == 0 || trCycles.fullPath[m + m_groupSize - 1] == 0)
			{
				memcpy(path, &trCycles.fullPath[m], pathLength - m);
				memcpy(path + pathLength - m, trCycles.fullPath, m);
			}
		}
		nce++;
		if (addRow(allPaths, MaxPaths, &np, path, pathLength))
		{
			printf("new path from common elemets\n");
			printTable("", path, 1, pathLength, 0);
			added = true;
		}
	}
	if (!added)
		printf("All 'common elements' based paths(%d) already present in %d n&d paths\n", nce, np);
	if (m_groupSize <= 3)
	{
		tchar sp[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
		memset(sp, 0, sizeof(sp));
		//printf("   ");
		for (int j = 0; j < m_numPlayers * nr; j += m_groupSize)
		{
			int a = result()[j], b = result()[j + 1];
			if (m_groupSize == 2)
			{
				tchar c = j / m_numPlayers;
				sp[a * m_numPlayers + b] = c;
				sp[b * m_numPlayers + a] = c;
			}
			else 
			{
				tchar c = result()[j + 2];
				sp[a * m_numPlayers + b] = c;
				sp[b * m_numPlayers + a] = c;
				sp[c * m_numPlayers + b] = a;
				sp[b * m_numPlayers + c] = a;
				sp[a * m_numPlayers + c] = b;
				sp[c * m_numPlayers + a] = b;
			}
		}
		for (int j = 0; j < m_numPlayers; j++)
		{
			sp[j * m_numPlayers + j] = j;
		}
		printTable("\nclassic", sp, m_numPlayers, m_numPlayers, m_groupSize);
	}
	cnvCheckNew(0, nr, false);
	exit(0);
}
void alldata::TestkmProcessMatrix(int nrows, unsigned char n, const tchar* tr, const tchar* ttr, int icmp) const
{
//  test for kmProcessMatrix
	auto* res = result();
	int playerIndex  = m_playerIndex;
	const int icmp2 = kmProcessMatrix(res, ttr, nrows);
	const auto flag = playerIndex != m_playerIndex;
	if (icmp != icmp2 || (icmp == -1 && flag))
	{
		printTransformed(nrows, m_numPlayers, m_groupSize, tr, ttr, res, m_Km, n, nLoops, m_finalKMindex);
		if (icmp != icmp2)
			printfRed(" TestkmProcessMatrix: ic=%d must be %d\n", icmp, icmp2);
		if (flag)
			printfRed("TestkmProcessMatrix: m_playerIndex %d must be %d\n", playerIndex, m_playerIndex);
		myExit(1);
		icmp = (this->*m_pProcessMatrix)(res, ttr, nrows, n, NULL);
	}
}
void alldata::testCanonizatorSpeed()
{
#if 0
	tchar tr[] = { 0,2,1,6,8,7,3,5,4,18,20,19,24,26,25,21,23,22,9,11,10,15,17,16,12,14,13 };
	kmProcessMatrix(result(), tr, iDay);
	printTable("ktmp", m_Ktmp, iDay, 27, 3);
#endif
	m_precalcMode = eCalcResult; // 
	// sort matrix from data.h
	kmProcessMatrix(result(), NULL, iDay);
	auto* pRes1 = m_Km;
	int iret = memcmp(pRes1, result(), iDay * m_numPlayers);
	if (iret)
	{
		printfYellow("Sorted matrix different than original\n");
		memcpy(result(), pRes1, iDay * m_numPlayers);
		printTable("Sorted", result(), iDay, m_numPlayers, m_groupSize);
	}
	int errLine = 0, errGroup = 0, dubLine = 0;
	if (!CheckMatrix(result(0), iDay, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
	{
		printf("Duplicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
		abort();
	}
	char stat[256];
	clock_t tTime = clock();
	const auto pProcessMatrix = m_pProcessMatrix;
	m_pProcessMatrix = &alldata::kmProcessMatrix;
	int i, nTests = 100;
	int order = 0;
	bool bRet = false;
	m_lastRowWithTestedTrs = 0;
	int iTest = 0;
	for (int j = 0; j < iDay; j++)
		u1fSetTableRow(neighbors(j), result(j));
	for (i = 1; i <= nTests; i++)
	{
		iTest++;
		bRet = cnvCheckNew(0, iDay);
		if (!bRet)
		{
			order = orderOfGroup();
			printf("Group Order=%d\n", order);
			iret = memcmp(pRes1, result(), iDay * m_numPlayers);
			if (iret >= 0)
			{
				printfYellow("Start matrix rejected by cnvCheckNew\n");
				break;
			}
			printTable("Result ", result(), iDay, m_numPlayers, m_groupSize, 0, true);
			printf("Result improved (%d):\n", i);
			printTable("", pRes1, iDay, m_numPlayers, m_groupSize, 0, true);
#if 1
			memcpy(result(0), pRes1, m_nLenResults);
			m_lastRowWithTestedTrs = 0;
			errLine = errGroup = dubLine = 0;
			if (!CheckMatrix(result(0), iDay, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
			{
				printf("Duplicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
				abort();
			}
			for (int j = 0; j < iDay; j++)
				u1fSetTableRow(neighbors(j), result(j));
			//printTableColor("Links improved", links(0), m_numPlayers, m_numPlayers, 0);
#endif
		}
		else
		{
			order = orderOfGroup();
			printf("Group Order=%d\n", order);
			//StatReportAfterEachResult(ResetStat, "Canonization time", (int)((clock() - tTime)) / i, true);
			//printf("order=%d\n",  order);
			break;
		}
	}

	//printTableColor("Links improved", links(0), m_numPlayers, m_numPlayers, 0);
	printf("End of Improvement (%d):\n", i);
	//time2 = __rdtscp(&junk);
	//GetProcessTimes(hwnd, &dum1, &dum2, &dum3, &time2);
	//QueryPerformanceCounter(&time2);
	//GetThreadTimes(hwnd, &dum1, &dum2, &dum3, &time2);
	if (nTests > 0)
		printf("+++ %.1f ms needed per one improvement check\n", (double(clock() - tTime)) / iTest);
	//printf(" %.1f ms (%.1f cpu) needed per one improvement check\n", (double(clock() - tTime)) / nTests);
	//	(time2.QuadPart - time1.QuadPart) / nTests / 1000.0);
	//	(double(time2.dwLowDateTime - time1.dwLowDateTime)) / nTests / 1000.0);

	for (int j = 0; j < iDay; j++)
		u1fSetTableRow(neighbors(j), result(j));

	bool needOutput = false;
	matrixStat(neighbors(), iDay, &needOutput);
	if (needOutput) {
		matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll);
		printf("%d rows: %s, AUT=%d, %s\n", iDay,
			bRet ? "Canonical" : "Not canonical", orderOfGroup(), stat);
	}
}

#if !USE_CUDA 
bool alldata::testGroupOrderEachSubmatrix(int iPrintMatrices, eThreadStartMode iCalcMode)
{
	if (!(iPrintMatrices & 128) || iCalcMode == eCalcSecondRow)
		return false;
	else {
		printf("Submatrices Automorphism and Cycles:\n");
		printTable("The Matrix", result(), iDay, m_numPlayers, m_groupSize);
		m_numDaysResult = iDay;
		for (int i = 2; i <= iDay; i++) {
			char stat[256];
			bool bRet = cnvCheckNew(0, i, false);
			bool bRet2 = matrixStat(neighbors(), i, NULL);
			bool needOutput = false;
			matrixStat(neighbors(), i, &needOutput);
			if (needOutput) {
				matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll);
				printf("%d rows: %s, AUT=%d, %s\n", i,
					bRet ? "Canonical" : "Not canonical", orderOfGroup(), stat);
			}
		}
		printf("\n");
		exit(1);
	}
	return true;
}

/*
void alldata::testPrintGroupRows()
{
	tchar trm[MAX_PLAYER_NUMBER];
	for (int i = 0; i < m_numDays; i++) {

	}
		
	for (int i = 0; i < numObjects(); i++) {
		tchar* tr;


	}
}
*/
#endif
