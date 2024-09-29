#include "TripleSys.h"

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
	const int nv1 = getAllV(v1, MAX_3PF_SETS, iRow0, iRow1);
	const auto* pV1 = v1;
	int nce = 0;
	for (int iv1 = 0; iv1 < nv1; iv1++, pV1 += m_nGroups)  // Andrei nv1 is equal to 1
	{
		TrCycles trCycles;

		if (!getCyclesAndPath3(&trCycles, pV1, p1ftable(iRow1), p1ftable(iRow0), result(iRow1), result(iRow0)))
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
	//const int nv1 = getAllV(v1, maxv1, indRow0, indRow1);
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
		icmp = (this->*m_pProcessMatrix)(res, ttr, nrows, n);
	}
}
void alldata::testCanonizatorSpeed()
{
#if 0
	tchar tr[] = { 0,2,1,6,8,7,3,5,4,18,20,19,24,26,25,21,23,22,9,11,10,15,17,16,12,14,13 };
	kmProcessMatrix(result(), tr, iDay);
	printTable("ktmp", m_Ktmp, iDay, 27, 3);
#endif

	// sort matrix from data.h
	kmProcessMatrix(result(), NULL, iDay);
	auto* bRes1 = m_Km;
	int iret = memcmp(bRes1, result(), iDay * m_numPlayers);
	if (iret)
	{
		printfYellow("Sorted matrix different than original\n");
		printTable("Sorted", m_Km, iDay, m_numPlayers, m_groupSize);
	}
	int errLine = 0, errGroup = 0, dubLine = 0;
	if (!CheckMatrix(result(0), iDay, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
	{
		printf("Duplicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
		abort();
	}

	tchar v0[MAX_GROUP_NUMBER * MAX_3PF_SETS];
	int nv0;

	memset(&m_TrCycles, 0, sizeof(m_TrCycles));

	if (m_groupSize == 3)
	{
		nv0 = getAllV(v0, MAX_3PF_SETS, 0, 1);
	}
	clock_t tTime = clock();
	int improveResult = m_improveResult;
	const auto pProcessMatrix = m_pProcessMatrix;
	m_pProcessMatrix = &alldata::kmProcessMatrix;
	m_improveResult = 2;
	int i, nTests = 1000;
	for (i = 1; i <= nTests; i++)
	{
		if (!cnvCheckNew(0, iDay))
		{
			int order = groupOrder();
			printf("Group Order=%d\n", order);
			iret = memcmp(bRes1, result(), iDay * m_numPlayers);
			if (iret >= 0)
			{
				printfYellow("Start matrix rejected by cnvCheckNew\n");
				break;
			}
			printTable("Result ", result(), iDay, m_numPlayers, m_groupSize, 0, true);
			printf("Result improved (%d):\n", i);
			printTable("", bRes1, iDay, m_numPlayers, m_groupSize, 0, true);
#if 1
			memcpy(result(0), bRes1, m_nLenResults);
			errLine = errGroup = dubLine = 0;
			if (!CheckMatrix(result(0), iDay, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
			{
				printf("Duplicate pair in group %d on line %d (already present in line %d)\n", errGroup, errLine, dubLine);
				abort();
			}
			for (int j = 0; j < iDay; j++)
				p1fSetTableRow(p1ftable(j), result(j));
			printTableColor("Links improved", links(0), m_numPlayers, m_numPlayers, m_groupSize);
#endif
		}
		else
		{
			iret = iret;
			int order = groupOrder();
			printf("Group Order=%d\n", order);
			//StatReportAfterEachResult(ResetStat, "Canonizator time", (int)((clock() - tTime)) / i, true);
			//printf("order=%d\n",  groupOrder());
			break;
		}
	}

	printTableColor("Links improved", links(0), m_numPlayers, m_numPlayers, m_groupSize);
	printf("End of Improved (%d):\n", i);
	m_improveResult = improveResult;
	m_pProcessMatrix = pProcessMatrix;
	//time2 = __rdtscp(&junk);
	//GetProcessTimes(hwnd, &dum1, &dum2, &dum3, &time2);
	//QueryPerformanceCounter(&time2);
	//GetThreadTimes(hwnd, &dum1, &dum2, &dum3, &time2);
	if (nTests > 0)
		printf("+++ %.1f ms needed per one improvement check\n", (double(clock() - tTime)) / nTests);
	//printf(" %.1f ms (%.1f cpu) needed per one improvement check\n", (double(clock() - tTime)) / nTests);
	//	(time2.QuadPart - time1.QuadPart) / nTests / 1000.0);
	//	(double(time2.dwLowDateTime - time1.dwLowDateTime)) / nTests / 1000.0);
}
