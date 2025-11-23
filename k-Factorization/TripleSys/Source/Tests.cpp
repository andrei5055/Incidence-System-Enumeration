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
bool alldata::canonizeMatrix(int nRows)
{
	auto precalcMode = m_precalcMode;
	m_precalcMode = eCalcResult;
	//printTable("Input matrix", result(), nRows, m_numPlayers, m_groupSize, 0, true);

	tchar tm[MAX_PLAYER_NUMBER];
	memset(tm, 0, sizeof(tm));
	memcpy(m_Km, result(), nRows * m_numPlayers);
	(this->*m_pSortGroups)(m_Km, nRows);
	auto* coi = m_Ktmp;
	auto* cii = m_Km;
	for (int i = 0; i < nRows; i++, coi += m_numPlayers, cii += m_numPlayers)
		kmSortGroupsByFirstValue(cii, coi);
	// Result of the loop above is in m_Ktmp, sort and send it to m_Km.
	kmSortRowsBy2ndValue(nRows, tm);
	memcpy(result(), m_Km, nRows * m_numPlayers);
	// sort matrix 
	kmProcessMatrix(result(), NULL, nRows);
	auto* pRes1 = m_Km;
	int iret = memcmp(pRes1, result(), nRows * m_numPlayers);
	if (iret)
	{
		if (m_printMatrices & 64)
			printfYellow("Sorted matrix different than original\n");
		memcpy(result(), pRes1, nRows * m_numPlayers);
		if (m_printMatrices & 64)
			printTable("Sorted", result(), nRows, m_numPlayers, m_groupSize);
	}
	char stat[1024];
	const auto pProcessMatrix = m_pProcessMatrix;
	m_pProcessMatrix = &alldata::kmProcessMatrix;
	int order = 0;
	bool bRet = false;
	m_lastRowWithTestedTrs = 0;
	int i = 0;
	while (1)
	{
		i++;

		m_lastRowWithTestedTrs = 0;

		for (int j = 0; j < nRows; j++) {
			setArraysForLastRow(nRows);
			u1fSetTableRow(neighbors(j), result(j));
		}
		if (!linksFromMatrix(links(), result(), nRows, false)) {

			if (m_printMatrices & 64)
				printTable("Input matrix", result(), nRows, m_numPlayers, m_groupSize, 0, true);
			return false;
		}

		int errLine = 0, errGroup = 0, dubLine = 0;
		if (!CheckMatrix(result(), nRows, m_numPlayers, m_groupSize, true, &errLine, &errGroup, &dubLine))
		{
			if (m_printMatrices & 64)
				printTable("Input matrix", result(), nRows, m_numPlayers, m_groupSize, 0, true);
			printfYellow("Duplicate pair in group=%d row=%d (already present in row=%d)\n", errGroup, errLine, dubLine);
			return false;
		}
		m_lastRowWithTestedTrs = 0;
		bRet = cnvCheckNew(0, nRows);
		if (bRet)
			break;

		if (i == 1) {
			if (m_printMatrices & 64)
				printTable("Initial matrix", result(), nRows, m_numPlayers, m_groupSize, 0, true);
		}
		else {
			if (m_printMatrices & 64) {
				printf("Result of canonization (iteration %d):\n", i);
				printTable("", result(), nRows, m_numPlayers, m_groupSize, 0, true);
			}
		}
		iret = memcmp(pRes1, result(), nRows * m_numPlayers);
		if (iret >= 0)
		{
			printfRed("Start matrix rejected. Check for conflict(s) with your params (AllowUndefinedCycles, Any2RowsConvertToFirst2)\n");
			return false;
		}
		memcpy(result(0), pRes1, m_nLenResults);
	}

	memcpy(m_pResultsPrev, result(), nRows * m_numPlayers);
	memcpy(m_pResultsPrev2, result(), nRows * m_numPlayers);

	bool needOutput = false;
	getAllCycles(neighbors(), nRows);
	matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll);
	if (m_printMatrices & 64) {
		printf("%d rows: AUT=%d, %s\n", nRows, orderOfGroup(), stat);
		if (i == 1)
			printTable("Input matrix (canonical)", result(), nRows, m_numPlayers, m_groupSize, 0, true);
		else {
			printf("Input matrix (after %d canonization iterations):\n", i);
			printTable("", result(), nRows, m_numPlayers, m_groupSize, 0, true);
		}
	}
	m_precalcMode = precalcMode;
	return true;
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
			//bool bRet2 = checkNewRow(neighbors(), i);
			getAllCycles(neighbors(), i);
			matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll);
			printf("%d rows: %s, AUT=%d, %s\n", i,
				bRet ? "Canonical" : "Not canonical", orderOfGroup(), stat);
		}
		printf("\n");
		exit(1);
	}
	return true;
}

bool alldata::generateMatrixExample() {
	if (m_groupSize == 3 && !(m_nGroups & 1)) {
		printfYellow("*** Warning: can't generate matrix (%dx%dx%d): with GroupSize=%d 'number of groups'(%d) must be odd number\n",
			m_numPlayers, m_numDays, m_groupSize, m_groupSize, m_nGroups);
		return false;
	}

	int is = param(t_generateMatrixExample);
	if (is > 1 && m_groupSize <= 2) {
		auto* ls = links(0); // use links as a temporary storage for latin square, we wll recreate links later in this function
		getLS(ls, m_nGroups, is);
		for (int i = 0; i < m_numDays; i++) {
			auto* r = result(i);
			auto* lsRow = ls + i * m_nGroups;
			for (int j = 0, k = 0; j < m_numPlayers; j += m_groupSize, k++, r += m_groupSize) {
				r[0] = k * 2;
				r[1] = lsRow[k] * 2 + 1;
			}
		}
		//printTableColor("r", result(), m_numDays, m_numPlayers, m_groupSize);
		tchar tr[MAX_PLAYER_NUMBER];
		memset(tr, 0, m_numPlayers);
		for (int i = 0;i < m_numPlayers; i++)
			tr[result(0)[i]] = i;
		kmTranslate(result(0), result(0), tr, m_numDays * m_numPlayers);
		//printTableColor("r", result(), m_numDays, m_numPlayers, m_groupSize);
	}
	else {

		for (int i = 0; i < m_numDays; i++) {
			auto* r = result(i);
			for (int j = 0; j < m_numPlayers; j += m_groupSize, r += m_groupSize) {
				for (int k = 0; k < m_groupSize; k++) {
					int ip = ((j / m_groupSize + i * k) % m_numDays);
					ip = ip * m_groupSize + k;
					r[k] = ip;
				}
			}
		}
	}

	if (!canonizeMatrix(m_numDays))
		return false;
	if (!linksFromMatrix(links(), result(), iDay, false)) {
		printfRed("*** Error: can't generate matrix for number of groups=%d and group size=%d\n", m_nGroups, m_groupSize);
		printfRed("***        number of groups (NPlayers / GroupSize) with group size > 3 must be prime number\n");
		return false;
	}
	return true;
}
#endif
