#include "TripleSys.h"

#if !USE_CUDA && USE_BINARY_CANONIZER
#include "k-SysSupport.h"
#include "CDTools.h"
#endif

CC void alldata::goBack()
{
	ASSERT(m_playerIndex >= iDay * m_numPlayers,
		printfRed("*** Request to go 'back in future' request (from player %d to player %d)\n", iDay * m_numPlayers - 1, m_playerIndex);
		printTable("Input matrix", result(), iDay, m_numPlayers, m_groupSize, 0, true);
		abort();
	)

	if (iDay == m_numDays)
	{
		if (m_playerIndex > (m_numDays - 1) * m_numPlayers - m_groupSize - 1)
			m_playerIndex = (m_numDays - 1) * m_numPlayers - m_groupSize - 1;
	}
	iDay--;
	while (iDay > m_playerIndex / m_numPlayers)
	{
		while (iPlayer >= 0)
		{
			getPrevPlayer();
		}
		initPrevDay();
	}
	while (iDay * m_numPlayers + iPlayer > m_playerIndex)
	{
		getPrevPlayer();
	}
	m_playerIndex = 0;
}
CC int alldata::checkCurrentResult(int iPrintMatrices, void* pIS_Canonizer)
{
	// function returns : -1 - prev result, 0 - continue, 1 - eoj
	m_playerIndex = iDay * m_numPlayers - m_groupSize - 1;
	if (iDay > 1)
	{
		//if (!checkLinksV2(links(), iDay))
		//	return -1;
		if (m_use2RowsCanonization) {
			if (iDay == 2)
			{
				switch (p1fCheck2ndRow()) {
					case -1: return -1;
					case  1: return 1;
					case  0: break;
				}
			}
		}

#if 1 // set to 0 to disable all improvements
#if 0
		if (m_bCheckLinkT && !checkLinksT(links(), iDay))
			return -1;
#endif
		if (param(t_useImproveMatrix) && improveMatrix(m_improveResult, NULL, 0/*, bResults, lenResult()*/))
			return -1;
#if 1
		if ((iDay == numDaysResult())
			|| checkCanonicity()
			|| (param(t_submatrixGroupOrderMin) > 0)
			|| (param(t_nestedGroups) > 1)
			|| (m_use2RowsCanonization && m_groupSize == 3) // slightly faster (for 15,7,3)
			|| (m_precalcMode == eCalculateRows && m_groupSize == 2) // significantly faster for gs=2
			)
		{
			bool bPrev = true;
#if !USE_CUDA && USE_BINARY_CANONIZER
			if (m_ppBinMatrStorage) {
				if (pIS_Canonizer && m_ppBinMatrStorage[iDay]) {
					const auto* pCanonBinaryMatr = runCanonizer(pIS_Canonizer, result(0), m_groupSize, iDay < m_numDays ? iDay : 0);
					if (m_ppBinMatrStorage[iDay]->updateRepo(pCanonBinaryMatr) < 0) {
						bPrev = false;
					}
				}
			}
			else
#endif
			{
#define LOOP_LENGTH1		0
#if LOOP_LENGTH1
				for (int i = 0; i < LOOP_LENGTH1; i++)
					bPrev = !cnvCheckNew(0, iDay);
#else
				bPrev = cnvCheckNew(0, iDay);
				if (iPrintMatrices & 8) {
					if (!bPrev)
						printf(" %d", m_playerIndex);
					else
						printf(".");
				}
#if !USE_CUDA
				if ((iPrintMatrices & 4) && iDay > 2)
					printPermutationMatrices(3);
#endif
#endif
			}
#endif
			// print stat result after all transitions applied
			StatReportAfterAllTr(ResetStat, "Stat for one improvement. iDay", iDay, bPrint);
			if (!bPrev || orderOfGroup() < param(t_submatrixGroupOrderMin))
				return -1;
			if (!semiCheck())
				return -1;
		}
#endif
	}
	m_playerIndex = 0;
	return 0;
}


#include <mutex>

extern std::mutex mtxLinks; // The mutex to protect the shared resource
extern CStorageIdx<tchar>** mpLinks;
extern CStorageIdx<tchar>** mShLinks;
extern CStorageIdx<tchar>** mShLinks2;
extern int SemiPhase;
tchar testL[] = {
 0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0
,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0
,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0
,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0
,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0
,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0
,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0
,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0
,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0
,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1
,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1
,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0
,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0
,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1
,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0
,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1
,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0
};
void L2L(tchar* out, tchar* in, int n, tchar notSet)
{
	tchar* outi = out;
	tchar* ini = in;
	memset(out, 0, n * n / 4);
	for (int i = 0; i < n; i++, ini += n) {
		for (int j = i + 1; j < n; j++) {
			if (ini[j] != notSet) {
				int iout = (i & 1) ? j / 2 : i / 2;
				int jout = (i & 1) ? i / 2 : j / 2;
				out[iout * (n / 2) + jout] = 1;
			}
		}
	}
}
bool add2DBs(tchar* lnkT) {
	// Phase 1: create DB
	/** check for trs with odd to even swap
	for (int i = 0; i < numObjects(); i++) {
		const auto tr = getObject(i);
		if ((tr[0] & 1) == 1 && (tr[1] & 1) == 0)
			return false;
	}**/
	int itr = mpLinks[0]->numObjects();
	if (mpLinks[0]->isProcessed(lnkT)) {
		int itr3 = mShLinks[0]->updateRepo(lnkT);
		mShLinks2[0]->updateRepo(lnkT);
		printf("Same %d(%d,%d)\n", itr, mShLinks[0]->numObjects(), itr3);
	}
	return true;
}
int nPairsLeft(tchar* lnkT, int nnd4) {
	// Phase 2 or 3: DBs with links already created
	int n = 0;
	for (int i = 0; i < nnd4; i++) {
		if (lnkT[i])
			n++;
	}
	return n;
}
int adjustPairsCount(tchar* tr, int np, tchar* lnkT, int nGroupsToCheck) {
	tchar gr0 = 0, gr1 = 0;
	for (int j = 1; j < np; j++) {
		if (tr[j] == 0)
			gr0 = j;
		else if (tr[j] == 1)
			gr1 = j;
	}
	if (gr0 & 1)
		SWAP(gr0, gr1);
	int i0 = gr0 / 2 * (np / 2) + gr1 / 2;
	ASSERT(i0 < 0 || i0 >= np * np / 4);
	if (lnkT[i0]) {
		lnkT[i0] = 0;
		nGroupsToCheck -= 1;
		ASSERT(nGroupsToCheck < 0);
	}
	return nGroupsToCheck;
}
bool alldata::semiCheck()
{
	tchar lnk[16 * 16];
	if (param(t_semiSymmetricGraphs) != iDay)
		return true;
	if (orderOfGroup() < 2)
		return false;
	
	int nn = m_numPlayers * m_numPlayers;
	if (orderOfGroup() == 400)
		nn = nn;
	int nnd4 = nn / 4;
	int nGroupsToCheck = 0;
	tchar* lnkT = links(m_numPlayers);
	bool bSharedLink = false;
	L2L(lnk, links(), m_numPlayers, unset);
	if (m_test & 0x4) {
		L2L(lnkT, testL, m_numPlayers, 0);
		if (MEMCMP(lnkT, lnk, nnd4))
			return false;
	}
	else
		memcpy(lnkT, lnk, nnd4);
	std::lock_guard<std::mutex> lock(mtxLinks);
	if (mpLinks == NULL) {
		SemiPhase = 1;
		mpLinks = new CStorageIdx<tchar>*[1];
		mpLinks[0] = new CStorageIdx<tchar>(1000, nnd4);
		mShLinks = new CStorageIdx<tchar>*[1];
		mShLinks[0] = new CStorageIdx<tchar>(100, nnd4);
		mShLinks2 = new CStorageIdx<tchar>*[1];
		mShLinks2[0] = new CStorageIdx<tchar>(100, nnd4);
	}
	switch (SemiPhase) {
	case 1: return add2DBs(lnkT);
	case 2:
	case 3: {
		// Phase 2 or 3: DBs with links already created 
		int idx = -1;
		if (!(m_test & 0x8))
			idx = mShLinks[0]->findObject(lnk, 0, mShLinks[0]->numObjects());
		if (idx >= 0) {
			bSharedLink = true;
			lnkT = mShLinks2[0]->getObject(idx);
			if (0 && idx == 17) {
				printf("\nidx = %d aut=%d\n", idx, orderOfGroup());
				printTableColor("R", result(), iDay, m_numPlayers, m_groupSize);
				printTableColor("links (rows:0,2,4..., columns:1,3,5,...", lnk, m_numPlayers / 2, m_numPlayers / 2, 0);
				printTableColor("lnkT", lnkT, m_numPlayers / 2, m_numPlayers / 2, 0);
			}
		}
		nGroupsToCheck = nPairsLeft(lnkT, nnd4);
		if (SemiPhase == 3) {
			if (bSharedLink) {
				static int cc = 999; if (cc > nGroupsToCheck) printf(" %d.%d ", idx, cc = nGroupsToCheck);
				if (!nGroupsToCheck) {
					printf("\nidx = %d\n", idx);
					printTableColor("Multi matrix result", result(), iDay, m_numPlayers, m_groupSize);
					printTableColor("links (rows:0,2,4..., columns:1,3,5,...", lnk, m_numPlayers / 2, m_numPlayers / 2, 0);
					return true;
				}
			}
			return false;
		}
		// phase 2
		int i;
		for (i = 0; i < numObjects(); i++) {
			if (!nGroupsToCheck)
				break;
			const auto tr = getObjAddr(i);
			nGroupsToCheck = adjustPairsCount(tr, m_numPlayers, lnkT, nGroupsToCheck);
		}
		if (0 && idx == 17) {
			printf("\nidx = %d aut=%d nch=%d\n", idx, orderOfGroup(), nGroupsToCheck);
			printTableColor("R", result(), iDay, m_numPlayers, m_groupSize);
			printTableColor("links (rows:0,2,4..., columns:1,3,5,...", lnk, m_numPlayers / 2, m_numPlayers / 2, 0);
			printTableColor("lnk", lnk, m_numPlayers / 2, m_numPlayers / 2, 0);
			printTableColor("lnkT", lnkT, m_numPlayers / 2, m_numPlayers / 2, 0);
		}
		if (!bSharedLink) {
			static int ccc = 999; if (ccc > nGroupsToCheck) printf(" .%d ", ccc = nGroupsToCheck);
			if (!nGroupsToCheck) {
				printTableColor("One matrix result", result(), iDay, m_numPlayers, m_groupSize);
				printTableColor("links", lnk, m_numPlayers / 2, m_numPlayers / 2, 0);
				return true;
			}
		}
	}
	}
	return false;
}