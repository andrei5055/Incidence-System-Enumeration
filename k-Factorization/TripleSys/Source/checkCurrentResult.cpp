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
				case  1: noMoreResults = true; return 1;
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
			|| (m_useRowsPrecalculation == eCalculateRows && m_groupSize == 2) // significantly faster for gs=2
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
		}
#endif
	}
	m_playerIndex = 0;
	return 0;
}
