#include "TripleSys.h"

CC bool alldata::cnvCheck2U1F(int nrows)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	TrCycles trCycles01;
	bool ok = getCyclesAndPath(&trCycles01, 1, neighbors(0), neighbors(1)) > 0;
	ASSERT(!ok);

	bool bUseTestedTrs = param(t_autSaveTestedTrs) > 0;

	// get first row
	for (int indRow0 = 0; indRow0 < nrows; indRow0++)
	{
		// get second row
		for (int indRow1 = 0; indRow1 < nrows; indRow1++)
		{
			if (indRow0 == indRow1)
				continue;
			TrCycles trCycles;
			CGroupInfo* pTestedTRs = bUseTestedTrs ? testedTRs(indRow1 * m_numDays + indRow0) : NULL;
			if (bUseTestedTrs) {
				if (indRow0 < nrows - 1 && indRow1 < nrows - 1) {
					int nTrs = pTestedTRs->numObjects();
					for (int itr = 0; itr < nTrs; itr++) {
						tchar* trt = pTestedTRs->getObject(itr);
						m_TrInd++;
						if (pTestedTRs && pTestedTRs->isProcessed(trt))
							continue;
						const int icmp = kmProcessMatrix(result(), trt, nrows);
						if (icmp == 0)
							updateGroup(tr);
						else if (icmp < 0)
						{
							bRet = false;
							goto ret;
						}
					}
					continue;
				}
				else
					pTestedTRs->resetGroupOrder();
			}
			bool ok = getCyclesAndPath(&trCycles, 1, neighbors(indRow0), neighbors(indRow1)) > 0 &&
				!MEMCMP(trCycles01.length, trCycles.length, MAX_CYCLES_PER_SET);
			if (!ok) {
				printTable("result", result(), nrows, m_numPlayers, 2);
				printTable("resi", result(indRow0), 1, m_numPlayers, 2);
				printTable("resj", result(indRow1), 1, m_numPlayers, 2);
				printTable("neii", neighbors(indRow0), 1, m_numPlayers, 2);
				printTable("neij", neighbors(indRow1), 1, m_numPlayers, 2);
			}

			ASSERT(!ok);
			ctchar *pDir, *pStartOut;
			auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, 2, &pDir, &pStartOut);

			do {
				const bool btr = createU1FTr(tr, &trCycles01, &trCycles, pDir, pIdx, pStartOut);

				ASSERT(!btr);

				m_TrInd++;
#if !USE_CUDA
				if (m_cnvMode) {
					cnvPrintAuto(tr, nrows);
					continue; // print only
				}
#endif
				const int icmp = kmProcessMatrix(result(), tr, nrows);
				//TestkmProcessMatrix(nrows, nrows, tr, tr, icmp);

				if (icmp == 1 || (bUseTestedTrs && pTestedTRs->isProcessed(tr)))
					continue;

				if (icmp == 0)
				{
					updateGroup(tr);
				}
				else if (icmp < 0)
				{
					bRet = false;
					goto ret;
				}

			} while (ProceedToNextMapping());
		}
	}
ret:
	return bRet;
}
