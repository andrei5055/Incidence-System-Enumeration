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
	for (int iRowLast = 1; iRowLast < nrows; iRowLast++) {
		bool bSaveTestedTrs = bUseTestedTrs && (iRowLast < m_numDaysResult - 1);
		for (int indRow = 0; indRow < iRowLast; indRow++)
		{
			// get second row
			for (int iRowSwap = 0; iRowSwap < 2; iRowSwap++)
			{
				const int indRow1 = iRowSwap ? indRow : iRowLast;
				const int indRow0 = iRowSwap ? iRowLast : indRow;
				TrCycles trCycles;
				auto *pTestedTRs = bUseTestedTrs ? m_pTrRepo->getTrSet(indRow1, indRow0) : NULL;
				if (bUseTestedTrs) {
					if (iRowLast <= m_lastRowWithTestedTrs) {
						int nTrs = pTestedTRs->numObjects();
						for (int itr = 0; itr < nTrs; itr++) {
							const auto* trt = pTestedTRs->getObject(itr);
							const int icmp = kmProcessMatrix(result(), trt, nrows);
							m_TrInd++;
							if (icmp < 0)
							{
								bRet = false;
								goto ret;
							}
							if (icmp == 0)
								updateGroup(tr);
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
				ctchar* pDir, * pStartOut;
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

					if (icmp < 0)
					{
						bRet = false;
						goto ret;
					}
					// save Tr if icmp is equal 0, 1, 2, or 3; continue if it was already processed
					if (bSaveTestedTrs && pTestedTRs->isProcessed(tr))
						continue;

					if (icmp == 0)
						updateGroup(tr);

				} while (ProceedToNextMapping());
			}
		}
		if (bSaveTestedTrs && m_lastRowWithTestedTrs < iRowLast)
			m_lastRowWithTestedTrs = iRowLast;
	}
ret:
	return bRet;
}
