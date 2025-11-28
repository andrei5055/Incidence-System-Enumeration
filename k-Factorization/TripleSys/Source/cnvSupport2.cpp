#include "TripleSys.h"

CC bool alldata::cnvCheck2U1F(int nrows, int nrowsToUseForTrs)
{
	if (nrows < 2)
		return true;

	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	TrCycles trCycles;
	TrCycles trCycles01;
#if 0
	static tchar a[] = { 0,3,1,4,2,7,5,8,6,9 };
	if (memcmp(a, result(1), sizeof(a)) == 0)
		bRet = bRet;
#endif
	bool bCBMP = !completeGraph();
	bool bok;
	if (bCBMP)
		bok = getCyclesAndPathCBMP(&trCycles01, neighbors(0), neighbors(1), result(0), result(1), 0, eCheckErrors) > 0;
	else
		bok = getCyclesAndPathFromNeighbors(&trCycles01, neighbors(0), neighbors(1), NULL, NULL, eCheckErrors) > 0;
	ASSERT_IF(!bok);
	bool bUseTestedTrs = param(t_autSaveTestedTrs) > 0;
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);

	// get first row
	for (int iRowLast = 1; iRowLast < nrowsToUseForTrs; iRowLast++) {
		bool bSaveTestedTrs = bUseTestedTrs && (iRowLast < m_numDaysResult - 1);
		for (int indRow = 0; indRow < iRowLast; indRow++)
		{
			// get second row
			for (int iRowSwap = 0; iRowSwap < 2; iRowSwap++)
			{
				const int indRow1 = iRowSwap ? indRow : iRowLast;
				const int indRow0 = iRowSwap ? iRowLast : indRow;
				auto* pTestedTRs = bUseTestedTrs ? m_pTrRepo->getTrSet(indRow1, indRow0) : NULL;
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
								if (m_doNotExitEarlyIfNotCanonical)
									continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
								goto ret;
							}
							if (icmp == 0)
								updateGroup(trt);
						}
						continue;
					}
					else if (bSaveTestedTrs)
						pTestedTRs->resetGroupOrder();
				}
				bool bok;
				int nitr = 1; //bCBMP ? m_groupSizeFactorial : 1;
				for (int itr = 0; itr < nitr; itr++) {
					if (bCBMP) {
						bok = getCyclesAndPathCBMP(&trCycles, neighbors(indRow0), neighbors(indRow1),
							result(indRow0), result(indRow1), itr, eNoErrorCheck) > 0;
					}
					else {
						bok = getCyclesAndPathFromNeighbors(&trCycles, neighbors(indRow0), neighbors(indRow1), NULL, NULL, eNoErrorCheck) > 0;
					}
					if (!bok || MEMCMP(trCycles01.length, trCycles.length, MAX_CYCLES_PER_SET)) {/**
						printTable("result", result(), nrows, m_numPlayers, 2);
						printTable("resi", result(indRow0), 1, m_numPlayers, 2);
						printTable("resj", result(indRow1), 1, m_numPlayers, 2);
						printTable("neii", neighbors(indRow0), 1, m_numPlayers, 2);
						printTable("neij", neighbors(indRow1), 1, m_numPlayers, 2);*/
						if (any2RowsConvertToFirst2) {
							bRet = false;
							if (m_doNotExitEarlyIfNotCanonical)
								continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
							goto ret;
						}
						continue;
					}

					ASSERT_IF(!bok);
					ctchar* pDir, * pStartOut;
					auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, 2, &pDir, &pStartOut);

					do {
#if 0
						printTable("pDir", pDir, 1, 2, 2);
						printTable("pIdx", pIdx, 1, 2, 2);
						printTable("pStartOut", pStartOut, 1, 2, 2);
#endif
						//if (pDir[0])
							//continue;
						const bool btr = createU1FTr(tr, &trCycles01, &trCycles, pDir, pIdx, pStartOut);

						ASSERT_IF(!btr);

						m_TrInd++;
#if !USE_CUDA
						if (m_cnvMode) {
							cnvPrintAuto(tr, nrows);
							continue; // print only
						}
#endif
						if (bCBMP) {
							tchar g0 = (tr[0] & 1), g1 = (tr[1] & 1);
							if (g0 == g1)
								continue;
							auto i = 2;
							for (; i < m_numPlayers; i += 2) {
								if ((tr[i] & 1) != g0 || (tr[i + 1] & 1) != g1)
									break;
							}
							if (i != m_numPlayers)
								continue;
						}
						const int icmp = kmProcessMatrix(result(), tr, nrows);

						if (icmp < 0)
						{
							//TestkmProcessMatrix(nrows, nrows, tr, tr, 0);//icmp);
							bRet = false;
							if (m_doNotExitEarlyIfNotCanonical)
								continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
							goto ret;
						}
						// save Tr (if icmp not -1 or 1), continue if Tr was already processed
						//if (bSaveTestedTrs && icmp != 1 && pTestedTRs->isProcessed(tr))
						if (bSaveTestedTrs && pTestedTRs->isProcessed(tr))
							continue;

						if (icmp == 0)
							updateGroup(tr);

					} while (ProceedToNextMapping());
				}
			}
		}
		if (bSaveTestedTrs && m_lastRowWithTestedTrs < iRowLast)
			m_lastRowWithTestedTrs = iRowLast;
	}
ret:
	return bRet;
}
