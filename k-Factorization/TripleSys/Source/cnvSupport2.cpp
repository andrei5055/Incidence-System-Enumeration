#include "TripleSys.h"

CC bool alldata::cnvCheck2U1F(int nrows)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	TrCycles trCycles01;
	bool ok = getCyclesAndPath(&trCycles01, 1, p1ftable(0), p1ftable(1)) > 0;
	ASSERT(!ok);

#define _StatAdd(x, y, z)  // StatAdd(x, y, z)
	_StatAdd("AllcnvCheck2u1F", 10, true);
	CGroupInfo* pTestedTRs = param(t_autSaveTestedTrs) >= 2 ? testedTRs() : NULL;

	// get first row
	for (int indRow0 = 0; indRow0 < nrows; indRow0++)
	{
		// get second row
		for (int indRow1 = 0; indRow1 < nrows; indRow1++)
		{
			if (indRow0 == indRow1)
				continue;
			TrCycles trCycles;
			bool ok = getCyclesAndPath(&trCycles, 1, p1ftable(indRow0), p1ftable(indRow1)) > 0 &&
				!MEMCMP(trCycles01.length, trCycles.length, MAX_CYCLES_PER_SET);
			ASSERT(!ok);
			ctchar *pDir, *pStartOut;
			auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, 2, &pDir, &pStartOut);

			do {
				_StatAdd("create2U1FTr", 11, true);

				const bool btr = createU1FTr(tr, &trCycles01, &trCycles, pDir, pIdx, pStartOut);

#if 0 && !USE_CUDA 			// if btr == false, print tr, cycles and full pathes for rows (0, 1) and (indRow0, indRow1)
				if (!btr)//&& indRow0 >= 1 && indRow1 >= 2)
				{
					printf("btr=%d\n", btr);
					int numCycles1 = trCycles01.ncycles, numCycles2 = trCycles.ncycles;
					printf("\nnCycles=(%d,%d) indRow0=%d indRow1=%d\nrows01:",
						numCycles1, numCycles2, indRow0, indRow1);
					for (int jl = 0; jl < numCycles1; jl++)
						printf(" cycle%d=%d(starts at %d)", jl, trCycles01.length[jl], trCycles01.start[jl]);
					printf("\nrows%d%d:", indRow0, indRow1);
					for (int jl = 0; jl < numCycles2; jl++)
						printf(" cycle%d=%d(starts at %d)", jl, trCycles.length[jl], m_TrCycles.start[jl]);
					printf("\ndir/idx/start:");
					for (int jl = 0; jl < numCycles2; jl++)
						printf(" %d/%d/%d", pDir[jl], pIdx[jl], pStartOut[jl]);
					printf("\n");

					printTable("tr  ", tr, 1, m_numPlayers, 2);
					printTable("res0", result(0), 1, m_numPlayers, 2);
					printTable("res1", result(1), 1, m_numPlayers, 2);
					printTable("resi", result(indRow0), 1, m_numPlayers, 2);
					printTable("resj", result(indRow1), 1, m_numPlayers, 2);
					printTable("fp01", trCycles01.fullPath, 1, m_numPlayers * 2, 2);
					printf("fp%d%d:", indRow0, indRow1);
					printTable("", trCycles.fullPath, 1, m_numPlayers * 2, 2);
				}
#endif
				ASSERT(!btr);

				_StatAdd("TR_created", 12, true);
				m_TrInd++;
				if (pTestedTRs && pTestedTRs->isProcessed(tr))
				{
					continue;
				}
#if !USE_CUDA
				if (m_cnvMode) {
					cnvPrintAuto(tr, nrows);
					continue; // print only
				}
#endif
				const int icmp = kmProcessMatrix(result(), tr, nrows);
				//TestkmProcessMatrix(nrows, nrows, tr, tr, icmp);

				_StatAdd("kmProcessMatrix", 13, bCurrentSet);
				Stat_cnvCheckKm1("cmp(2)", 2, icmp == 2);
				Stat_cnvCheckKm1("cmp(all)", 3, true);

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
	Stat_cnvCheckKm1("can(all)", 0, true);
	Stat_cnvCheckKm1("(-1)", 1, !bRet);
	return bRet;
}
