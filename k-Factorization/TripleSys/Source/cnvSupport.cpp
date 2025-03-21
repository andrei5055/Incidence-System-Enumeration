#include "TripleSys.h"
CC bool alldata::cnvCheck2P1F(int nrows)
{
	if (nrows < 2)
		return true;

	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;

	auto* neigbors0 = neighbors(0);
	auto* neigbors1 = neighbors(1);
	// get first row
	//printTable("tr", result(), nrows, m_numPlayers);
	for (int indRow0 = 0; indRow0 < nrows; indRow0++)
	{
		auto* neigborsi = neighbors(indRow0);
		// get second row
		for (int indRow1 = 0; indRow1 < nrows; indRow1++)
		{
			if (indRow0 == indRow1)
				continue;
			auto* neigborsj = neighbors(indRow1);
			//printTable("matrix", neigbors0, nrows, m_numPlayers);
			// get value of 0
			int nk = m_numPlayers;
			for (int k = 0; k < nk; k++)
			{
				if (!create2P1FTr(tr, k, neigbors0, neigbors1, neigborsi, neigborsj))
					continue;
				m_TrInd++;
#if !USE_CUDA
				if (m_cnvMode) {
					cnvPrintAuto(tr, nrows);
					continue; // print only
				}
#endif
				int icmp = nrows < 3 ? 0 : kmProcessMatrix2p1f(tr, nrows, indRow0, indRow1);
				//TestkmProcessMatrix(nrows, 0, tr, tr, icmp);
				Stat_cnvCheckKm1("cmp(2)", 2, icmp == 2);
				Stat_cnvCheckKm1("cmp(all)", 3, true);
				if (icmp == 0)
				{
					updateGroup(tr);
				}
				else if (icmp < 0)
				{
#if PRINT_TRANSFORMED
					printTransformed(nrows, m_numPlayers, m_groupSize, tr, tr, result(), m_Ktmp, 0, 0, 0);
#endif
					bRet = false;
					goto ret;
				}
			}
		}
	}
ret:
	Stat_cnvCheckKm1("can(all)", 0, true);
	Stat_cnvCheckKm1("(-1)", 1, !bRet);
	return bRet;
}
CC bool alldata::cnvCheck3U1F(int nrows)
{
	if (nrows < 2)
		return true;
	tchar trLocal[MAX_PLAYER_NUMBER];
	tchar* tr = trLocal;
	bool bRet = true;

	auto v0 = getV0();
	auto v1 = getV1();
	int ip1 = 0;
	auto* neigbors0 = neighbors(0);
	tchar* p1 = NULL;
	bool bCurrentSet = false;
	const auto nRowStart = param(t_nRowsInStartMatrix);
	//if (result(1)[5] == 9 && iDay == 4)
	//	printf("%d ", iDay);
	const auto u1fCycles = sysParam()->u1fCycles[0];
	int maxv0 = (u1fCycles && (u1fCycles[0] != 1 || u1fCycles[1] != m_numPlayers)) ? MAX_3PF_SETS : 1;
	if (u1fCycles && u1fCycles[0] == 1 && u1fCycles[1] == 9 && u1fCycles[2] == 9 && u1fCycles[3] == 9)
		maxv0 = 1;
	//const int maxv0 = MAX_3PF_SETS;
	const int maxv1 = MAX_3PF_SETS;

#define _StatAdd(x, y, z)  // StatAdd(x, y, z)
	_StatAdd("AllcnvCheck3U1F", 10, true);
	CGroupInfo* pTestedTRs = param(t_autSaveTestedTrs) >= 2 ? testedTRs() : NULL;
	//static tchar a[] = {0,3,6, 1,9,12, 2,15,18, 4,10,16, 5,13,19, 7,11,20, 8,14,17, };
	//if (memcmp(a, result(1), sizeof(a)) == 0)
	//	ip1 = ip1;
	while (1)
	{
		if (m_pSecondRowsDB->numObjects() > ip1)
			p1 = m_pSecondRowsDB->getObject(ip1);
		else if (m_createSecondRow)
			p1 = result(1);
		else {
			bRet = false;
			break;
		}
		ip1++;
		tchar neigbors1[MAX_PLAYER_NUMBER];
		memset(&m_TrCycles, 0, sizeof(m_TrCycles));
		memset(&m_TrCyclesAll, 0, sizeof(m_TrCyclesAll));
		if (MEMCMP(p1, result(1), m_numPlayers) == 0)
		{
			bCurrentSet = true;  // ANDREI
			memcpy(neigbors1, neighbors(1), m_numPlayers);
		}
		else
			u1fSetTableRow(neigbors1, p1);

		const int nv0 = getAllV(v0, maxv0, 0, 1, neigbors1);
		ASSERT(!nv0);
		for (int i = 0; i < nv0; i++)
		{
			ctchar* vtr = v0 + m_nGroups * i; // Common Values ANDREI
			const auto ncycles = p3Cycles(&m_TrCycles, MAX_CYCLE_SETS, neighbors(0), neigbors1, vtr, result(0), p1);
			//collectCyclesAndPath(&m_TrCycles);
		}

		// get first row
		for (int indRow0 = 0; indRow0 < nrows; indRow0++)
		{
			// get second row
			for (int indRow11 = nRowStart; indRow11 < nrows + nRowStart; indRow11++) 
			{
				int indRow1 = indRow11 % nrows;
				if (indRow0 == indRow1)
					continue;
				const int nv1 = getAllV(v1, maxv1, indRow0, indRow1);
				//if (nv1 > 53)
				//	indRow1 = indRow1;
				bool bPair = false;
				//bool bLastRow = indRow0 == nrows - 1 || indRow1 == nrows - 1;
				const auto* pV1 = v1;
				for (int iv1 = 0; iv1 < nv1; iv1++, pV1 += m_nGroups)  // Andrei nv1 is equal to 1
				{
					TrCycles trCycles;

					if (!getCyclesAndPath3(&trCycles, pV1, neighbors(indRow0), neighbors(indRow1), result(indRow0), result(indRow1)))
						continue;

					for (int itr0 = 0; itr0 < MAX_CYCLE_SETS; itr0++)
					{
						if (m_TrCyclesAll[itr0].counter == 0)
							break;
						if (MEMCMP(m_TrCyclesAll[itr0].length, trCycles.length, MAX_CYCLES_PER_SET))
							continue;

						ctchar *pDir, *pStartOut;
						auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, 3, &pDir, &pStartOut);

						do {
							_StatAdd("create3U1FTr", 11, true);
							
							const bool btr = createU1FTr(tr, &m_TrCyclesAll[itr0], &trCycles, pDir, pIdx, pStartOut);
							/**if (btr)
								bPair = true;
							continue;**/
#if 0 && !USE_CUDA 			// if btr == false, print tr, cycles and full pathes for rows (0, 1) and (indRow0, indRow1)
							if (btr)// && indRow0 >= 1 && indRow1 >= 2)
							{
								int numCycles1 = m_TrCyclesAll[itr0].ncycles, numCycles2 = trCycles.ncycles;
								printf("\nnCycles=(%d,%d) indRow0=%d indRow1=%d iv1=%d\nrows01:",
									numCycles1, numCycles2, indRow0, indRow1, iv1);
								for (int jl = 0; jl < numCycles1; jl++)
									printf(" cycle%d=%d(starts at %d)", jl, m_TrCyclesAll[itr0].length[jl], m_TrCyclesAll[itr0].start[jl]);
								printf("\nrows%d%d:", indRow0, indRow1);
								for (int jl = 0; jl < numCycles2; jl++)
									printf(" cycle%d=%d(starts at %d)", jl, trCycles.length[jl], m_TrCycles.start[jl]);
								printf("\ndir/idx/start:");
								for (int jl = 0; jl < numCycles2; jl++)
									printf(" %d/%d/%d", pDir[jl], pIdx[jl], pStartOut[jl]);
								printf("\n");
								
								printTable("tr  ", tr, 1, m_numPlayers, 3);
								printTable("res0", result(0), 1, m_numPlayers, 3);
								printTable("res1", result(1), 1, m_numPlayers, 3);
								printTable("resi", result(indRow0), 1, m_numPlayers, 3);
								printTable("resj", result(indRow1), 1, m_numPlayers, 3);
								printf("v%d%d:", indRow0, indRow1);
								printTable("", pV1, 1, m_nGroups, 0);
								printTable("fp01", m_TrCyclesAll[itr0].fullPath, 1, m_numPlayers * 2, 3);
								printf("fp%d%d:", indRow0, indRow1);
								printTable("", trCycles.fullPath, 1, m_numPlayers * 2, 3);
							}
#endif
							if (btr) {
								if (trCycles.ncycles != 1)
									bPair = bPair;
								_StatAdd("TR_created", 12, true);
								bPair = true;
								m_TrInd++;
								if (pTestedTRs && pTestedTRs->isProcessed(tr))
								{
#if !USE_CUDA 
									if (0)//indRow0 + indRow1 >= 2)
									{
										int k = 0, numCycles = trCycles.ncycles;
										printf("%4d:", ++k);
										for (int j = 0; j < numCycles; j++)
											printf("  %2d %2d %2d", pDir[j], pIdx[j], pStartOut[j]);
										printf("\n");
									}
#endif
									continue;
								}
#if !USE_CUDA
								if (m_cnvMode) {
									cnvPrintAuto(tr, nrows);
									continue; // print only
								}
#endif
								const int icmp = !bCurrentSet ? -1 : kmProcessMatrix(result(), tr, nrows);
								//TestkmProcessMatrix(nrows, nrows, tr, tr, icmp);

								_StatAdd("kmProcessMatrix", 13, bCurrentSet);
								Stat_cnvCheckKm1("cmp(2)", 2, icmp == 2);
								Stat_cnvCheckKm1("cmp(all)", 3, true);

								if (icmp == 0)
								{
#if 1
									updateGroup(tr);
#elif !USE_CUDA
									int mm = groupOrder();
									updateGroupOrder(tr);
									if (mm == groupOrder())
									{
										if (0)//indRow0 + indRow1 >= 2)
										{
											int k = 0, numCycles = trCycles.ncycles;
											printf("%4d:", ++k);
											for (int j = 0; j < numCycles; j++)
												printf("  %2d %2d %2d", pDir[j], pIdx[j], pStartOut[j]);
											printf("\n");
										}
									}
#endif
								}
								else if (icmp < 0)
								{
									bRet = false;
									goto ret;
								}
							}

						} while (ProceedToNextMapping());
					}
				}
				if (bCurrentSet)
				{
					//if (bPair)
					//	printf("indRow0=%d indRow1=%d can be converted to 0,1\n", indRow0, indRow1);
					if (!bPair)
					{
#if Any2RowsConvertToFirst2 == 1
						bRet = false;
						goto ret;
#endif
					}
				}
			}
		}
		if (bCurrentSet)
			break;
	}
ret:
	Stat_cnvCheckKm1("can(all)", 0, true);
	Stat_cnvCheckKm1("(-1)", 1, !bRet);
	if (bRet)
	{
		//printf("Trs total=%d\n", m_TrInd); Andrei 
		if (m_createSecondRow && nrows == numDaysResult() && p1 == result(1))
			memcpy(m_pSecondRowsDB->getNextObject(), p1, m_numPlayers);
	}
	return bRet;
}
#if !USE_CUDA
void alldata::cnvPrintAuto(ctchar* tr, int nrows)
{
	ASSERT(nrows < 2);
	tchar ttr[MAX_PLAYER_NUMBER];
	memset(ttr, 0, m_numPlayers);
	for (int iSwap = 0; iSwap < 2; iSwap++)
	{
		const auto imode = abs(m_cnvMode);
		const auto id = imode > 1 ? imode - 2 : nrows - 2;
		const auto* resn = result(id + iSwap);
		for (int i = 0; i < m_numPlayers; i++)
			ttr[resn[i]] = tr[i];
		auto icmp = kmProcessMatrix(result(id), ttr, 2, NULL);
		if ((icmp = MEMCMP(m_Km, result(), m_numPlayers * 2)) == 0)
		{
			printf("%c-TrInd=%d (transition rows %d,%d to 0,1)\n", m_cnvMode < 0 ? 'K' : 'P', m_TrInd, id + iSwap, id + 1 - iSwap);
			printTable("input", result(), id + 2, m_numPlayers, m_groupSize);
			printTable("ttr", ttr, 1, m_numPlayers, m_groupSize);
		}
	}
}
#endif

