#include "TripleSys.h"
#define FastMode 0
CC bool alldata::cnvCheck2P1F(int nrows)
{
	if (nrows < 2) // || (m_useRowsPrecalculation == eCalculateRows && nrows != m_numDaysResult))
	//if (nrows < 2)
		return true;

	if (nrows == 3)
		nrows = nrows;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;

	auto* neighbors0 = neighbors(0);
	auto* neighbors1 = neighbors(1);
	int iRowStart = (FastMode && nrows != m_numDaysResult) ? nrows - 2 : 0;
	if ((param(t_p1f_counter) && !(m_p1f_counter % param(t_p1f_counter))))
		iRowStart = 0;
	// get first row
	//printTable("tr", result(), nrows, m_numPlayers);
	for (int indRow0 = iRowStart; indRow0 < nrows; indRow0++)
	{
		auto* neighborsi = neighbors(indRow0);
		// get second row
		for (int indRow1 = iRowStart; indRow1 < nrows; indRow1++)
		{
			if (indRow0 == indRow1)
				continue;
			auto* neighborsj = neighbors(indRow1);
			//printTable("matrix", neighbors0, nrows, m_numPlayers);
			// get value of 0
			int nk = m_numPlayers;
			for (int k = 0; k < nk; k++)
			{
				if (!create2P1FTr(tr, k, neighbors0, neighbors1, neighborsi, neighborsj))
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
#define Test1 0
	if (nrows < 2)// || (m_useRowsPrecalculation == eCalculateRows && nrows != m_numDaysResult))
		return true;
	tchar trLocal[MAX_PLAYER_NUMBER];
	tchar* tr = trLocal;
	bool bRet = true;

	auto v1 = getV1();
	int ip1 = 0;
	auto* neighbors0 = neighbors(0);
	tchar* p1 = NULL;
	bool bCurrentSet = false;
	//if (result(1)[5] == 9 && iDay == 4)
	//	printf("%d ", iDay);
	const int maxv1 = MAX_3PF_SETS;
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);

#define _StatAdd(x, y, z)  // StatAdd(x, y, z)
	_StatAdd("AllcnvCheck3U1F", 10, true);
	CGroupInfo* pTestedTRs = param(t_autSaveTestedTrs) >= 2 ? testedTRs() : NULL;
	if (pTestedTRs)
		pTestedTRs->resetGroupOrder();
	//static tchar a[] = {0,3,6, 1,9,12, 2,15,18, 4,10,16, 5,13,19, 7,11,20, 8,14,17, };
	//if (memcmp(a, result(1), sizeof(a)) == 0)
	//	ip1 = ip1;
	int iRowStart = (FastMode && nrows != m_numDaysResult) ? nrows - 2 : 0;
	if ((param(t_p1f_counter) && !(m_p1f_counter % param(t_p1f_counter))))
		iRowStart = 0;
	while (1)
	{
		//if (!m_pSecondRowsDB)
		if ((any2RowsConvertToFirst2 && nrows != 2) || !m_pSecondRowsDB)
			p1 = result(1);
		else if (m_pSecondRowsDB->numObjects() > ip1)
			p1 = m_pSecondRowsDB->getObject(ip1);
		else if (m_createSecondRow) // do not merge with first if 
			p1 = result(1);
		else {
			bRet = false;
			break;
		}
		ip1++;
		memset(&m_TrCycles, 0, sizeof(m_TrCycles));
		cyclesFor2Rows(p1); // result is in m_TrCyclesAll
		// check that all and only requested cycles are present
		if (nrows == 2) {
			bool bAllCyclesOk = false;
			auto u1fPntr = sysParam()->u1fCycles[0];
			int itr0 = 0;
			if (!u1fPntr) {
				bAllCyclesOk = m_TrCyclesAll[0].length[0] == m_numPlayers && !m_TrCyclesAll[1].counter;
				itr0 = 1;
			}
			else {
				for (; itr0 < MAX_3PF_SETS; itr0++)
				{
					if (m_TrCyclesAll[itr0].counter == 0) {
						bAllCyclesOk = itr0 == *u1fPntr;
						break;
					}
					if (itr0 >= *u1fPntr || 
						MEMCMP(m_TrCyclesAll[itr0].length, u1fPntr + 1 + itr0 * MAX_CYCLES_PER_SET, MAX_CYCLES_PER_SET))
						break;
				}
			}
			if (itr0 != MAX_3PF_SETS && !bAllCyclesOk) {
				bRet = false;
				break;
			}
		}
		if (MEMCMP(p1, result(1), m_numPlayers) == 0)
			bCurrentSet = true; 
		//for (int isw = 1; isw >= 0; isw--)
		{
			// get first row
			for (int indRow0 = iRowStart; indRow0 < nrows; indRow0++)
			{
				// get second row
				for (int indRow1 = iRowStart; indRow1 < nrows; indRow1++)
				{
					//if (isw == 0 && indRow0 != nrows - 1 && indRow1 != nrows - 1)
					//	continue;
					if (indRow0 == indRow1)
						continue;
					const int nv1 = getAllV(v1, maxv1, indRow0, indRow1);
					//if (nv1 > 53)
					//	indRow1 = indRow1;
					bool bPair = false;
					//bool bLastRow = indRow0 == nrows - 1 || indRow1 == nrows - 1;
					auto* pV1 = v1;
					if (any2RowsConvertToFirst2 && nrows != 2 && (indRow0 > 2 || indRow1 > 2)) {
						TrCycles trCycles;
						memset(m_TrCyclesPair, 0, sizeof(m_TrCyclesPair));
						for (int iv1 = 0; iv1 < nv1; iv1++, pV1 += m_nGroups)
						{
							if (!getCyclesAndPath3(&trCycles, pV1, neighbors(indRow0), neighbors(indRow1), result(indRow0), result(indRow1))) {
								ASSERT(1);
								continue;
							}
							collectCyclesAndPath(m_TrCyclesPair, &trCycles, true);
						}

						for (int itr0 = 0; itr0 < MAX_3PF_SETS; itr0++)
						{
							if (m_TrCyclesAll[itr0].counter != m_TrCyclesPair[itr0].counter ||
								MEMCMP(m_TrCyclesAll[itr0].length, m_TrCyclesPair[itr0].length, sizeof(m_TrCyclesAll[0].length))) {
								bRet = false;
								//printf("+");
								goto ret;
							}
							if (m_TrCyclesAll[itr0].counter == 0)
								break;
						}
					}
					pV1 = v1;
					for (int iv1 = 0; iv1 < nv1; iv1++, pV1 += m_nGroups)
					{
						TrCycles trCycles;
						if (!getCyclesAndPath3(&trCycles, pV1, neighbors(indRow0), neighbors(indRow1), result(indRow0), result(indRow1))) {
							ASSERT(1);
							continue;
							//bRet = false;
							//goto ret;
						}
						bool bCycleSelected = false;
						for (int itr0 = 0; itr0 < MAX_3PF_SETS; itr0++)
						{
							if (m_TrCyclesAll[itr0].counter == 0)
								break;
							if (MEMCMP(m_TrCyclesAll[itr0].length, trCycles.length, MAX_CYCLES_PER_SET))
								continue;
							bCycleSelected = true;
							ctchar* pDir, * pStartOut;
							auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, 3, &pDir, &pStartOut);

							do {
								_StatAdd("create3U1FTr", 11, true);

								const bool btr = createU1FTr(tr, &m_TrCyclesAll[itr0], &trCycles, pDir, pIdx, pStartOut);
								if (btr)
									bPair = bPair;
								//continue;
#if 0 && !USE_CUDA 			// if btr == false, print tr, cycles and full pathes for rows (0, 1) and (indRow0, indRow1)
								if (btr && indRow0 == 2 && indRow1 == 3)
								{
									int numCycles1 = m_TrCyclesAll[itr0].ncycles, numCycles2 = trCycles.ncycles;
									printf("\nRows(%d,%d) indRow0=%d indRow1=%d iv1=%d\nrows01:",
										numCycles1, numCycles2, indRow0, indRow1, iv1);
									for (int jl = 0; jl < numCycles1; jl++)
										printf(" cycle%d=%d(starts at %d)", jl, m_TrCyclesAll[itr0].length[jl], m_TrCyclesAll[itr0].start[jl]);
									printf("\nrows(%d%d):", indRow0, indRow1);
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
									_StatAdd("TR_created", 12, true);
									if (Test1) {
										int imax = MAX2(indRow0, indRow1);
										int imin = MIN2(indRow0, indRow1);
										if (bCurrentSet && (imax != 1 && imax != 6) || imin > 1) {
											bRet = false;
											//goto ret;
										}
									}
									bPair = true;
 									m_TrInd++;
#if !USE_CUDA 
									if (0)//indRow0 == 2 && indRow1 == 3)
									{
										int k = 0, numCycles = trCycles.ncycles;
										printf(" %d:%d:%d", itr0, m_TrInd, pTestedTRs->numObjects());
										//for (int j = 0; j < numCycles; j++)
										//	printf("  %2d %2d %2d", pDir[j], pIdx[j], pStartOut[j]);
										//printf("\n");
									}
#endif
									if (pTestedTRs && pTestedTRs->isProcessed(tr))
									{
										continue;
									}
#if !USE_CUDA
									//if (itr0 > 0)
									//	printfRed("%d ", itr0);
									if (m_cnvMode) {
										cnvPrintAuto(tr, nrows);
										continue; // print only
									}
#endif
									if (!bCurrentSet)
										bCurrentSet = bCurrentSet;


									const int icmp = kmProcessMatrix(result(), tr, nrows);
									//const int icmp = !bCurrentSet ? -1 : kmProcessMatrix(result(), tr, nrows);
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
#if 0 //!USE_CUDA
										static int a;
										if (nrows > 2)
											printf(" %d:%d", a++, m_playerIndex);
										if (a > 100)
											exit(0);
#endif
										bRet = false;
										goto ret;
									}
								}
							} while (ProceedToNextMapping());
						}
						if (!bCycleSelected && !allowNotSelectedCycles()) {
							bRet = false;
							goto ret;
						}
					}
					if (bCurrentSet)
					{
						//if (bPair)
						//	printf("indRow0=%d indRow1=%d can be converted to 0,1\n", indRow0, indRow1);
						if (!bPair)
						{
							//if (nrows == 10)
							//printf("indRow0=%d indRow1=%d can't be converted to 0,1\n", indRow0, indRow1);
							//bPair = bPair;
							if (any2RowsConvertToFirst2) {
								bRet = false;
								goto ret;
							}
						}
					}
				}
			}
		}
		if (bCurrentSet)
			break;
	}
ret:
	//if (bRet && nrows == 4)
	//	exit(1);
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

