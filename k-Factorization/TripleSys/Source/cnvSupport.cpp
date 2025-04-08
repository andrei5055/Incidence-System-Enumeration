#include "TripleSys.h"

CC bool alldata::cnvCheck2P1F(int nrows)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;

	auto* neighbors0 = neighbors(0);
	auto* neighbors1 = neighbors(1);
	// get first row
	for (int indRow0 = 0; indRow0 < nrows; indRow0++)
	{
		auto* neighborsi = neighbors(indRow0);
		// get second row
		for (int indRow1 = 0; indRow1 < nrows; indRow1++)
		{
			if (indRow0 == indRow1)
				continue;
			auto* neighborsj = neighbors(indRow1);
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
	return bRet;
}
CC bool alldata::cnvCheck3U1F(int nrows)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	bool bCollectInfo = param(t_printMatrices) & 16;
	bool bPrintInfo = (nrows == numDaysResult()) && bCollectInfo;
	auto v1 = getV1();
	int ip1 = 0;
	auto* neighbors0 = neighbors(0);
	tchar* p1 = NULL;
	bool bCurrentSet = false;
	const int maxv1 = MAX_3PF_SETS;
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);
	setAllowNotSelectedCycles(nrows);

	int nrr = param(t_useRowsPrecalculation);
	bool bPrecalcRow = m_useRowsPrecalculation == eCalculateRows && nrows > nrr;
	bool bUseTestedTrs = param(t_autSaveTestedTrs) > 0 && !bCollectInfo;

	//static tchar a[] = {0,3,6, 1,9,12, 2,15,18, 4,10,16, 5,13,19, 7,11,20, 8,14,17, };
	//if (memcmp(a, result(1), sizeof(a)) == 0)
	//	ip1 = ip1;
	while (1)
	{
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
		// for second row check that all and only requested cycles are present
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
					// warning! cycles defined in params must be sorted
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
			TrCycles trCycles;
			// get first row
			for (int indRow0 = 0; indRow0 < nrows; indRow0++)
			{
				// get second row
				for (int indRow1 = 0; indRow1 < nrows; indRow1++)
				{
					if (indRow0 == indRow1)
						continue;
					int nv1 = 0;
					int nTrsForPair = 0;
					CGroupInfo* pTestedTRs = testedTRs(indRow0 * m_numDays + indRow1);
					int nTrs = pTestedTRs->numObjects();
					if (indRow0 < nrows - 1 && indRow1 < nrows - 1 && bUseTestedTrs && !bPrecalcRow && nrows > 2) {
						for (int itr = 0; itr < nTrs; itr++) {
							tchar* trt = tr;
							trt = pTestedTRs->getObject(itr);
#if 0    // To take tr's in order of their 
							trt = pTestedTRs->CStorage<tchar>::getObject(itr);
#endif
							m_TrInd++;
							nTrsForPair++;
							const int icmp = kmProcessMatrix(result(), trt, nrows);
							if (icmp == 0)
								updateGroup(trt);
							else if (icmp < 0) {
#if !USE_CUDA
								if (bCollectInfo) {
									printTable("Tr", trt, 1, m_numPlayers, m_groupSize);
									printTable("Matrix with Tr applied", m_Ktmp, nrows, m_numPlayers, m_groupSize);
								}
#endif
								bRet = false;
								goto ret;
							}
						}
					}
					else {
						//if (nv1 > 53)
						//	indRow1 = indRow1
						bool bSaveTestedTr = (indRow0 < m_numDays - 1) && (indRow1 < m_numDays - 1) && !bPrecalcRow && bUseTestedTrs;
						auto* pV1 = v1;
						pTestedTRs->resetGroupOrder();
						nv1 = getAllV(v1, maxv1, indRow0, indRow1);
#if !USE_CUDA
						if (bPrintInfo) {
							printCyclesInfo(v1, nv1, indRow0, indRow1);
						}
#endif
						for (int itr = 0; itr < nv1; itr++, pV1 += m_nGroups) {
							if (!getCyclesAndPath3(&trCycles, pV1, neighbors(indRow0), neighbors(indRow1), result(indRow0), result(indRow1))) {
								//ASSERT(1);
								continue;
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
									const bool btr = createU1FTr(tr, &m_TrCyclesAll[itr0], &trCycles, pDir, pIdx, pStartOut);
									if (btr) {
										m_TrInd++;
										nTrsForPair++;
										const int icmp = kmProcessMatrix(result(), tr, nrows);

										switch (icmp) {
										case -1:
#if !USE_CUDA
											if (bCollectInfo) {
												printCyclesInfoSummary(&trCycles, tr, indRow0, indRow1, nrows);
												printfRed("Not canonical matrix. See Tr below from Cycles(%d:%d:%d), rows %d,%d\n",
													trCycles.length[0], trCycles.length[1], trCycles.length[2], indRow0, indRow1);
												printTable(" Tr", tr, 1, m_numPlayers, m_groupSize);
												printTable(" Matrix with Tr applied", m_Ktmp, nrows, m_numPlayers, m_groupSize);
											}
#endif
											bRet = false;
											goto ret;
										case 1: continue;
										default: break;
										}
										if (bSaveTestedTr && pTestedTRs->isProcessed(tr))
											continue;
										if (icmp == 0)
											updateGroup(tr);
									}
								} while (ProceedToNextMapping());
								break;
							}
							if (!bCycleSelected) {
								if (!allowNotSelectedCycles()) {
#if !USE_CUDA
									if (bCollectInfo)
										printfRed("Cycles (%d:%d:%d) for rows %d,%d are present, but not selected\n",
											trCycles.length[0], trCycles.length[1], trCycles.length[2], indRow0, indRow1);
#endif
									bRet = false;
									goto ret;
								}
								break;
							}
						}
					}
#if !USE_CUDA
					if (bPrintInfo) {
						char stat[128];
						matrixStatOutput(stat, sizeof(stat), m_TrCyclesPair);
						printf("nTr=(%3d,%-3d) GroupOrder=%-2d Cycles for rows %d,%d: %s\n",
							nTrsForPair, pTestedTRs->numObjects(), groupOrder(), indRow0, indRow1, stat);
					}
#endif
					if (bCurrentSet)
					{
						if (!nTrsForPair)
						{
							if (any2RowsConvertToFirst2) {
#if !USE_CUDA
								if (bCollectInfo)
									printfRed("Rows %d,%d with Cycles(%d:%d:%d) cant be converted to rows 0,1\n",
										indRow0, indRow1, trCycles.length[0], trCycles.length[1], trCycles.length[2]);
#endif
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
void alldata::printCyclesInfo(tchar* pV1, int nv1, int indRow0, int indRow1) 
{
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
}
void alldata::printCyclesInfoSummary(TrCycles* trCycles, tchar* tr, int indRow0, int indRow1, int nrows)
{
	printfRed("Not canonical matrix. See Tr below from Cycles(%d:%d:%d), rows %d,%d\n",
		trCycles->length[0], trCycles->length[1], trCycles->length[2], indRow0, indRow1);
	printTable(" Tr", tr, 1, m_numPlayers, m_groupSize);
	printTable(" Matrix with Tr applied", m_Ktmp, nrows, m_numPlayers, m_groupSize);
}
#endif

