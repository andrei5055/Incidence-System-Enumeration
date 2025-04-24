#include "TripleSys.h"

CC bool alldata::cnvCheck2P1F(int nrows, int nrowsToUseForTrs)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;

	const bool bUseTestedTrs = param(t_autSaveTestedTrs) > 0;
	const auto* neighbors0 = neighbors(0);
	const auto* neighbors1 = neighbors(1);
	// get first row
	const auto cbmpGraph = !completeGraph();
	for (int iRowLast = 1; iRowLast < nrowsToUseForTrs; iRowLast++) {
		const bool bSaveTestedTrs = bUseTestedTrs && (iRowLast < m_numDaysResult - 1);
		for (int indRow = 0; indRow < iRowLast; indRow++)
		{
			// get second row
			for (int iRowSwap = 0; iRowSwap < 2; iRowSwap++)
			{
				const int indRow1 = iRowSwap ? indRow : iRowLast;
				const int indRow0 = iRowSwap ? iRowLast : indRow;
				auto* pTestedTRs = bUseTestedTrs ? m_pTrRepo->getTrSet(indRow0, indRow1) : NULL;
				if (bUseTestedTrs) {
					if (iRowLast <= m_lastRowWithTestedTrs) {
						const int nTrs = pTestedTRs->numObjects();
						for (int itr = 0; itr < nTrs; itr++) {
							tchar* trt = pTestedTRs->getObject(itr);

							const int icmp = nrows < 3 ? 0 : cbmpGraph ? kmProcessMatrix(result(), trt, nrows) : kmProcessMatrix2p1f(trt, nrows, indRow0, indRow1);
							//const int icmp = nrows < 3 ? 0 : kmProcessMatrix2p1f(trt, nrows, indRow0, indRow1);
							m_TrInd++;
							if (icmp < 0)
							{
								bRet = false;
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
				const auto* neighborsi = neighbors(indRow0);
				const auto* neighborsj = neighbors(indRow1);
				for (int k = 0; k < m_numPlayers; k++)
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
					const int icmp = nrows < 3 ? 0 : cbmpGraph ? kmProcessMatrix(result(), tr, nrows) : kmProcessMatrix2p1f(tr, nrows, indRow0, indRow1);
					//TestkmProcessMatrix(nrows, 0, tr, tr, icmp);
					// save Tr if icmp not -1, and not 1; continue if it was already processed
					if (icmp < 0) {
#if PRINT_TRANSFORMED
						printTransformed(nrows, m_numPlayers, m_groupSize, tr, tr, result(), cmpGraph? m_Km:m_Ktmp, 0, 0, 0);
#endif
						bRet = false;
						goto ret;
					}
					if (bSaveTestedTrs && pTestedTRs->isProcessed(tr))
						continue;
					if (icmp == 0) {
						updateGroup(tr);
					}
				}
			}
		}
		if (bSaveTestedTrs && m_lastRowWithTestedTrs < iRowLast)
			m_lastRowWithTestedTrs = iRowLast;
	}
ret:
	return bRet;
}
CC bool alldata::cnvCheck3U1F(int nrows, int nrowsToUseForTrs)
{
	if (nrows < 2)
		return true;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	auto v1 = getV1();
	int ip1 = 0;
	const auto* neighbors0 = neighbors(0);
	tchar* p1 = NULL;
	bool bCurrentSet = false;
	const int maxv1 = MAX_3PF_SETS;
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);
	setAllowUndefinedCycles(nrows);
	/*
	static tchar a[] = { 0, 4, 8,  1, 5, 6,  2, 9,13,  3,10,14,  7,11,12 };
	//static tchar a[] = { 0 , 4 , 8,   1,  5,  9,   2, 10, 12,   3,  7, 14,   6, 11, 13 };
	if (memcmp(a, result(1), sizeof(a)) == 0)
		ip1 = ip1;*/
	while (1)
	{
		//if ((any2RowsConvertToFirst2 && nrows != 2) || !m_pSecondRowsDB)
		if (!m_pSecondRowsDB)
			p1 = result(1);
		else if (m_pSecondRowsDB->numObjects() > ip1)
			p1 = m_pSecondRowsDB->getObject(ip1);
		else if (m_createSecondRow) // do not merge with "first if" above
			p1 = result(1);
		else {
			bRet = false;
			break;
		}
		ip1++;
		memset(&m_TrCycles, 0, sizeof(m_TrCycles));
		if (MEMCMP(p1, result(1), m_numPlayers) == 0)
			bCurrentSet = true;
		cyclesFor2Rows(p1); // result is in m_TrCyclesAll 
		if (nrows == 2) {
			if (!cyclesOfTwoRowsOk(m_TrCyclesAll)) {
				if (bCurrentSet) {
					bRet = false;
					break;
				}
				ASSERT(1);
			}
		}
		bool bCollectInfo = bCurrentSet && (param(t_printMatrices) & 16);
		bool bPrintInfo = bCollectInfo && (nrows == numDaysResult());
		bool bUseTestedTrs = bCurrentSet && ((param(t_autSaveTestedTrs) > 0) || bCollectInfo);
		TrCycles trCycles;
		// get first row
		for (int iRowLast = 1; iRowLast < nrowsToUseForTrs; iRowLast++) {
			bool bSaveTestedTrs = bUseTestedTrs && ((iRowLast < m_numDaysResult - 1) || bCollectInfo);
			for (int indRow = 0; indRow < iRowLast; indRow++)
			{
				// get second row
				for (int iRowSwap = 0; iRowSwap < 2; iRowSwap++)
				{
					const int indRow1 = iRowSwap ? indRow : iRowLast;
					const int indRow0 = iRowSwap ? iRowLast : indRow;
					int nv1 = 0;
					int nTrsForPair = 0;
					auto* pTestedTRs = bUseTestedTrs ? m_pTrRepo->getTrSet(indRow1, indRow0) : NULL;
					if (bUseTestedTrs && (iRowLast <= m_lastRowWithTestedTrs && !bCollectInfo)) {
						const auto nTrs = pTestedTRs->numObjects();
						for (int itr = 0; itr < nTrs; itr++) {
#if 1
							const auto * trt = pTestedTRs->getObject(itr);
#else    
							// Take tr in the order they were added to the database
							const auto* trt = pTestedTRs->CStorage<tchar>::getObject(itr);
#endif
							m_TrInd++;
							nTrsForPair++;
							int icmp = kmProcessMatrix(result(), trt, nrows);
							if (!bCurrentSet && icmp != -1) {
								icmp = -1;
								m_playerIndex = (iRowLast + 1) * m_numPlayers - m_groupSize - 1;
							}
							if (icmp == 0)
								updateGroup(trt);
							else if (icmp < 0) {
								bRet = false;
								goto ret;
							}
						}
					}
					else {
						//if (nv1 > 53)
						//	indRow1 = indRow1
						auto* pV1 = v1;
						if (bSaveTestedTrs)
							pTestedTRs->resetGroupOrder();
						nv1 = getAllV(v1, maxv1, indRow0, indRow1);
#if !USE_CUDA
						if (bPrintInfo) {
							collectCyclesInfo(v1, nv1, indRow0, indRow1);
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
										int icmp;
										if (1) {//bCurrentSet) {

											if (bPrintInfo) {
												tchar ts[MAX_PLAYER_NUMBER];
												icmp = kmProcessMatrix(result(), tr, nrows, 0, ts);
#if !USE_CUDA
												if (!icmp)
													printTable("Aut", ts, 1, nrows, 1);
#endif
											}
											else {
												icmp = kmProcessMatrix(result(), tr, nrows);
												if (!bCurrentSet && icmp != -1) {
													icmp = -1;
													m_playerIndex = (iRowLast + 1) * m_numPlayers - m_groupSize - 1;
												}
											}
										}
										if (icmp == -1) {
											//if (nrows>2)
											//printTransformed(nrows, m_numPlayers, m_groupSize, tr, tr, result(), m_Km, nrows, nLoops, m_finalKMindex);
#if !USE_CUDA
											if (bCollectInfo) {
												printCyclesInfoNotCanonical(&trCycles, tr, indRow0, indRow1, nrows);
											}
#endif
											bRet = false;
											goto ret;
										}
										// save Tr if icmp is not -1; continue if it was already processed
										if (bSaveTestedTrs && pTestedTRs->isProcessed(tr))
											continue;

										if (icmp == 0)
											updateGroup(tr);
									}
								} while (ProceedToNextMapping());
								break;
							}
							if (!bCycleSelected) {
								if (!allowUndefinedCycles()) {
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
						printf("nTr(generated=%-3d, new=%-3d) GroupOrder(accumulated)=%-2d Cycles for rows %d,%d: %s\n",
							nTrsForPair, pTestedTRs->numObjects(), orderOfGroup(), indRow0, indRow1, stat);
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
			if (bSaveTestedTrs && m_lastRowWithTestedTrs < iRowLast)
				m_lastRowWithTestedTrs = iRowLast;
		}
		if (bCurrentSet)
			break;
	}
ret:
	if (bRet)
	{
		//printf("Trs total=%d\n", m_TrInd); Andrei 
		if (m_createSecondRow && nrows == numDaysResult() && p1 == result(1) && nrows == nrowsToUseForTrs)
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
void alldata::collectCyclesInfo(tchar* pV1, int nv1, int indRow0, int indRow1) 
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
void alldata::printCyclesInfoNotCanonical(TrCycles* trCycles, tchar* tr, int indRow0, int indRow1, int nrows)
{
	printfRed("Not canonical matrix. See Tr below from Cycles(%d:%d:%d), rows %d,%d\n",
		trCycles->length[0], trCycles->length[1], trCycles->length[2], indRow0, indRow1);
	printTable(" Tr", tr, 1, m_numPlayers, m_groupSize);
	printTable(" Matrix with Tr applied", m_Ktmp, nrows, m_numPlayers, m_groupSize);
}
#endif

