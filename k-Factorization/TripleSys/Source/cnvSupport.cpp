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
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);
	// get first row
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
							const int icmp = nrows < 3 ? 0 : kmProcessMatrix2p1f(trt, nrows, indRow0, indRow1);
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
				const auto* neighborsi = neighbors(indRow0);
				const auto* neighborsj = neighbors(indRow1);
				for (int k = 0; k < m_numPlayers; k++)
				{
					if (!create2P1FTr(tr, k, neighbors0, neighbors1, neighborsi, neighborsj)) {
						if (any2RowsConvertToFirst2) {
							bRet = false;
							if (m_doNotExitEarlyIfNotCanonical)
								break; // Calculate |Aut| and minimum player index to comeback for all such tr's 
							goto ret;
						}
						break;
						// ??? continue;
					}
					m_TrInd++;
#if !USE_CUDA
					if (m_cnvMode) {
						cnvPrintAuto(tr, nrows);
						continue; // print only
					}
#endif
					const int icmp = nrows < 3 ? 0 : kmProcessMatrix2p1f(tr, nrows, indRow0, indRow1);
					//TestkmProcessMatrix(nrows, 0, tr, tr, icmp);
					// save Tr if icmp not -1, and not 1; continue if it was already processed
					if (icmp < 0) {
						bRet = false;
						if (m_doNotExitEarlyIfNotCanonical)
							continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
#if PRINT_TRANSFORMED
						printTransformed(nrows, m_numPlayers, m_groupSize, tr, tr, result(), cmpGraph? m_Km:m_Ktmp, 0, 0, 0);
#endif
						goto ret;
					}
					// save Tr (if icmp not -1 or 1), continue if Tr was already processed
					//if (bSaveTestedTrs && icmp != 1 && pTestedTRs->isProcessed(tr))
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
	TrCycles trCycles;
	tchar tr[MAX_PLAYER_NUMBER];
	bool bRet = true;
	auto v1 = getV1();
	int ip1 = 0;
	const auto* neighbors0 = neighbors(0);
	tchar* p1 = NULL;
	bool bCurrentSet = false;
	bool bCBMP = !completeGraph();
	const auto any2RowsConvertToFirst2 = param(t_any2RowsConvertToFirst2);
#if 0
	static tchar a[] = {0,3,1,4,2,7,5,8,6,9};
	if (memcmp(a, result(1), sizeof(a)) == 0)
		ip1 = ip1;
#endif
	while (1)
	{
		//if ((any2RowsConvertToFirst2 && nrows != 2) || !m_pSecondRowsDB)
		if (!m_pSecondRowsDB || (param(t_autSaveTestedTrs) > 0 && nrows < m_numDaysResult))
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
		if (MEMCMP(p1, result(1), m_numPlayers) == 0) {
			bCurrentSet = true;
			cyclesFor2Rows(m_TrCyclesAll, &trCycles, neighbors(0), neighbors(1), result(0), result(1));
			if (m_createSecondRow && nrows == numDaysResult() && nrows == nrowsToUseForTrs) {
				if (param(t_rejectCycleLength)) {
					for (int i = 1; i < m_TrCyclesAll[0].ncycles; i++) {
						if (m_TrCyclesAll[0].length[0] != m_TrCyclesAll[0].length[i]) {
							bRet = false;
							goto ret;
						}
					}
				}
			}
		}
		else {
			tchar neighbors1[MAX_PLAYER_NUMBER];
			u1fSetTableRow(neighbors1, p1);
			cyclesFor2Rows(m_TrCyclesAll, &trCycles, neighbors(0), neighbors1, result(0), p1);
		}
		if (nrows == 2) {
			if (!cyclesOfTwoRowsOk(m_TrCyclesAll)) {
				if (bCurrentSet) {
					bRet = false;
					break;
				}
				ASSERT_IF(1);
			}
		}
		bool bCollectInfo = bCurrentSet && (param(t_printMatrices) & 16);
		bool bPrintInfo = bCollectInfo && (nrows == numDaysResult());
		bool bUseTestedTrs = bCurrentSet && ((param(t_autSaveTestedTrs) > 0) || bCollectInfo);
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
								if (m_playerIndex > (iRowLast + 1) * m_numPlayers - m_groupSize - 1)
									m_playerIndex = (iRowLast + 1) * m_numPlayers - m_groupSize - 1;
							}
							if (icmp == 0)
								updateGroup(trt);
							else if (icmp < 0) {
								bRet = false;
								if (m_doNotExitEarlyIfNotCanonical)
									continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
								goto ret;
							}
						}
					}
					else {
						auto* pV1 = v1;
						if (bSaveTestedTrs)
							pTestedTRs->resetGroupOrder();
						nv1 = bCBMP ? m_maxCommonVSets : ((m_groupSize == 3) ? getAllV(v1, m_maxCommonVSets, indRow0, indRow1) :
							MAX_CYCLE_SETS);
#if !USE_CUDA
						if (bPrintInfo) {
							cyclesFor2Rows(m_TrCyclesPair, &m_TrCycles, neighbors(indRow0), neighbors(indRow1), 
								result(indRow0), result(indRow1));
						}
#endif
						for (int itr = 0; itr < nv1; itr++, pV1 += m_nGroups) {
							if (collectOneCyclesSet(&trCycles, pV1, itr, indRow0, indRow1, eNoErrorCheck) <= 0)
								continue;
							bool bCycleSelected = false;
							for (int itr0 = 0; itr0 < MAX_CYCLE_SETS; itr0++)
							{
								if (m_TrCyclesAll[itr0].counter == 0)
									break;
								if (MEMCMP(m_TrCyclesAll[itr0].length, trCycles.length, MAX_CYCLES_PER_SET))
									continue;
								bCycleSelected = true;
								ctchar* pDir, * pStartOut;
								auto pIdx = InitCycleMapping(trCycles.length, trCycles.start, trCycles.ncycles, m_groupSize, &pDir, &pStartOut);
								do {
									const bool btr = createU1FTr(tr, &m_TrCyclesAll[itr0], &trCycles, pDir, pIdx, pStartOut);
									if (btr) {
										if (bCBMP) {
											if (!checkCBMPtr(tr))
												continue;
										}
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
													if (m_playerIndex > (iRowLast + 1) * m_numPlayers - m_groupSize - 1)
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
											if (m_doNotExitEarlyIfNotCanonical)
												continue; // Calculate |Aut| and/or minimum player to return
											goto ret;
										}
										// save Tr (if icmp not -1 or 1), continue if Tr was already processed
										//if (bSaveTestedTrs && icmp != 1 && pTestedTRs->isProcessed(tr))
										if (bSaveTestedTrs && pTestedTRs->isProcessed(tr))
											continue;

										if (icmp == 0)
											updateGroup(tr);
									}
								} while (ProceedToNextMapping());
								//break; // if not p1f we can't break;
							}
							if (!bCurrentSet)
								continue;
							if (!bCycleSelected) {
								if (!m_allowUndefinedCycles) {
#if !USE_CUDA
									if (bCollectInfo)
										printfRed("Cycles (%d:%d:%d) for rows %d,%d are present, but not selected\n",
											trCycles.length[0], trCycles.length[1], trCycles.length[2], indRow0, indRow1);
#endif
									bRet = false;
									if (m_doNotExitEarlyIfNotCanonical)
										continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
									goto ret;
								}
								break;
							}
						}
					}
#if !USE_CUDA
					if (bPrintInfo) {
						char stat[256];
						matrixStatOutput(stat, sizeof(stat), m_TrCyclesPair);
						printf("nTr(generated=%-3d, new=%-3d) GroupOrder(accumulated)=%-2d Cycles for rows %d,%d: %s\n",
							nTrsForPair, pTestedTRs->numObjects(), orderOfGroup(), indRow0, indRow1, stat);
					}
#endif
					if (bCurrentSet) {
						if (!nTrsForPair) {/**
							if (nrows == 7 && trCycles.length[0] == 21)//bCollectInfo)
								printfRed("nrows=%d: Rows %d,%d with Cycles(%d:%d:%d) cannot be converted to rows 0,1\n",
									nrows, indRow0, indRow1, trCycles.length[0], trCycles.length[1], trCycles.length[2]);**/
							if (any2RowsConvertToFirst2) {
#if !USE_CUDA
								if (bCollectInfo)
									printfRed("nrows=%d: Rows %d,%d with Cycles(%d:%d:%d) cannot be converted to rows 0,1\n",
										nrows, indRow0, indRow1, trCycles.length[0], trCycles.length[1], trCycles.length[2]);
#endif
								bRet = false;
								if (m_doNotExitEarlyIfNotCanonical)
									continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
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
	return bRet;
}
CC int alldata::collectOneCyclesSet(TrCycles* trc, tchar* pV1, int ind, int indRow0, int indRow1, eCheckForErrors checkErrors)
{
	int iRet = -1;
	if (!completeGraph())
		iRet = getCyclesAndPathCBMP(trc, neighbors(indRow0), neighbors(indRow1),
			result(indRow0), result(indRow1), ind, checkErrors);
	else if (m_groupSize == 2)
		iRet = getCyclesAndPathFromNeighbors(trc, neighbors(indRow0), neighbors(indRow1), NULL, NULL, checkErrors);
	else
		iRet = getCyclesAndPath3(trc, pV1, neighbors(indRow0), neighbors(indRow1), result(indRow0), result(indRow1), checkErrors);
	return iRet;
}

bool alldata::checkCBMPtr(tchar* tr) {
	tchar gtest[MAX_PLAYER_NUMBER];
	int is = 0;
	ASSERT_IF(m_groupSize >= sizeof(is) * 8);// the code below will not work if group size > 31
	for (int i = 0; i < m_groupSize; i++) {
		auto gind = m_groupSizeRemainder[tr[i]];
		gtest[i] = gind;
		is |= (1 << gind);
	}
	if (is != (1 << m_groupSize) - 1)
		return false;
	auto i = m_groupSize;
	for (; i < m_numPlayers; i += m_groupSize) {
		for (int j = 0; j < m_groupSize; j++) {
			if (m_groupSizeRemainder[tr[i + j]] != gtest[j])
				return false;
		}
	}
	return i == m_numPlayers;
}

#if !USE_CUDA
void alldata::cnvPrintAuto(ctchar* tr, int nrows)
{
	ASSERT_IF(nrows < 2);
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
void alldata::printCyclesInfoNotCanonical(TrCycles* trCycles, tchar* tr, int indRow0, int indRow1, int nrows)
{
	printfRed("Not canonical matrix. See Tr below from Cycles(%d:%d:%d), rows %d,%d\n",
		trCycles->length[0], trCycles->length[1], trCycles->length[2], indRow0, indRow1);
	printTable(" Tr", tr, 1, m_numPlayers, m_groupSize);
	printTable(" Matrix with Tr applied", m_Ktmp, nrows, m_numPlayers, m_groupSize);
}
#endif

