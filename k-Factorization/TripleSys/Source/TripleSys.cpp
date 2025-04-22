#include <iostream>
#include "TripleSys.h"
#include "Table.h"

#if !USE_CUDA && USE_BINARY_CANONIZER
#include "k-SysSupport.h"
#include "CDTools.h"
#endif

#if 0
#ifdef _MSC_VER
#include <intrin.h> /* for rdtscp and clflush */
#pragma optimize("gt",on)
#else
#include <x86intrin.h> /* for rdtscp and clflush */
#endif
#endif

using namespace std;

#if !USE_CUDA
void _printf(FILE* f, bool toScreen, const char* format, const char* pStr) {
	if (f)
		fprintf(f, format, pStr);

	if (toScreen)
		printf(format, pStr);
}

#if PrintImprovedResults
void alldata::outputResults(int iDay, const unsigned char *pResult, int cntr) const
{
	char buffer[256];
	const bool toScreen = PrintImprovedResults > 1;
	FOPEN_W(f, ImprovedResultFile, canonCalls(1) || cntr? "a" : "w", m_file);

	const unsigned char* pDayPerm = NULL;
	auto flag = true;
	if (cntr) {	
		flag = m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainMatr);
		if (flag) {
			pDayPerm = pResult + iDay * numPlayers();
			sprintf_s(buffer, "Improved Result #%d for %d days  m_groupIndex = %d:\n", cntr, iDay, m_groupIndex);
			_printf(f, toScreen, buffer);
		}

		if (m_pCheckCanon->improvedResultIsReady(t_bResultFlags::t_readyToExplainTxt))
			_printf(f, toScreen, "%s\n", m_pCheckCanon->comment());
	}
	else {
		sprintf_s(buffer, "Initial Result #%zd (%zd):\n", canonCalls(0), ((Table<char>*) m_pRes)->counter());
		_printf(f, toScreen, buffer);
	}

	if (flag)
		outMatrix(pResult, iDay, numPlayers(), m_groupSize, 0, f, false, toScreen, cntr, pDayPerm);

	FCLOSE_W(f, m_file);
}
#endif

#if CHECK_PERMUTATIONS
void alldata::outputError() const {
	extern char lastError[];
	FOPEN_W(f, ImprovedResultFile, "a", m_file);
	_printf(f, false, lastError);
	FCLOSE_W(f, m_file);
}
#endif
#endif

CC sLongLong alldata::Run(int threadNumber, eThreadStartMode iCalcMode, CStorageSet<tchar>* secondRowsDB,
	tchar* mStart0, ctchar* mfirst, int nrowsStart, sLongLong* pcnt, string* pOutResult, int iThread) {
	// Input parameters:
	int iPrintMatrices = 0;
	m_lastRowWithTestedTrs = 0;
#if !USE_CUDA
	iPrintMatrices = iThread == 0 ? param(t_printMatrices) : 0;
	int* edges = NULL;
	if (iPrintMatrices & 32)
		edges = new int[numDaysResult() * m_numPlayers * m_numPlayers]();
	const auto iTime = clock();
	const char* fHdr = getFileNameAttr(sysParam());
	auto rTime = iTime;
	auto cTime = iTime;
	const auto iCalcModeOrg = iCalcMode;
	const auto bSavingMatricesToDisk = iCalcModeOrg != eCalcSecondRow? param(t_savingMatricesToDisk) : false;
	int nMatricesMax = 0;
	int startMatrixCount = 0;
#endif
	int iDaySaved = 0;
	auto nPrecalcRows = param(t_useRowsPrecalculation);
	m_pSecondRowsDB = secondRowsDB;
	if (iCalcMode == eCalcSecondRow) {
		iCalcMode = eCalcResult;
		m_createSecondRow = 1;
		m_numDaysResult = 2;
		nrowsStart = nPrecalcRows = 0;
	}
	else
		m_createSecondRow = 0;

	const auto groupSize_2 = m_groupSize == 2;
	if (iCalcMode == eCalcResult)
		m_useRowsPrecalculation = (nPrecalcRows && nPrecalcRows <= 3 && m_groupSize <= 3 && nrowsStart <= nPrecalcRows) ? eCalculateRows : eDisabled;
	else
		m_useRowsPrecalculation = iCalcMode;
	tchar secondPlayerMax = m_numPlayers - (m_groupSize == 2 ? 1 : m_groupSize + 1 + m_numDays - m_numDaysResult);
	tchar secondPlayerInRow4First = groupSize_2 ? nPrecalcRows + 1 : m_groupSize + 2;
	tchar secondPlayerInRow4Last = MIN2(param(t_lastRowSecondPlayer), secondPlayerMax);
	if (!secondPlayerInRow4Last)
		secondPlayerInRow4Last = groupSize_2 ? m_numDaysResult : secondPlayerMax;

	m_secondPlayerInRow4 = nPrecalcRows ? secondPlayerInRow4First : 0; // used only if UseRowsPrecalculation not 0
	int nRows4 = 0;
	int nRows4Day = 0;
	const auto bPrint = iPrintMatrices != 0;
	int minRows = nrowsStart;

#if !USE_CUDA
	//extern void aq();
	//aq();
	unsigned char* bResults = NULL;

	COutGroupHandle* pAutGroup[3] = { NULL, NULL, NULL };
	TableAut Result("\n\n|Aut(M)|", m_numDays, m_numPlayers, 0, m_groupSize, true, true);

	if (bSavingMatricesToDisk) {
		string fName = format("{:0>10}.txt", threadNumber);
		if (m_improveResult) {
			const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
			bResults = new unsigned char[(m_improveResult > 1 ? 2 : 1) * lenResult];
			createFolderAndFileName(ImprovedResultFile, sysParam(), t_ImprovedResultFolder, numDaysResult(), fName);
		}

		createFolderAndFileName(ResultFile, sysParam(), t_ResultFolder, numDaysResult(), fName);

		Result.allocateBuffer(32);
		const auto pResFile = ResultFile.c_str();
		Result.setOutFileName(pResFile);
		m_pRes = &Result;

		const auto outAutGroup = param(t_outAutomorphismGroup);
		if (outAutGroup) {
			if (outAutGroup & 1) {
				pAutGroup[0] = new Generators(outAutGroup, "\nOrbits and generators of the Aut(M) acting on elements", m_numPlayers, m_groupSize);
				pAutGroup[0]->setOutFileName(pResFile, false);
			}
			if (outAutGroup & 2) {
				pAutGroup[1] = new COutGroupHandle(outAutGroup, "\nAut(M) acting on elements", m_numPlayers, m_groupSize);
				pAutGroup[1]->setOutFileName(pResFile, false);
			}

			if (outAutGroup & 12) {
				pAutGroup[2] = new RowGenerators(outAutGroup, "", numDaysResult());
				pAutGroup[2]->setOutFileName(pResFile, false);
			}
		}
	}

#endif
	CUDA_PRINTF("*** mStart0 = %p\n", mStart0);
#if 0 // make 27 player matrix on the basis of first 5 rows
	if (m_numPlayers != 27 || iDay != 13) {
		printfRed("Can't create 27x13 matrix. Exit"); myExit(1);
	}
	for (int i = 4; i < 13; i++) {
		tchar* r = result(i);
		for (int j = 0; j < 27; j += 3, r += 3) {
			r[0] = j / 3; r[1] = ((j / 3 + i + 5) % 9) + 9; r[2] = ((j / 3 + (i - 4) * 2) % 9) + 18;
		}
	}
	printTable("circular matrix", result(), iDay, m_numPlayers, m_groupSize);
#endif
	if (mStart0)
	{
		iDay = nrowsStart;
		memcpy(result(), mStart0, nrowsStart * m_numPlayers);
		linksFromMatrix(links(), mStart0, nrowsStart);
	}

	if (iDay > numDaysResult() && !(iPrintMatrices & 16))
		iDay = numDaysResult(); // warning?

	memset(m_rowTime, 0, m_numDays * sizeof(m_rowTime[0]));
	for (int i = 0; i < iDay; i++)
		u1fSetTableRow(neighbors(i), result(i));

	if (iDay >= 2 && (m_use2RowsCanonization || param(t_u1f)))
	{
		if (groupSize_2)		// need to be implemented for 3?
			p1fCheckStartMatrix(iDay);
	}
#if 1 // preset automorphism groups
	const auto iCalc = m_useRowsPrecalculation;
	m_useRowsPrecalculation = eCalcResult;
	if (!param(t_autGroupNumb)) {
		if (param(t_useAutForPrecRows) > 1 && iDay >= param(t_useAutForPrecRows)) {
			auto* pGroupInfo = groupInfo(param(t_useAutForPrecRows));
			if (pGroupInfo) {
				cnvCheckNew(0, param(t_useAutForPrecRows), false); // create initial set of tr for first i rows
				pGroupInfo->copyIndex(*this);
				resetGroupOrder();
			}
		}
	}
	else {
		for (int i = firstGroupIdx(); i <= lastGroupIdx(); i++) {
			auto* pGroupInfo = groupInfo(i);
			if (!pGroupInfo)
				break;
			if (iDay < i)
				break;
			cnvCheckNew(0, i, false); // create initial set of tr for first i rows
			pGroupInfo->copyIndex(*this);
			resetGroupOrder();
		}
	}
	m_useRowsPrecalculation = iCalc;
#endif
	m_playerIndex = 0;
#if !USE_CUDA
	sLongLong nMCreated = 0;
	auto mTime = clock();
#if 0
	if (!mStart0 && iDay > 0)
	{
		testCanonizatorSpeed();
		exit(0);
	}
#endif
#endif
#if 0
	const auto checkFlag = iDay == m_numDays || iDay == numDaysResult();
	if (checkFlag)
	{
		const auto retVal = improveMatrix(m_improveResult, NULL, 0/*, bResults, lenResult()*/);
#if 0
		if (retVal) {
			noMoreResults = true;
			goto noResult;
		} 
#endif
		int ir = checkCurrentResult(bPrint);
		//if (improveMatrix(m_improveResult, NULL, 0/*, bResults, lenResult()*/))
		//	return -1;
		// special case if input is full result matrix or UseTwoLastDaysG2 mode
		switch (ir) {
		case -1:
		case  1: noMoreResults = true; goto noResult;
		default: break;
		}
	}
#endif
	void* pIS_Canonizer = NULL;
#if !USE_CUDA && USE_BINARY_CANONIZER
	if (m_ppBinMatrStorage)
		pIS_Canonizer = createCanonizer(numPlayers(), m_groupSize);
#endif

	bPrevResult = false;
	const auto p1f_counter = param(t_p1f_counter);

#if 1 && !USE_CUDA 
	if (testGroupOrderEachSubmatrix(iPrintMatrices, iCalcModeOrg)) {
		noMoreResults = true;
		goto noResult;
	}
#endif

	if (iCalcMode == eCalculateMatrices)
		m_pRowUsage->init(iThread, param(t_numThreads));

	else if (iDay > 0) {
		if (m_useRowsPrecalculation == eCalculateRows)
			m_pRowStorage->initPlayerMask(mfirst);
		m_lastRowWithTestedTrs = 0;
		setArraysForLastRow(iDay);
		//printTable("p1f", neighbors(), iDay, m_numPlayers, 0);
		iDay--;
		//updateIndexPlayerMinMax();
		//testRightNeighbor(iDay + 1);
		goto checkCurrentMatrix;
	}

	while (nLoops < LoopsMax)
	{
#if !USE_CUDA
		mTime = clock();
#endif
		CUDA_PRINTF(" *** nLoops = %lld\n", nLoops);
		while (iDay < numDaysResult() || bPrevResult)
		{
			if (nPrecalcRows && m_useRowsPrecalculation == eCalculateMatrices) {

			ProcessPrecalculatedRow:
#if 0
				static int a, b, c;
				if (m_playerIndex)
					a++;
				if (checkCanonicity())
					c++;
				b++;
				if ((b%1000000)==0)
					printf("%d:%d:%d ", a, c, b);
#endif
				if (iDaySaved){
					if (m_pRowUsage->getMatrix2(result(), neighbors(), numDaysResult(), iDaySaved)) {
#if !USE_CUDA
						cTime = clock();
						for (int i = nPrecalcRows + 3; i < numDaysResult(); i++)
							m_rowTime[i] = cTime - iTime;
#endif
						iDay = numDaysResult() - 1;
						goto checkCurrentMatrix;
					}
					iDay = iDaySaved;
					iDaySaved = 0;
					bPrevResult = false;
					m_playerIndex = 0;
				}
				if (bPrevResult) {
					iDay--;
					bPrevResult = false;
				}
				if (iDay >= nPrecalcRows) {
					int ipx = 0;
					if (m_playerIndex)
					{
						iDay = m_playerIndex / m_numPlayers;
						ipx = m_playerIndex % m_numPlayers;
						m_playerIndex = 0;
						if (iDay < nPrecalcRows || iDay >= numDaysResult())
						{
							//ASSERT(1);
							noMoreResults = true;
							goto noResult;
						}
					}

					if (m_lastRowWithTestedTrs >= iDay)
						m_lastRowWithTestedTrs = iDay - 1;
					const auto retVal = m_pRowUsage->getRow(iDay, ipx);
					switch (retVal) {
					case 2:
						iDaySaved = iDay;
						break;
					case 1:
						///m_pRowUsage->getMatrix(result(), neighbors(), iDay + 1);
						//printTable("tbl", result(), iDay + 1, m_numPlayers, groupSize());

						m_playerIndex = 0;
#if !USE_CUDA
						if (bPrint && iDay < nPrecalcRows + 3) {
							cTime = clock();
							m_rowTime[iDay] = cTime - iTime;
						}
#endif
						if (++iDay < numDaysResult() && !checkCanonicity()) {

#if 0   // Temporary
							if (!p1f_counter || ((++m_p1f_counter) % p1f_counter))
#endif
#if !USE_CUDA
								if (!bPrint || cTime - rTime < ReportInterval)
#endif
									goto ProcessPrecalculatedRow;
						}

						m_pRowUsage->getMatrix(result(), neighbors(), iDay);
						if (m_bCheckLinkV && iDay < m_numDays) {
							linksFromMatrix(links(), result(), iDay);
							if (!checkLinks(links(), iDay)) {
								iDay--;
								goto ProcessPrecalculatedRow;
							}
						}
						iDay--;
#if !USE_CUDA
						cTime = clock();
						for (int i = nPrecalcRows + 3; i <= iDay; i++)
							m_rowTime[i] = cTime - iTime;
#endif
						goto checkCurrentMatrix;
					case -1: // reported if requested row not in solution (can happen only if LastRowSecondPlayer is incorrect)
					{
						if (secondPlayerInRow4Last < secondPlayerMax) {
							m_pRowUsage->getMatrix(result(), neighbors(), iDay);
							linksFromMatrix(links(), result(), iDay);
							if (cnvCheckNew(0, iDay, false) && (m_groupSize != 3 || checkLinks(links(), iDay))) {
								// matrix ok, error cannot be ignored
#if !USE_CUDA
								printfRed("*** Failing getRow(%d) = %d ***\n", iDay, retVal);
								printTable("Current result", result(), iDay, m_numPlayers, groupSize());
								exit(1);
#endif					
							}
						}
						// continue to case 0:
					}
					case 0:
						iDay--;
						break;
					}
					goto ProcessPrecalculatedRow;
				}
				if (iCalcMode == eCalcResult) {
					m_useRowsPrecalculation = eCalculateRows;
					m_pRowStorage->init();
					iDay++;
					bPrevResult = true;
					continue;
				}
				noMoreResults = true;
				goto noResult;
			}

			CUDA_PRINTF("   *** iDay = %d  bPrevResult = %d\n", iDay, bPrevResult);
			if (iDay < 0)
			{
				noMoreResults = true;
				goto noResult;
			}
			if (bPrevResult)
			{
				if (iDay < minRows)
				{
					noMoreResults = true;
					goto noResult;
				}
				if (iDay == 2 && groupSize_2 && m_use2RowsCanonization)
				{
					noMoreResults = true;
					goto noResult;
				}
				if (!initPrevDay())
					continue;
			}
			else if (!initCurrentDay())
				continue;

		ProcessOneDay:
			if (m_lastRowWithTestedTrs >= iDay)
				m_lastRowWithTestedTrs = iDay - 1;
			if (iDay < minRows)
			{
				noMoreResults = true;
				goto noResult;
			}
			// temporary
#if 0 //!USE_CUDA
			if (!iPlayer && iDay){
				tchar* l = links(0);
				for (int i = 0;i < m_numPlayers * m_numPlayers; i++) {
					if (l[i] == iDay) {
						printfRed("Bad link(%d) at %dx%d\n", l[i], i / m_numPlayers, i % m_numPlayers);
						//exit(1);
					}
				}
			}
#endif
			if (!processOneDay())
			{
				if (nPrecalcRows && nPrecalcRows == iDay && m_useRowsPrecalculation == eCalculateRows) {
					m_secondPlayerInRow4++;
					if (!nRows4Day && groupSize_2)
					{
						m_secondPlayerInRow4 = secondPlayerInRow4First;
						bPrevResult = true;
						if (nRows4) {
							nRows4 = 0;
							if (iCalcMode == eCalculateRows) {
								nLoops = nRows4;
								noMoreResults = true;
								goto noResult;
							}
							m_pRowUsage->init();
							m_pRowStorage->reset();
						}
						continue;
						//noMoreResults = true;
						//goto noResult;
					}
					nRows4Day = 0;
					if (m_secondPlayerInRow4 <= secondPlayerInRow4Last)
						continue;
					m_secondPlayerInRow4 = secondPlayerInRow4First;
					if (nRows4) {
						iDay = nPrecalcRows;
						if (bPrint) {
							printf("Total number of precalculated row solutions = %5d\n", nRows4);
							m_lastRowWithTestedTrs = 0;
#if !USE_CUDA
							if (iPrintMatrices & 32) {
								for (int j = 0; j < m_numPlayers; j++) {
									bool bnl = true;
									for (int k = j + 1; k < m_numPlayers; k++) {
										if (links(j)[k] == unset) {
											if (bnl) {
												bnl = false;
												printf("\n\n     ");
												for (int i = nPrecalcRows; i < numDaysResult(); i++)
													printf(" R%-2d ", i);
											}
											printf("\n%2d:%-2d", j, k);
											for (int i = nPrecalcRows; i < numDaysResult(); i++)
												printf(" %-4d", edges[(j * m_numPlayers + k) * numDaysResult() + i]);
										}
									}
								}
								printf("\n");
								delete[] edges;
								exit(0);
							}
#endif
						}
						m_useRowsPrecalculation = eCalculateMatrices;
						m_playerIndex = 0;

						if (!m_pRowStorage->initCompatibilityMasks()) {
#if !USE_CUDA
							//printfRed("*** Unexpected error returned by initCompatibilityMask()\n");
							//ASSERT(1);
#endif
							noMoreResults = true;
							goto noResult;
						}

						if (iCalcMode == eCalculateRows) {
							nLoops = nRows4;
							noMoreResults = true;
							goto noResult;
						}
						m_pRowUsage->init();
						m_pRowStorage->reset();
						nRows4 = 0;
						continue;
					}
				}
				bPrevResult = true;
				continue;
			}
			/*
			static tchar a[] = { 0, 4, 8,  1, 5, 6,  2, 9,13,  3,10,14,  7,11,12 };
			//static tchar a[] = { 0 , 4 , 8,   1,  5,  9,   2, 10, 12,   3,  7, 14,   6, 11, 13 };
			if (memcmp(a, result(1), sizeof(a)) == 0)
				iDay = iDay;*/
			if (m_useRowsPrecalculation == eCalculateMatrices) {
				ASSERT(1);
				noMoreResults = true;
				goto noResult;
			}
			memcpy(result(iDay), tmpPlayers, m_numPlayers);

		checkCurrentMatrix:

			if (m_lastRowWithTestedTrs >= iDay)
				m_lastRowWithTestedTrs = iDay - 1;

#if ReportPeriodically && !USE_CUDA
			cTime = clock();
			m_rowTime[iDay] = cTime - iTime;
			if (maxDays < iDay || cTime - rTime > ReportInterval)
			{
				if (bPrint)
				{
					printf("Thread %d: Current data for %s-matrix %zd: rows=%d, build time=%d, time since start=%d\n",
						threadNumber, fHdr, nLoops + 1, iDay + 1, cTime - rTime, cTime - iTime);
					printResultWithHistory("Current matrix", iDay + 1);
				}
#if ReportCheckLinksData
				if (bFirstThread)
					reportCheckLinksData();
#endif
				StatReportPeriodically(ResetStat, "Stat current. iDay", iDay, iThread == 0);
				rTime = cTime;
				maxDays = iDay;
			}
#endif
			iDay++;
#if 0
			Stat("<10", 1, iDay < 10);
			Stat("=10", 2, iDay == 10);
			Stat("=11", 3, iDay == 11);
			Stat("=12", 4, iDay == 12);
			Stat(">12", 5, iDay > 12);
#endif
#if !USE_CUDA
			if (bPrint && (int)((++nMCreated) % 10000000) == 0)
				printf(" %zdM calls to checkCurrentResult\n", nMCreated / 1000000);
#endif

#if 1
			switch (checkCurrentResult(iPrintMatrices, pIS_Canonizer)) {
			case -1:  
				if (nPrecalcRows && m_useRowsPrecalculation == eCalculateMatrices) {
					goto ProcessPrecalculatedRow;
				}
				else {
					goBack();
					goto ProcessOneDay;
				}
			case  1: noMoreResults = true; goto noResult;
			default: break;
			}
#endif
			if (nPrecalcRows && m_useRowsPrecalculation == eCalculateRows) {
				if (iDay == nPrecalcRows + 1) {
					if (!nRows4++ && iDay >= 2) {
						cyclesFor2Rows(result(1));
						memcpy(m_TrCyclesFirst2Rows, m_TrCyclesAll, sizeof(m_TrCyclesFirst2Rows));
					}

					nRows4Day++;
#if 0
					printf("%6d:", nRows4);
					printTable("", neighbors(3), 1, m_numPlayers, 2);
					printTable("r", result(3), 1, m_numPlayers, 2);
#endif
#if 0
					bool bP1F = p1fCheck3(result(0), result(nPrecalcRows), neighbors(0), neighbors(nPrecalcRows));
					if (!bP1F)
						printf("not p1f\n");
					bP1F = p1fCheck3(result(2), result(nPrecalcRows), neighbors(2), neighbors(nPrecalcRows));
					if (!bP1F)
						printf("not p1f\n");
#endif
					const auto retVal = m_pRowStorage->addRow(result(nPrecalcRows), neighbors(nPrecalcRows));

#if !USE_CUDA
					if (iPrintMatrices & 32) {
						tchar* c = result(nPrecalcRows);
						for (int i = 0;i < m_numPlayers;i += m_groupSize) {
							int j = c[1] - 1 + (c[i] * m_numPlayers + c[i + 1]) * numDaysResult();
							edges[j]++;
							if (m_groupSize > 2) {
								j = c[1] - 1 + (c[i + 1] * m_numPlayers + c[i + 2]) * numDaysResult();
								edges[j]++;
							}
						}
					}
#endif
					if (!retVal){/**
						if (param(t_MultiThreading))
						{
							noMoreResults = true;
							goto noResult;
						}**/
						m_playerIndex = nPrecalcRows * m_numPlayers;
						goBack();
						m_secondPlayerInRow4 = secondPlayerInRow4First;
						m_playerIndex = 0;
						nRows4 = nRows4Day = 0;
					}
					bPrevResult = true;
					continue;
				}
			}
		}
		ASSERT(iDay < numDaysResult());
		if (orderOfGroup() >= param(t_resultGroupOrderMin))
		{
#if !USE_CUDA && USE_BINARY_CANONIZER && 0
			if (pIS_Canonizer) {
				const auto* pCanonBinaryMatr = runCanonizer(pIS_Canonizer, result(0), m_groupSize, iDay);
				if (m_ppBinMatrStorage[iDay]->updateRepo(pCanonBinaryMatr) < 0) {
					// The binary matrix has already been encountered
					printfRed("Error in canonizator\n");
					myExit(1);
				}
			}
#endif
			nLoops++;
			m_finalKMindex++;
#if !USE_CUDA
			if (m_createSecondRow) {
				if (groupSize_2)
					m_pSecondRowsDB->addObject(result(1));
				if (bPrint)
					printTableColor("Second Row", result(1), 1, m_numPlayers, m_groupSize);
				goto cont1;
			}

			if (bPrint)
			{
				rTime = cTime = clock();
				setConsoleOutputMode();
				//report result
#if !DEBUG_NextPermut
				printf("%5zd: %s-Matrix, build time=%d, time since start=%d\n", nLoops, fHdr, cTime - mTime, cTime - iTime);
#else
				extern int matr_cntr;
				printf("%5zd: %s-Matrix, matr_cntr = %d\n", nLoops, fHdr, matr_cntr);
#endif
			}

			if (bPrint || bSavingMatricesToDisk)
			{
				char stat[128];
				bool needOutput = false;
				matrixStat(neighbors(), iDay, &needOutput);
				if (needOutput)
					m_matrixDB.addMatrix(orderOfGroup(), matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll));
				else
					stat[0] = '\0';

				Result.setInfo(stat);
				Result.setGroupOrder(orderOfGroup());
#if 0			// record result and print on screen (if bPrint==true)
				Result.printTable(result(), true, bPrint, numDaysResult());
#else			// record result without print on screen
				Result.printTable(result(), true, false, numDaysResult());
				if (bPrint) {
					printf("%5zd: " AUT "%d, % s\n", nLoops, orderOfGroup(), stat);
					// print on screen result with highlighted differences from prev result
					printResultWithHistory("", iDay);
					if (iPrintMatrices & 2)
						printPermutationMatrices(2);
				}

				if (orderOfGroup() > 1) {
					for (int i = 0; i < countof(pAutGroup); i++) {
						if (pAutGroup[i])
							pAutGroup[i]->makeGroupOutput(this, bPrint);
					}
				}
#endif
				//cnvCheckNew(2, iDay);
			}

			//checkCommonValues();

			//Result.printTable(neighbors(), true, ResultFile, bPrint, numDaysResult());
			//reportCheckLinksData();
			//printTable("p1f", neighbors(), iDay, m_numPlayers, 2);
cont1:
			StatReportAfterEachResult(ResetStat, "Stat for matrix result. iDay", iDay, bPrint); // see stat.h to activate
			if (pcnt) {
				*pcnt = -m_finalKMindex - 1;
				*(pcnt + 1) = nLoops;
			}

			const auto pMatrTest = sysParam()->strVal[t_matrTest];
			if (pMatrTest) {
				auto testDescr = string(*pMatrTest);
				string testParam;
				const auto pos = testDescr.find(':');
				if (pos != string::npos) {
					testParam = testDescr.substr(pos + 1);
					testDescr = testDescr.substr(0, pos);
				}

				if (testDescr == "FindIsomorphicBaseElements") {
					FindIsomorphicBaseElements(testParam);
				}
				else {
					printfRed("Test \"%s\" is not implemented\n", testDescr.c_str());
					myExit(1);
				}
			}
#endif
		}
		// 60 sec
#define LOOP_LENGTH	0//100 //10000
#if LOOP_LENGTH
		for (int i = 0; i < LOOP_LENGTH; i++)
			cnvCheckNew(0, iDay);
		goto noResult;
#endif
		if (iDay <= minRows)
			break;
		bPrevResult = true;
	}

noResult:

#if !USE_CUDA
#if USE_BINARY_CANONIZER
	releaseCanonizer(pIS_Canonizer);
#endif
	if (pcnt)
	{
		*pcnt = m_finalKMindex;
		*(pcnt + 1) = nLoops;
	}

	else {
		if (param(t_MultiThreading) <= 1) {
			std::string str;
			if (m_createSecondRow) {
				str = format("{} row(s) calculated for '2nd rows DB' for matrices ({},{},{}), Execution time = {} ms\n",
					m_finalKMindex, m_numPlayers, numDaysResult(), groupSize(), clock() - iTime);
			}
			else {
				str = format("\nThread {}: {} non-isomorphic matrices ({},{},{}) created\n",
					threadNumber, m_finalKMindex, m_numPlayers, numDaysResult(), m_groupSize);

				str += format("Thread execution time = {} ms\n", clock() - iTime);
			}
			printf(str.c_str());
			if (pOutResult)
				*pOutResult += str;
		}
	}
	StatReportAfterThreadEnd(ResetStat, "Thread ended, processed", (int)nLoops, bFirstThread); // see stat.h to enable
	delete[] bResults;

	for (int i = 0; i < countof(pAutGroup); i++)
		delete pAutGroup[i];

#endif
#if COUNT_GET_ROW_CALLS && !USE_CUDA
	extern ll cntr;
	static int nMatr;
	static ll cntrTotal;
	printf("Matr# %4d: ********** cntr = %lld  = %lld\n", ++nMatr, cntr, cntrTotal += cntr);
	cntr = 0;
#endif
	return nLoops;
}
