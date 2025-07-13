#include <iostream>
#include "TripleSys.h"
#include "kOrbits.h"

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
	ctchar* mStart0, ctchar* mfirst, int nrowsStart, sLongLong* pcnt, string* pOutResult, int iThread) {
	// Input parameters:
	memset(m_rowTime, 0, m_numDaysResult * sizeof(m_rowTime[0]));
	m_allRowPairsSameCycles = param(t_allowUndefinedCycles) == 0 || param(t_any2RowsConvertToFirst2) > 0;
	int m_printMatrices = 0;
	m_lastRowWithTestedTrs = 0;
	m_threadNumber = threadNumber;
	const auto p1f_counter = param(t_p1f_counter);
	void* pIS_Canonizer = NULL;
#if !USE_CUDA && USE_BINARY_CANONIZER
	if (m_ppBinMatrStorage)
		pIS_Canonizer = createCanonizer(numPlayers(), m_groupSize);
#endif
	//if (iThread == 0 && iCalcMode == eCalculateMatrices)
	//return 0;
#if !USE_CUDA
	///m_printMatrices = (iThread == 1 || param(t_numThreads) < 2) ? param(t_printMatrices) : 0;
	m_printMatrices = (iThread == 0 || param(t_numThreads) < 2) ? param(t_printMatrices) : 0;
	int* edges = NULL;
	if (m_printMatrices & 32)
		edges = new int[numDaysResult() * m_numPlayers * m_numPlayers]();
	m_iTime = clock();
	m_fHdr = getFileNameAttr(sysParam());
	m_rTime = m_iTime;
	m_cTime = m_iTime;
	const auto iCalcModeOrg = iCalcMode;
	const auto bSavingMatricesToDisk = (iCalcModeOrg != eCalcSecondRow && iCalcModeOrg != eCalculateRows) ?
		param(t_savingMatricesToDisk) : false;
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
	const bool bCBMP = !completeGraph();
	const auto groupSize_2 = m_groupSize == 2;
	if (iCalcMode == eCalcResult)
		m_useRowsPrecalculation = (nPrecalcRows && nPrecalcRows <= 3 && m_groupSize <= 3 && nrowsStart <= nPrecalcRows) ? eCalculateRows : eDisabled;
	else
		m_useRowsPrecalculation = iCalcMode;

	tchar secondPlayerMax = m_numPlayers - (m_groupSize == 2 ? 1 : m_groupSize + 1 + m_numDays - m_numDaysResult);
	tchar m_secondPlayerInRow4First = 0;
	if (bCBMP) {
		secondPlayerMax = m_numPlayers - 1 - (m_groupSize - 2) * m_groupSize;
	}
	tchar secondPlayerInRow4Last = MIN2(param(t_lastRowSecondPlayer), secondPlayerMax);
	if (!secondPlayerInRow4Last)
		secondPlayerInRow4Last = groupSize_2 ? (bCBMP ? m_numDaysResult * 2 - 1 : m_numDaysResult) : secondPlayerMax;
	m_secondPlayerInRow4 = 0;

	int nRows4 = 0;
	int nRows4Day = 0;
	m_bPrint = m_printMatrices != 0;
	int minRows = nrowsStart;

	const auto semiSymGraph = !m_createSecondRow && numDaysResult() > 2 && param(t_semiSymmetricGraphs) == 1;
	const auto minGroupSize = semiSymGraph ? m_numDaysResult * m_numPlayers / 2 : 0;
	const auto outAutGroup = param(t_outAutomorphismGroup);
	IOutGroupHandle* pAutGroup[4] = { NULL };
	if (semiSymGraph || outAutGroup) {
		if (bSavingMatricesToDisk) {
			if (outAutGroup & 1)
				pAutGroup[0] = new Generators<tchar>(outAutGroup, "\nOrbits and generators of the Aut(M) acting on elements", m_numPlayers);

			if (outAutGroup & 2)
				pAutGroup[1] = new COutGroupHandle<tchar>(outAutGroup, "\nAut(M) acting on elements", m_numPlayers);
		}

		if (semiSymGraph || bSavingMatricesToDisk && (outAutGroup & 12))
			pAutGroup[2] = new RowGenerators(outAutGroup, numDaysResult());

		if (semiSymGraph || bSavingMatricesToDisk && (outAutGroup & 48))
			pAutGroup[3] = new CKOrbits(outAutGroup, m_numPlayers, m_groupSize, numDaysResult());
	}
#if !USE_CUDA
	//extern void aq();
	//aq();
	unsigned char* bResults = NULL;

	TableAut Result(MATR_ATTR, m_numDays, m_numPlayers, 0, m_groupSize, true, true);

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

		if (outAutGroup) {
			for (int i = 0; i < countof(pAutGroup); i++) {
				if (pAutGroup[i])
					pAutGroup[i]->setOutFileName(pResFile, false);
			}
		}
	}
#endif
	CUDA_PRINTF("*** mStart0 = %p\n", mStart0);
#if 0 // Generate KC matrix
#define N m_numPlayers
#define G m_groupSize
#define M (N/G)
	if (m_numPlayers != N || m_numDaysResult != M || m_groupSize != G) {
		printfRed("Can't create %dx%dx%d matrix. Exit\n", N, M, G); 
		myExit(1);
	}
	iDay = M;
	nrowsStart = M;

	for (int i = 0; i < N; i++)
		result()[i] = i;

	for (int i = 1; i < M; i++) {
		tchar* r = result(i);
		for (int j = 0; j < N; j += G, r += G) {
			for (int k = 0; k < G; k++)
				r[k] = ((j / G + i * k) % M) + M * k;
		}
	}
	linksFromMatrix(links(), result(), nrowsStart);
	memcpy(m_pResultsPrev, result(), m_nLenResults);
	memcpy(m_pResultsPrev2, result(), m_nLenResults);
	printResultWithHistory("Generated KC-matrix", m_numDaysResult);
	//printTableColor("circular matrix", result(), iDay, m_numPlayers, m_groupSize);
#else
	if (mStart0)
	{
		iDay = nrowsStart;
		memcpy(result(), mStart0, nrowsStart * m_numPlayers);
		linksFromMatrix(links(), result(), nrowsStart);
	}
#endif


	if (iDay > numDaysResult() && !(m_printMatrices & 16))
		iDay = numDaysResult(); // warning?

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
		m_lastRowWithTestedTrs = 0;
		setArraysForLastRow(iDay);
		minRows = iDay;
		iDay--;
		goto checkCurrentMatrix;
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
		int ir = checkCurrentResult(m_bPrint);
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

	bPrevResult = false;

#if 1 && !USE_CUDA 
	if (testGroupOrderEachSubmatrix(m_printMatrices, iCalcModeOrg)) {
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
						m_cTime = clock();
						for (int i = nPrecalcRows + 3; i < numDaysResult(); i++)
							m_rowTime[i] = m_cTime - m_iTime;
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
				if (iDay < nPrecalcRows)
					m_secondPlayerInRow4First = 0;
				else {
					int ipx = 0;
					if (m_playerIndex)
					{
						iDay = m_playerIndex / m_numPlayers;
						ipx = m_playerIndex % m_numPlayers;
						m_playerIndex = 0;
						if (iDay < nPrecalcRows)
						{
							if (param(t_MultiThreading))
							{
								noMoreResults = true;
								goto noResult;
							}
							m_pRowStorage->init();
							m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
							m_playerIndex = 0;
							nRows4 = nRows4Day = 0;
							continue;
						}
						if (iDay >= numDaysResult())
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
						if (m_bPrint && iDay < nPrecalcRows + 3) {
							m_cTime = clock();
							m_rowTime[iDay] = m_cTime - m_iTime;
						}
#endif
						if (++iDay < numDaysResult() && !checkCanonicity()) {

#if 0   // Temporary
							if (!p1f_counter || ((++m_p1f_counter) % p1f_counter))
#endif
#if !USE_CUDA
								if (!m_bPrint || m_cTime - m_rTime < ReportInterval)
#endif
									goto ProcessPrecalculatedRow;
						}

						m_pRowUsage->getMatrix(result(), neighbors(), iDay);
						/**
						if (m_bCheckLinkV && iDay < m_numDays) {
							linksFromMatrix(links(), result(), iDay);
							if (!checkLinks(links(), iDay)) {
								iDay--;
								goto ProcessPrecalculatedRow;
							}
						}*/
						iDay--;
#if !USE_CUDA
						m_cTime = clock();
						for (int i = nPrecalcRows + 3; i <= iDay; i++)
							m_rowTime[i] = m_cTime - m_iTime;
#endif
						goto checkCurrentMatrix;
					case -1: // reported if requested row not in solution (can happen only if LastRowSecondPlayer is incorrect)
					{
						if (secondPlayerInRow4Last < secondPlayerMax) {/**
							m_pRowUsage->getMatrix(result(), neighbors(), iDay);
							linksFromMatrix(links(), result(), iDay);
							if (cnvCheckNew(0, iDay, false) && (m_groupSize != 3 || checkLinks(links(), iDay))) {
								// matrix ok, error cannot be ignored*/
#if !USE_CUDA
								printfRed("*** Failing getRow(%d) = %d ***\n", iDay, retVal);
								printTable("Current result", result(), iDay, m_numPlayers, groupSize());
								exit(1);
#endif					
							//}
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
				if (iDay <= minRows || iDay < 2)
				{
					noMoreResults = true;
					goto noResult;
				}
				if (0)//iDay == 2 && groupSize_2 && m_use2RowsCanonization)
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
				if (nPrecalcRows && m_useRowsPrecalculation == eCalculateRows && m_secondPlayerInRow4) {

					if (nPrecalcRows > iDay) {
						m_useRowsPrecalculation = eCalculateRows;
						m_pRowStorage->init();
					}
					else if (nPrecalcRows == iDay) {
						ASSERT(m_secondPlayerInRow4 == 0);
						m_secondPlayerInRow4++;
						if (bCBMP && !m_groupSizeRemainder[m_secondPlayerInRow4])
							m_secondPlayerInRow4++;
						if (!nRows4Day && groupSize_2)
						{
							m_secondPlayerInRow4 = m_secondPlayerInRow4First;
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
						m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
						if (nRows4) {
							iDay = nPrecalcRows;
							if (m_bPrint) {
								printf("Total number of precalculated row solutions = %5d\n", nRows4);
								m_lastRowWithTestedTrs = 0;
#if !USE_CUDA
								if (m_printMatrices & 32) {
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
								if (param(t_MultiThreading)) {
									nLoops = 0;
									goto noResult;
								}
								iDay = nPrecalcRows;
								m_useRowsPrecalculation = eCalculateRows;
								m_pRowStorage->init();
								nRows4 = 0;
								iDay++;
								//linksFromMatrix(links(), result(), iDay);
								//iPlayer = m_numPlayers;
								bPrevResult = true;
								continue;
							}

							if (iCalcMode == eCalculateRows) {
								nLoops = nRows4;
								noMoreResults = true;
								goto noResult;
							}
							m_pRowUsage->init();
//						m_pRowStorage->reset();
							nRows4 = 0;
							continue;
						}
					}
				}
				bPrevResult = true;
				continue;
			}
			if (m_useRowsPrecalculation == eCalculateMatrices) {
				ASSERT(1);
				noMoreResults = true;
				goto noResult;
			}
			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			if (!m_secondPlayerInRow4First && nPrecalcRows && m_useRowsPrecalculation == eCalculateRows) {
				if (nPrecalcRows - 1 == iDay) {
					m_pRowStorage->initPlayerMask(NULL);
				//	setArraysForLastRow(nPrecalcRows);
				}
				else if (nPrecalcRows == iDay) {
					m_secondPlayerInRow4 = m_secondPlayerInRow4First = result(nPrecalcRows)[1];
					m_lastRowWithTestedTrs = 0;
				}
			}

		checkCurrentMatrix:

			if (m_lastRowWithTestedTrs >= iDay)
				m_lastRowWithTestedTrs = iDay - 1;

#if !USE_CUDA
			if (m_bPrint)
				reportCurrentMatrix();
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
			if (m_bPrint && (int)((++nMCreated) % 10000000) == 0)
				printf(" %zdM calls to checkCurrentResult\n", nMCreated / 1000000);
#endif

#if 1
			switch (checkCurrentResult(m_printMatrices, pIS_Canonizer)) {
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
						TrCycles trc;
						cyclesFor2Rows(m_TrCyclesFirst2Rows, &trc, neighbors(0), neighbors(1), result(0), result(1));
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
					ASSERT(!m_secondPlayerInRow4First);
#if 0
					if (iTest) {
						printTable("r3", result(nPrecalcRows), 1, m_numPlayers, m_groupSize);
						if (result(nPrecalcRows)[1] == 8)
							iTest = iTest;
					}
#endif
					tchar nb[MAX_PLAYER_NUMBER]; 
					u1fSetTableRow(nb, result(nPrecalcRows), true);
					const auto retVal = m_pRowStorage->addRow(result(nPrecalcRows), neighbors(nPrecalcRows), nb);

#if !USE_CUDA
					if (m_printMatrices & 32) {
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
						m_pRowStorage->init();
						m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
						m_playerIndex = 0;
						nRows4 = nRows4Day = 0;
					}
					bPrevResult = true;
					continue;
				}
			}
		}
		ASSERT(iDay < numDaysResult());

		auto flag = true;
		if (semiSymGraph && (flag = orderOfGroup() >= minGroupSize)) {
			int i = 2;
			for (; i <= 3; i++) {
				auto* pGroup = static_cast<RowGenerators*>(pAutGroup[i]);
				pGroup->createGroupAndOrbits(this);
				const auto *pObj = pGroup->getObject(0);
				int j = pGroup->lenObject();
				while (j-- && !pObj[j]);
				if (j >= 0)
					break;
			}

			flag = i > 3;
		}

		if (flag && orderOfGroup() >= param(t_resultGroupOrderMin)) {
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
				if (m_bPrint) {
					const auto nr = secondRowsDB->numObjects();
					if (nr) {
						printf("Second row %04d:", nr);
						printTableColor("", secondRowsDB->getObject(nr - 1), 1, m_numPlayers, m_groupSize);
					}
				}
				if (groupSize_2 && sysParam()->u1fCycles[0] && sysParam()->u1fCycles[0][0] <= secondRowsDB->numObjects())
					goto noResult;
				goto cont1;
			}

			if (m_bPrint)
			{
				m_rTime = m_cTime = clock();
				setConsoleOutputMode();
				//report result
#if !DEBUG_NextPermut
				printf("\n%d(%zd): %s-Matrix, build time=%d, time since start=%d\n", threadNumber, nLoops, m_fHdr, m_cTime - mTime, m_cTime - m_iTime);
#else
				extern int matr_cntr;
				printf("\n%5zd: %s-Matrix, matr_cntr = %d\n", nLoops, m_fHdr, matr_cntr);
#endif
			}

			if (m_bPrint || bSavingMatricesToDisk)
			{
				char stat[256];
				bool needOutput = false;
				matrixStat(neighbors(), iDay, &needOutput);
				if (needOutput)
					m_matrixDB.addObjDescriptor(orderOfGroup(), matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll));
				else
					stat[0] = '\0';

				Result.setInfo(stat);
				Result.setGroupOrder(orderOfGroup());
#if 0			// record result and print on screen (if m_bPrint==true)
				Result.printTable(result(), true, m_bPrint, numDaysResult());
#else			// record result without print on screen
				Result.printTable(result(), true, false, numDaysResult());
				if (m_bPrint) {
					printf("%d(%zd): " AUT "%d, % s\n", threadNumber, nLoops, orderOfGroup(), stat);
					// print on screen result with highlighted differences from prev result
					printResultWithHistory("", iDay);
					if (m_printMatrices & 2)
						printPermutationMatrices(2);
					//printTable("links", links(), m_numPlayers, m_numPlayers, 0);
				}

				if (orderOfGroup() > 1) {
					for (int i = 0; i < countof(pAutGroup); i++) {
						if (pAutGroup[i])
							pAutGroup[i]->makeGroupOutput(this, m_bPrint);
					}
				}
#endif
				//cnvCheckNew(2, iDay);
			}

			//checkCommonValues();

			//Result.printTable(neighbors(), true, ResultFile, m_bPrint, numDaysResult());
			//reportCheckLinksData();
			//printTable("p1f", neighbors(), iDay, m_numPlayers, 2);
cont1:
			StatReportAfterEachResult(ResetStat, "Stat for matrix result. iDay", iDay, m_bPrint); // see stat.h to activate
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
					m_finalKMindex, m_numPlayers, numDaysResult(), groupSize(), clock() - m_iTime);
			}
			else {
				str = format("\nThread {}: {} non-isomorphic matrices ({},{},{}) created\n",
					threadNumber, m_finalKMindex, m_numPlayers, numDaysResult(), m_groupSize);

				str += format("Thread execution time = {} ms\n", clock() - m_iTime);
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
