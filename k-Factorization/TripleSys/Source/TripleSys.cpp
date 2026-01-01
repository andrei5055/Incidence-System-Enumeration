#include <iostream>
#include "TableLS.h"
#include "TripleSys.h"
#include "kOrbits.h"

#if !USE_CUDA && USE_BINARY_CANONIZER
#include "k-SysSupport.h"
#include "CDTools.h"
#endif
#include <filesystem>

using namespace std;

CC sLongLong alldata::Run(int threadNumber, eThreadStartMode iCalcMode, CStorageSet<tchar>* secondRowsDB,
	ctchar* mStart0, ctchar* mfirst, int nrowsStart, sLongLong* pcnt, string* pOutResult, int iThread) {
	// Input parameters:
	const auto iCalcModeOrg = iCalcMode;
#if !USE_CUDA
	m_printMatrices = (iThread == 0 || param(t_numThreads) < 2) ? param(t_printMatrices) : 0;
	m_bPrint = (m_printMatrices & 1) != 0;
	bool bPrintPeriodic = (m_printMatrices & 32) != 0;
	m_bPrintAll = m_bPrint || bPrintPeriodic;
	m_fHdr = getFileNameAttr(sysParam());
	m_cTime = m_rTime = m_iTime = clock();

	const bool bNotSpecialMode = iCalcModeOrg != eCalcSecondRow && iCalcModeOrg != eCalculateRows;
	const auto bSavingMatricesToDisk = bNotSpecialMode ? param(t_savingMatricesToDisk) : false;
#endif
	int minRows = nrowsStart;
	m_doNotExitEarlyIfNotCanonical = (m_test & 2) != 0; // || param(t_generateMatrixExample) != 0;
	memset(m_rowTime, 0, m_numDaysResult * sizeof(m_rowTime[0]));
	m_allRowPairsSameCycles = param(t_any2RowsConvertToFirst2) != 0 || param(t_allowUndefinedCycles) == 0;
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
	m_pSecondRowsDB = secondRowsDB;
	if (iCalcMode == eCalcSecondRow) {
		iCalcMode = eCalcResult;
		m_createSecondRow = 1;
		m_numDaysResult = 2;
		nrowsStart = 0;
	} 
	else
		m_createSecondRow = 0;
	const bool bCBMP = !completeGraph();
	const auto groupSize_2 = m_groupSize == 2;

	if (iDay == 0)
		iDay = nrowsStart;
	if (iDay > numDaysResult() && !(m_printMatrices & 16) && !(m_test & 1))
		iDay = numDaysResult(); // warning?

	if (mStart0)
	{
		memcpy(result(), mStart0, nrowsStart * m_numPlayers);
		if (m_printMatrices & 64)
			printTable("Input matrix", result(), iDay, m_numPlayers, m_groupSize, 0, true);
		linksFromMatrix(links(), result(), nrowsStart);
	}

	for (int i = 0; i < iDay; i++)
		u1fSetTableRow(neighbors(i), result(i));

	initPrecalculationData(iCalcModeOrg, nrowsStart);

	int maxPlayerInFirstGroup = groupSize_2 ? m_secondPlayerInRow4Last : m_numPlayers;

	const auto semiSymGraph = !m_createSecondRow && numDaysResult() > 2 && param(t_semiSymmetricGraphs) == 1;
	const auto minGroupSize = semiSymGraph ? m_numDaysResult * m_numPlayers / 2 : 0;
	const auto outAutGroup = param(t_outAutomorphismGroup);
	IOutGroupHandle<tchar>* pAutGroup[4] = { NULL };
	if (semiSymGraph || outAutGroup) {
		if (bSavingMatricesToDisk) {
			if (outAutGroup & 1)
				pAutGroup[0] = new Generators<tchar>(outAutGroup, "\nOrbits and generators of the Aut(M) acting on elements", m_numPlayers);

			if (outAutGroup & 2)
				pAutGroup[1] = new COutGroupHandle<tchar>(outAutGroup, "\nAut(M) acting on elements", m_numPlayers);
		}

		if (semiSymGraph || bSavingMatricesToDisk && (outAutGroup & 12))
			pAutGroup[2] = new RowGenerators<tchar>(outAutGroup, numDaysResult());

		if (semiSymGraph || bSavingMatricesToDisk && (outAutGroup & 48))
			pAutGroup[3] = new CKOrbits(outAutGroup, m_numPlayers, m_groupSize, numDaysResult());
	}
#if !USE_CUDA
	sLongLong nMCreated = 0;
	auto mTime = clock();
	unsigned char* bResults = NULL;

	
	TableAut* pResult = NULL;
	int iSaveLS = m_groupSize <= 3 ? param(t_saveLatinSquareType) : 0;
	if (iSaveLS)
		pResult = new TableLS(MATR_ATTR, m_numDays, m_numPlayers, 0, m_groupSize, true, true, iSaveLS, bCBMP);
	else
		pResult = new TableAut(MATR_ATTR, m_numDays, m_numPlayers, 0, m_groupSize, true, true);

	if (bSavingMatricesToDisk) {
		string fName = format("{:0>10}.txt", threadNumber);
		if (m_improveResult) {
			const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
			bResults = new unsigned char[(m_improveResult > 1 ? 2 : 1) * lenResult];
			createFolderAndFileName(ImprovedResultFile, sysParam(), t_ImprovedResultFolder, numDaysResult(), fName);
		}

		createFolderAndFileName(ResultFile, sysParam(), t_ResultFolder, numDaysResult(), fName);

		pResult->allocateBuffer(32);
		const auto pResFile = ResultFile.c_str();
		if (param(t_keepPrevResult) && std::filesystem::exists(pResFile)) {
			printfRed("*** Error: Request to override existing file rejected (KeepPrevResult=%d, file name=%s)\n", 
				param(t_keepPrevResult), pResFile);
			myExit(1);
		}
		pResult->setOutFileName(pResFile);
		m_pRes = pResult;

		if (outAutGroup) {
			for (int i = 0; i < countof(pAutGroup); i++) {
				if (pAutGroup[i])
					pAutGroup[i]->setOutFileName(pResFile, false);
			}
		}
	}
#endif
	CUDA_PRINTF("*** mStart0 = %p\n", mStart0);
#if !USE_CUDA
	if (param(t_generateMatrixExample)) {
		iDay = m_numDaysResult;
		nrowsStart = m_numDaysResult;
		if (!generateMatrixExample())
			goto noResult;
		if (m_createSecondRow) {
			memcpy(m_pSecondRowsDB->getNextObject(), result(1), m_numPlayers);
			m_finalKMindex = 1;
			nLoops = 1;
			goto noResult;
		}
	}
	else if (sysParam()->strVal[t_InputDataFileName]->length() && !m_createSecondRow) {

		tchar tr[MAX_PLAYER_NUMBER];
		memset(tr, 0, m_numPlayers);
		for (int i = 0;i < m_numPlayers; i++)
			tr[result(0)[i]] = i;
		kmTranslate(result(0), result(0), tr, m_numDays * m_numPlayers);
		if (!canonizeMatrix(m_numDays))
			goto noResult;
		if (!linksFromMatrix(links(), result(), m_numDays, false)) {
			printfRed("*** Internal error\n");
			goto noResult;
		}
	}
#endif

	if (iDay >= 2 && (m_use2RowsCanonization || param(t_u1f)))
	{
		if (groupSize_2 && !param(t_generateMatrixExample))		// need to be implemented for 3?
			p1fCheckStartMatrix(iDay);
	}
	m_playerIndex = 0;
#if !USE_CUDA
#if 1   // Matrix from data.h
	if ((m_test & 1) && (iCalcModeOrg != eCalcSecondRow)){
		if (!iDay) {
			printfRed("*** Error: parameter Test=%d, but matrix not defined in data.h\n", m_test);
			myExit(1);
		}
		printfGreen(" SW uses input matrix from data.h (Test=%d)\n", m_test);
		canonizeMatrix(iDay);
#if 0
		minRows = iDay--;
		goto checkCurrentMatrix;
		//exit(0);
#endif
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
			goto noResult;
		} 
#endif
		int ir = checkCurrentResult(m_bPrint);
		//if (improveMatrix(m_improveResult, NULL, 0/*, bResults, lenResult()*/))
		//	return -1;
		// special case if input is full result matrix or UseTwoLastDaysG2 mode
		switch (ir) {
		case -1:
		case  1: goto noResult;
		default: break;
		}
	}
#endif

	bPrevResult = false;

#if 1 && !USE_CUDA 
	if (testGroupOrderEachSubmatrix(m_printMatrices, iCalcModeOrg)) {
		goto noResult;
	}
#endif

	if (iCalcMode == eCalculateMatrices)
		m_pRowUsage->init(iThread, param(t_numThreads));

	else if (iDay > 0) {
		if (m_precalcMode == eCalculateRows)
			m_pRowStorage->initPlayerMask(mfirst, maxPlayerInFirstGroup);
		setArraysForLastRow(iDay);
		//printTable("p1f", neighbors(), iDay, m_numPlayers, m_groupSize);
		minRows = iDay--;

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
			if (m_nPrecalcRows && m_precalcMode == eCalculateMatrices) {
			ProcessPrecalculatedRow:
				const auto iRet = precalculatedSolutions(iCalcMode);
				switch (iRet) {
				case eCheckCurrentMatrix: goto checkCurrentMatrix;
				case eContinue: continue;
				case eNoResult: goto noResult;
				}
			}
			else {
				CUDA_PRINTF("   *** iDay = %d  bPrevResult = %d\n", iDay, bPrevResult);
				if (iDay < 0)
				{
					goto noResult;
				}
				if (bPrevResult)
				{
					if (iDay <= minRows || iDay < 2)
						goto noResult;
					if (!initPrevDay())
						continue;
				}
				else if (!initCurrentDay())
					continue;
			}
		ProcessOneDay:
			if (m_lastRowWithTestedTrs >= iDay)
				m_lastRowWithTestedTrs = iDay - 1;
			if (iDay < minRows)
				goto noResult;
			// temporary check that link table is ok
#if 0 //!USE_CUDA
			if (!iPlayer && iDay) {
				tchar* l = links(0);
				for (int i = 0;i < m_numPlayers * m_numPlayers; i++) {
					if (l[i] == iDay) {
						printfRed("Bad link(%d) at %dx%d\n", l[i], i / m_numPlayers, i % m_numPlayers);
						//exit(1);
					}
				}
			}
#endif
			if (!processOneDay()) {
				if (m_nPrecalcRows && m_precalcMode == eCalculateRows && m_secondPlayerInRow4) {
					if (m_bPrint)
						printf("m_secondPlayerInRow4=%d m_nRows4Day=%d\n", m_secondPlayerInRow4, m_nRows4Day);
					if (!m_nRows4Day && m_groupSize == 2 && param(t_MultiThreading))
						goto noResult;
					const auto iRet = endOfRowPrecalculation(iCalcMode);
					switch (iRet) {
					case eCheckCurrentMatrix: goto checkCurrentMatrix;
					case eContinue: continue;
					case eNoResult: goto noResult;
					case eOk: break;
					}
				}
				bPrevResult = true;
				continue;
			}
			if (m_precalcMode == eCalculateMatrices) {
				ASSERT_IF(1); // error in SW
				goto noResult;
			}
			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			if (!m_secondPlayerInRow4First && m_nPrecalcRows && m_precalcMode == eCalculateRows) {
				if (m_nPrecalcRows - 1 == iDay) {
					m_pRowStorage->initPlayerMask(NULL, maxPlayerInFirstGroup);
				}
				else if (m_nPrecalcRows == iDay) {
					m_secondPlayerInRow4 = m_secondPlayerInRow4First = result(m_nPrecalcRows)[1];
					m_lastRowWithTestedTrs = 0;
				}
			}

		checkCurrentMatrix:

			if (m_lastRowWithTestedTrs >= iDay)
				m_lastRowWithTestedTrs = iDay - 1;

#if !USE_CUDA
			if (bPrintPeriodic)
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
			if (m_bPrintAll && (int)((++nMCreated) % 10000000) == 0)
				printf(" %zdM calls to checkCurrentResult\n", nMCreated / 1000000);
#endif

#if 1
			switch (checkCurrentResult(m_printMatrices, pIS_Canonizer)) {
			case -1:
				if (param(t_generateMatrixExample)) {
					printfYellow("\n*** GenerateMatrixExample=%d: generated matrix below processed but it is not canonical", param(t_generateMatrixExample));
					break;
				}
				goBack();
				if (m_nPrecalcRows && m_precalcMode == eCalculateMatrices)
					goto ProcessPrecalculatedRow;
				goto ProcessOneDay;
			case 1: goto noResult;
			default: break;
			}
#endif
			if (m_nPrecalcRows && m_precalcMode == eCalculateRows) {
				if (iDay == m_nPrecalcRows + 1) {
					addPrecalculatedRow();
					bPrevResult = true;
					continue;
				}
			}
		}
		ASSERT_IF(iDay < numDaysResult());

		auto flag = true;
		if (semiSymGraph && (flag = orderOfGroup() >= minGroupSize)) {
			int i = 2;
			for (; i <= 3; i++) {
				auto* pGroup = static_cast<RowGenerators<tchar>*>(pAutGroup[i]);
				pGroup->createGroupAndOrbits(this);
				const auto *pObj = pGroup->getObject(0);
				int j = pGroup->lenObject();
				while (j-- && !pObj[j]);
				if (j >= 0)
					break;
			}

			flag = i > 3;
		}

		if (flag && (m_createSecondRow || orderOfGroup() >= param(t_resultGroupOrderMin))) {
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
				memcpy(m_pSecondRowsDB->getNextObject(), result(1), m_numPlayers);
				if (m_bPrintAll) {
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
				char stat[1024];
				getAllCycles(neighbors(), iDay);
				m_matrixDB.addObjDescriptor(orderOfGroup(), matrixStatOutput(stat, sizeof(stat), m_TrCyclesAll));
				pResult->setInfo(stat);
				pResult->setGroupOrder(orderOfGroup());
#if 0			// record result and print on screen (if m_bPrint==true)
				pResult->printTable(result(), true, m_bPrint, numDaysResult());
#else			// record result without print on screen
				if (bSavingMatricesToDisk) {
					pResult->printTable(result(), true, false, numDaysResult());
					if (pResult->isLSCreated()) {
						m_nLS++;
						auto* pResultLS = static_cast<const TableLS*>(pResult);
						if (pResultLS->isAtomicLS())
							m_atomicLS++;
						if (pResultLS->isSymmetricLS())
							m_symmetricLS++;
						if (m_bPrint)
							printfYellow("Latin Square (%sAtomic, %sTotally symmetric) saved with the matrix below\n",
								pResultLS->isAtomicLS() ? "" : "Not ", pResultLS->isSymmetricLS() ? "" : "Not ");
					}
				}

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
			//printTable("Result", result(), iDay, m_numPlayers, m_groupSize, -1); // print comma separated values

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
					//myExit(1);
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

	if (m_bPrint && param(t_saveLatinSquareType) && m_nLS)
		printfYellow("Thread %d: Latin Squares=%d, Atomic=%d, Totally symmetric=%d\n", threadNumber, m_nLS, m_atomicLS, m_symmetricLS);
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
				/**
				if (param(t_saveLatinSquareType))
					str += format("Latin Squares {}, Atomic: {}, Totally symmetric: {}\n", m_nLS, m_atomicLS, m_symmetricLS);
				**/
				str += format("Thread execution time = {} ms\n", clock() - m_iTime);
			}
			printf(str.c_str());
			if (pOutResult)
				*pOutResult += str;
		}
	}
	StatReportAfterThreadEnd(ResetStat, "Thread ended, processed", (int)nLoops, bFirstThread); // see stat.h to enable

	delete pResult;
	delete[] bResults;

	for (int i = 0; i < countof(pAutGroup); i++)
		delete pAutGroup[i];

#endif
#if COUNT_GET_ROW_CALLS && !USE_CUDA
	extern ll getRowCallsCalls;
	extern size_t totalWeighChange;
	static int nMatr;
	static ll cntrTotal;
	printf("Matr# %4d: ********** cntr = %lld  = %lld  totalWeighChange = %zd\n", ++nMatr, getRowCallsCalls, cntrTotal += getRowCallsCalls, totalWeighChange);
	getRowCallsCalls = 0;
#endif
	return nLoops;
}
