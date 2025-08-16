#include <iostream>
#include "TripleSys.h"
#include "kOrbits.h"

#if !USE_CUDA && USE_BINARY_CANONIZER
#include "k-SysSupport.h"
#include "CDTools.h"
#endif

using namespace std;

CC sLongLong alldata::Run(int threadNumber, eThreadStartMode iCalcMode, CStorageSet<tchar>* secondRowsDB,
	ctchar* mStart0, ctchar* mfirst, int nrowsStart, sLongLong* pcnt, string* pOutResult, int iThread) {
	// Input parameters:
	m_test = param(t_test);
	m_ignoreCanonizationMinus1 = (m_test & 0x22) != 0;
	const auto iCalcModeOrg = iCalcMode;
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
	m_iTime = clock();
	m_fHdr = getFileNameAttr(sysParam());
	m_rTime = m_iTime;
	m_cTime = m_iTime;
	const auto bSavingMatricesToDisk = (iCalcModeOrg != eCalcSecondRow && iCalcModeOrg != eCalculateRows) ?
		param(t_savingMatricesToDisk) : false;
#endif
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

	initPrecalculationData(iCalcModeOrg, nrowsStart);

	int maxPlayerInFirstGroup = groupSize_2 ? m_secondPlayerInRow4Last : m_numPlayers;

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
	sLongLong nMCreated = 0;
	auto mTime = clock();
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
#if !USE_CUDA
	// Generate KC matrix
	if (m_test & 0x20) {
		if (!m_createSecondRow && m_groupSize == 3 && !(m_nGroups & 1)) {
			printfRed("*** Error: can't generate matrix. Number of groups(%d) must be odd number with group size=%d\n", m_nGroups, m_groupSize);
			goto noResult; 
		}
		iDay = m_numDaysResult;
		nrowsStart = m_numDaysResult;

		for (int i = 0; i < m_numDays; i++) {
			tchar* r = m_Km + i * m_numPlayers;
			for (int j = 0; j < m_numPlayers; j += m_groupSize, r += m_groupSize) {
				for (int k = 0; k < m_groupSize; k++) {
					int ip = ((j / m_groupSize + i * k) % m_numDays);
					ip = ip * m_groupSize + k;
					r[k] = ip;
				}
			}
		}
		tchar tm[MAX_PLAYER_NUMBER];
		memset(tm, 0, sizeof(tm));
		(this->*m_pSortGroups)(m_Km, m_numDays);
		auto* coi = m_Ktmp;
		auto* cii = m_Km;
		for (int i = 0; i < m_numDays; i++, coi += m_numPlayers, cii += m_numPlayers)
			kmSortGroupsByFirstValue(cii, coi);
		// Result of the loop above is in m_Ktmp, sort and send it to m_Km.
		kmSortRowsBy2ndValue(m_numDays, tm);
		memcpy(result(), m_Km, m_nLenResults);
		if (m_createSecondRow) {
			memcpy(m_pSecondRowsDB->getNextObject(), result(1), m_numPlayers);
			m_finalKMindex = 1;
			nLoops = 1;
			goto noResult;
		}
		if (!linksFromMatrix(links(), result(), nrowsStart, false)) {
			printfRed("*** Error: can't generate matrix for number of groups=%d and group size=%d\n", m_nGroups, m_groupSize);
			printfRed("***        number of groups (NPlayers / GroupSize) with group size > 3 must be prime number\n");
			goto noResult;
		}
		memcpy(m_pResultsPrev, result(), m_nLenResults);
		memcpy(m_pResultsPrev2, result(), m_nLenResults);
		printResultWithHistory("Generated KC-matrix", m_numDaysResult);
		iDay = m_numDaysResult;
	}
#endif
	if (mStart0)
	{
		iDay = nrowsStart;
		memcpy(result(), mStart0, nrowsStart * m_numPlayers);
		linksFromMatrix(links(), result(), nrowsStart);
	}
	if (iDay > numDaysResult() && !(m_printMatrices & 16))
		iDay = numDaysResult(); // warning?

	for (int i = 0; i < iDay; i++)
		u1fSetTableRow(neighbors(i), result(i));

	if (iDay >= 2 && (m_use2RowsCanonization || param(t_u1f)))
	{
		if (groupSize_2 && !(m_test & 0x20))		// need to be implemented for 3?
			p1fCheckStartMatrix(iDay);
	}
	m_playerIndex = 0;
#if !USE_CUDA
#if 1   // Matrix from data.h
	if (m_test & 1) {
		if (!iDay || !mStart0) {
			printfRed("*** Error: parameter Test=%d, but matrix not defined in data.h\n", m_test);
			//testCanonizatorSpeed(););
			myExit(1);
		}
		printfGreen(" SW uses input matrix from data.h (Test=%d)\n", m_test);
		//testCanonizatorSpeed();
		setArraysForLastRow(iDay);
		minRows = iDay--;
		goto checkCurrentMatrix;
		//exit(0);
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
		if (m_useRowsPrecalculation == eCalculateRows)
			m_pRowStorage->initPlayerMask(mfirst, maxPlayerInFirstGroup);
		setArraysForLastRow(iDay);
		//printTable("p1f", neighbors(), iDay, m_numPlayers, m_groupSize);
		minRows = iDay--;
		//testRightNeighbor(iDay + 1);
		if ((m_test & 0x10) && !p1f16())
			goto noResult;
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
			if (m_nPrecalcRows && m_useRowsPrecalculation == eCalculateMatrices) {
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
			if (!processOneDay())
			{
				if (m_nPrecalcRows && m_useRowsPrecalculation == eCalculateRows && m_secondPlayerInRow4) {
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
			if (m_useRowsPrecalculation == eCalculateMatrices) {
				ASSERT(1); // error in SW
				goto noResult;
			}
			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			if (!m_secondPlayerInRow4First && m_nPrecalcRows && m_useRowsPrecalculation == eCalculateRows) {
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
				if (m_test & 0x20) {
					printfYellow("\n*** Matrix below processed and accepted but it is not canonical (Test=%d)", m_test);
					break;
				}
				goBack();
				if (m_nPrecalcRows && m_useRowsPrecalculation == eCalculateMatrices)
					goto ProcessPrecalculatedRow;
				goto ProcessOneDay;
			case 1: goto noResult;
			default: break;
			}
#endif
			if (m_nPrecalcRows && m_useRowsPrecalculation == eCalculateRows) {
				if (iDay == m_nPrecalcRows + 1) {
					addPrecalculatedRow();
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
				memcpy(m_pSecondRowsDB->getNextObject(), result(1), m_numPlayers);
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
				if (bSavingMatricesToDisk)
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
