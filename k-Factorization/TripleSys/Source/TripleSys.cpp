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

CC sLongLong alldata::Run(int threadNumber, int iCalcMode,
	tchar* mStart0, tchar* mStart, int nrowsStart, int nrowsOut, sLongLong* pcnt, string *pOutResult, bool bFirstThread) {
	// Input parameter:
#if !USE_CUDA
	const auto iTime = clock();
	const char* fHdr = getFileNameAttr(sysParam());
	auto rTime = iTime;
	auto cTime = iTime;
	const auto bSavingMatricesToDisk = param(t_savingMatricesToDisk);
	int nMatricesMax = 0;
	int startMatrixCount = 0;
	int nBytesInStartMatrix = nrowsStart * m_numPlayers;
#endif
	const auto bPrint = bFirstThread && param(t_printMatrices);
	int minRows = nrowsStart;
	if (!nrowsOut)
		nrowsOut = m_numDays;
	m_numDaysResult = nrowsOut;

	CUDA_PRINTF("*** threadNumber = %d nrowsOut = %d, numDaysAdj = %d \n", threadNumber, nrowsOut, numDaysAdj);

#if !USE_CUDA
	unsigned char* bResults = NULL;
	string fName = format("{:0>10}.txt", threadNumber);
	if (m_improveResult) {
		const auto lenResult = (m_numDays + 1) * (m_numPlayers + m_numDays);
		bResults = new unsigned char[(m_improveResult > 1 ? 2 : 1) * lenResult];
		createFolderAndFileName(ImprovedResultFile, sysParam(), t_ImprovedResultFolder, m_numDaysResult, &fName);
	}

	createFolderAndFileName(ResultFile, sysParam(), t_ResultFolder, m_numDaysResult, &fName);

	TableAut Result("|Aut(M)|", m_numDays, m_numPlayers, 0, m_groupSize, true, true);
	Result.allocateBuffer(32);
	m_pRes = &Result;

	TableAut* pAutGroup = NULL;
	if (param(t_outAutomorphismGroup)) {
		pAutGroup = new TableAut("\nAutomorphisms of matrix with |Aut(M)|", groupOrder(), m_numPlayers, -1, m_groupSize, false, true);
		pAutGroup->allocateBuffer(64);
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
	CUDA_PRINTF("*** iDay = %d  m_numDaysResult = %d\n", iDay, m_numDaysResult);
	if (iDay > m_numDaysResult)
		iDay = m_numDaysResult; // warning?

	memset(m_rowTime, 0, m_numDays * sizeof(m_rowTime[0]));
	for (int i = 0; i < iDay; i++)
		p1fSetTableRow(p1ftable(i), result(i));

	if (m_p1f || param(t_u1f))
	{
		if (m_groupSize == 2) // need to be implemented for 3?
			p1fCheckStartMatrix(iDay);
	}
#if 0 // print group order for each submatrix
	printf("Rows-Group Order: ");
	for (int i = 2; i <= iDay; i++) {
		cnvCheckNew(0, i, false);
		printf("%d-%d ", i, groupOrder());
	}
	printf("\n");
#endif
	for (int i = firstGroupIdx(); i <= lastGroupIdx(); i++) {
		auto* pGroupInfo = groupInfo(i);
		if (!pGroupInfo)
			break;

		if (iDay < i)
			break;

		cnvCheckNew(0, i, false); // create initial set of tr for first i rows
		pGroupInfo->copyIndex(*this);
		resetGroupOrder();
#if 0
		int grOrder = 0;
		if (i > 2) {
			auto* pPrevGroup = groupInfo(i - 1);
			auto* pTmp = pPrevGroup;
			pPrevGroup = pGroupInfo;
			pGroupInfo = pTmp;
			grOrder = pGroupInfo->groupOrder();
			auto* cmpTr = pPrevGroup->getObject();
			for (int j = pPrevGroup->groupOrder(); --j;) {
				cmpTr += m_numPlayers;
				pGroupInfo->updateGroupOrder(cmpTr);
			}
			grOrder = pGroupInfo->groupOrder() - grOrder;
		}
#endif
	}

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
	const auto checkFlag = iDay == m_numDays || iDay == m_numDaysResult;
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
	if (iDay > 0) {
		setArraysForLastRow(iDay);
		//printTable("p1f", p1ftable(), iDay, m_numPlayers, 0);
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
		CUDA_PRINTF(" **** nLoops = %lld, iDay = %d  m_numDaysResult = %d  bPrevResult = %d\n", nLoops, iDay, m_numDaysResult, bPrevResult);
		CUDA_PRINTF(" *** nLoops = %lld\n", nLoops);
		while (iDay < m_numDaysResult || bPrevResult)
		{
			CUDA_PRINTF("   *** iDay = %d  bPrevResult = %d\n", iDay, bPrevResult);
			if (iDay < 0)
			{
				noMoreResults = true;
				goto noResult;
			}
			if (bPrevResult)
			{
				if (iDay == 2 && m_groupSize == 2 && (m_p1f || param(t_u1f)))
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
			//CUDA_PRINTF("   *** ProcessOneDay\n");
			if (iDay < minRows)
			{
				noMoreResults = true;
				goto noResult;
			}

			if (!processOneDay())
			{
				bPrevResult = true;
				continue;
			}

			memcpy(result(iDay), tmpPlayers, m_numPlayers);
		checkCurrentMatrix:

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
				StatReportPeriodically(ResetStat, "Stat current. iDay", iDay, bFirstThread);
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
			if (bPrint && (int)((++nMCreated) % 1000000) == 0)
				printf(" %zdM", nMCreated / 1000000);
#endif

#if 1
			switch (checkCurrentResult(bPrint, pIS_Canonizer)) {
				case -1: goBack(); goto ProcessOneDay;
				case  1: noMoreResults = true; goto noResult;
				default: break;
			}
#endif
		}
		ASSERT(iDay < m_numDaysResult);
		if (groupOrder() >= param(t_resultGroupOrderMin))
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
#if !USE_CUDA
			m_finalKMindex++;
			if (bPrint)
			{
				cTime = clock();
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
				matrixStat(p1ftable(), iDay, &needOutput);
				if (needOutput)
					m_matrixDB.addMatrix(groupOrder(), matrixStatOutput(stat, sizeof(stat)));
				else
					stat[0] = '\0';

				Result.setInfo(stat);
				Result.setGroupOrder(groupOrder());
#if 0			// record result and print on screen (if bPrint==true)
				Result.printTable(result(), true, ResultFile.c_str(), bPrint, m_numDaysResult);
#else			// record result without print on screen
				Result.printTable(result(), true, ResultFile.c_str(), false, m_numDaysResult);
				if (bPrint) {
					printf("%5zd: |Aut(M)| = %d, %s\n", nLoops, groupOrder(), stat);
					// print on screen result with highlighted differences from prev result
					printResultWithHistory("", iDay);
				}
#endif
				//cnvCheckNew(2, iDay);

				if (pAutGroup) {
					pAutGroup->setGroupOrder(groupOrder());
					string outFile, fName("_AutPermut.txt");
					createFolderAndFileName(outFile, sysParam(), t_ResultFolder, m_numDaysResult, &fName);
					pAutGroup->printTable(getObject(), true, outFile.c_str(), false, groupOrder(), getIndices());
				}
			}

			//checkCommonValues();

			//Result.printTable(p1ftable(), true, ResultFile, bPrint, m_numDaysResult);
			//reportCheckLinksData();
			//printTable("p1f", p1ftable(), iDay, m_numPlayers, 2);

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
		auto str = format("\nThread {}: {} non-isomorphic matrices ({},{},{}) created\n",
			threadNumber, m_finalKMindex, m_numPlayers, m_numDaysResult, m_groupSize);

#if GenerateSecondRowsFor3U1F
		printTable("Calculated second rows set", m_p3fSecondRows, m_p3fNumSecondRowsAct, m_numPlayers, m_groupSize);
#endif
		str += format("Thread execution time = {} ms\n", clock() - iTime);
		printf(str.c_str());
		if (pOutResult)
			*pOutResult += str;
	}
	StatReportAfterThreadEnd(ResetStat, "Thread ended, processed", (int)nLoops, bFirstThread); // see stat.h to enable
	delete[] bResults;
	delete pAutGroup;
#endif

	return nLoops;
}
