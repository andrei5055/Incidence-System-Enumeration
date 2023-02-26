#include "stdafx.h"
#include "DataTypes.h"
#include "EnumInfo.h"
#include "Enumerator.h"
#include "ThreadEnumerator.h"

template class CEnumInfo<TDATA_TYPES>;
template class CInsSysEnumInfo<TDATA_TYPES>;

static const int outDiv[] = { ':', ':', ':', '.' };

FClass2(CEnumInfo, size_t)::convertTime(float time, char *buffer, size_t lenBuf, bool alignment) const
{
	const int secTime = (int)time;
	const int minTime = secTime / 60;
	const int hourTime = minTime / 60;
	const int out[] = { hourTime / 24, hourTime % 24, minTime % 60, secTime % 60 };

	bool flag = false;
	char *pBuf = buffer;
	for (int i = 0; i < countof(outDiv); i++) {
		if (out[i] > 0 || flag || outDiv[i] == '.') {
			pBuf += SNPRINTF(pBuf, lenBuf, (flag? "%02d%c" : "%2d%c"), out[i], outDiv[i]);
			flag = true;
		} else
		if (alignment)
			pBuf += SNPRINTF(pBuf, lenBuf, "   ");
		else
			continue;

		lenBuf -= 3;
	}

	pBuf += SNPRINTF(pBuf, lenBuf, "%02d sec%s", (int)(100 * (time - secTime)), alignment? "" : ",\n");
	return pBuf - buffer;
}

FClass2(CEnumInfo, double)::stringToTime(char *pTime)
{
	const double mult[] = { 3600 * 24, 3600, 60, 1, 0.01 };

	double time = 0;
	size_t i = countof(outDiv);
	char* pDivPlace[countof(outDiv)] = {};
	char *pDiv = NULL;
	while (i-- && (pDiv = strrchr(pTime, outDiv[i])) != NULL) {
		time += atoi(pDiv + 1) * mult[i + 1];
		*(pDivPlace[i] = pDiv) = '\0';
	}

	pDiv = strrchr(pTime, ' ');
	if (!pDiv)
		pDiv = pTime;
	else
		pDiv++;

	time += atoi(pDiv) * mult[i + 1];

	// Recreating input value
	while (++i < countof(pDivPlace))
		*pDivPlace[i] = outDiv[i];

	return time;
}

FClass2(CEnumInfo, bool)::compareTime(char *pTime1, char *pTime2)
{
	const double time1 = stringToTime(pTime1);
	const double time2 = stringToTime(pTime2);
	return time1 < time2;
}

FClass2(CEnumInfo, void)::updateEnumInfo(const CEnumInfo *pInfo)
{
	incrConstrCanonical(pInfo->numMatrOfType(t_design_type::t_canonical));
	incrConstrTotal(pInfo->numMatrOfType(t_design_type::t_totalConstr));
}

FClass2(CEnumInfo, void)::setReportFileName(const char *pntr)
{ 
	delete [] reportFileName();
	if (pntr) {
		const auto len = strlen(pntr) + 1;
		m_pReportFileName = new char[len];
		memcpy_s(m_pReportFileName, len, pntr, len);
	} else
		m_pReportFileName = NULL;
}

#define ARROW   "==>"

FClass2(CEnumInfo, void)::reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo, const WorkingInfo* pNumLevels)
{
	// Only master will report the progress and
	// only when there are no output of the solutions on the screen i.e. PRINT_SOLUTIONS = 0
	if (!m_bReportProgress)
		return;

	const auto nCanon = numMatrOfType(t_design_type::t_canonical);
	const ulonglong *pTestNumber;
	bool reportNeeded = false;
	const auto currClock = clock();
	switch (reportType) {
	case t_reportCriteria::t_reportByTime:
								if (currClock - prevClockReport() < CLOCKS_PER_SEC * 30)
									return;

								setPrevClockReport(currClock);
								reportNeeded = true;
		case t_reportCriteria::t_reportNow:
								pTestNumber = NULL;
								break;
		case t_reportCriteria::t_matrConstructed:
								pTestNumber = &nCanon;
								break;
		case t_reportCriteria::t_treadEnded:
								if (currClock - prevClockReport() < CLOCKS_PER_SEC)
									return;
								pTestNumber = &m_nCounter;
	}

	if (pTestNumber && *pTestNumber % reportInt())
		return;

	const auto numMatrTotalConstructed = numMatrOfType(t_design_type::t_totalConstr);
	double recentCanonProportion = -1;
	if ((reportNeeded || nCanon >= reportBound()) && nCanon >= prevReportCounter() || nCanon > prevReportCounter()) {
		if (nCanon > prevReportCounter() && numMatrTotalConstructed > prevReportCounter(1))
			recentCanonProportion = (static_cast<double>(numMatrTotalConstructed - prevReportCounter(1))) / (nCanon - prevReportCounter());

		setPrevReportCounter(nCanon);
		setPrevReportCounter(numMatrTotalConstructed, 1);
		switch (reportType) {
			case t_reportCriteria::t_matrConstructed:
			case t_reportCriteria::t_reportNow:
			case t_reportCriteria::t_treadEnded:
				pGroupInfo = this;
			case t_reportCriteria::t_reportByTime: {
					FOPEN(file, reportFileName(), "w");
					if (file) {
						fprintf(file, "%s\n", strToScreen());
						outEnumInfo(&file, false, pGroupInfo);
						setReportBound(REPORT_INTERVAL_OBJ_NUMB * ((nCanon + REPORT_INTERVAL_OBJ_NUMB - 1) / REPORT_INTERVAL_OBJ_NUMB));
					}
				}
		}
	}

	char buffer[512];
//	std::stringstream buffer;
	const float runTime = (float)(currClock - startTime()) / (60 * CLOCKS_PER_SEC);
	auto len = SPRINTF(buffer, "  Canon: %lld  NRB: %s%lld  Total: %lld  RunTime: %.2f min.",
		nCanon, constructedAllNoReplBlockMatrix() ? "=" : "", numbSimpleDesign(), numMatrTotalConstructed, runTime);

	if (nCanon) {
		len += snprintf(buffer + len, countof(buffer) - len, "  Efficiency: %.1f", static_cast<double>(numMatrTotalConstructed) / nCanon);

		if (recentCanonProportion > 0)
			len += snprintf(buffer + len, countof(buffer) - len, "  latest: %.1f", recentCanonProportion);
	}

	if (pNumLevels) {
		len += snprintf(buffer + len, countof(buffer), "  Threads:");
		for (int i = pNumLevels->minRow; i <= pNumLevels->maxRow; i++) {
			const auto numThreads = pNumLevels->pNumThreadsOnRow[i];
			if (numThreads)
				len += snprintf(buffer + len, countof(buffer) - len, " %d:%d", i, numThreads);
		}
	}

	cleanEndOfLine(buffer);
	std::cout << '\r' << strToScreen() << (reportType == t_reportCriteria::t_reportByTime ? ARROW : "   ") << buffer;
	fflush(stdout);

#if TEST
	startPrinting = nCanon >= START_PRINTING_AFTER;
#endif
	// Adjust reportInt to report approx. 1 time per 60 seconds
	if (!pTestNumber || currClock == prevClock())
		return;

	const float averPermMin = 60 * CLOCKS_PER_SEC * (*pTestNumber - prevCounter()) / (float)(currClock - prevClock());
	setPrevCounter(*pTestNumber);
	setPrevClock(currClock);
	size_t repInt = 1;
	while (averPermMin > 5 * repInt)
		repInt *= 10;

	setReportInt(repInt);
}

FClass2(CEnumInfo, void)::reportProgress(Class2(CThreadEnumerator) **ppThreadEnum, size_t nThread)
{
	if (nThread >= 1) {
		// Save already collected information 
		const auto nCanon = numMatrOfType(t_design_type::t_canonical);
		const auto nTotal = numMatrOfType(t_design_type::t_totalConstr);
		const auto nrbTotal = numbSimpleDesign();
		CGroupsInfo groupsInfo;
		groupsInfo.updateGroupInfo(this);

		int nTreadsOnRow[64];  // Number of threads currently working on each row
		int nRowMin, nRowMax = 0;
		const auto nRows = nRowMin = designInfo()->v + 2;
		int* pTreadsOnRow = nRows < countof(nTreadsOnRow) ? nTreadsOnRow : new int[nRows];
		memset(pTreadsOnRow, 0, nRows * sizeof(pTreadsOnRow[0]));
		for (size_t i = 0; i < nThread; i++, ppThreadEnum++) {
			const auto* pThreadEnum = *ppThreadEnum;
			if (pThreadEnum->code() != t_threadRunning)
				continue;

			updateEnumInfo(pThreadEnum->enumInfo());
			groupsInfo.updateGroupInfo(pThreadEnum->enumInfo());
			const auto nRow = pThreadEnum->enumerator()->numRow();
			if (nRowMin > nRow) nRowMin = nRow;
			if (nRowMax < nRow) nRowMax = nRow;
			++pTreadsOnRow[nRow];
		}

		WorkingInfo WorkInfo(pTreadsOnRow, nRowMin, nRowMax);
		reportProgress(t_reportCriteria::t_reportByTime, &groupsInfo, nRowMin <= nRowMax? &WorkInfo : NULL);

		// Restore information
		setNumMatrOfType(nCanon, t_design_type::t_canonical);
		setNumMatrOfType(nTotal, t_design_type::t_totalConstr);
		setNumbSimpleDesign(nrbTotal);
		if (pTreadsOnRow != nTreadsOnRow)
			delete[] pTreadsOnRow;
	} else {
		updateEnumInfo(ppThreadEnum[0]->enumInfo());
		reportThreadProgress();
	}
}

FClass2(CEnumInfo, void)::outAdditionalInfo(ulonglong nMatr, FILE *outFile, char *buff, size_t lenBuff) const {
	SNPRINTF(buff, lenBuff, "%10llu " OF_THEM " %s transitive on the rows.\n", nMatr, nMatr == 1 ? "is" : "are");
	outString(buff, outFile);
}

FClass2(CEnumInfo, void)::outEnumInfo(FILE **pOutFile, bool removeReportFile, const CGroupsInfo *pGroupInfo, const char* pComment)
{
	setRunTime();
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	if (!pGroupInfo)
		pGroupInfo = this;

	COrderInfo total(0, 1, 0);
	pGroupInfo->printGroupInfo(outFile, total);
	if (!(m_pParam->outType & t_Summary))
		return;

	auto nMatr = numMatrOfType(t_design_type::t_canonical);
	char buff[256];
	SPRINTF(buff, "\n%10llu matri%s" CONSTRUCTED_IN " ", nMatr, nMatr == 1 ? "x was" : "ces were");
	const auto len = strlen(buff);
	convertTime(runTime(), buff + len, countof(buff) - len, false);
	outString(buff, outFile);

	if (nMatr) {
		outAdditionalInfo(total.numMatrOfType(t_design_type::t_transitive), outFile, buff, sizeof(buff));
		nMatr = total.numMatrOfType(t_design_type::t_simple);
		SPRINTF(buff, "%10llu matri%s no replicated blocks,\n", nMatr, nMatr == 1 ? "x has" : "ces have");
		outString(buff, outFile);

		outAdditionalInfo(total.numMatrOfType(t_design_type::t_simpleTrans), outFile, buff, sizeof(buff));
	}

	nMatr = numMatrOfType(t_design_type::t_totalConstr);
	SPRINTF(buff, "%10llu matri%s fully constructed.\n", nMatr, nMatr == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile, true, pComment);
	if (removeReportFile) // Remove temporary file with the intermediate results	
		remove(reportFileName());
}

FClass2(CEnumInfo, void)::outEnumAdditionalInfo(FILE **pOutFile) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	char buff[256];
	const auto nTotal = numMatrOfType(t_design_type::t_totalConstr);
	SPRINTF(buff, "%10llu matri%s fully constructed\n", nTotal, nTotal == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile);
}

FClass2(CEnumInfo, void)::outEnumInformation(FILE **pOutFile, bool printMTlevel, const char *pComment) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

#define SPACE "        "
	char buff[256];
	SPRINTF(buff, "\n" SPACE "Using %zd-bit program with%s Assembly\n", sizeof(size_t) << 3, USE_ASM? "" : " no");
	auto outLen = outString(buff, outFile);

	const auto dLen = sizeof(SIZE_TYPE);
	char* pDataType = dLen == 1 ? "char" : (dLen == 2 ? "int16" : "int32");
	SPRINTF(buff, SPACE "Main data storage format: \"unsigned %s\"\n", pDataType);
	outLen += outString(buff, outFile);

	const auto nThreads = designInfo()->threadNumb;
	if (nThreads > 1) {
		char buffer[32];
		if (printMTlevel)
			SPRINTF(buffer, "row %d", multiThreadLevel());
		else
			SPRINTF(buffer, "different rows");

		SPRINTF(buff, "%10zd threads launched on %s (%swaiting to finish mode)\n", nThreads, buffer, WAIT_THREADS ? "" : "not ");
	} else
		SPRINTF(buff, SPACE "Single thread mode\n");

	outLen += outString(buff, outFile);

	if (nThreads >= 1) {
		SPRINTF(buff, SPACE "Solutions obtained by master are %s by the thread%s\n", designInfo()->use_master_sol ? "used" : "copied", nThreads>1? "s" : "");
		outLen += outString(buff, outFile);
		if (CANON_ON_GPU) {
			SPRINTF(buff, SPACE"        Canonicity was tested on GPU by %zd checkers (%d for each thread)\n", NUM_GPU_WORKERS * nThreads, NUM_GPU_WORKERS);
			outLen += outString(buff, outFile);
		}
	}

	if (designInfo()->find_master_design)
		SPRINTF(buff, "        Thread master DB was %sused\n", designInfo()->thread_master_DB ? "" : "not ");

	SPRINTF(buff, "        Canonicity of partial constructed matrix was %sused\n", USE_CANON_GROUP ? "" : "not ");
	outLen += outString(buff, outFile);

	SPRINTF(buff, "        Strong canonicity was %sused\n", USE_STRONG_CANONICITY ? "" : "not ");
	outLen += outString(buff, outFile);

	SPRINTF(buff, "        Super strong canonicity was %sused\n\n", USE_STRONG_CANONICITY_A ? "" : "not ");
	outLen += outString(buff, outFile);
	if (pComment) {
		SPRINTF(buff, "%s\n", pComment);
		outLen += outString(buff, outFile);
	}

	FCLOSE(outFile);
	*pOutFile = NULL;
	designInfo()->rewindLen = outLen;
}

FClass2(CEnumInfo, void)::updateConstrCounters(int matrFlags, const EnumeratorPntr pEnum) {
	incrConstrCanonical();
	const ulonglong simpleMatrFlag = matrFlags & t_noReplicatedBlock ? 1 : 0;
	if (simpleMatrFlag)
		incNumbSimpleDesign();

	auto *pOrderInfo = addGroupOrder(pEnum->groupOrder(), pEnum->extraGroupOrder()->groupOrder(), 1, simpleMatrFlag);
	if (matrFlags & t_transitiveGroup || pEnum->groupIsTransitive())
		pOrderInfo->addMatrixTrans(1, simpleMatrFlag);
}

FClass2(CEnumInfo, size_t)::cleanEndOfLine(char* pBuffer) const {
	const auto lenOut = strlen(pBuffer);
	auto retVal = lenOut;
	if (m_lenPrev > lenOut) {
		// To cover the previous output with the spaces.
		memset(pBuffer + lenOut, ' ', m_lenPrev - lenOut);
		pBuffer[retVal = m_lenPrev] = '\0';
	}

	m_lenPrev = lenOut;
	return retVal;
}

FClass2(CEnumInfo, void)::updateCounters(EnumInfoPntr pNumbInfo) {
	updateGroupInfo(pNumbInfo);
	addMatrices(pNumbInfo);
	pNumbInfo->resetNumbInfo();
}

#if CANON_ON_GPU
FClass2(CEnumInfo, void)::RecalcCountersByGroupOrders(const COrderInfo* pOrderInfo, size_t nElem) {
	// Recalculating counters by the group order info
	updateGroupInfo(pOrderInfo, nElem);

	COrderInfo total;
	calcCountersTotal(&total);
	init();
	resetNumbInfo(1);
	const ulonglong nSimple = total.numMatrOfType(t_simple);
	addMatrix(total.numMatrOfType(t_design_type::t_canonical), nSimple);
	incNumbSimpleDesign(nSimple);
}
#endif

FClass2(CInsSysEnumInfo, void)::updateEnumInfo(const EnumInfoPntr pInfo)
{
	if (!pInfo) return;
	Class2(CEnumInfo)::updateEnumInfo(pInfo);
	incNumbSimpleDesign(pInfo->numbSimpleDesign());
}


FClass2(CInsSysEnumInfo, void)::reportResult(char *buffer, int lenBuffer) const
{
	size_t len = SNPRINTF(buffer, lenBuffer, "%s       %9llu       %9llu",
						  this->strToScreen(), this->numMatrOfType(t_design_type::t_canonical), numbSimpleDesign());
	len += SNPRINTF(buffer + len, lenBuffer - len, "       %8llu    ", this->numMatrOfType(t_design_type::t_totalConstr));
	len += this->convertTime(this->runTime(), buffer + len, lenBuffer - len);

	// Prepare the comments regarding the results
	const char *pResComment = NULL;
	switch (this->getResType()) {
	case t_resType::t_resNew:			pResComment = "N"; break;
	case t_resType::t_resBetter:		pResComment = "B"; break;
	case t_resType::t_resWorse:			pResComment = "W"; break;
	case t_resType::t_resPostponed:     pResComment = "P"; break;
	case t_resType::t_resInconsistent:	pResComment = "???";
	}

	SNPRINTF(buffer + len, lenBuffer - len, "  %s", pResComment);
	len = strlen(this->strToScreen()) + strlen(ARROW);
	len += this->cleanEndOfLine(buffer + len);
	SNPRINTF(buffer + len, lenBuffer - len, "\n");
}

