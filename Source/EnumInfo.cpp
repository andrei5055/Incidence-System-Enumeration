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

	pBuf += SNPRINTF(pBuf, lenBuf, "%02d%s", (int)(100 * (time - secTime)), alignment? "" : "\n");
	return pBuf - buffer;
}

FClass2(CEnumInfo, double)::stringToTime(char *pTime)
{
	const double mult[] = { 3600 * 24, 3600, 60, 1, 0.01 };

	double time = 0;
	size_t i = countof(outDiv);
	char *pDivPlace[countof(outDiv)];
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
	incrConstrCanonical(pInfo->constrCanonical());
	incrConstrTotal(pInfo->numMatrOfType(t_totalConstr));
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

FClass2(CEnumInfo, void)::reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo)
{
	// Only master will report the progress
	if (PRINT_SOLUTIONS || !strToScreen() || !reportFileName())
		return;

	const ulonglong nCanon = constrCanonical();
	const ulonglong *pTestNumber;
	bool reportNeeded = false;
	clock_t currClock = clock();
	switch (reportType) {
		case t_reportByTime:	if (currClock - prevClockReport() < CLOCKS_PER_SEC * 30)
									return;

								setPrevClockReport(currClock);
								reportNeeded = true;
		case t_reportNow:		pTestNumber = NULL; break;
		case t_matrConstructed:	pTestNumber = &nCanon;
								break;
		case t_treadEnded:		if (currClock - prevClockReport() < CLOCKS_PER_SEC)
									return;
								pTestNumber = &m_nCounter;
	}

	if (pTestNumber && *pTestNumber % reportInt())
		return;

	if (reportNeeded || nCanon >= reportBound()) {
		switch (reportType) {
			case t_matrConstructed:
			case t_reportNow:
			case t_treadEnded:	  pGroupInfo = this;
			case t_reportByTime:
				{
					FOPEN(file, reportFileName(), "w");
					if (file) {
						fprintf(file, "%s\n", strToScreen());
						outEnumInfo(&file, false, pGroupInfo);
						setReportBound(REPORT_INTERVAL_OBJ_NUMB * ((nCanon + REPORT_INTERVAL_OBJ_NUMB - 1) / REPORT_INTERVAL_OBJ_NUMB));
					}
				}
		}
	}

	const float runTime = (float)(currClock - startTime()) / (60 * CLOCKS_PER_SEC);
	std::cout << '\r' << strToScreen() << (reportType == t_reportByTime ? "==>" : "   ")
			  << "  Canon: " << nCanon 
			  << "  NRB: "   << (constructedAllNoReplBlockMatrix() ? "=" : "") << numbSimpleDesign() 
			  << "  Total: " << numMatrOfType(t_totalConstr)
			  << "  RunTime: " << runTime << " min.";
	fflush(stdout);

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

FClass2(CEnumInfo, void)::reportProgress(const Class2(CThreadEnumerator) *pThreadEnum, size_t nThread)
{
	if (nThread >= 1) {
		// Save already collected information 
		const ulonglong nCanon = constrCanonical();
		const ulonglong nTotal = numMatrOfType(t_totalConstr);
		const ulonglong nrbTotal = numbSimpleDesign();
		CGroupsInfo groupsInfo;
		groupsInfo.updateGroupInfo(this);
		for (size_t i = 0; i < nThread; i++, pThreadEnum++) {
			switch (pThreadEnum->code()) {
				case t_threadRunning:
				case t_threadFinished: break;
				default:	continue;
			}

			updateEnumInfo(pThreadEnum->enumInfo());
			groupsInfo.updateGroupInfo(pThreadEnum->enumInfo());
		}
			
		reportProgress(t_reportByTime, &groupsInfo);
		// Restore information
		setNumMatrOfType(nCanon, t_canonical);
		setNumMatrOfType(nTotal, t_totalConstr);
		setNumbSimpleDesign(nrbTotal);
	} else {
		updateEnumInfo(pThreadEnum->enumInfo());
		reportThreadProgress();
	}
}

FClass2(CEnumInfo, void)::outEnumInfo(FILE **pOutFile, bool removeReportFile, const CGroupsInfo *pGroupInfo)
{
	setRunTime();
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	if (!pGroupInfo)
		pGroupInfo = this;

	pGroupInfo->printGroupInfo(outFile);
	if (!(m_pParam->outType & t_Summary))
		return;

	const ulonglong nConstrMatr = constrCanonical();
	char buff[256];
	SPRINTF(buff, "\n%10llu matri%s" CONSTRUCTED_IN " ", nConstrMatr, nConstrMatr == 1 ? "x" : "ces");
	const size_t len = strlen(buff);
	convertTime(runTime(), buff + len, countof(buff) - len, false);
	outString(buff, outFile);

	const ulonglong nMatr = numbSimpleDesign();
	if (nConstrMatr > 0) {
		SPRINTF(buff, "%10llu matri%s ha%s no replicated blocks\n", nMatr, nMatr == 1 ? "x" : "ces", nMatr == 1 ? "s" : "ve");
		outString(buff, outFile);
	}

	const auto nTotal = numMatrOfType(t_totalConstr);
	SPRINTF(buff, "%10llu matri%s fully constructed\n", nTotal, nTotal == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile);
	if (removeReportFile) // Remove temporary file with the intermediate results	
		remove(reportFileName());
}

FClass2(CEnumInfo, void)::outEnumAdditionalInfo(FILE **pOutFile) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	char buff[256];
	const auto nTotal = numMatrOfType(t_totalConstr);
	SPRINTF(buff, "%10llu matri%s fully constructed\n", nTotal, nTotal == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile);
}

FClass2(CEnumInfo, void)::outEnumInformation(FILE **pOutFile, bool printMTlevel) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	char buff[256];
	SPRINTF(buff, "\n        Using %zd-bit program with%s Assembly\n", sizeof(size_t) << 3, USE_ASM? "" : " no");
	auto outLen = outString(buff, outFile);

	const auto dLen = sizeof(SIZE_TYPE);
	char* pDataType = dLen == 1 ? "char" : (dLen == 2 ? "int16" : "int32");
	SPRINTF(buff, "        Main data storage format: \"unsigned %s\"\n", pDataType);
	outLen += outString(buff, outFile);

	const auto nThreads = designInfo()->threadNumb;
	if (USE_THREADS >= 1) {
		char buffer[32];
		if (printMTlevel)
			SPRINTF(buffer, "row %d", multiThreadLevel());
		else
			SPRINTF(buffer, "different rows");

		SPRINTF(buff, "%10zd threads launched on %s (%swaiting to finish mode)\n", nThreads, buffer, WAIT_THREADS ? "" : "not ");
	} else
		SPRINTF(buff, "        Single thread mode\n");

	outLen += outString(buff, outFile);

	if (nThreads >= 1 && CANON_ON_GPU) {
		SPRINTF(buff, "        Canonicity was tested on GPU by %zd checkers (%d for each thread)\n", NUM_GPU_WORKERS * nThreads, NUM_GPU_WORKERS);
		outLen += outString(buff, outFile);
	}

	SPRINTF(buff, "        Canonicity of partial constructed matrix was %sused\n", USE_CANON_GROUP ? "" : "not ");
	outLen += outString(buff, outFile);

	SPRINTF(buff, "        Strong canonicity was %sused\n", USE_STRONG_CANONICITY ? "" : "not ");
	outLen += outString(buff, outFile);

	SPRINTF(buff, "        Super strong canonicity was %sused\n\n", USE_STRONG_CANONICITY_A ? "" : "not ");
	outLen += outString(buff, outFile);
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

#if CANON_ON_GPU
FClass2(CEnumInfo, void)::RecalcCountersByGroupOrders(const COrderInfo* pOrderInfo, size_t nElem) {
	// Recalculating counters by the group order info
	updateGroupInfo(pOrderInfo, nElem);

	COrderInfo total;
	calcCountersTotal(&total);
	init();
	resetNumbInfo(1);
	const ulonglong nSimple = total.numMatrOfType(t_simple);
	addMatrix(total.numMatrOfType(t_canonical), nSimple);
	incNumbSimpleDesign(nSimple);
}
#endif

FClass2(CInsSysEnumInfo, void)::updateEnumInfo(const EnumInfoPntr pInfo)
{
	Class2(CEnumInfo)::updateEnumInfo(pInfo);
	incNumbSimpleDesign(pInfo->numbSimpleDesign());
}

FClass2(CInsSysEnumInfo, void)::reportResult(char *buffer, int lenBuffer) const
{
	size_t len = SNPRINTF(buffer, lenBuffer, "%s       %9llu       %9llu",
						   this->strToScreen(), this->constrCanonical(), numbSimpleDesign());
	len += SNPRINTF(buffer + len, lenBuffer - len, "       %8llu    ", this->numMatrOfType(t_totalConstr));
	len += this->convertTime(this->runTime(), buffer + len, lenBuffer - len);

	// Prepare the comments regarding the results
	const char *pResComment = NULL;
	switch (this->getResType()) {
	case t_resNew:			pResComment = "N"; break;
	case t_resBetter:		pResComment = "B"; break;
	case t_resWorse:		pResComment = "W"; break;
	case t_resInconsistent:	pResComment = "???";
	}

	SNPRINTF(buffer + len, lenBuffer - len, "  %s       \n", pResComment);
}

