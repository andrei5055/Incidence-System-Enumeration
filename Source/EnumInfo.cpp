#include "stdafx.h"
#include "DataTypes.h"
#include "EnumInfo.h"
#include "Enumerator.h"
#include "ThreadEnumerator.h"

template class CEnumInfo<MATRIX_ELEMENT_TYPE>;
template class CInsSysEnumInfo<MATRIX_ELEMENT_TYPE>;

static const int outDiv[] = { ':', ':', ':', '.' };

template<class T>
size_t CEnumInfo<T>::convertTime(float time, char *buffer, size_t lenBuf, bool alignment) const
{
	const int secTime = (int)time;
	const int minTime = secTime / 60;
	const int hourTime = minTime / 60;
	const int out[] = { hourTime / 24, hourTime % 24, minTime % 60, secTime % 60 };

	bool flag = false;
	char *pBuf = buffer;
	for (int i = 0; i < countof(outDiv); i++) {
		if (out[i] > 0 || flag || outDiv[i] == '.') {
			pBuf += sprintf_s(pBuf, lenBuf, (flag? "%02d%c" : "%2d%c"), out[i], outDiv[i]);
			flag = true;
		} else
		if (alignment)
			pBuf += sprintf_s(pBuf, lenBuf, "   ");
		else
			continue;

		lenBuf -= 3;
	}

	pBuf += sprintf_s(pBuf, lenBuf, "%02d%s", (int)(100 * (time - secTime)), alignment? "" : "\n");
	return pBuf - buffer;
}

template<class T>
double CEnumInfo<T>::stringToTime(char *pTime)
{
	const double mult[] = { 3600 * 24, 3600, 60, 1, 0.01 };

	double time = 0;
	size_t i = countof(outDiv);
	char *pDiv = NULL;
	while (i-- && (pDiv = strrchr(pTime, outDiv[i])) != NULL) {
		time += atoi(pDiv + 1) * mult[i + 1];
		*pDiv = '\0';
	}

	pDiv = strrchr(pTime, ' ');
	if (!pDiv)
		pDiv = pTime;
	else
		pDiv++;

	time += atoi(pDiv) * mult[i + 1];
	return time;
}

template<class T>
bool CEnumInfo<T>::compareTime(char *pTime1, char *pTime2)
{
	const double time1 = stringToTime(pTime1);
	const double time2 = stringToTime(pTime2);
	return time1 > time2;
}

template<class T>
void CEnumInfo<T>::updateEnumInfo(const CEnumInfo<T> *pInfo)
{
	incrConstrCanonical(pInfo->constrCanonical());
	incrConstrTotal(pInfo->constrTotal());
}

template<class T>
void CEnumInfo<T>::setReportFileName(const char *pntr)
{ 
	delete [] reportFileName();
	if (pntr) {
		const size_t len = strlen(pntr) + 1;
		m_pReportFileName = new char[len];
		memcpy_s(m_pReportFileName, len, pntr, len);
	} else
		m_pReportFileName = NULL;
}

template<class T>
void CEnumInfo<T>::reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo)
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
	std::cout << '\xd' << strToScreen() << (reportType == t_reportByTime ? "==>" : "   ") 
			  << "  Canon: " << nCanon 
			  << "  NRB: "   << (constructedAllNoReplBlockMatrix() ? "=" : "") << numbSimpleDesign() 
			  << "  Total: " << constrTotal()
			  << "  RunTime: " << runTime << " min.                 ";
			

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

template<class T>
void CEnumInfo<T>::reportProgress(const CThreadEnumerator<T> *pThreadEnum, size_t nThread)
{
	if (nThread >= 1) {
		// Save already collected information 
		const ulonglong nCanon = constrCanonical();
		const ulonglong nTotal = constrTotal();
		const ulonglong nrbTotal = numbSimpleDesign();
		CGroupsInfo groupsInfo;
		groupsInfo.updateGroupInfo(this);
		for (size_t i = 0; i < nThread; i++, pThreadEnum++) {
			if (pThreadEnum->code() != t_threadRunning)
				continue;

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

template<class T>
void CEnumInfo<T>::outEnumInfo(FILE **pOutFile, bool removeReportFile, const CGroupsInfo *pGroupInfo)
{
	setRunTime();
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	if (!pGroupInfo)
		pGroupInfo = this;

	pGroupInfo->printGroupInfo(outFile);
	const ulonglong nConstrMatr = constrCanonical();
	char buff[256];
	SPRINTF(buff, "\n%10llu matri%s"CONSTRUCTED_IN" ", nConstrMatr, nConstrMatr == 1 ? "x" : "ces");
	const size_t len = strlen(buff);
	convertTime(runTime(), buff + len, countof(buff) - len, false);
	outString(buff, outFile);

	const ulonglong nMatr = numbSimpleDesign();
	if (nConstrMatr > 0) {
		SPRINTF(buff, "%10llu matri%s ha%s no replicated blocks\n", nMatr, nMatr == 1 ? "x" : "ces", nMatr == 1 ? "s" : "ve");
		outString(buff, outFile);
	}

	SPRINTF(buff, "%10llu matri%s fully constructed\n", constrTotal(), constrTotal() == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile);
	if (removeReportFile) // Remove temporary file with the intermediate results	
		remove(reportFileName());
}

template<class T>
void CEnumInfo<T>::outEnumAdditionalInfo(FILE **pOutFile) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	char buff[256];
	SPRINTF(buff, "%10llu matri%s fully constructed\n", constrTotal(), constrTotal() == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	outEnumInformation(pOutFile);
}

template<class T>
void CEnumInfo<T>::outEnumInformation(FILE **pOutFile, bool printMTlevel) const
{
	FILE *outFile = pOutFile ? *pOutFile : NULL;
	if (!outFile)
		return;

	char buff[256];
	SPRINTF(buff, "\nUsing %Id-bit program, Assembly flag: %d\n", sizeof(size_t) << 3, USE_ASM);
	outString(buff, outFile);

	const size_t nThreads = designInfo()->threadNumb;
	if (USE_THREADS >= 1) {
		if (printMTlevel)
			SPRINTF(buff, "%10zu threads launched on level %d (%swaiting to finish mode)\n", nThreads, multiThreadLevel(), WAIT_THREADS ? "" : "not ");
		else
			SPRINTF(buff, "%10zu threads launched on difefrent levels (%swaiting to finish mode)\n", nThreads, WAIT_THREADS ? "" : "not ");
	} else
		SPRINTF(buff, "        Single thread mode\n");

	outString(buff, outFile);

	if (nThreads >= 1 && CANON_ON_GPU) {
		SPRINTF(buff, "        Canonicity was tested on GPU by %zd checkers (%d for each thread) (\n", NUM_GPU_WORKERS * nThreads, NUM_GPU_WORKERS);
		outString(buff, outFile);
	}

	SPRINTF(buff, "        Canonicity of partial constructed matrix was %sused\n", USE_CANON_GROUP ? "" : "not ");
	outString(buff, outFile);

	SPRINTF(buff, "        Strong canonicity was %sused\n", USE_STRONG_CANONICITY ? "" : "not ");
	outString(buff, outFile);

	SPRINTF(buff, "        Super strong canonicity was %sused\n\n", USE_STRONG_CANONICITY_A ? "" : "not ");
	outString(buff, outFile);
	FCLOSE(outFile);
	*pOutFile = NULL;
}

template<class T>
void CInsSysEnumInfo<T>::updateEnumInfo(const CEnumInfo *pInfo)
{
	CEnumInfo<T>::updateEnumInfo(pInfo);
	incNumbSimpleDesign(pInfo->numbSimpleDesign());
}

template<class T>
void CInsSysEnumInfo<T>::reportResult(char *buffer, int lenBuffer) const
{
	size_t len = sprintf_s(buffer, lenBuffer, "%s       %9llu       %9llu", strToScreen(), constrCanonical(), numbSimpleDesign());
	len += sprintf_s(buffer + len, lenBuffer - len, "       %8llu    ", constrTotal());
	len += convertTime(runTime(), buffer + len, lenBuffer - len);

	// Prepare the comments regarding the results
	char *pResComment = NULL;
	switch (getResType()) {
	case t_resNew:			pResComment = "N"; break;
	case t_resBetter:		pResComment = "B"; break;
	case t_resWorse:		pResComment = "W"; break;
	case t_resInconsistent:	pResComment = "???";
	}

	sprintf_s(buffer + len, lenBuffer - len, "  %s       \n", pResComment);
}

