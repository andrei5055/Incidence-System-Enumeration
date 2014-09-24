#include "stdafx.h"
#include "DataTypes.h"
#include "EnumInfo.h"
#include "Enumerator.h"

CEnumInfo::CEnumInfo(const char *pStrToScreen) : m_pStrToScreen(pStrToScreen)
{
	resetCounters(); 
	setReportBound(REPORT_INTERVAL_OBJ_NUMB);
	m_pReportFileName = NULL;
}

void CEnumInfo::convertTime(float time, char *buffer, size_t lenBuf, bool alignment) const
{
	const int secTime = (int)time;
	const int minTime = secTime / 60;
	const int hourTime = minTime / 60;
	const int out[] = { hourTime / 24, hourTime % 24, minTime % 60, secTime % 60 };
	const int outDiv[] = { ':', ':', ':', '.' };

	bool flag = false;
	char *pBuf = buffer;
	for (int i = 0; i < sizeof(out) / sizeof(out[0]); i++) {
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

	sprintf_s(pBuf, lenBuf, "%02d\n", (int)(100 * (time - secTime)));
}

void CEnumInfo::updateEnumInfo(const CEnumInfo *pInfo)
{
	incrConstrCanonical(pInfo->constrCanonical());
	incrConstrTotal(pInfo->constrTotal());
}

void CEnumInfo::updateConstrCounters(const CEnumerator *pInum)
{
	incrConstrCanonical();
	if (simpleMatrFlag())
		incNumbSimpleDesign();

	const CCanonicityChecker *pCanon = pInum->canonChecker();
	const size_t groupOrder = pCanon->groupOrder();
	COrderInfo *pOrderInfo = addGroupOrder(groupOrder, 1, simpleMatrFlag() ? 1 : 0);
	if (pCanon->groupIsTransitive())
		pOrderInfo->addMatrixTrans(1, simpleMatrFlag() ? 1 : 0);
}

void CEnumInfo::setReportFileName(const char *pntr)
{ 
	delete [] reportFileName();
	if (pntr) {
		const size_t len = strlen(pntr) + 1;
		m_pReportFileName = new char[len];
		memcpy_s(m_pReportFileName, len, pntr, len);
	} else
		m_pReportFileName = NULL;
}

void CEnumInfo::reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo)
{
	ulonglong nMatr = constrCanonical();
	if (!strToScreen() || PRINT_SOLUTIONS)
		return;

	ulonglong *pTestNumber;
	switch (reportType) {
		case t_reportNow:		
		case t_reportByTime:	pTestNumber = NULL; break;
		case t_matrConstructed:	pTestNumber = &nMatr;
								break;
		case t_treadEnded:		pTestNumber = &m_nCounter;
	}

	if (!strToScreen() || (pTestNumber && *pTestNumber % reportInt()) || PRINT_SOLUTIONS)
		return;

	std::cout << '\xd' << strToScreen() << (reportType == t_reportByTime ? "==>" : "   ") << nMatr << "  NRB: " << (constructedAllNoReplBlockMatrix() ? "=" : "") << numbSimpleDesign() << "  Total: " << constrTotal();

	switch (reportType) {
		case t_matrConstructed:
		case t_reportNow:
		case t_treadEnded:	  pGroupInfo = this;
		case t_reportByTime:
			if (constrCanonical() >= reportBound()) {				
				FOPEN(file, reportFileName(), "w");
                if (file) {
                    fprintf(file, "%s\n", strToScreen());
                    outEnumInfo(file, false, pGroupInfo);
                    incReportBound(REPORT_INTERVAL_OBJ_NUMB * (constrCanonical() - reportBound() + REPORT_INTERVAL_OBJ_NUMB - 1) / REPORT_INTERVAL_OBJ_NUMB);
                }
			}
	}

	// Adjust reportInt to report approx. 1 time per 60 seconds
	clock_t currClock;
	if (!pTestNumber || (currClock = clock()) == prevClock())
		return;

	const float timeInterval = (float)(currClock - prevClock()) / CLOCKS_PER_SEC;
	const float averPermMin = 60 * (*pTestNumber - prevCounter()) / timeInterval;
	setPrevCounter(*pTestNumber);
	setPrevClock(currClock);
	size_t repInt = 1;
	while (averPermMin > 5 * repInt)
		repInt *= 10;

	setReportInt(repInt);
}

void CEnumInfo::reportProgress(const CThreadEnumerator *pThreadEnum, int nThread)
{
	if (nThread >= 1) {
		// Save already collected information 
		const ulonglong nCanon = constrCanonical();
		const ulonglong nTotal = constrTotal();
		const ulonglong nrbTotal = numbSimpleDesign();
		CGroupsInfo groupsInfo;
		groupsInfo.updateGroupInfo(this);
		for (int i = 0; i < nThread; i++, pThreadEnum++) {
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

void CEnumInfo::outEnumInfo(FILE *outFile, bool removeReportFile, const CGroupsInfo *pGroupInfo)
{
	if (!pGroupInfo)
		pGroupInfo = this;

	pGroupInfo->printGroupInfo(outFile);
	setRunTime();
	const ulonglong nConstrMatr = constrCanonical();
	char buff[256];
	SPRINTF(buff, "\n%10llu matri%s constructed in ", nConstrMatr, nConstrMatr == 1 ? "x" : "ces");
	const size_t len = strlen(buff);
	convertTime(runTime(), buff + len, countof(buff) - len, false);
	outString(buff, outFile);

	ulonglong nMatr = numbSimpleDesign();
	if (nConstrMatr > 0) {
		SPRINTF(buff, "%10llu matri%s ha%s no replicated blocks\n", nMatr, nMatr == 1 ? "x" : "ces", nConstrMatr == 1 ? "s" : "ve");
		outString(buff, outFile);
	}

	SPRINTF(buff, "%10llu matri%s fully constructed\n", constrTotal(), constrTotal() == 1 ? "x was" : "ces were");
	outString(buff, outFile);

	if (USE_THREADS >= 1)
		SPRINTF(buff, "%10d threads launched on level %d \n", USE_THREADS, MT_LEVEL);
	else
		SPRINTF(buff, "        Single thread mode\n");

	outString(buff, outFile);

	SPRINTF(buff, "        Strong canonisity was %sused\n", USE_STRONG_CANONICITY_A? "" : "not ");
	outString(buff, outFile);
	FCLOSE(outFile);

	if (removeReportFile) // Remove temporary file with the intermediate results	
		remove(reportFileName());
}

void CInsSysEnumInfo::updateEnumInfo(const CEnumInfo *pInfo)
{
	CEnumInfo::updateEnumInfo(pInfo);
	incNumbSimpleDesign(pInfo->numbSimpleDesign());
}

void CInsSysEnumInfo::reportResult(char *buffer, int lenBuffer) const
{
	int len = sprintf_s(buffer, lenBuffer, "%s       %9llu       %9llu", strToScreen(), constrCanonical(), numbSimpleDesign());
	len += sprintf_s(buffer + len, lenBuffer - len, "       %8llu    ", constrTotal());
	convertTime(runTime(), buffer + len, lenBuffer - len);
}

