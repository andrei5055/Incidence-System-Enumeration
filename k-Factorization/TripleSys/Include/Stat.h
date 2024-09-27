#pragma once

#define ReportCheckLinksData	0  && !USE_CUDA // print information about CheckLinks execution
#define PrintNVminmax 0 && !USE_CUDA
#define ReportInterval			120000
#define Report					1			// to enable and disable reporting
#define ReportPeriodically		Report && 1 // print current matrix if "current number of rows" changed, or after ReportInterval ms timeout
#define ReportStatPeriodically	Report && 0 // print stat when "ReportPeriodically" triggered
#define ReportAfterAllTr		Report && 0 // print stat after canonization finished
#define ReportAfterEachResult	Report && 0 // print stat after each result
#define ReportAfterThreadEnd	Report && 0 // print stat just before thread exit

#define GetStat_cnvCheckKm1		0 // get stat for result of matrix comparison in canonization 
#define GetStat_checkLinksT		0
#define GetStat_p1fCheck		0 // implemented only for group size = 2
//-------  all user settings are above this line -------

#define ResetStat true
#define KeepStat false
#if !USE_CUDA
void StatAdd(const char* t, int ind, bool bAdd);
void StatEnd(bool bReset, const char* cHdr, int iHdr, bool bPrint);
#else
#define StatAdd(...)
#define StatEnd(...)
#endif

#if GetStat_cnvCheckKm1
#define Stat_cnvCheckKm1 StatAdd
#else
#define Stat_cnvCheckKm1(a, b, c) 
#endif
#if GetStat_checkLinksT
#define Stat_checkLinksT StatAdd
#else
#define Stat_checkLinksT(a, b, c) 
#endif
#if GetStat_p1fCheck
#define Stat_p1fCheck StatAdd
#else
#define Stat_p1fCheck(a, b, c) 
#endif

#define StatReport(a, b, c, d) StatEnd(a, b, c, d)

#if ReportStatPeriodically
#define StatReportPeriodically(a, b, c, d) StatReport(a, b, c, d)
#else
#define StatReportPeriodically(a, b, c, d)
#endif
#if ReportAfterAllTr
#define StatReportAfterAllTr(a, b, c, d) StatReport(a, b, c, d)
#else
#define StatReportAfterAllTr(a, b, c, d)
#endif
#if ReportAfterEachResult
#define StatReportAfterEachResult(a, b, c, d) StatReport(a, b, c, d)
#else
#define StatReportAfterEachResult(a, b, c, d)
#endif
#if ReportAfterThreadEnd
#define StatReportAfterThreadEnd(a, b, c, d) StatReport(a, b, c, d)
#else
#define StatReportAfterThreadEnd(a, b, c, d)
#endif


