#pragma once
#include <time.h>
#include "CudaInfo.h"

class CTimerInfo
{
public:
	CC inline CTimerInfo() {}
	CC inline ~CTimerInfo() {}
	CK void startClock();
	CK void setRunTime();
	inline auto runTime() const { return m_fRunTime; }
	inline auto procTime(int idx) const { return m_ProcTime[idx]; }
	size_t convertTime(double time, char* buffer, size_t lenBuf, bool alignment = true) const;
	static bool compareTime(char* pTime1, char* pTime2);
protected:
	inline clock_t prevClock() const { return m_prevClock; }
	CK inline void setPrevClock(const clock_t& val) { m_prevClock = val; }
	inline clock_t prevClockReport() const { return m_prevClockReport; }
	CK inline void setPrevClockReport(const clock_t& val) { m_prevClockReport = val; }
	inline size_t reportInt() const { return m_reportInt; }
	CK inline void setReportInt(size_t val) { m_reportInt = val; }
	inline ulonglong prevCounter() const { return m_prevCounter; }
	CK inline void setPrevCounter(ulonglong value) { m_prevCounter = value; }
	inline clock_t startTime() const { return m_startClock; }
private:
	CK void get_cpu_usage(double& user, double& system);
	static double stringToTime(char* pTime);

	clock_t m_startClock = 0;
	clock_t m_prevClock = 0;
	clock_t m_prevClockReport = 0;
	size_t m_reportInt = 0;
	ulonglong m_prevCounter = 0;
	double m_fRunTime = 0.0f;
	double m_ProcTime[4] = { 0,0,0,0 };
};
