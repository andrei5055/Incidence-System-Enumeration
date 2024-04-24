#include <float.h>
#include "TimerInfo.h"

#undef MAC
#undef _MAC
#include <rpc.h>

static const int outDiv[] = { ':', ':', ':', '.' };

#if defined(_WIN32)
void CTimerInfo::get_cpu_usage(double& user, double& system) {
	FILETIME creation_time, exit_time, kernel_time, user_time;
	if (GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time, &kernel_time, &user_time)) {
		ULARGE_INTEGER kernel, user_t;
		kernel.HighPart = kernel_time.dwHighDateTime;
		kernel.LowPart = kernel_time.dwLowDateTime;
		user_t.HighPart = user_time.dwHighDateTime;
		user_t.LowPart = user_time.dwLowDateTime;
		user = (static_cast<double>(kernel.QuadPart) +
			static_cast<double>(user_t.QuadPart)) * 1e-7;
		system = static_cast<double>(kernel.QuadPart) * 1e-7;
	}
	else
		user = system = DBL_MAX;
}
#else
double tv_to_double(struct timeval& tv) {
	return tv.tv_sec + tv.tv_usec / 1e6;
}

void CTimerInfo::get_cpu_usage(double& user, double& system) {
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	user = tv_to_double(usage.ru_utime);
	system = tv_to_double(usage.ru_stime);
}
#endif

void CTimerInfo::startClock() {
	setPrevClock(m_prevClockReport = m_startClock = clock());
	setReportInt(1);
	setPrevCounter(0);

	// Save initial values for user and system times
	get_cpu_usage(m_ProcTime[2], m_ProcTime[3]);
}

void CTimerInfo::setRunTime() {
	m_fRunTime = (double)(clock() - startTime()) / CLOCKS_PER_SEC;
	double user_time, system_time;
	get_cpu_usage(user_time, system_time);
	m_ProcTime[0] = user_time - m_ProcTime[2];
	m_ProcTime[1] = system_time - m_ProcTime[3];
}

size_t CTimerInfo::convertTime(double time, char* buffer, size_t lenBuf, bool alignment) const
{
	const int secTime = (int)time;
	const int minTime = secTime / 60;
	const int hourTime = minTime / 60;
	const int out[] = { hourTime / 24, hourTime % 24, minTime % 60, secTime % 60 };

	bool flag = false;
	char* pBuf = buffer;
	for (int i = 0; i < countof(outDiv); i++) {
		if (out[i] > 0 || flag || outDiv[i] == '.') {
			pBuf += SNPRINTF(pBuf, lenBuf, (flag ? "%02d%c" : "%2d%c"), out[i], outDiv[i]);
			flag = true;
			continue;
		}

		if (alignment) {
			pBuf += SNPRINTF(pBuf, lenBuf, "   ");
			lenBuf -= 3;
		}
	}

	pBuf += SNPRINTF(pBuf, lenBuf, "%02d sec%s", (int)(100 * (time - secTime)), alignment ? "" : ",\n");
	return pBuf - buffer;
}

double CTimerInfo::stringToTime(char* pTime)
{
	const double mult[] = { 3600 * 24, 3600, 60, 1, 0.01 };

	double time = 0;
	size_t i = countof(outDiv);
	char* pDivPlace[countof(outDiv)] = {};
	char* pDiv = NULL;
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

bool CTimerInfo::compareTime(char* pTime1, char* pTime2)
{
	const double time1 = stringToTime(pTime1);
	const double time2 = stringToTime(pTime2);
	return time1 < time2;
}

