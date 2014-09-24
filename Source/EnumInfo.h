#include <time.h>
#include "GroupsInfo.h"

class CTimerInfo
{
public:
	inline void startClock()								{ setPrevClock(m_startClock = clock());  setReportInt(1);  setPrevCounter(0); }
	inline void setRunTime()								{ m_fRunTime = (float)(clock() - startTime()) / CLOCKS_PER_SEC; }
	inline float runTime() const							{ return m_fRunTime; }
protected:
	inline clock_t prevClock() const						{ return m_prevClock; }
	inline void setPrevClock(clock_t val)					{ m_prevClock = val; }
	inline size_t reportInt() const							{ return m_reportInt; }
	inline void setReportInt(size_t val)					{ m_reportInt = val; }
	inline ulonglong prevCounter() const					{ return m_prevCounter; }
	inline void setPrevCounter(ulonglong value)				{ m_prevCounter = value; }
private:
	inline clock_t startTime() const						{ return m_startClock; }

	clock_t m_startClock;
	clock_t m_prevClock;
	size_t m_reportInt;
	ulonglong m_prevCounter;
	float m_fRunTime;
};

typedef enum {
	t_reportNow,
	t_reportByTime,
	t_matrConstructed,
	t_treadEnded
} t_reportCriteria;

class CThreadEnumerator;
class CEnumerator;

class CEnumInfo : public CTimerInfo, public CGroupsInfo, public CNumbInfo
{
public:
	CEnumInfo(const char *pStrToScreen = NULL);
	virtual ~CEnumInfo()									{ delete[] reportFileName(); }
public:
	#define constrCanonical()	numMatrices()
	inline void incrConstrTotal(ulonglong val = 1)			{ addMatrOfType(val, t_totalConstr); }
	#define constrTotal()		numMatrOfType(t_totalConstr)
	inline const char *strToScreen() const					{ return m_pStrToScreen; }
	void reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo = NULL);
	void reportProgress(const CThreadEnumerator *pTreadEnum, int nThread = 0);
	virtual ulonglong numbSimpleDesign() const				{ return 0; }
	virtual void reportResult(char *buffer, int lenBuffer) const {}
	virtual void init()										{ resetCounters(); }
	virtual void updateEnumInfo(const CEnumInfo *pInfo);
	virtual void setNoReplBlockFlag(bool val)				{}
	virtual bool constructedAllNoReplBlockMatrix() const	{ return false; }
	virtual void setSimpleMatrFlag(bool val)				{}
	virtual bool simpleMatrFlag() const						{ return false; }
	void updateConstrCounters(const CEnumerator *pInum);
	void setReportFileName(const char *pntr);
	void outEnumInfo(FILE *outFile, bool removeReportFile = true, const CGroupsInfo *pGroupInfo = NULL);
	void convertTime(float time, char *buffer, size_t lenBuf, bool alignment = true) const;
private:
	inline void incrConstrCanonical(ulonglong val = 1)		{ addMatrOfType(val, t_canonical); }
	virtual void incNumbSimpleDesign(ulonglong v = 1)		{}
	inline void resetCounters()								{ m_nCounter = 0; }
	inline void incCounter()								{ m_nCounter++; }
	inline void reportThreadProgress()						{ incCounter(); reportProgress(t_treadEnded); }
	virtual void setNumbSimpleDesign(ulonglong v)			{}
	inline void setReportBound(ulonglong val)				{ m_nReportBound = val; }
	inline void incReportBound(ulonglong val)				{ m_nReportBound += val; }
	inline ulonglong reportBound() const					{ return m_nReportBound; }
	inline char *reportFileName() const						{ return m_pReportFileName; }

	ulonglong m_nCounter;
	const char *m_pStrToScreen;
	ulonglong m_nReportBound;
	char *m_pReportFileName;
}; 

class CInsSysEnumInfo : public CEnumInfo
{
public:
	CInsSysEnumInfo(const char *pStrToScreen = NULL) : CEnumInfo(pStrToScreen)	{ setNoReplBlockFlag(false); resetCounter(); }
public:
	virtual ulonglong numbSimpleDesign() const						{ return m_nSimpleDesign; }
	virtual void reportResult(char *buffer, int lenBuffer) const;
	virtual void init()												{ resetCounter();  CEnumInfo::init(); }
	virtual void updateEnumInfo(const CEnumInfo *pInfo);
	virtual void setNoReplBlockFlag(bool val)						{ m_bNoReplBlockFlag = val; }
	virtual bool constructedAllNoReplBlockMatrix() const			{ return m_bNoReplBlockFlag; }
	virtual void setSimpleMatrFlag(bool val)						{ m_bSimpleMatrFlag = val; }
	virtual bool simpleMatrFlag() const								{ return m_bSimpleMatrFlag; }
private:
	virtual void incNumbSimpleDesign(ulonglong v = 1)				{ m_nSimpleDesign += v; }
	inline void resetCounter()										{ setNumbSimpleDesign(0); }
	virtual void setNumbSimpleDesign(ulonglong v)					{ m_nSimpleDesign = v; }

	ulonglong m_nSimpleDesign;
	bool m_bNoReplBlockFlag;
	bool m_bSimpleMatrFlag;
};