#pragma once
#include <time.h>
#include "GroupsInfo.h"

class CTimerInfo
{
public:
	CC inline CTimerInfo()									{}
	CC inline ~CTimerInfo()									{}
	CK inline void startClock()								{ setPrevClock(m_prevClockReport = m_startClock = clock());  setReportInt(1);  setPrevCounter(0); }
	inline void setRunTime()								{ m_fRunTime = (float)(clock() - startTime()) / CLOCKS_PER_SEC; }
	inline float runTime() const							{ return m_fRunTime; }
protected:
	inline clock_t prevClock() const						{ return m_prevClock; }
	CK inline void setPrevClock(clock_t &val)				{ m_prevClock = val; }
	inline clock_t prevClockReport() const					{ return m_prevClockReport; }
	CK inline void setPrevClockReport(clock_t &val)			{ m_prevClockReport = val; }
	inline size_t reportInt() const							{ return m_reportInt; }
	CK inline void setReportInt(size_t val)					{ m_reportInt = val; }
	inline ulonglong prevCounter() const					{ return m_prevCounter; }
	CK inline void setPrevCounter(ulonglong value)			{ m_prevCounter = value; }
	inline clock_t startTime() const						{ return m_startClock; }
private:

	clock_t m_startClock;
	clock_t m_prevClock;
	clock_t m_prevClockReport;
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

typedef enum {
	t_resNew,
	t_resBetter,
	t_resWorse,
	t_resInconsistent
} t_resType;

typedef enum {
	t_Summary =		  1<<0,		// default
	t_AllObject =	  1<<1,
	t_GroupOrbits =			1<<2,
	t_GroupGeneratingSet = t_GroupOrbits + (1<<3),
	t_Transitive =			1<<4,
	t_GroupOrderEQ =		1<<5,
	t_GroupOrderGT =		1<<6,
	t_GroupOrderLT =		1<<7,
} t_outputType;

typedef enum {
	t_noReplicatedBlock = 1 << 1,
	t_transitiveGroup   = 1 << 2,
} t_MatrixFlags;

Class2Def(CThreadEnumerator);

Class2Def(CEnumInfo) : public CTimerInfo, public CGroupsInfo, public CNumbInfo
{
public:
	CC CEnumInfo(const char *pStrToScreen = NULL) : m_pStrToScreen(pStrToScreen) {
        resetCounters();
		setReportBound(REPORT_INTERVAL_OBJ_NUMB);
		m_pReportFileName = NULL;
	}
	CC virtual ~CEnumInfo()									{ delete[] reportFileName(); }
public:
	#define constrCanonical()	numMatrices()
	CK inline void incrConstrTotal(ulonglong val = 1)		{ addMatrOfType(val, t_totalConstr); }
	#define constrTotal()		numMatrOfType(t_totalConstr)
	inline const char *strToScreen() const					{ return m_pStrToScreen; }
#if !CONSTR_ON_GPU
	void reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo = NULL);
	void reportProgress(const ThreadEnumeratorPntr pTreadEnum, size_t nThread = 0);
#endif
	virtual ulonglong numbSimpleDesign() const				{ return 0; }
	CC virtual void incNumbSimpleDesign(ulonglong v = 1)	{}
	virtual void reportResult(char *buffer, int lenBuffer) const {}
	CC virtual void init()									{ resetCounters(); }
	virtual void updateEnumInfo(const CEnumInfo * pInfo);
	CC virtual void setNoReplBlockFlag(bool val)			{}
	CC virtual bool constructedAllNoReplBlockMatrix() const	{ return false; }
	CC void updateConstrCounters(int matrFlag, size_t groupOrder, bool groupIsTransitive) {
		incrConstrCanonical();
		const ulonglong simpleMatrFlag = matrFlag & t_noReplicatedBlock ? 1 : 0;
		if (simpleMatrFlag)
			incNumbSimpleDesign();

		COrderInfo *pOrderInfo = addGroupOrder(groupOrder, 1, simpleMatrFlag);
		if (groupIsTransitive)
			pOrderInfo->addMatrixTrans(1, simpleMatrFlag);
	}
	void RecalcCountersByGroupOrders(const COrderInfo *pOrderInfo, size_t nElem) {
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
	void setReportFileName(const char *pntr);
	CK void outEnumInfo(FILE **pOutFile, bool removeReportFile = true, const CGroupsInfo *pGroupInfo = NULL);
	void outEnumAdditionalInfo(FILE **pOutFile) const;
	size_t convertTime(float time, char *buffer, size_t lenBuf, bool alignment = true) const;
	CK void setResType(t_resType resType)					{ m_nResType = resType; }
	CC void resetEnumInfo()									{ init(); resetGroupsInfo(); }
	static bool compareTime(char *time1, char *time2);
	void outEnumInformation(FILE **pOutFile, bool printMTlevel = true) const;
	inline void setDesignInfo(designParam *pParam)			{ m_pParam = pParam; }
protected:
	t_resType getResType() const							{ return m_nResType; }
private:
	static double stringToTime(char *pTime);
	CC inline void incrConstrCanonical(ulonglong val = 1)	{ addMatrOfType(val, t_canonical); }
	CC inline void resetCounters()							{ m_nCounter = 0; }
	inline void incCounter()								{ m_nCounter++; }
	inline void reportThreadProgress()						{ incCounter(); reportProgress(t_treadEnded); }
	CC virtual void setNumbSimpleDesign(ulonglong v)		{}
	CC inline void setReportBound(ulonglong val)			{ m_nReportBound = val; }
	inline ulonglong reportBound() const					{ return m_nReportBound; }
	CC inline char *reportFileName() const					{ return m_pReportFileName; }
	inline int multiThreadLevel() const 					{ return m_pParam->mt_level; }
	inline designParam *designInfo() const					{ return m_pParam; }

	ulonglong m_nCounter;
	const char *m_pStrToScreen;
	ulonglong m_nReportBound;
	char *m_pReportFileName;
	t_resType m_nResType;
	int m_mtlevel;
	designParam *m_pParam;
}; 

Class2Def(CInsSysEnumInfo) : public Class2(CEnumInfo)
{
public:
	CC CInsSysEnumInfo(const char *pStrToScreen = NULL) : Class2(CEnumInfo)(pStrToScreen)
																	{ setNoReplBlockFlag(false); resetCounter(); }
	CC ~CInsSysEnumInfo()											{}
public:
	virtual ulonglong numbSimpleDesign() const						{ return m_nSimpleDesign; }
	CC virtual void incNumbSimpleDesign(ulonglong v = 1)			{ m_nSimpleDesign += v; }
	virtual void reportResult(char *buffer, int lenBuffer) const;
	CC virtual void init()											{ resetCounter();  Class2(CEnumInfo)::init(); }
	virtual void updateEnumInfo(const CEnumInfo * pInfo);
	CC virtual void setNoReplBlockFlag(bool val)					{ m_bNoReplBlockFlag = val; }
	CC virtual bool constructedAllNoReplBlockMatrix() const			{ return m_bNoReplBlockFlag; }
private:
	CC inline void resetCounter()									{ setNumbSimpleDesign(0); }
	CC virtual void setNumbSimpleDesign(ulonglong v)				{ m_nSimpleDesign = v; }

	ulonglong m_nSimpleDesign;
	bool m_bNoReplBlockFlag;
};