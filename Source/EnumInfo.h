#pragma once
#include <time.h>
#include "GroupsInfo.h"
#include "TimerInfo.h"

#define BEG_OUT_BLOCK			"<<<< "		// Marks for beginning and end of the output info, which
#define END_OUT_BLOCK			">>>> "		// will be skipped during the comparison of output results
#define ONE_LINE_BLOCK			"#### "		// To skip just one line
#define OF_THEM		            "of them"   // Marker of supporting information to skip

enum class t_reportCriteria {
	t_reportNow,
	t_reportByTime,
	t_matrConstructed,
	t_treadEnded
};

enum class t_resType {
	t_resNew,
	t_resBetter,
	t_resWorse,
	t_resInconsistent,
	t_resPostponed
};


Class2Def(CThreadEnumerator);
Class2Def(CEnumerator);

struct WorkingInfo {
	WorkingInfo(int* pntr, int min, int max) : pNumThreadsOnRow(pntr), minRow(min), maxRow(max) {}
	int* pNumThreadsOnRow;
	const int minRow;
	const int maxRow;
};

Class2Def(CEnumInfo) : public CTimerInfo, public CGroupsInfo, public CNumbInfo
{
public:
	CC CEnumInfo(const char *pStrToScreen = NULL) {
		m_pStrToScreen = pStrToScreen? new std::string(pStrToScreen) : NULL;
        resetCounters();
		setReportBound(REPORT_INTERVAL_OBJ_NUMB);
		m_pReportFileName = NULL;
		m_lenPrev = 0;
		m_bReportProgress = !PRINT_SOLUTIONS && m_pStrToScreen;
	}
	CC virtual ~CEnumInfo()									{ delete[] reportFileName();
															  delete m_pStrToScreen;
															}
public:
	CK inline void incrConstrTotal(ulonglong val = 1)		{ addMatrOfType(val, t_design_type::t_totalConstr); }
	inline const char *strToScreen() const					{ return m_pStrToScreen? m_pStrToScreen->c_str() : NULL; }
#if !CONSTR_ON_GPU
	void reportProgress(t_reportCriteria reportType, const CGroupsInfo *pGroupInfo = NULL, const WorkingInfo* pWorkInfo = NULL);
	void reportProgress(Class2(CThreadEnumerator)** ppTreadEnum, size_t nThread = 0);
#endif
	virtual ulonglong numbSimpleDesign() const				{ return 0; }
	CC virtual void incNumbSimpleDesign(ulonglong v = 1)	{}
	virtual void reportResult(char *buffer, int lenBuffer) const {}
	CC virtual void init()									{ resetCounters(); }
	virtual void updateEnumInfo(const CEnumInfo * pInfo);
	CC virtual void setNoReplBlockFlag(bool val)			{}
	CC virtual bool constructedAllNoReplBlockMatrix() const	{ return false; }
	CC void updateConstrCounters(int matrFlag, const EnumeratorPntr pEnum);
#if CANON_ON_GPU
	void RecalcCountersByGroupOrders(const COrderInfo* pOrderInfo, size_t nElem);
#endif
	void setReportFileName(const char *pntr);
	CK void outEnumInfo(FILE **pOutFile, bool removeReportFile = true, const CGroupsInfo *pGroupInfo = NULL, const char* pComment = NULL);
	CK inline void setResType(t_resType resType)			{ m_nResType = resType; }
	CC inline void resetEnumInfo()							{ init(); resetGroupsInfo(); }
	void outEnumInformation(FILE** pOutFile, const uint enumInfo, bool printMTlevel = true, const char* pComment = NULL) const;
	inline void setDesignInfo(designParam *pParam)			{ m_pParam = pParam; }
	void updateCounters(CEnumInfo* p);
	void outRunTimeInfo(FILE* outFile, const char* pOutString = NULL) const;
protected:
	inline t_resType getResType() const						{ return m_nResType; }
	size_t cleanEndOfLine(char* pBuffer) const;
private:
	CC inline void incrConstrCanonical(ulonglong val = 1)	{ addMatrOfType(val, t_design_type::t_canonical); }
	void outAdditionalInfo(ulonglong nMatr, FILE* outFile, char* buff, size_t lenBuf) const;
	CC inline void resetCounters()							{ m_prevReportCounter[0] = m_prevReportCounter[1] = m_nCounter = 0; }
	inline void incCounter()								{ m_nCounter++; }
	inline void reportThreadProgress()						{ incCounter(); reportProgress(t_reportCriteria::t_treadEnded); }
	inline auto prevReportCounter(int idx = 0) const		{ return m_prevReportCounter[idx]; }
	inline void setPrevReportCounter(ulonglong nCanon, int idx = 0)	{ m_prevReportCounter[idx] = nCanon; }
	CC virtual void setNumbSimpleDesign(ulonglong v)		{}
	CC inline void setReportBound(ulonglong val)			{ m_nReportBound = val; }
	inline ulonglong reportBound() const					{ return m_nReportBound; }
	CC inline char *reportFileName() const					{ return m_pReportFileName; }
	inline int multiThreadLevel() const 					{ return m_pParam->MT_level(); }
	inline designParam *designInfo() const					{ return m_pParam; }

	ulonglong m_nCounter;
	ulonglong m_prevReportCounter[2];
	std::string *m_pStrToScreen;
	ulonglong m_nReportBound;
	char *m_pReportFileName;
	t_resType m_nResType;
	int m_mtlevel;
	mutable size_t m_lenPrev;			// Length of previous output
	designParam *m_pParam;
	bool m_bReportProgress;
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