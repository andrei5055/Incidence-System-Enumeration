#pragma once
#include <string>
#include "Global.h"
#include "Stat.h"
#include "GroupOrder.h"
#include "Storage.h"
#include "TrRepo.h"

typedef enum {
	eNoErrorCheck,
	eCheckErrors,
} eCheckForErrors;

class SizeParam {
public:
	CC SizeParam(const kSysParam& p) :
		m_numPlayers(p.val[t_numPlayers]),
		m_groupSize(p.val[t_groupSize]),
		m_use2RowsCanonization(p.val[t_use2RowsCanonization]),
		m_groupSizeFactorial(p.groupSizeFactorial),
		m_numDays(p.numFactors()),
		m_allowUndefinedCycles(p.val[t_allowUndefinedCycles]) {}
protected:
	CC void linksFromMatrix(tchar* lnk, ctchar* iv, int nr) const;
	CC bool setLinksForOnePlayer(tchar id, tchar* lnk, tchar* p, int ip, tchar v) const;
	void createFolderAndFileName(std::string& fnn, const kSysParam* param, int tFolder, int nr, const std::string& fName = "") const;
	const int m_numPlayers;
	const int m_groupSize;
	const int m_use2RowsCanonization;
	const int m_groupSizeFactorial;
	const int m_numDays;
	const int m_allowUndefinedCycles;
};

#if ReportCheckLinksData
#define CChecklLinkBase		CReportCheckLinksData
#define PrepareIDs()		m_cnt++;	\
							const auto idx = id * m_numPlayers;	\
							auto* counts_id = m_counts + idx;	\
							auto* faults_id = m_faults + idx
#define UpdateCntOK(ret)	if (ret) m_cntOk++;
#define SetLinksFault(i)	faults_id[i] |= 1; counts_id[i]++; m_cntErr++
#define UpdateFaultsID(i)	faults_id[i] |= 2;
#else
#define CChecklLinkBase		SizeParam
#define PrepareIDs()
#define UpdateCntOK(ret)
#define SetLinksFault(i)
#define UpdateFaultsID(i)
#endif

#if ReportCheckLinksData
class CReportCheckLinksData : protected SizeParam {
protected:
	CReportCheckLinksData(const SizeParam& sizeParam);
	~CReportCheckLinksData();
	void reportCheckLinksData();
	bool checkLinksT(ctchar* pLinks, int id, bool printLinksStatTime = false);
private:
	bool checkLinksTR(ctchar* v, int nvAll, int nv, int ind) const;
	int getAllUnlinked(int ic, tchar* v, int nvMax) const;
protected:
	sLongLong m_cnt = 0;
	sLongLong m_tmtotal = 0;
	sLongLong m_tmtotalFalse = 0;
	sLongLong m_tmtotalOk = 0;
	sLongLong m_cntErr = 0;
	sLongLong m_cntOk = 0;
	sLongLong* m_tmok = NULL;
	sLongLong* m_tmfalse = NULL;
	sLongLong* m_counts = NULL;
	tchar* m_pLinksCopy = NULL;
	tchar* m_faults = NULL;
};
#endif

class CChecklLink : protected CChecklLinkBase {
public:
	CC CChecklLink(const SizeParam& sizeParam, const kSysParam* p);
	CC ~CChecklLink();
	CC bool checkLinks(tchar* c, int id, bool printLinksStatTime = false);
	CC bool cycleLengthOk(tchar length) const;
	CC bool cyclesNotOk(int ncycles, tchar* cycles, eCheckForErrors checkErrors) const;
	CC bool checkLinksV2(ctchar* lnks, int nr) const;
	CC inline int param(paramID id)	const				{ return m_param->val[id]; }
	CC inline auto completeGraph() const				{ return m_param->completeGraph();}
protected:
	CC void u1fSetTableRow(tchar* ro, ctchar* ri, bool bNeighbors = false) const;
	CC inline auto* neighbors(int nDay = 0) const		{ return m_pU1Ftable + nDay * m_numPlayers; }
	CC inline void set_kSysParam(const kSysParam* p)	{ m_param = p; }
	CC inline const kSysParam* sysParam() const			{ return m_param; }
private:
	CC bool checkLinksV(ctchar* links, ctchar* v, int nv, int ind, tchar* vo) const;

	tchar* m_remainder3;
	tchar* m_pU1Ftable;
	tchar* m_pLinksCopy = NULL;
	tchar* m_v = NULL;
	tchar* m_vo = NULL;
	const kSysParam* m_param;
#if PrintNVminmax
	char* m_nvmn = NULL;
	char* m_nvmx = NULL;
	void setNV_MinMax(int id, int idx, char nv);
#else
#define setNV_MinMax(id, idx, nv)   // empty macros
#endif
};

template<typename T>
class CRepository : public CStorageIdx<T> {
public:
	CC CRepository(int lenObj, int numObjects = 0) :
		CStorageIdx<T>(numObjects, lenObj) {}
};

typedef CRepository<tchar> CGroupInfo;
#define incGroupOrder()		m_numObjects++
#define resetGroupOrder()	m_numObjects = 0
#define orderOfGroup()		numObjects()
#define updateGroup(x)		updateRepo(x)

class CGroupUtilisation {
public:
	CC inline CGroupInfo* groupInfo(int nRow) const {
		const int idx = nRow - m_autLevelDef[0];
		if (idx < 0 || nRow > m_autLevelDef[1])
			return NULL;

		return m_ppAutGroups[idx];
	}
protected:
	CC CGroupUtilisation(const kSysParam* p) : 
		m_autLevelDef{ p->val[t_autLevelMinDef], p->val[t_autLevelMaxDef] },
		m_autLevel{ p->val[t_autLevelMin], p->val[t_autLevelMax] },
		m_bDirection(p->val[t_autDirection] == 0) {
		const auto& val = p->val;
		const auto numPlayers = val[t_numPlayers];
		const auto nRows = p->numFactors();
		if ((val[t_printMatrices] & 16) || val[t_autSaveTestedTrs])
			m_pTrRepo = new CTrRepo(nRows, numPlayers);

		m_numLevels = m_autLevelDef[1] - m_autLevelDef[0] + 1;
		if (m_numLevels <= 0)
			m_numLevels = 1;

		m_numTrGroups = nRows * nRows;
		m_ppAutGroups = new CGroupInfo* [m_numLevels];

		for (int i = 0; i < m_numLevels; i++)
			m_ppAutGroups[i] = new CGroupInfo(numPlayers);
	}
	CC ~CGroupUtilisation() {
		if (m_ppAutGroups) {
			for (int i = 0; i < m_numLevels; i++)
				delete m_ppAutGroups[i];

			delete[] m_ppAutGroups;
		}
		delete m_pTrRepo;
	}

	CC void saveGroup(const CGroupInfo& grInfo, int nRows) {
		auto* pRowGroup = groupInfo(nRows);
		if (pRowGroup)
			*pRowGroup = grInfo;
	}
	
	CC inline int setupIteratorByGroups(int* pLastVal, int* pStep) {
		if (m_bDirection) {
			*pStep = 1; *pLastVal = m_autLevelDef[1];
			return m_autLevelDef[0] - 1;
		}
		*pStep = -1; *pLastVal = m_autLevelDef[0];
		return m_autLevelDef[1] + 1;
	}
	CC inline auto rowGroup(int nRow) const			{ return m_ppAutGroups[nRow - m_autLevelDef[0]]; }
	CC inline auto utilizeGroups(int nRow) const	{ return m_autLevel[0] <= nRow && nRow <= m_autLevel[1]; }
	CC inline auto firstGroupIdx() const			{ return m_autLevelDef[0]; }
	CC inline auto lastGroupIdx() const				{ return m_autLevelDef[1]; }

	CTrRepo* m_pTrRepo = NULL;
private:
	CGroupInfo** m_ppAutGroups = NULL;
	const int m_autLevelDef[2];
	const int m_autLevel[2];
	const bool m_bDirection;
	int m_numLevels;
	int m_numTrGroups;
};


