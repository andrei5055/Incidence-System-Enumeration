#pragma once
#include <string>
#include "Global.h"
#include "Stat.h"
#include "GroupOrder.h"

class SizeParam {
public:
	CC SizeParam(const kSysParam& p) :
		m_numPlayers(p.val[t_numPlayers]),
		m_groupSize(p.val[t_groupSize]),
		m_p1f(p.val[t_p1f]),
		m_groupSizeFactorial(p.groupSizeFactorial),
		m_numDays((m_numPlayers - 1) / (m_groupSize - 1)) {}
protected:
	CC void linksFromMatrix(tchar* lnk, tchar* iv, int nr) const;
	CC bool setLinksForOnePlayer(tchar id, tchar* lnk, tchar* p, int ip, tchar v) const;
	void createFolderAndFileName(std::string& fnn, const kSysParam* param, int tFolder, int nr, const std::string* fName = NULL) const;
	const int m_numPlayers;
	const int m_groupSize;
	const int m_p1f;
	const int m_groupSizeFactorial;
	const int m_numDays;
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
	CC bool cycleLengthOk(tchar length);
	CC bool cyclesNotOk(int ncr, int ncycles, tchar* cycles);
protected:
	CC int p1fCheck(int nr, ctchar* row) const;
	CC void p1fSetTableRow(tchar* ro, ctchar* ri) const;
	CC inline auto* p1ftable(int nDay = 0) const { return m_pP1Ftable + nDay * m_numPlayers; }
	CC inline void set_kSysParam(const kSysParam* p) { m_param = p; }
	CC inline int param(paramID id)	const { return m_param->val[id]; }
	CC inline const kSysParam* sysParam() const { return m_param; }
private:
	CC int u1fCheckFunc(const int nr, ctchar* rowm) const;
	CC bool checkLinksV(ctchar* links, ctchar* v, int nv, int ind, tchar* vo) const;

	tchar* m_pP1Ftable;
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
class CStorage {
public:
	CC CStorage(int numObjects, int lenObj = 1) : m_lenObj(lenObj) {
		if (numObjects)
			m_pObjects = new T[numObjects * m_lenObj]; 
	}
	CC ~CStorage()				{ delete[] m_pObjects; }
	CC T *reallocStorageMemory(int numObjects) {
		return ::reallocStorageMemory(&m_pObjects, m_lenObj * numObjects);
	}
	CC inline T* getObject(int idx = 0) const	{ return m_pObjects + idx * m_lenObj; }
	CC inline T** getObjectsPntr()				{ return &m_pObjects; }
protected:
	const int m_lenObj;
private:
	T* m_pObjects = NULL;
};

template<typename T>
class CStorageIdx : public CStorage<T> {
public:
	CC CStorageIdx(int numObjects, int lenObj = 1) : CStorage<T>(numObjects, lenObj) {
		if (numObjects)
			m_pIdx = new int[m_numObjectsMax = numObjects];
	}
	CC ~CStorageIdx()					{ delete [] m_pIdx; }
	CC inline T* getObjPntr(int idx) const { return CStorage<T>::getObject(m_pIdx[idx]); }
	CC inline auto numObjects() const	{ return m_numObjects; }
	CC T* reallocStorageMemory() {
		m_numObjectsMax <<= 1;
		::reallocStorageMemory(&m_pIdx, m_numObjectsMax * sizeof(m_pIdx[0]));
		return CStorage<T>::reallocStorageMemory(m_numObjectsMax);
	}
	CC void push_back(int idx)			{ m_pIdx[m_numObjects++] = idx; }
	CC void insert(int idx, int value)  {
		const auto len = (m_numObjects++ - idx) * sizeof(m_pIdx[0]);
		const auto pSrc = m_pIdx + idx;
		MEMMOVE((void *)(pSrc + 1), pSrc, len);
		m_pIdx[idx] = value;
	}
	CC T *getObjAddr(int idx)			{ return idx < m_numObjectsMax ? CStorage<T>::getObject(idx) : reallocStorageMemory(); }
	CC inline int* getIndices() const	{ return m_pIdx; }
	CC CStorageIdx& operator = (const CStorageIdx& other) {
		auto** ppStorage = CStorage<T>::getObjectsPntr();
		delete[] *ppStorage;
		delete[] m_pIdx;
		m_pIdx = NULL;	// We don't need indices for now 	
		m_numObjectsMax = m_numObjects = other.numObjects();
		const auto len = m_numObjects * CStorage<T>::m_lenObj;
		*ppStorage = new T[len];
		memcpy(*ppStorage, other.getObject(), len);
		return *this;
	}
	CC void copyIndex(const CStorageIdx& other) {
		m_pIdx = new int[other.numObjects()];
		memcpy(m_pIdx, other.getIndices(), m_numObjects * sizeof(m_pIdx[0]));
	}

	int m_numObjects = 0;
protected:
	int m_numObjectsMax = 0;
private:
	int* m_pIdx = NULL;
};

class CGroupInfo : public CStorageIdx<tchar> {
public:
	CC CGroupInfo(int numPlayers, int nPerm = 0) : 
		CStorageIdx<tchar>(nPerm, numPlayers) {}

	CC int updateGroupOrder(ctchar* tr);
	CC inline auto groupOrder()		{ return numObjects(); }
	CC inline auto isProcessed(ctchar* tr) {
		const auto numRegisteredTrs = numObjects();
		updateGroupOrder(tr);
		return numRegisteredTrs == numObjects();
	}

#define incGroupOrder()		m_numObjects++
#define resetGroupOrder()	m_numObjects = 0
	int m_numDaysResult;
};

class CGroupUtilisation {
protected:
	CC CGroupUtilisation(const kSysParam* p) : 
		m_autLevelDef{ p->val[t_autLevelMinDef], p->val[t_autLevelMaxDef] },
		m_autLevel{ p->val[t_autLevelMin], p->val[t_autLevelMax] },
		m_bDirection(p->val[t_autDirection] == 0) {
		m_numLevels = m_autLevelDef[1] - m_autLevelDef[0] + 1;
		if (m_numLevels <= 0)
			return;

		const auto numPlayers = p->val[t_numPlayers];
		m_ppAutGroups = new CGroupInfo* [m_numLevels];
		for (int i = 0; i < m_numLevels; i++)
			m_ppAutGroups[i] = new CGroupInfo(numPlayers);

		if (p->val[t_autSaveTestedTrs])
			m_pTestedTrs = new CGroupInfo(numPlayers, 32000);
	}
	CC ~CGroupUtilisation() {
		for (int i = 0; i < m_numLevels; i++)
			delete m_ppAutGroups[i];

		delete[] m_ppAutGroups;
		delete m_pTestedTrs;
	}

	CC inline CGroupInfo* groupInfo(int nRow) const {
		const int idx = nRow - m_autLevelDef[0];
		if (idx < 0 || nRow > m_autLevelDef[1])
			return NULL;

		return m_ppAutGroups[idx];
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
	CC inline auto testedTRs() const				{ return m_pTestedTrs; }
private:
	CGroupInfo* m_pTestedTrs = NULL;
	CGroupInfo** m_ppAutGroups = NULL;
	const int m_autLevelDef[2];
	const int m_autLevel[2];
	const bool m_bDirection;
	int m_numLevels;
};
