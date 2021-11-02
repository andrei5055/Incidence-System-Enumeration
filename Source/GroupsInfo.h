#pragma once
#include "DataTypes.h"

typedef enum {
	t_totalConstr,
	t_canonical,
	t_simple,
	t_transitive,
	t_simpleTrans,
	t_design_type_total
} t_design_type;

class CNumbInfo
{
public:
	CC CNumbInfo()												{ resetNumbInfo(); }
	CC ~CNumbInfo()												{}
	CK inline ulonglong numMatrOfType(t_design_type t)	const	{ return m_nCounter[t]; }
	inline void setNumMatrOfType(ulonglong val, t_design_type t){ m_nCounter[t] = val; }
	CC inline CNumbInfo *addMatrix(ulonglong numb, ulonglong nSimple) {
		addMatrOfType(numb, t_canonical);
		addMatrOfType(nSimple, t_simple);
		return this;
	}
	CC void addMatrixTrans(ulonglong numb, ulonglong nSimple) {
		addMatrOfType(numb, t_transitive);
		addMatrOfType(nSimple, t_simpleTrans);
	}

	void addMatrices(const CNumbInfo *pNumbInfo);
	void outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const;
	CC inline void resetNumbInfo()								{ memset(m_nCounter, 0, sizeof(m_nCounter)); }
protected:
	CC inline void addMatrOfType(ulonglong v, t_design_type t)	{ m_nCounter[t] += v; }
private:
	ulonglong m_nCounter[t_design_type_total];
};

class COrderNumb : public CNumbInfo
{
public:
	COrderNumb(size_t groupOrder=1) : m_groupOrder(groupOrder)	{}
	CC inline size_t groupOrder() const							{ return m_groupOrder; }
private:
	const size_t m_groupOrder;
};

#define NUMB_INFO_TYPE	  			COrderNumb *
#define NUMB_INFO_ACCESS_TYPE		NUMB_INFO_TYPE
typedef CArray<NUMB_INFO_TYPE, NUMB_INFO_ACCESS_TYPE> COrderNumbArray;

class COrderInfo
{
public:
	CC COrderInfo()												{ m_cNumbInfo.Add(new COrderNumb()); }
	CC COrderInfo(size_t order, size_t extraOrder, ulonglong num, ulonglong numSimple=0) : m_groupOrder(order) {
		m_cNumbInfo.Add(new COrderNumb(extraOrder));
		// Only one (just added) element of array exists, so no need to do search by extraOrder
		m_cNumbInfo[0]->addMatrix(num, numSimple);
	}
	CC ~COrderInfo() {
		for (auto i = m_cNumbInfo.GetSize(); i--;)
			delete m_cNumbInfo[i];
	}
	CC inline size_t groupOrder() const							{ return m_groupOrder; }
	CC inline CNumbInfo *addMatrix(size_t extraGroupOrder, ulonglong num, ulonglong nSimple) {
		return GetByKey(extraGroupOrder)->addMatrix(num, nSimple);
	}
	CC inline void resetNumbInfo() {
		for (auto i = m_cNumbInfo.GetSize(); i--;)
			m_cNumbInfo[i]->resetNumbInfo();
	}
	CC inline auto nonEmptyInfo() const	{
		return m_cNumbInfo.GetSize() > 1 || numMatrOfType(t_canonical);
	}
	CK inline ulonglong numMatrOfType(t_design_type t)	const	{ return m_cNumbInfo[0]->numMatrOfType(t); }
	CC inline void addMatrixTrans(ulonglong n, ulonglong nS)	{ m_cNumbInfo[0]->addMatrixTrans(n, nS); }
#if CANON_ON_GPU
	// Andrei (10/25/2021) Following function was NOT tested after COrderInfo refactoring
	inline void addMatrix(const COrderInfo *pOrderInfo)			{
		m_cNumbInfo[0]->addMatrices(pOrderInfo->getNumInfoPtr());
	}

	inline CNumbInfo *getNumInfoPtr() const						{ return m_cNumbInfo[0]; }
#endif
	COrderNumb *GetByKey(size_t extraOrder);
	inline size_t numOrderNumbers() const						{ return m_cNumbInfo.GetSize(); }
	inline COrderNumb *getOrderNumbers(size_t idx=0) const		{ return m_cNumbInfo[idx]; }
	inline CNumbInfo *getCombinedNumbInfo() const				{ return m_cNumbInfo[0]; }
private:
	size_t m_groupOrder;
	COrderNumbArray m_cNumbInfo;
};

#define GROUP_INFO_TYPE	  			COrderInfo *
#define GROUP_INFO_ACCESS_TYPE		GROUP_INFO_TYPE
typedef CArray<GROUP_INFO_TYPE, GROUP_INFO_ACCESS_TYPE> CArrayGroupInfo;

class CGroupsInfo : public CArrayGroupInfo
{
public:
	CC CGroupsInfo() {
		// Since most of the matrices of the big sets will have trivial 
		// automorphism group, it makes sence to deal with them separately.
		// Add the group of order 1 with 0 counter
		Add(new COrderInfo(1, 1, 0));
	}

	CC ~CGroupsInfo() { 
		for (auto i = GetSize(); i--;)
			delete GetAt(i);
	}

	CC CNumbInfo *addGroupOrder(size_t groupOrder, size_t extraGroupOrder = 1, ulonglong numb = 1, ulonglong numSimple = 0) {
		if (groupOrder == 1)
			return GetAt(0)->addMatrix(extraGroupOrder, numb, numSimple);

		size_t left = 1;   // Starting at 1 because groups of order 1 are processed separately
		size_t right = GetSize() - 1;
		while (left <= right) {
			const auto i = (right + left) >> 1;
			const auto grOrder = GetAt(i)->groupOrder();
			if (grOrder == groupOrder)
				return GetAt(i)->addMatrix(extraGroupOrder, numb, numSimple);

			if (grOrder < groupOrder)
				left = i + 1;
			else
				right = i - 1;
		}

		auto *pOrderInfo = new COrderInfo(groupOrder, extraGroupOrder, numb, numSimple);
		InsertAt(left, pOrderInfo);
		return pOrderInfo->getOrderNumbers(0);
	}
	CC void resetGroupsInfo() {
		for (auto i = GetSize(); i--;)  // Loop over all group numbers
			GetAt(i)->resetNumbInfo();
	}

	void printGroupInfo(FILE *file) const;
	void calcCountersTotal(COrderInfo *pTotal);
	CK void updateGroupInfo(const CGroupsInfo *pGroupInfo);
	inline auto GetStartIdx() const				{ return GetAt(0)->nonEmptyInfo()? 0 : 1; }
protected:
	CK void updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem);
private:
	void addGroupOrders(const COrderInfo* pOrderInfo);
};

