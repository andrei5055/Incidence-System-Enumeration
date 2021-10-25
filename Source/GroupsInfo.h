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
#define numMatrices()		numMatrOfType(t_canonical)
#define numSimpleMatrices() numMatrOfType(t_simple)
	inline void setNumMatrOfType(ulonglong val, t_design_type t){ m_nCounter[t] = val; }
	CC inline void addMatrix(ulonglong numb, ulonglong nSimple)	{
		addMatrOfType(numb, t_canonical);
		addMatrOfType(nSimple, t_simple);
	}
	CC void addMatrixTrans(ulonglong numb, ulonglong nSimple)	{
		addMatrOfType(numb, t_transitive);
		addMatrOfType(nSimple, t_simpleTrans);
	}

	void addMatrix(const CNumbInfo *pNumbInfo);
	void outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const;
	CC inline void resetNumbInfo(int idx = 0)					{ memset(m_nCounter + idx, 0, sizeof(m_nCounter) - idx * sizeof(m_nCounter[0])); }
protected:
	CC inline void addMatrOfType(ulonglong v, t_design_type t)	{ m_nCounter[t] += v; }
private:
	ulonglong m_nCounter[t_design_type_total];
};

class COrderNumb : public CNumbInfo
{
public:
	COrderNumb(size_t groupOrder=1) : m_groupOrder(groupOrder) {}
	CC inline size_t groupOrder() const { return m_groupOrder; }
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
		m_cNumbInfo[0]->addMatrix(num, numSimple);
	}
	CC ~COrderInfo() {
		for (auto i = m_cNumbInfo.GetSize(); i--;)
			delete m_cNumbInfo[i];
	}
	CC inline size_t groupOrder() const							{ return m_groupOrder; }
	CC inline void addMatrix(size_t extraGroupOrder, ulonglong num, ulonglong nSimple) {
		m_cNumbInfo[0]->addMatrix(num, nSimple);
	}
	CC inline void resetNumbInfo(int idx = 0)					{ m_cNumbInfo[0]->resetNumbInfo(idx); }
	CK inline ulonglong numMatrOfType(t_design_type t)	const	{ return m_cNumbInfo[0]->numMatrOfType(t); }
	CC inline void addMatrixTrans(ulonglong n, ulonglong nS)	{ m_cNumbInfo[0]->addMatrixTrans(n, nS); }
	inline void addMatrix(const COrderInfo *pOrderInfo)			{ m_cNumbInfo[0]->addMatrix(pOrderInfo->getNumInfoPtr()); }
	inline CNumbInfo *getNumInfoPtr() const						{ return m_cNumbInfo[0]; }
	inline void outNumbInfo(char* buffer, const size_t lenBuf, size_t poz) const {
		m_cNumbInfo[0]->outNumbInfo(buffer, lenBuf, poz);
	}
	COrderNumb* GetByKey(size_t extraOrder) const				{ return m_cNumbInfo[0]; }
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

	CC COrderInfo/*COrderNumb*/ * addGroupOrder(size_t groupOrder, size_t extraGroupOrder = 1, ulonglong numb = 1, ulonglong numSimple = 0) {
		if (groupOrder == 1) {
			GetAt(0)->addMatrix(extraGroupOrder, numb, numSimple);
			return GetAt(0); // ->GetByKey(extraGroupOrder);
		}

		size_t left = 1;
		size_t right = GetSize() - left;
		while (left <= right) {
			const size_t i = (right + left) >> 1;

			const size_t grOrder = GetAt(i)->groupOrder();
			if (grOrder == groupOrder) {
				GetAt(i)->addMatrix(extraGroupOrder, numb, numSimple);
				return GetAt(i); // ->GetByKey(extraGroupOrder);
			}

			if (grOrder < groupOrder)
				left = i + 1;
			else
				right = i - 1;
		}

		COrderInfo *pOrderInfo = new COrderInfo(groupOrder, extraGroupOrder, numb, numSimple);
		InsertAt(left, pOrderInfo);
		return pOrderInfo; // pOrderInfo->GetByKey(extraGroupOrder);
	}
	CC void resetGroupsInfo() {
		for (size_t i = GetSize(); i--;)
			GetAt(i)->resetNumbInfo();
	}

	void printGroupInfo(FILE *file) const;
	void calcCountersTotal(COrderInfo *pTotal);
	CK void updateGroupInfo(const CGroupsInfo *pGroupInfo);
	auto GetStartIdx() const				{ return GetAt(0)->numMatrices() ? 0 : 1; }
protected:
	CK void updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem);
};

