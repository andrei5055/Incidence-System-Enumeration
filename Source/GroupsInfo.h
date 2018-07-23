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

class COrderInfo : public CNumbInfo
{
public:
	CC COrderInfo()												{}
	CC COrderInfo(size_t order, ulonglong num = 1, ulonglong numSimple = 0) : m_groupOrder(order)
																{ addMatrix(num, numSimple); }
	CC inline size_t groupOrder() const							{ return m_groupOrder; }
private:
	size_t m_groupOrder;
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
		Add(new COrderInfo(1, 0));
	}

	CC ~CGroupsInfo() { 
		for (size_t i = GetSize(); i--;)
			delete GetAt(i);
	}

	CC COrderInfo *addGroupOrder(size_t groupOrder, ulonglong numb = 1, ulonglong numSimple = 0) {
		if (groupOrder == 1) {
			GetAt(0)->addMatrix(numb, numSimple);
			return GetAt(0);
		}

		size_t left = 1;
		size_t right = GetSize() - left;
		while (left <= right) {
			const size_t i = (right + left) >> 1;

			const size_t grOrder = GetAt(i)->groupOrder();
			if (grOrder == groupOrder) {
				GetAt(i)->addMatrix(numb, numSimple);
				return GetAt(i);
			}

			if (grOrder < groupOrder)
				left = i + 1;
			else
				right = i - 1;
		}

		COrderInfo *pOrderInfo = new COrderInfo(groupOrder, numb, numSimple);
		InsertAt(left, pOrderInfo);
		return pOrderInfo;
	}
	CC void resetGroupsInfo() {
		for (size_t i = GetSize(); i--;)
			GetAt(i)->resetNumbInfo();
	}

	void printGroupInfo(FILE *file) const;
	void calcCountersTotal(COrderInfo *pTotal);
	CK void updateGroupInfo(const CGroupsInfo *pGroupInfo);
protected:
	CK void updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem);
};

