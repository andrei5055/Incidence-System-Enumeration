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
	CNumbInfo()													{ memset(m_nCounter, 0, sizeof(m_nCounter)); }
	inline ulonglong numMatrOfType(t_design_type t)	const		{ return m_nCounter[t]; }
#define numMatrices()		numMatrOfType(t_canonical)
#define numSimpleMatrices() numMatrOfType(t_simple)
	inline void setNumMatrOfType(ulonglong val, t_design_type t){ m_nCounter[t] = val; }
	void addMatrix(ulonglong numb, ulonglong numSimple);
	void addMatrixTrans(ulonglong numb, ulonglong numSimple);
	void addMatrix(const CNumbInfo *pNumbInfo);
	void outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const;
protected:
	inline void addMatrOfType(ulonglong val, t_design_type t)	{ m_nCounter[t] += val; }
private:
	ulonglong m_nCounter[t_design_type_total];
};

class COrderInfo : public CNumbInfo
{
public:
	COrderInfo(int order, ulonglong num = 1, ulonglong numSimple = 0) : m_groupOrder(order)
																{ addMatrix(num, numSimple); }
	inline uint groupOrder() const								{ return m_groupOrder; }
private:
	uint const m_groupOrder;
};

#define GROUP_INFO_TYPE	  			COrderInfo *
#define GROUP_INFO_ACCESS_TYPE		GROUP_INFO_TYPE
typedef CArray<GROUP_INFO_TYPE, GROUP_INFO_ACCESS_TYPE> CArrayGroupInfo;

class CGroupsInfo : public CArrayGroupInfo
{
public:
	CGroupsInfo();
	~CGroupsInfo();
	COrderInfo *addGroupOrder(uint groupOrder, ulonglong numb = 1, ulonglong numSimple = 0);
	void printGroupInfo(FILE *file) const;
	void updateGroupInfo(const CGroupsInfo *pGroupInfo);
};

