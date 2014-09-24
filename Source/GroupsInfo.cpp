#include "GroupsInfo.h"

void CNumbInfo::addMatrix(ulonglong numb, ulonglong numSimple)
{
	addMatrOfType(numb, t_canonical);
	addMatrOfType(numSimple, t_simple);
}

void CNumbInfo::addMatrix(const CNumbInfo *pNumbInfo)
{
	for (auto i = t_totalConstr; i < t_design_type_total; i = (t_design_type)(i + 1))
		addMatrOfType(pNumbInfo->numMatrOfType(i), i);
}

void CNumbInfo::addMatrixTrans(ulonglong numb, ulonglong numSimple)
{
	addMatrOfType(numb, t_transitive);
	addMatrOfType(numSimple, t_simpleTrans);
}

void CNumbInfo::outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const
{
	for (auto j = t_canonical; j < t_design_type_total; j = (t_design_type)(j + 1))
		poz += sprintf_s(buffer + poz, lenBuf - poz, "      %10llu", numMatrOfType(j));

	sprintf_s(buffer + poz, lenBuf - poz, "\n");
}

CGroupsInfo::CGroupsInfo()
{
	// Since most of the matrices of the big sets will have trivial 
	// automorphism group, it makes sence to deal with hem separately.
	// Add the group of order 1 with 0 counter
	Add(new COrderInfo(1, 0));
}


CGroupsInfo::~CGroupsInfo()
{
	for (size_t i = GetSize(); i--;)
		delete GetAt(i);
}

COrderInfo *CGroupsInfo::addGroupOrder(uint groupOrder, ulonglong numb, ulonglong numSimple)
{
	if (groupOrder == 1) {
		GetAt(0)->addMatrix(numb, numSimple);
		return GetAt(0);
	}

 	size_t left = 1;
	size_t right = GetSize() - left;
	int diffOrder = 0;
	while (left <= right) {
		const size_t i = (right + left) >> 1;

		diffOrder = GetAt(i)->groupOrder() - groupOrder;
		if (!diffOrder) {
			GetAt(i)->addMatrix(numb, numSimple);
			return GetAt(i);
		}

		if (diffOrder < 0)
			left = i + 1;
		else
			right = i - 1;
	}

	COrderInfo *pOrderInfo = new COrderInfo(groupOrder, numb, numSimple);
	InsertAt(left, pOrderInfo);
	return pOrderInfo;
}

void CGroupsInfo::updateGroupInfo(const CGroupsInfo *pGroupInfo)
{
	size_t i = pGroupInfo->GetAt(0)->numMatrices() ? 0 : 1;
	const size_t iMax = pGroupInfo->GetSize();
	for (; i < iMax; i++) {
		const COrderInfo *pOrderInfo = pGroupInfo->GetAt(i);
		COrderInfo *pInfo = addGroupOrder(pOrderInfo->groupOrder(), pOrderInfo->numMatrices(), pOrderInfo->numSimpleMatrices());
		pInfo->addMatrixTrans(pOrderInfo->numMatrOfType(t_transitive), pOrderInfo->numMatrOfType(t_simpleTrans));
	}
}

void CGroupsInfo::printGroupInfo(FILE *file) const
{
#define SHIFT "    "
	char buffer[256], line[256];
	size_t len = SPRINTF(buffer, "\n" SHIFT "    |Aut(D)|          Nd:             Ns:            Ndt:            Nst:\n");
	outString(buffer, file);

	strcpy_s(line, countof(line), SHIFT);
	const size_t l_Shift = strlen(SHIFT);
	memset(line + l_Shift, '_', len);
	len += l_Shift;
	strcpy_s(line + len, countof(line) - len, "\n");
	outString(line, file);

	COrderInfo total(0, 0);
	size_t i = GetAt(0)->numMatrices()? 0 : 1;
	const size_t iMax = GetSize();
	for (; i < iMax; i++) {
		const COrderInfo *pInfo = GetAt(i);
		total.addMatrix(pInfo);
		len = SPRINTF(buffer, SHIFT"%10d", pInfo->groupOrder());
		pInfo->outNumbInfo(buffer, countof(buffer) - len, len);
		outString(buffer, file);
	}

	outString(line, file);
	len = SPRINTF(buffer, "        Total:");
	total.outNumbInfo(buffer, countof(buffer) - len, len);
	outString(buffer, file);
}