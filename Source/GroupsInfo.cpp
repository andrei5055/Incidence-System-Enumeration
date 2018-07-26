#include "GroupsInfo.h"

void CNumbInfo::addMatrix(const CNumbInfo *pNumbInfo)
{
	for (auto i = t_totalConstr; i < t_design_type_total; i = (t_design_type)(i + 1))
		addMatrOfType(pNumbInfo->numMatrOfType(i), i);
}

void CNumbInfo::outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const
{
	for (auto j = t_canonical; j < t_design_type_total; j = (t_design_type)(j + 1))
		poz += SNPRINTF(buffer + poz, lenBuf - poz, "      %10llu", numMatrOfType(j));

	SNPRINTF(buffer + poz, lenBuf - poz, "\n");
}

void CGroupsInfo::updateGroupInfo(const CGroupsInfo *pGroupInfo)
{
	size_t i = pGroupInfo->GetAt(0)->numMatrices() ? 0 : 1;
	const size_t nElem = pGroupInfo->GetSize();
	for (; i < nElem; i++) {
		const COrderInfo *pOrderInfo = pGroupInfo->GetAt(i);
		COrderInfo *pInfo = addGroupOrder(pOrderInfo->groupOrder(), pOrderInfo->numMatrices(), pOrderInfo->numSimpleMatrices());
		pInfo->addMatrixTrans(pOrderInfo->numMatrOfType(t_transitive), pOrderInfo->numMatrOfType(t_simpleTrans));
	}
}

void CGroupsInfo::updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem)
{
	size_t i = pOrderInfoBase->numMatrices() ? 0 : 1;
	for (; i < nElem; i++) {
		const COrderInfo *pOrderInfo = pOrderInfoBase + i;
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
		len = SPRINTF(buffer, SHIFT"%10zd", pInfo->groupOrder());
		pInfo->outNumbInfo(buffer, countof(buffer) - len, len);
		outString(buffer, file);
	}

	outString(line, file);
	len = SPRINTF(buffer, "        Total:");
	total.outNumbInfo(buffer, countof(buffer) - len, len);
	outString(buffer, file);
}

void CGroupsInfo::calcCountersTotal(COrderInfo *pTotal)
{
	const size_t iMax = GetSize();
	for (size_t i = GetAt(0)->numMatrices()? 0 : 1; i < iMax; i++)
		pTotal->addMatrix(GetAt(i));
}