#include "GroupsInfo.h"

void CNumbInfo::addMatrices(const CNumbInfo *pNumbInfo)
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
	const auto nElem = pGroupInfo->GetSize();
	for (auto i = pGroupInfo->GetStartIdx(); i < nElem; i++) {
		const auto *pOrderInfo = pGroupInfo->GetAt(i);
		const auto numCanons = pOrderInfo->numMatrOfType(t_canonical);
		const auto numSimple = pOrderInfo->numMatrOfType(t_simple);
		auto *pInfo = addGroupOrder(pOrderInfo->groupOrder(), 1, numCanons, numSimple);
		pInfo->addMatrixTrans(pOrderInfo->numMatrOfType(t_transitive), pOrderInfo->numMatrOfType(t_simpleTrans));
	}
}

void CGroupsInfo::updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem)
{
	size_t i = pOrderInfoBase->numMatrOfType(t_canonical) ? 0 : 1;
	for (; i < nElem; i++) {
		const COrderInfo *pOrderInfo = pOrderInfoBase + i;
		const auto numCanons = pOrderInfo->numMatrOfType(t_canonical);
		const auto numSimple = pOrderInfo->numMatrOfType(t_simple);
		auto *pInfo = addGroupOrder(pOrderInfo->groupOrder(), 1, numCanons, numSimple);
		pInfo->addMatrixTrans(pOrderInfo->numMatrOfType(t_transitive), pOrderInfo->numMatrOfType(t_simpleTrans));
	}
}

void CGroupsInfo::printGroupInfo(FILE *file) const
{
	auto i = GetStartIdx();
	const auto iMax = GetSize();
	if (i == iMax)
		return;			// Nothing was constructed

	char buffer[256], line[256];
	// Loop to find biggest additional space, used by the group on parts
	size_t maxLen = 0;
	for (; i < iMax; i++) {  // Loop over the orders of the groups
		const auto* pInfo = GetAt(i);
		for (size_t j = 0; j < pInfo->numOrderNumbers(); j++) {
			const auto* p = pInfo->getOrderNumbers(j);
			if (p->groupOrder() <= 1)
				continue;

			const size_t len = SPRINTF(buffer, "*%zd", p->groupOrder());
			if (maxLen < len)
				maxLen = len;
		}
	}

#define SHIFT "    "
	size_t len = SPRINTF(buffer, "\n" SHIFT "    |Aut(D)|");
	if (maxLen) {
		// At least one group order correspond to the design with nontrivial group on parts
		// Adjusting the title of the table
		memset(buffer + len, ' ', maxLen);
		len += maxLen;
	}

	len += SNPRINTF(buffer+len, countof(line)-len, "          Nd:             Ns:            Ndt:            Nst:\n");
	outString(buffer, file);

	strcpy_s(line, countof(line), SHIFT);
	const auto l_Shift = strlen(SHIFT);
	memset(line + l_Shift, '_', len);
	len += l_Shift;
	strcpy_s(line + len, countof(line) - len, "\n");
	outString(line, file);

	COrderInfo total(0, 1, 0);
	auto *pCombinedNumbInfo = total.getCombinedNumbInfo();
	for (i = GetStartIdx(); i < iMax; i++) {
		const auto *pInfo = GetAt(i);
		for (size_t j = 0; j < pInfo->numOrderNumbers(); j++) {
			const auto *p = pInfo->getOrderNumbers(j);
			pCombinedNumbInfo->addMatrices(p);
			len = SPRINTF(buffer, SHIFT"%10zd", pInfo->groupOrder());
			if (maxLen) {
				if (p->groupOrder() <= 1) {
					memset(buffer + len, ' ', maxLen);
					len += maxLen;
				} else
					len += SNPRINTF(buffer + len, countof(buffer) - len, "*%zd", p->groupOrder());
			}

			pInfo->outNumbInfo(buffer, countof(buffer) - len, len);
		}

		outString(buffer, file);
	}

	outString(line, file);
	len = SPRINTF(buffer, "        Total:");
	if (maxLen) {
		memset(buffer + len, ' ', maxLen);
		len += maxLen;
	}
	total.outNumbInfo(buffer, countof(buffer) - len, len);
	outString(buffer, file);
}

#if CANON_ON_GPU
void CGroupsInfo::calcCountersTotal(COrderInfo *pTotal)
{
	const auto iMax = GetSize();
	for (auto i = GetStartIdx(); i < iMax; i++)
		pTotal->addMatrix(GetAt(i));
}
#endif
