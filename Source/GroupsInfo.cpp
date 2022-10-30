#include "GroupsInfo.h"

void CNumbInfo::addMatrices(const CNumbInfo *pNumbInfo)
{
	for (auto i = t_design_type::t_totalConstr; i < t_design_type::t_design_type_total; ++i)
		addMatrOfType(pNumbInfo->numMatrOfType(i), i);
}

void CNumbInfo::outNumbInfo(char *buffer, const size_t lenBuf, size_t poz) const
{
	const char *formats[] = { "", "   %12llu", "   %10llu", "   %4llu", "   %4llu" };
	for (auto j = t_design_type::t_canonical; j < t_design_type::t_design_type_total; ++j)
		poz += SNPRINTF(buffer + poz, lenBuf - poz, formats[+j], numMatrOfType(j));

	SNPRINTF(buffer + poz, lenBuf - poz, "\n");
}

COrderNumb *COrderInfo::GetByKey(size_t groupOrder) {
	size_t left = 0;
	if (groupOrder == 1) {
		if (m_cNumbInfo.GetAt(0)->groupOrder() == 1)
			return m_cNumbInfo.GetAt(0);
	}
	else {
		size_t right = m_cNumbInfo.GetSize() - 1;
		while (left <= right) {
			const auto i = (right + left) >> 1;
			const auto grOrder = m_cNumbInfo.GetAt(i)->groupOrder();
			if (grOrder == groupOrder) {
				return m_cNumbInfo.GetAt(i);
			}

			if (grOrder < groupOrder)
				left = i + 1;
			else {
				if (i)
					right = i - 1;
				else
					break;  // index cannot be negative
			}
		}
	}

	auto *pOrderInfo = new COrderNumb(groupOrder);
	m_cNumbInfo.InsertAt(left, pOrderInfo);
	return pOrderInfo;
}

void CGroupsInfo::addGroupOrders(const COrderInfo *pOrderInfo) {
	const auto groupOrder = pOrderInfo->groupOrder();
	for (size_t j = 0; j < pOrderInfo->numOrderNumbers(); j++) {
		const auto* p = pOrderInfo->getOrderNumbers(j);
		const auto numCanons = p->numMatrOfType(t_design_type::t_canonical);
		const auto numSimple = p->numMatrOfType(t_design_type::t_simple);
		auto* pInfo = addGroupOrder(groupOrder, p->groupOrder(), numCanons, numSimple);
		pInfo->addMatrixTrans(p->numMatrOfType(t_design_type::t_transitive), p->numMatrOfType(t_design_type::t_simpleTrans));
	}
}

void CGroupsInfo::updateGroupInfo(const CGroupsInfo *pGroupInfo) {
	const auto nElem = pGroupInfo->GetSize();
	for (auto i = pGroupInfo->GetStartIdx(); i < nElem; i++)
		addGroupOrders(pGroupInfo->GetAt(i));
}

void CGroupsInfo::updateGroupInfo(const COrderInfo *pOrderInfoBase, size_t nElem) {
	auto i = pOrderInfoBase->numMatrOfType(t_design_type::t_canonical) ? 0 : 1;
	for (; i < nElem; i++)
		addGroupOrders(pOrderInfoBase + i);
}

void CGroupsInfo::printGroupInfo(FILE *file, COrderInfo& total) const
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
		const auto a = pInfo->groupOrder();
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

	len += SNPRINTF(buffer+len, countof(line)-len, "          Nd:           Ns:   Ndt:   Nst:\n");
	outString(buffer, file);

	strcpy_s(line, countof(line), SHIFT);
	const auto l_Shift = strlen(SHIFT);
	memset(line + l_Shift, '_', len);
	len += l_Shift;
	strcpy_s(line + len, countof(line) - len, "\n");
	outString(line, file);

	const auto lenToCount = maxLen + SPRINTF(buffer, SHIFT"%10zd", (size_t)1);
	auto *pCombinedNumbInfo = total.getCombinedNumbInfo();
	for (i = GetStartIdx(); i < iMax; i++) {
		const auto *pInfo = GetAt(i);
		// For i == 0, first element, which corresponds to the group of order 1
		// of the COrderNumbArray  could be empty
		size_t j = i || pInfo->numMatrOfType(t_design_type::t_canonical) > 0? 0 : 1;
		for (; j < pInfo->numOrderNumbers(); j++) {
			const auto *p = pInfo->getOrderNumbers(j);
			pCombinedNumbInfo->addMatrices(p);
			len = SPRINTF(buffer, SHIFT"%10zd", pInfo->groupOrder());
			if (maxLen) {
				if (p->groupOrder() > 1)
					len += SNPRINTF(buffer + len, countof(buffer) - len, "*%zd", p->groupOrder());

				if (len < lenToCount) {
					memset(buffer + len, ' ', lenToCount - len);
					len = lenToCount;
				}
			}

			p->outNumbInfo(buffer, countof(buffer) - len, len);
			outString(buffer, file);
		}
	}

	outString(line, file);
	len = SPRINTF(buffer, "        Total:");
	if (maxLen) {
		memset(buffer + len, ' ', maxLen);
		len += maxLen;
	}
	total.getOrderNumbers()->outNumbInfo(buffer, countof(buffer) - len, len);
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
