#include "GroupInfo.h"

CC int CGroupInfo::updateGroupOrder(ctchar* tr) {
	// search for tr 	
	int itr;
	int low = 0;
	const auto grOrder = groupOrder();
	auto high = itr = grOrder - 1;
	int cmp = -1;
	while (low <= high) {
		itr = low + ((high - low) >> 1);
		cmp = MEMCMP(getObjPntr(itr), tr, m_lenObj);
		if (!cmp)
			return itr;

		if (cmp < 0)
			low = itr + 1;  // ignore left half
		else
			high = itr - 1; // ignore right half
	}

	if (cmp < 0)
		itr++;

	auto* cmpTr = getObjAddr(grOrder);

	if (itr < grOrder)
		insert(itr, grOrder);
	else
		push_back(grOrder);

	memcpy(cmpTr, tr, m_lenObj);
	return itr;
}
