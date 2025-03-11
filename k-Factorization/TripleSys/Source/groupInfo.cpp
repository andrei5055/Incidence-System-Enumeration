#include "GroupInfo.h"

CC int CRepository::getElementIndex(ctchar* tr) const {
	// search for element 	
	int itr;
	int low = 0;
	const auto nElem = numObjects();
	auto high = itr = nElem - 1;
	int cmp = -1;
	while (low <= high) {
		itr = low + ((high - low) >> 1);
		cmp = MEMCMP(getObjPntr(itr), tr, m_lenObj);
		if (!cmp)
			return -itr - 1;

		if (cmp < 0)
			low = itr + 1;  // ignore left half
		else
			high = itr - 1; // ignore right half
	}

	if (cmp < 0)
		itr++;

	return itr;
}

CC int CRepository::updateRepo(ctchar* tr) {
	// search for element 	
	const auto itr = getElementIndex(tr);
	if (itr < 0)
		return itr;

	const auto nElem = numObjects();
	auto* cmpTr = getObjAddr(nElem);

	if (itr < nElem)
		insert(itr, nElem);
	else
		push_back(nElem);

	memcpy(cmpTr, tr, m_lenObj);
	return itr;
}
