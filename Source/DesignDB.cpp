#include <cstring>
#include "DesignDB.h"
#include <assert.h>


void CDesignDB::AddRecord(recPtr pRecord, size_t groupOrder) {
	int cmpRes;
	const auto iMin = FindRecord(pRecord, &cmpRes);
	if (!cmpRes) {
		// Increase counter for the record just found
		auto* pMasterInfo = (masterInfo *)getRecord(m_pRecPermutation[iMin]);
		pMasterInfo->numbDecomp++;
		assert(pMasterInfo->groupOrder == groupOrder);
		return;
	}

	if (m_nRecNumb == m_nRecNumbMax) {
		// All previously allocated memory were used - need to reallocate
		reallocateMemory();
	}

	auto* pntr = (unsigned char *)getRecord(m_nRecNumb);
	memcpy(pntr + LEN_HEADER, pRecord, recordLength() - LEN_HEADER);
	auto* pMasterInfo = (masterInfo*)(pntr);
	pMasterInfo->numbDecomp = 1;
	pMasterInfo->groupOrder = groupOrder;

	size_t i = m_nRecNumb;
	for (; i > iMin; i--)
		m_pRecPermutation[i] = m_pRecPermutation[i - 1];

	m_pRecPermutation[i] = m_nRecNumb++;
}

size_t CDesignDB::FindRecord(recPtr pRecord, int *pResCmp) {
	*pResCmp = 1;
	size_t mid(0), low(0), high(m_nRecNumb);
	// Repeat until the pointers low and high meet each other
	while (low < high) {
		mid = low + (high - low) / 2;

		recPtr pRecordMid = getRecord(m_pRecPermutation[mid]) + LEN_HEADER;
		if (!(*pResCmp = memcmp(pRecordMid, pRecord, recordLength() - LEN_HEADER)))
			return mid;

		if (*pResCmp < 0)
			low = mid + 1;
		else
			high = mid;
	}

	return *pResCmp > 0? mid : mid + 1;
}

bool CDesignDB::reallocateMemory() {
	const auto newRecNumber = 2 * m_nRecNumbMax;
	auto *pRecPerm = new size_t[newRecNumber];
	memcpy(pRecPerm, m_pRecPermutation, m_nRecNumbMax * sizeof(pRecPerm[0]));
	delete[] m_pRecPermutation;
	m_pRecPermutation = pRecPerm;

	auto* pNewRecStorage = new unsigned char[newRecNumber * recordLength()];
	memcpy(pNewRecStorage, firstRecord(), m_nRecNumbMax * recordLength() * sizeof(pNewRecStorage[0]));
	setRecordStorage(pNewRecStorage);
	m_nRecNumbMax = newRecNumber;
	return true;
}

int compareRecordsA(recPtr pRec1, recPtr pRec2) {
	const auto* pMaster_1 = (masterInfo*)pRec1;
	const auto* pMaster_2 = (masterInfo*)pRec2;
	const auto decompNumber_1 = pMaster_1->numbDecomp;
	const auto decompNumber_2 = pMaster_2->numbDecomp;
	if (decompNumber_1 == decompNumber_2) {
		return (pMaster_1->groupOrder == pMaster_2->groupOrder) ? 0 :
			pMaster_1->groupOrder > pMaster_2->groupOrder ? 1 : -1;
	}

	return decompNumber_1 > decompNumber_2 ? 1 : -1;
}

int compareRecords(recPtr pRec1, recPtr pRec2) {
	const auto* pMaster_1 = (masterInfo*)pRec1;
	const auto* pMaster_2 = (masterInfo*)pRec2;
	if (pMaster_1->groupOrder == pMaster_2->groupOrder) {
		return (pMaster_1->numbDecomp == pMaster_2->numbDecomp) ? 0 :
			pMaster_1->numbDecomp > pMaster_2->numbDecomp ? 1 : -1;
	}

	return pMaster_1->groupOrder > pMaster_2->groupOrder ? 1 : -1;
}

void CDesignDB::SortRecods(FILE *file) {
	if (!recNumb())
		return;   // nothing to sort;

	setCompareFunc(compareRecords);
	const auto *pSortedRecords = Sort(recNumb());
	if (!file)
		return;

	masterInfo* pRec;
#if 0
	fprintf(file, "\nMaster #:    Number of Decomp:   |Aut(M)|:\n");
	for (size_t i = 0; i < recNumb(); i++) {
		pRec = (masterInfo*)getRecord(pSortedRecords[i]);
		fprintf(file, "%4zd:          %6zd           %5zd\n", i+1, pRec->numbDecomp, pRec->groupOrder);
	}
#else
#define SHIFT "    "
	char buff[256];
	const auto lenStr = sprintf_s(buff, SHIFT " |Aut(M)|:    Decomp:    Masters: ");
	fprintf(file, "\n%s\n", buff);
	memset(buff, '_', lenStr);
	fprintf(file, SHIFT "%s\n", buff);

	const auto last = recNumb() - 1;
	size_t i = 0;
	pRec = (masterInfo*)getRecord(pSortedRecords[0]);
	unsigned long long total = 0, totalDIBD = 0;
	while (i != recNumb()) {
		bool newGroupOrder = true;
		auto groupOrder = pRec->groupOrder;
		auto numbDecomp = pRec->numbDecomp;
		size_t numbMasters = 1;
		bool sameGroupOrder = true;
		bool sameDecomp = true;
		while (sameGroupOrder && (!sameDecomp || ++i < recNumb())) {
			pRec = (masterInfo*)getRecord(pSortedRecords[i]);
			sameGroupOrder = groupOrder == pRec->groupOrder;
			sameDecomp = numbDecomp == pRec->numbDecomp;
			if (sameGroupOrder && sameDecomp && i != last) {
				numbMasters++;
				continue;
			}

			fprintf(file, SHIFT "%7zd   %7zd    %7zd\n", groupOrder, numbDecomp, numbMasters);
			totalDIBD += numbMasters;
			total += numbMasters * numbDecomp;
			numbDecomp = pRec->numbDecomp;
			numbMasters = 1;
		}
	}
	fprintf(file, SHIFT "%s\n", buff);
	fprintf(file, SHIFT "  Total:  %7zd      %5zd\n", total, totalDIBD);
#endif
}


