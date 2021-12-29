#include <cstring>
#include "DesignDB.h"


bool CDesignDB::AddRecord(recPtr pRecord) {
	int cmpRes;
	size_t iMin = FindRecord(pRecord, &cmpRes);
	if (!cmpRes) {
		// Increase counter for the record just found
		++*(unsigned int *)(m_pRecStorage + iMin * m_nRecLen);
		return false;
	}

	if (m_nRecNumb == m_nRecNumbMax) {
		// All previously allocated memory were used - need to reallocate
		reallocateMemory();
	}

	auto* pntr = m_pRecStorage + m_nRecNumb * m_nRecLen;
	memcpy(pntr + LEN_HEADER, pRecord, m_nRecLen - LEN_HEADER);
	*(unsigned int *)(pntr) = 1;

	size_t i = m_nRecNumb;
	for (; i > iMin; i--)
		m_pRecPermutation[i] = m_pRecPermutation[i - 1];

	m_pRecPermutation[i] = m_nRecNumb++;
	return true;
}

size_t CDesignDB::FindRecord(recPtr pRecord, int *pResCmp) {
	*pResCmp = 1;
	size_t mid(0), low(0), high(m_nRecNumb);
	// Repeat until the pointers low and high meet each other
	while (low < high) {
		mid = low + (high - low) / 2;

		recPtr pRecordMid = m_pRecStorage + m_pRecPermutation[mid] * recordLength() + LEN_HEADER;
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
	memcpy(pNewRecStorage, m_pRecStorage, m_nRecNumbMax * recordLength() * sizeof(pNewRecStorage[0]));
	delete[] m_pRecStorage;
	m_pRecStorage = pNewRecStorage;
	m_nRecNumbMax = newRecNumber;
	return true;
}

void CDesignDB::SortRecods(FILE *file) {
	if (recNumb() <= 1)
		return;   // nothing to sort;

	delete[] m_pSortedRecords;
	m_pSortedRecords = new size_t[recNumb()];
	for (auto i = recNumb(); i--;)
		m_pSortedRecords[i] = i;

	quickSort(m_pSortedRecords, 0, recNumb()-1);
	if (!file)
		return;

	fprintf(file, "\nMaster Design:   Numb of Decompositions:\n");
	for (size_t i = 0; i < recNumb(); i++) {
		auto* pRec = firstRecord() + m_pSortedRecords[i] * recordLength();
		fprintf(file, "   %6zd:                 %5d\n", i+1, *(int*)pRec);
	}
}

int CDesignDB::compareRecords(const size_t idx, recPtr pSecnd) const {
	static int cntr = 0; cntr++;
	auto *pFirst = firstRecord() + idx * recordLength();
	const auto decompNumber_1 = *(unsigned int *)pFirst;
	const auto decompNumber_2 = *(unsigned int *)pSecnd;
	if (decompNumber_1 == decompNumber_2)
		return 0;

	return decompNumber_1 > decompNumber_2 ? 1 : -1;
}

void CDesignDB::quickSort(size_t *arr, size_t left, size_t right) const {
	size_t i = left, j = right;
	const auto pivotIdx = (left + right) >> 1;
	auto pivot = firstRecord() + arr[pivotIdx] * recordLength();

	/* partition */
	while (i <= j) {
		while (i != pivotIdx && compareRecords(arr[i], pivot) == -1)
			i++;

		while (j != pivotIdx && compareRecords(arr[j], pivot) == 1)
			j--;

		if (i < j) {
			const auto tmp = arr[i];
			arr[i++] = arr[j];
			arr[j--] = tmp;
		} else {
			i++;
			if (j) j--;
			else break;
		}
	}

	/* recursion */
	if (left < j)
		quickSort(arr, left, j);

	if (i < right)
		quickSort(arr, i, right);
}
