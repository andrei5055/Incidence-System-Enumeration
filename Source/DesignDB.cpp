#include <cstring>
#include "DesignDB.h"


bool CDesignDB::AddRecord(const unsigned char* pRecord) {
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

size_t CDesignDB::FindRecord(const unsigned char* pRecord, int *pResCmp) {
	*pResCmp = 1;
	size_t mid(0), low(0), high(m_nRecNumb);
	// Repeat until the pointers low and high meet each other
	while (low < high) {
		mid = low + (high - low) / 2;

		const unsigned char* pRecordMid = m_pRecStorage + m_pRecPermutation[mid] * m_nRecLen + LEN_HEADER;
		if (!(*pResCmp = memcmp(pRecordMid, pRecord, m_nRecLen - LEN_HEADER)))
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

	auto* pNewRecStorage = new unsigned char[newRecNumber * m_nRecLen];
	memcpy(pNewRecStorage, m_pRecStorage, m_nRecNumbMax * m_nRecLen * sizeof(pNewRecStorage[0]));
	delete[] m_pRecStorage;
	m_pRecStorage = pNewRecStorage;
	m_nRecNumbMax = newRecNumber;
	return true;
}