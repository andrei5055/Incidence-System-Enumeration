#include <cstring>
#include "DesignDB.h"
#include <assert.h>


size_t CDesignDB::AddRecord(recPtr pRecord, size_t groupOrder, size_t numbDecomp) {
	int cmpRes;
	const auto iMin = FindRecord(pRecord, &cmpRes);
	if (!cmpRes) {
		// Increase counter for the record just found
		const auto idx = m_pRecPermutation[iMin];
		auto* pMasterInfo = (masterInfo *)getRecord(idx);
		pMasterInfo->numbDecomp += numbDecomp;
		assert(pMasterInfo->groupOrder == groupOrder);
		return idx;
	}

	if (m_nRecNumb == m_nRecNumbMax) {
		// All previously allocated memory were used - need to reallocate
		reallocateMemory();
	}

	auto* pntr = (unsigned char *)getRecord(m_nRecNumb);
	memcpy(pntr + LEN_HEADER, pRecord, recordLength() - LEN_HEADER);
	auto* pMasterInfo = (masterInfo*)(pntr);
	pMasterInfo->numbDecomp = numbDecomp;
	pMasterInfo->groupOrder = groupOrder;

	size_t i = m_nRecNumb;
	for (; i > iMin; i--)
		m_pRecPermutation[i] = m_pRecPermutation[i - 1];

	return m_pRecPermutation[i] = m_nRecNumb++;
}

size_t CDesignDB::FindRecord(recPtr pRecord, int *pResCmp) const {
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

void CDesignDB::mergeDesignDB(const CDesignDB* pDB) {
	for (size_t i = 0; i < pDB->recNumb(); i++) {
		auto *pRec = (const masterInfo*)pDB->getRecord(i);
		AddRecord((unsigned char*)pRec + LEN_HEADER, pRec->groupOrder, pRec->numbDecomp);
	}
}

void CDesignDB::mergeDesignDBs(const CDesignDB* pDB_A, const CDesignDB* pDB_B) {
	const auto len = recordLength() - LEN_HEADER;
	size_t ind_A = 0, ind_B = 0;
	const auto* perm_A = pDB_A->getPermut();
	const auto* perm_B = pDB_B->getPermut();
	const unsigned char* pRec_A, *pRec_B = NULL;
	int state = 3;
	while (true) {
		if (state & 1) {
			if (ind_A >= pDB_A->recNumb()) {
				perm_A = perm_B;
				pDB_A = pDB_B;
				ind_A = ind_B;
				if (state & 2)
					pRec_A = ind_A < pDB_B->recNumb() ? (const unsigned char*)(pDB_A)->getRecord(perm_A[ind_A++]) : NULL;
				else
					pRec_A = pRec_B;
				break;
			}

			pRec_A = (const unsigned char*)pDB_A->getRecord(perm_A[ind_A++]);
		}

		if (state & 2) {
			if (ind_B >= pDB_B->recNumb())
				break;

			pRec_B = (const unsigned char*)pDB_B->getRecord(perm_B[ind_B++]);
		}

		if (m_nRecNumb == m_nRecNumbMax) {
			// All previously allocated memory were used - need to reallocate
			reallocateMemory();
		}
		auto* pntr = (unsigned char*)getRecord(m_nRecNumb);
		getPermut()[m_nRecNumb] = m_nRecNumb;
		m_nRecNumb++;
		const auto cmpResult = memcmp(pRec_A + LEN_HEADER, pRec_B + LEN_HEADER, len);
		if (cmpResult <= 0) {
			memcpy(pntr, pRec_A, recordLength());
			if (!cmpResult) {
				((masterInfo*)pntr)->numbDecomp += ((const masterInfo*)pRec_B)->numbDecomp;
				state = 3;
			}
			else
				state = 1;
		}
		else {
			memcpy(pntr, pRec_B, recordLength());
			state = 2;
		}
	}

	// Copying remaining records
	while (pRec_A) {
		if (m_nRecNumb == m_nRecNumbMax) {
			// All previously allocated memory were used - need to reallocate
			reallocateMemory();
		}

		auto* pntr = (unsigned char*)getRecord(m_nRecNumb);
		getPermut()[m_nRecNumb] = m_nRecNumb;
		m_nRecNumb++;
		memcpy(pntr, pRec_A, recordLength());
		pRec_A = ind_A < pDB_A->recNumb() ? (const unsigned char*)pDB_A->getRecord(perm_A[ind_A++]) : NULL;
	}
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

void CDesignDB::SortRecods(FILE* file, int formatID) {
	if (!recNumb())
		return;   // nothing to sort;

	setCompareFunc(compareRecords);
	const auto* pSortedRecords = Sort(recNumb());
	if (!file)
		return;

	switch (formatID) {
	case 0: return outWithFormat_0(pSortedRecords, file);
	case 1: return outWithFormat_1(pSortedRecords, file);
	case 2: return outWithFormat_2(pSortedRecords, file);
	}
}

#define SHIFT			""
#define LEN_DECOMP		70
#define SHIFT_TO_DECOMP	26

void CDesignDB::outWithFormat_0(const size_t * pSortedRecords, FILE * file) const {
	char buff[256], *pBuff = buff;
	const size_t len_buff = sizeof(buff);
	sprintf_s(buff, SHIFT " |Aut(M)|: Masters: Decomp.      Decomposition Look:    ");
	const auto lenStr = 91;
	fprintf(file, "\n%s\n", buff);
	memset(buff, '_', lenStr);
	buff[lenStr] = '\0';
	fprintf(file, SHIFT "%s\n", buff);

	const auto last = recNumb() - 1;
	size_t i = 0;
	auto* pRec = (const masterInfo*)getRecord(pSortedRecords[i]);
	unsigned long long totalCombined = 0, totalMasters = 0;

	auto groupOrder = pRec->groupOrder;
	auto numbDecomp = pRec->numbDecomp;
	size_t numbDecompMaxGlobal, numbDecompMax, jMax, numbMasters, numbMastersGroup, numDecompGroup;
	numDecompGroup = numbMastersGroup = 0;
	numbDecompMaxGlobal = numbDecompMax = numbMasters = jMax = 1;
	while (++i < recNumb()) {
		pRec = (const masterInfo*)getRecord(pSortedRecords[i]);
		if (groupOrder == pRec->groupOrder && numbDecomp == pRec->numbDecomp) {
			numbMasters++;
			if (i != last)
				continue;
		}
		else {
			if (i == last)
				jMax = 2;
		}

		for (auto j = jMax; j--;) {
			numDecompGroup += numbMasters * numbDecomp;
			numbMastersGroup += numbMasters;
			if (numbMasters > 1)
				pBuff += sprintf_s(pBuff, len_buff - (pBuff - buff), pBuff == buff? "%zd*" : " + %zd*", numbMasters);

			if (numbDecompMax < numbDecomp) {
				if (numbDecompMaxGlobal < (numbDecompMax = numbDecomp))
					numbDecompMaxGlobal = numbDecompMax;
			}

			pBuff += sprintf_s(pBuff, len_buff - (pBuff - buff), numbMasters > 1 || pBuff == buff? "%zd" : " + %zd", numbDecomp);

			if (groupOrder != pRec->groupOrder || !j && i == last) {
				char saved;
				int j = LEN_DECOMP;
				int k = 0;
				while (true) {
					const bool flg = strlen(buff + k) > j;
					if (flg) {
						// Need to split decomposition information for given groupOrder into two strings
						while (buff[k+j] != '+') j--;
						saved = buff[k + ++j];
						buff[k + j] = '\0';
					}

					if (!k)
						fprintf(file, SHIFT "%7zd  %7zd %7zd     %-70s\n", groupOrder, numbMastersGroup, numbMastersGroup, buff);
					else
						fprintf(file, SHIFT "%-s\n", buff+k);

					if (!flg)
						break;

					memset(buff + k + j - SHIFT_TO_DECOMP, ' ', SHIFT_TO_DECOMP);
					buff[k + j] = saved;
					k += j - SHIFT_TO_DECOMP;
				}

				groupOrder = pRec->groupOrder;
				totalMasters += numbMastersGroup;
				totalCombined += numDecompGroup;
				numDecompGroup = numbMastersGroup = 0;
				pBuff = buff;
			}

			numbDecomp = pRec->numbDecomp;
			numbDecompMax = numbMasters = 1;
		}
	}

	memset(buff, '_', lenStr);
	buff[lenStr] = '\0';
	fprintf(file, SHIFT "%s\n", buff);
	fprintf(file, SHIFT "  Total: %7zd %7zd    MaxDecomp for master: %zd\n", totalMasters, totalCombined, numbDecompMaxGlobal);
}

void CDesignDB::outWithFormat_1(const size_t* pSortedRecords, FILE* file) const {
	char buff[256];
	const auto lenStr = sprintf_s(buff, SHIFT " |Aut(M)|:    Masters:    Decomp: ");
	fprintf(file, "\n%s\n", buff);
	memset(buff, '_', lenStr);
	buff[lenStr] = '\0';
	fprintf(file, SHIFT "%s\n", buff);

	const auto last = recNumb() - 1;
	size_t i = 0;
	auto* pRec = (const masterInfo*)getRecord(pSortedRecords[0]);
	unsigned long long totalCombined = 0, totalMasters = 0;

	auto groupOrder = pRec->groupOrder;
	auto numbDecomp = pRec->numbDecomp;
	size_t jMax = 1, numbMasters = 1;
	while (++i < recNumb()) {
		pRec = (masterInfo*)getRecord(pSortedRecords[i]);
		if (groupOrder == pRec->groupOrder && numbDecomp == pRec->numbDecomp) {
			numbMasters++;
			if (i != last)
				continue;
		}
		else {
			if (i == last)
				jMax = 2;
		}

		for (auto j = jMax; j--;) {
			fprintf(file, SHIFT "%7zd   %7zd    %7zd\n", groupOrder, numbMasters, numbDecomp);
			totalMasters += numbMasters;
			totalCombined += numbMasters * numbDecomp;
			groupOrder = pRec->groupOrder;
			numbDecomp = pRec->numbDecomp;
			numbMasters = 1;
		}
	}

	fprintf(file, SHIFT "%s\n", buff);
	fprintf(file, SHIFT "  Total:  %7zd    %7zd\n", totalMasters, totalCombined);
}

void CDesignDB::outWithFormat_2(const size_t * pSortedRecords, FILE * file) const {
	unsigned long long totalCombined = 0;
	fprintf(file, "\n" SHIFT "Master #:    Number of Decomp:   |Aut(M)|:\n");
	for (size_t i = 0; i < recNumb(); i++) {
		auto *pRec = (const masterInfo*)getRecord(pSortedRecords[i]);
		totalCombined += pRec->numbDecomp;
		fprintf(file, SHIFT "%4zd:          %6zd           %5zd\n", i + 1, pRec->numbDecomp, pRec->groupOrder);
	}
	fprintf(file, SHIFT "Total:        %7zd\n", totalCombined);
}



