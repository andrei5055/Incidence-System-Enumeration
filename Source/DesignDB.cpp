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

void CDesignDB::SortRecods(FILE* file) {
	if (!recNumb())
		return;   // nothing to sort;

	setCompareFunc(compareRecords);
	const auto* pSortedRecords = Sort(recNumb());
	if (file)
		outWithFormat(pSortedRecords, file);
}


void CDesignDB::outWithFormat(const size_t * pSortedRecords, FILE * file) const {
#define SHIFT			""
#define LEN_DECOMP		62
#define SHIFT_TO_DECOMP	30
#define MASTER			 0
#define DECOMP			 1
#define LEN_BUFFER    95

	size_t len_buff = LEN_BUFFER;
	char* buff = new char[len_buff];
	char *pBuff = buff;
	const auto lenStr = LEN_BUFFER-1;
	fprintf(file, "\n%s\n", SHIFT " |Aut(M)|: Masters:  Decomp:        Decomposition Look:");
	memset(buff, '_', lenStr);
	buff[lenStr] = '\0';
	fprintf(file, SHIFT "%s\n", buff);

	const auto last = recNumb() - 1;
	size_t i = 0;
	auto* pRec = (const masterInfo*)getRecord(pSortedRecords[i]);
	unsigned long long totalCombined = 0, totalMasters = 0;

	auto groupOrder = pRec->groupOrder;
	size_t numb[2] = { 1, pRec->numbDecomp };

	size_t numbDecompMaxGlobal, numbDecompMax, jMax, numbMastersGroup, numDecompGroup;
	numDecompGroup = numbMastersGroup = 0;
	numbDecompMaxGlobal = numbDecompMax = jMax = 1;
	while (++i < recNumb()) {
		pRec = (const masterInfo*)getRecord(pSortedRecords[i]);
		if (groupOrder == pRec->groupOrder && numb[DECOMP] == pRec->numbDecomp) {
			numb[MASTER]++;
			if (i != last)
				continue;
		}
		else {
			if (i == last)
				jMax = 2;
		}

		for (auto j = jMax; j--;) {
			numDecompGroup += numb[MASTER] * numb[DECOMP];
			numbMastersGroup += numb[MASTER];
			for (int i = numb[MASTER] == 1 ? 1 : 0; i < 2; i++) {
				char format[16], tmpBuff[32];
				const auto used = (pBuff != buff && (!i || numb[MASTER] == 1))? strcpy_s(format, " + "), 3 : 0;
				sprintf_s(format + used, sizeof(format) - used, "%%zd%s", i ? "" : "*");
				const auto len = sprintf_s(tmpBuff, format, numb[i]);
				if (len_buff - (pBuff - buff) <= len) {
					// Reallocating buffer
					len_buff = 2 * (len_buff + len);
					char* buffTmp = new char[len_buff];
					strcpy_s(buffTmp, len_buff, buff);
					pBuff = buffTmp + (pBuff - buff);
					delete[] buff;
					buff = buffTmp;

				}
				strcpy_s(pBuff, len_buff - (pBuff - buff), tmpBuff);
				pBuff += len;
			}

			if (numbDecompMax < numb[DECOMP]) {
				if (numbDecompMaxGlobal < (numbDecompMax = numb[DECOMP]))
					numbDecompMaxGlobal = numbDecompMax;
			}

			if (groupOrder != pRec->groupOrder || !j && i == last) {
				char saved;
				int k = 0;
				int j = LEN_DECOMP;
				while (true) {
					const bool flg = strlen(buff + k) > j;
					if (flg) {
						// Need to split decomposition information for given groupOrder into two strings
						while (buff[k+j] != '+') j--;
						saved = buff[k + ++j];
						buff[k + j] = '\0';
					}

					if (!k)
						fprintf(file, SHIFT "%7zd   %8zd %8zd    %-s\n", groupOrder, numbMastersGroup, numDecompGroup, buff);
					else 
						fprintf(file, SHIFT "%-s\n", buff+k);

					if (!flg)
						break;

					buff[k + j] = saved;
					k += j - SHIFT_TO_DECOMP;
					memset(buff + k, ' ', SHIFT_TO_DECOMP);
					j = LEN_DECOMP + SHIFT_TO_DECOMP;
				}

				groupOrder = pRec->groupOrder;
				totalMasters += numbMastersGroup;
				totalCombined += numDecompGroup;
				numDecompGroup = numbMastersGroup = 0;
				pBuff = buff;
			}

			numb[MASTER] = numbDecompMax = 1;
			numb[DECOMP] = pRec->numbDecomp;
		}
	}

	memset(buff, '_', lenStr);
	buff[lenStr] = '\0';
	fprintf(file, SHIFT "%s\n", buff);
	fprintf(file, SHIFT "  Total:  %8zd %8zd    MaxDecomp for master: %zd\n", totalMasters, totalCombined, numbDecompMaxGlobal);
	delete[] buff;
}


