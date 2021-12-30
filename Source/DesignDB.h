#pragma once
#include <stdio.h>
#include "Sorter.h"

#define STARTING_DB_VOLUME	100
typedef struct {
	size_t numbDecomp;
	size_t groupOrder;
} masterInfo;

#define LEN_HEADER			sizeof(masterInfo)

typedef const unsigned char * recPtr;

class CDesignDB : public CSorter<unsigned char>
{
public:
	CDesignDB(size_t len) : CSorter(len+LEN_HEADER) {
		m_nRecNumbMax = STARTING_DB_VOLUME;
		m_pRecPermutation = new size_t[m_nRecNumbMax];
		setRecordStorage(new unsigned char[recordLength() * m_nRecNumbMax]);
	}
	~CDesignDB() {
		delete[] m_pRecPermutation;
		delete[] m_pSortedRecords;
	}
	void AddRecord(recPtr pRecord, size_t groupOrder);
	void SortRecods(FILE* file = NULL);
private:
	size_t FindRecord(recPtr pRecord, int* pResCmp);
	bool reallocateMemory();
//	void quickSort(size_t * arr, size_t left, size_t right) const;
//	int compareRecords(const size_t idx, recPtr pSecnd) const;

	inline auto recNumb() const					{ return m_nRecNumb; }

	size_t m_nRecNumb = 0;					// number of DB records
	size_t m_nRecNumbMax = 0;				// max number of DB records
	size_t* m_pRecPermutation = NULL;		// permutation defining the order of records
	size_t *m_pSortedRecords = NULL;		// the array of indices of sorted records
};

