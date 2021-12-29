#pragma once
#include <stdio.h>

#define STARTING_DB_VOLUME	100
#define LEN_HEADER			sizeof(unsigned int)

typedef const unsigned char * recPtr;

class CDesignDB
{
public:
	CDesignDB(size_t len) : m_nRecLen(len+LEN_HEADER) {
		m_nRecNumbMax = STARTING_DB_VOLUME;
		m_pRecPermutation = new size_t[m_nRecNumbMax];
		m_pRecStorage = new unsigned char[recordLength() * m_nRecNumbMax];
	}
	~CDesignDB() {
		delete[] m_pRecPermutation;
		delete[] m_pRecStorage;
		delete[] m_pSortedRecords;
	}
	bool AddRecord(recPtr pRecord);
	void SortRecods(FILE* file = NULL);
private:
	size_t FindRecord(recPtr pRecord, int* pResCmp);
	bool reallocateMemory();
	void quickSort(size_t * arr, size_t left, size_t right) const;
	int compareRecords(const size_t idx, recPtr pSecnd) const;

	inline auto firstRecord() const				{ return m_pRecStorage; }
	inline size_t recordLength() const			{ return m_nRecLen; }
	inline auto recNumb() const					{ return m_nRecNumb; }

	const size_t m_nRecLen;					// length of each record
	size_t m_nRecNumb = 0;					// number of DB records
	size_t m_nRecNumbMax = 0;				// max number of DB records
	size_t* m_pRecPermutation = NULL;		// permutation defining the order of records
	unsigned char* m_pRecStorage = NULL;	// memory for record storing
	size_t *m_pSortedRecords = NULL;		// the array of indices of sorted records
};

