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
		setStorageOwner();
		setRecordStorage(new unsigned char[recordLength() * m_nRecNumbMax]);
	}
	~CDesignDB() {
		delete[] m_pRecPermutation;
	}
	void AddRecord(recPtr pRecord, size_t groupOrder, size_t numbDecomp = 1);
	void SortRecods(FILE* file = NULL, int formatID = 0);
	void mergeDesignDB(const CDesignDB* pDB);
	void mergeDesignDBs(const CDesignDB* pDB_A, const CDesignDB* pDB_B);
	inline auto recNumb() const					{ return m_nRecNumb; }
	inline size_t* getPermut() const			{ return m_pRecPermutation; }
private:
	size_t FindRecord(recPtr pRecord, int* pResCmp) const;
	bool reallocateMemory();
	void outWithFormat_0(const size_t * pSortedRecords, FILE * file) const;
	void outWithFormat_1(const size_t * pSortedRecords, FILE * file) const;
	void outWithFormat_2(const size_t * pSortedRecords, FILE * file) const;

	size_t m_nRecNumb = 0;					// number of DB records
	size_t m_nRecNumbMax = 0;				// max number of DB records
	size_t* m_pRecPermutation = NULL;		// permutation defining the order of records
};
