#pragma once
#include <stdio.h>
#include "Sorter.h"
#include "DataTypes.h"

#define STARTING_DB_VOLUME	100
template<typename T>
struct masterInfo {
	T numbDecomp;
	T groupOrder;
	inline void setDesignNumber(T number)	{ numbDecomp = number; }
	inline auto designNumber() const						{ return numbDecomp; }
};

#define LEN_HEADER			sizeof(masterInfo<DB_INFO_DATA_TYPE>)

typedef const unsigned char * recPtr;

class CDesignDB : public CSorter<unsigned char>
{
public:
	CDesignDB(size_t len) : CSorter(len) {
		m_nRecNumbMax = STARTING_DB_VOLUME;
		m_pRecPermutation = new size_t[m_nRecNumbMax];
		setStorageOwner();
		setRecordStorage(new unsigned char[recordLength() * m_nRecNumbMax]);
	}
	~CDesignDB()								{ delete[] m_pRecPermutation; }
	size_t AddRecord(recPtr pRecord, DB_INFO_DATA_TYPE groupOrder, DB_INFO_DATA_TYPE numbDecomp = 1);
	void SortRecods(FILE* file = NULL);
	void combineDesignDBs(const CDesignDB* pDB_A, const CDesignDB* pDB_B, bool complFlag = false, bool intersecFlag = false);
	inline void resetRecNumb()					{ m_nRecNumb = 0; }
	inline auto recNumb() const					{ return m_nRecNumb; }
	inline size_t* getPermut() const			{ return m_pRecPermutation; }
private:
	size_t FindRecord(recPtr pRecord, int* pResCmp) const;
	bool reallocateMemory();
	void outWithFormat(const size_t * pSortedRecords, FILE * file) const;

	size_t m_nRecNumb = 0;					// number of DB records
	size_t m_nRecNumbMax = 0;				// max number of DB records
	size_t* m_pRecPermutation = NULL;		// permutation defining the order of records
};
