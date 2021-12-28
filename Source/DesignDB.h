#pragma once

#define STARTING_DB_VOLUME	100
#define LEN_HEADER			sizeof(unsigned int)

class CDesignDB
{
public:
	CDesignDB(size_t len) : m_nRecLen(len+LEN_HEADER) {
		m_nRecNumb = 0;
		m_nRecNumbMax = STARTING_DB_VOLUME;
		m_pRecPermutation = new size_t[m_nRecNumbMax];
		m_pRecStorage = new unsigned char[m_nRecLen * m_nRecNumbMax];
	}
	~CDesignDB() {
		delete[] m_pRecPermutation;
		delete[] m_pRecStorage;
	}
	bool AddRecord(const unsigned char* pRecord);

private:
	size_t FindRecord(const unsigned char* pRecord, int* pResCmp);
	bool reallocateMemory();

	size_t m_nRecNumb;				// number of DB records
	size_t m_nRecNumbMax;			// max number of DB records
	const size_t m_nRecLen;			// length of each record
	size_t* m_pRecPermutation;		// permutation defining the order of records
	unsigned char* m_pRecStorage;	// memory for record storing 
};

