#pragma once

#include "k-SysSupport.h"
#include "CudaAttributes.h"
#include "Storage.h"

class CTrRepo {
public:
	CC CTrRepo(int nRows, int numPlayers);
	CC ~CTrRepo();
	CC void addTr(ctchar* pTr, int firstRow, int secondRow);
	CC bool initIterator(int firstRow, int secondRow);
	CC ctchar* getNextTr();
	CC void releaseTrs(int firstRow, int secondRow) const;
private:
	const int m_numRows;
	CStorageIdx<tchar>** m_pTrStorage = NULL;
	int m_iteratorIdx;
	CStorageIdx<tchar>* m_trStorage;
};
