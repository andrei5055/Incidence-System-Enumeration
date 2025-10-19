#pragma once

#include "k-SysSupport.h"
#include "CudaAttributes.h"
#include "Storage.h"

#define CHECK_ROW_INDICES(firstRow, secondRow)  ASSERT_IF(firstRow == secondRow); \
												ASSERT_IF(firstRow < 0 || firstRow >= m_numRows); \
												ASSERT_IF(secondRow < 0 || secondRow >= m_numRows);

class CTrRepo {
public:
	CC CTrRepo(int nRows, int numPlayers);
	CC ~CTrRepo();
	CC inline auto getTrSet(int firstRow, int secondRow) const {
		CHECK_ROW_INDICES(firstRow, secondRow);
		return m_pTrStorage[firstRow * m_numRows + secondRow];
	}
	CC inline void addTr(ctchar* pTr, int firstRow, int secondRow) {
		getTrSet(firstRow, secondRow)->updateRepo(pTr);
	}
	CC bool initIterator(int firstRow, int secondRow) {
		m_iteratorIdx = 0;
		m_trStorage = getTrSet(firstRow, secondRow);
		return m_trStorage->numObjects() > 0;
	}
	CC inline auto * getNextTr() {
		return m_trStorage->numObjects() <= m_iteratorIdx?
			NULL : m_trStorage->getObject(m_iteratorIdx++);
	}
	CC inline void releaseTrs(int firstRow, int secondRow) const {
		getTrSet(firstRow, secondRow)->releaseAllObjects();
	}
private:
	const int m_numRows;
	CStorageIdx<tchar>** m_pTrStorage = NULL;
	int m_iteratorIdx;
	CStorageIdx<tchar>* m_trStorage;
};
