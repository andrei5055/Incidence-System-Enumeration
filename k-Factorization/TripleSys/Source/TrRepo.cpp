#include "TrRepo.h"

#define CHECK_ROW_INDICES(firstRow, secondRow)  ASSERT(firstRow == secondRow); \
												ASSERT(firstRow < 0 || firstRow >= m_numRows); \
												ASSERT(secondRow < 0 || secondRow >= m_numRows);

CC CTrRepo::CTrRepo(int nRows, int numPlayers) : m_numRows(nRows) {
	const auto nDBs = nRows * nRows;
	m_pTrStorage = new CStorageIdx<tchar>* [nDBs];
	memset(m_pTrStorage, 0, nDBs * sizeof(m_pTrStorage[0]));
	int idx = 0;
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nRows; j++, idx++) {
			if (i == j)
				continue;

			m_pTrStorage[idx] = new CStorageIdx<tchar>(8, numPlayers);
		}
	}
}

CC CTrRepo::~CTrRepo() {
	for (int i = m_numRows * m_numRows; i--;)
		delete m_pTrStorage[i];

	delete[] m_pTrStorage;
}

CC void CTrRepo::addTr(ctchar* pTr, int firstRow, int secondRow) {
	CHECK_ROW_INDICES(firstRow, secondRow);
	m_pTrStorage[firstRow * m_numRows + secondRow]->updateRepo(pTr);
}

CC bool CTrRepo::initIterator(int firstRow, int secondRow) {
	CHECK_ROW_INDICES(firstRow, secondRow);
	m_iteratorIdx = 0;
	m_trStorage = m_pTrStorage[firstRow * m_numRows + secondRow];
	return m_trStorage->numObjects() > 0;
}

CC ctchar * CTrRepo::getNextTr() {
	if (m_trStorage->numObjects() <= m_iteratorIdx)
		return NULL;

	return m_trStorage->getObject(m_iteratorIdx++);
}

CC void CTrRepo::releaseTrs(int firstRow, int secondRow) const {
	CHECK_ROW_INDICES(firstRow, secondRow);
	m_pTrStorage[firstRow * m_numRows + secondRow]->releaseAllObjects();
}