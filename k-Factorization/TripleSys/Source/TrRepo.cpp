#include "TrRepo.h"

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

