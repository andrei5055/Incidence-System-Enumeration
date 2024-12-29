#include "TripleSys.h"

CC bool alldata::cnvCheck45(int nrows)
{
	if (nrows < 2)
		return true;
	tchar tgStartOffsets4[] = {
		0, 1, 2, 3,  1, 0, 3, 2,  2, 3, 0, 1,  3, 2, 1, 0 };
	tchar tgStartOffsets5[] = {
		0, 1, 2, 3, 4,  1, 0, 3, 4, 2,  2, 3, 4, 0, 1,   3, 4, 1, 2, 0,  4, 2, 0, 1, 3 };
	tchar tgStart4[] = {
		0, 0, 0, 0,  0, 1, 2, 3,  0, 2, 3, 1,  0, 3, 1, 2 };
	tchar tgStart5[] = {
		0, 0, 0, 0, 0,  0, 1, 2, 3, 4,  0, 2, 3, 4, 1,  0, 3, 4, 1, 2,  0, 4, 1, 2, 3 };
	
	int nTr0 = m_groupSize * (m_groupSize - 1);
	tchar* tr = new tchar[m_numPlayers * nTr0];
	tchar* tgStartOffsets = m_groupSize == 4 ? tgStartOffsets4 : tgStartOffsets5;
	tchar* tgStart = m_groupSize == 4 ? tgStart4 : tgStart5;
	for (int i = 0; i < m_groupSize - 1; i++) {
		for (int j = 0; j < m_groupSize; j++) {
			for (int k = 0; k < m_groupSize; k++) {
				for (int m = 0; m < m_groupSize; m++) {
					tr[(i * m_groupSize + j) * m_numPlayers + k * m_groupSize + m] = 
						m_groups[i * m_groupSize + k] * m_groupSize + tgStartOffsets[tgStart[j * m_groupSize + k] * m_groupSize + m];
				}
			}
		}
	}
	for (int n = 0; n < m_groupSizeFactorial; n++) {
		const auto* pGroups1 = m_groups + n * m_groupSize;
		for (int j = 0; j < m_groupSizeFactorial; j++) {
			const auto* pGroups = m_groups + j * m_groupSize;
			for (int i = 0; i < nTr0; i++) {
				tchar ttr[MAX_PLAYER_NUMBER];
				tchar* tri = tr + i * m_numPlayers;
				for (int k = 0; k < m_numPlayers; k += m_groupSize) {
					for (int m = 0; m < m_groupSize; m++)
						ttr[k + m] = tri[k + pGroups[m]];
				}
				tchar ttrr[MAX_PLAYER_NUMBER];
				for (int mm = 0; mm < m_numPlayers; mm += m_groupSize) {
					memcpy(ttrr + mm, ttr + pGroups1[mm / m_groupSize] * m_groupSize, m_groupSize);
				}

				int iRet = cnvCheckKm1(ttrr, nrows, NULL);

				//if (ttrr[0] == 4 && ttrr[1] == 5 && ttrr[2] == 6)
				//	printTable("ttrr", ttrr, 1, m_numPlayers, m_groupSize);
				if (iRet < 0)
				{
					delete[] tr;
					return false;
				}
				m_TrInd++;
			}
		}
	}
	delete[] tr;
	return true;
}