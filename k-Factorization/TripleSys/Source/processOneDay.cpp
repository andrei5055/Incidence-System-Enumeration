#include "TripleSys.h"
CC bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	//if (iPlayer == 0)
	//	updateIndexPlayerMinMax();
	unsigned int cnt = 0;
	while (iPlayer < m_numPlayers)
	{
		if (m_printMatrices && cnt++ > 300000000) {
#if !USE_CUDA
			printTableColor("processOneDay", tmpPlayers, 1, m_numPlayers, m_groupSize);
#endif
			cnt = 0;
		}
#if 0
		static tchar a[] = { 0, 4,10,  1, 8, 9,  2, 7,13,  3,11,12,  5, 6, };
		if (memcmp(a, tmpPlayers, 14) == 0 && iDay == 2)
			printTableColor("processOneDay", tmpPlayers, 1, m_numPlayers, m_groupSize);
		iDay = iDay;
#endif
		if (iPlayer < 0)
			return false;
		/* const */ auto iPlayerNumber = getNextPlayer();
		if ((m_test & 8) && iPlayer == 3 && iDay >= 2 && param(t_CBMP_Graph) == 2 && m_firstCycleSet && m_firstCycleSet[0] == 4) {
			if (iPlayerNumber >= tmpPlayers[1])
				iPlayerNumber = m_numPlayers;
			else iPlayerNumber = tmpPlayers[1] - 1;
		}

		if (iPlayerNumber >= m_numPlayers)
		{
			getPrevPlayer();
			continue;
		}
		else if (iPlayerNumber >= 0)
		{
			if (iPlayer + 1 < m_numPlayers)// && iPlayerNumber != indexPlayer[iPlayer])
				indexPlayer[iPlayer + 1] = m_indexPlayerMin[iPlayer + 1];
			tmpPlayers[iPlayer] = indexPlayer[iPlayer] = iPlayerNumber;
			selPlayers[iPlayerNumber] = iPlayer++; // check values of selPlayers only for equal or not to unset (-1)
		}
		else
		{
			iPlayer = m_numPlayers;
			break;
		}
		if (iPlayer == m_numPlayers && iDay > 0)
		{
			u1fSetTableRow(neighbors(iDay), tmpPlayers);
			if (m_use2RowsCanonization || param(t_u1f))
			{
				memcpy(result(iDay), tmpPlayers, m_numPlayers);
				//printTableColor("r2", result(1), 1, m_numPlayers, m_groupSize);
				m_playerIndex = m_numPlayers * (iDay + 1) - m_groupSize - 1;
				if (/**m_groupSize <= 3 && **/ !checkNewRow(neighbors(), iDay + 1))
				{
					while (iDay * m_numPlayers + iPlayer > m_playerIndex)
						getPrevPlayer();
				}
				m_playerIndex = 0;
			}
		}
	}
	return true;
}
