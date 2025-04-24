#include "TripleSys.h"
CC bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	while (iPlayer < m_numPlayers)
	{
#if 0
		if (iDay == 2) {
			tchar c[] = { 0,4,9, 1,3,15, 2,6,20, 5,12,16, 7,14,18, 8,10,19, 11,13,17 };
			if (memcmp(c, tmpPlayers, 15) == 0) {
				printf(" %2d:%2d ", iDay, iPlayer);
				printTable("p", tmpPlayers, 1, m_numPlayers, 3);
				if (iPlayer == 14 && memcmp(c, tmpPlayers, 14) == 0)
					iDay = iDay;
				if (iPlayer > 14 && memcmp(c, tmpPlayers, 15) <= 0)
					iDay = iDay;
			}
		}
#endif
		if (iPlayer < 0)
			return false;
		/**
		if (m_bCheckLinkT && !checkLinksT(links(), iDay))
		{
			getPrevPlayer();
			continue;
		}*/

		const auto iPlayerNumber = getNextPlayer();
		
		if (iPlayerNumber >= m_numPlayers)
		{
			//printf("d=%d p=%d go back\n", iDay, iPlayer);
			getPrevPlayer();
			continue;
		}
		else if (iPlayerNumber >= 0)
		{
			//if (iDay == 3)
			//	printf("%d p=%d v=%d\n", iDay, iPlayer, iPlayerNumber);
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
				m_playerIndex = m_numPlayers * (iDay + 1) - m_groupSize - 1;
#if 1
				if (m_groupSize <= 3 && !matrixStat(neighbors(), iDay + 1))
				{
					while (iDay * m_numPlayers + iPlayer > m_playerIndex)
						getPrevPlayer();
				}
				m_playerIndex = 0;
#else
				if ((param(t_u1f) && m_groupSize == 2 && p1fCheck(iDay + 1, tmpPlayers) >= 0) ||
					(!param(t_u1f) && m_groupSize <= 3 && !matrixStat(neighbors(), iDay + 1)))
				{
					while (iDay * m_numPlayers + iPlayer > m_playerIndex)
						getPrevPlayer();
				}
#endif
			}
		}
	}
	return true;
}
