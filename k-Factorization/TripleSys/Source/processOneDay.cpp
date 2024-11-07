#include "TripleSys.h"
CC bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	if (iPlayer >= m_numPlayers)
		return true;
	while (iPlayer < m_numPlayers)
	{
		/**
		if (iDay == 3) {
			tchar c[] = { 0,4,9,  1,5,11,  2,7,13, 3,8,12, 6,10,14 };
			tchar d[] = { 0,6,13, 1,10,12, 2,4,14, 3,7,11, 5,6,9 };
			printf(" %2d:%2d ", iDay, iPlayer);
			printTable("p", tmpPlayers, 1, m_numPlayers, 3);
			if (memcmp(c, result(2), m_numPlayers) == 0)
				iDay = iDay;
			if (memcmp(d, tmpPlayers, m_numPlayers - 13) == 0)
				iDay = iDay;
		}**/
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
			indexPlayer[iPlayer] = iPlayerNumber;
			tmpPlayers[iPlayer] = iPlayerNumber;
			selPlayers[iPlayerNumber] = iPlayer; // check values of selPlayers only for equal or not to unset (-1)
			iPlayer++;
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
