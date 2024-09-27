#include "TripleSys.h"
CC bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	if (iPlayer >= m_numPlayers)
		return true;
	while (iPlayer < m_numPlayers)
	{
		//printf(" %d:%d", iDay, iPlayer);
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
			if (iPlayer + 1 < m_numPlayers && iPlayerNumber != indexPlayer[iPlayer])
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
			p1fSetTableRow(p1ftable(iDay), tmpPlayers);
			if (m_p1f)
			{
				memcpy(result(iDay), tmpPlayers, m_numPlayers);
				m_playerIndex = m_numPlayers * (iDay + 1) - m_groupSize - 1;
#if 1
				if (m_groupSize <= 3 && !matrixStat(p1ftable(), iDay + 1))
				{
					while (iDay * m_numPlayers + iPlayer > m_playerIndex)
						getPrevPlayer();
				}
#else
				if ((param(t_u1f) && m_groupSize == 2 && p1fCheck(iDay + 1, tmpPlayers) >= 0) ||
					(!param(t_u1f) && m_groupSize <= 3 && !matrixStat(p1ftable(), iDay + 1)))
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
