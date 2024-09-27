#include "TripleSys.h"

CC int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int iRet;
	int m0 = iPlayer % m_groupSize;
	
	for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
	{
	checkPlayerNumber:
		if (iPlayerNumber > m_indexPlayerMax[iPlayer])
			return m_numPlayers;

		if (selPlayers[iPlayerNumber] != unset)
			continue;

		if ((iRet = checkPlayer1(iPlayerNumber)) >= m_numPlayers)
			return m_numPlayers;

		if (iPlayerNumber != iRet)
		{
			iPlayerNumber = iRet;
			goto checkPlayerNumber;
		}

		if (m0 != 0)
		{
			if (m_groupSize == 3)
			{
				if ((iRet = checkPlayer3(iPlayerNumber, m_numPlayers)) >= m_numPlayers)
					return m_numPlayers;
				if (iPlayerNumber == iRet)
					break;
				if (iRet < 0)
				{
					if (iPlayer < m_numPlayers)
					{
						m0 = iPlayer % m_groupSize;
						iPlayerNumber = indexPlayer[iPlayer];
						continue;
					}
					return iRet;
				}
				iPlayerNumber = iRet - 1;
				continue;
			}
			else
			{
				if (setLinksForOnePlayer(iDay, links(), tmpPlayers, iPlayer, (tchar)iPlayerNumber))
				{
					break;
				}
				continue;
			}
		}
		else
			break;
	}
	return iPlayerNumber;
}
