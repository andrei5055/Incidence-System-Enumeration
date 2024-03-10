#include <iostream>
#include "TripleSys.h"

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int iRet;
	int m0 = iPlayer % m_groupSize;
	for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
	{
		if (selPlayers[iPlayerNumber] != unset)
			continue;

		if ((iRet = checkPlayer1(iPlayerNumber)) >= m_numPlayers)
			return m_numPlayers;

		if (iPlayerNumber != iRet)
		{
			iPlayerNumber = iRet - 1;
			continue;
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
						m0 = iPlayer % GroupSize;
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
				if (setLinksForOnePlayer(tmpPlayers, iPlayer, (char)iPlayerNumber))
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