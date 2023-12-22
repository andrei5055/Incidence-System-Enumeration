#include <iostream>
#include "TripleSys.h"

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int iRet;
	int m0 = iPlayer % GroupSize;
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
				iPlayerNumber = iRet - 1;
				continue;
			}
			else
			{
				tmpPlayers[iPlayer] = iPlayerNumber;
				if (setLinksForOnePlayer(tmpPlayers, iPlayer, 1))
				{
					break;
				}
				tmpPlayers[iPlayer] = unset;
				continue;
			}
		}
		else
			break;
	}
	return iPlayerNumber;
}