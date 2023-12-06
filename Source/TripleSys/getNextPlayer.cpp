#include <iostream>
#include "TripleSys.h"

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int m0 = iPlayer % GroupSize;
	int m1 = m0 == 0 ? GroupSize : 1;
	if ((iPlayerNumber = checkPlayer1(iPlayerNumber)) >= m_numPlayers)
		return m_numPlayers;

	if (m0 != 0)
	{
		if (m_groupSize == 3)
		{
			if ((iPlayerNumber = checkPlayer3(iPlayerNumber, m_numPlayers)) >= m_numPlayers)
				return m_numPlayers;
		}
		else
		{
			for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
			{
				if (selPlayers[iPlayerNumber] != unset)
					continue;

				tmpPlayers[iPlayer] = iPlayerNumber;
				if (setLinksForOnePlayer(tmpPlayers, iPlayer, 1))
				{
					break;
				}
				tmpPlayers[iPlayer] = unset;
			}
		}
	}

	return iPlayerNumber;
}