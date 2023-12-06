#include <iostream>
#include "TripleSys.h"

int alldata::checkPlayer1(int iPlayerNumber)
{
	int m0 = iPlayer % GroupSize;
	int m1 = m0 == 0 ? GroupSize : 1;
	int ifixedPlayer = -1;
	if (iDay == 0)
	{
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
			iPlayerNumber = iPlayer;
	}
	else if (GroupSize == 3)
	{
		if (iPlayer == 1) // iDay != 0
		{
			if (iPlayerNumber <= result(iDay - 1)[iPlayer])
				iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
		}
		//	if (!m_bCheckLinkV) // with m_bCheckLinkV equal true we canot check days because they are not sorted
		{
			if (iDay == 1)
			{
				if (iPlayer > 0 && iPlayer < 3)
				{
					ifixedPlayer = iPlayer * 3;
					if (ifixedPlayer >= iPlayerNumber && selPlayers[ifixedPlayer] == unset)
					{
						if (links(0)[ifixedPlayer] == unset)
							return ifixedPlayer;
					}
					return m_numPlayers;
				}
			}/**
			else if (iDay == 3)
			{
				if (iPlayer == 2 && iPlayerNumber <= result(iDay - 1)[iPlayer])
					iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
			}*/

		}
		if (iPlayer < 7)
		{
			if (m0 == 0) // iDay != 0
			{
				ifixedPlayer = iPlayer / 3;
				if (ifixedPlayer >= iPlayerNumber && selPlayers[ifixedPlayer] == unset)
					return ifixedPlayer;
				return m_numPlayers;
			}
		}
	}
	if (iPlayer >= m1)
	{
		// the following 2 lines makes signinficant change in speed
		if (iPlayerNumber <= tmpPlayers[iPlayer - m1])
			iPlayerNumber = tmpPlayers[iPlayer - m1] + 1;
	}

	if (m0 == 0)
	{
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
		int firstNotSel = 0;
		for (int i = 0; i < m_numPlayers; i++)
		{
			if (selPlayers[i] == unset)
			{
				firstNotSel = i;
				break;
			}
		}
		if (iPlayerNumber > firstNotSel)
			return m_numPlayers;
		else
			iPlayerNumber = firstNotSel;
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
	}
	return iPlayerNumber;
}