#include <iostream>
#include "TripleSys.h"

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
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
		if (iDay == 1)
		{
			if (iPlayer > 0 && iPlayer < 3)
			{
				ifixedPlayer = iPlayer * 3;
				if (ifixedPlayer >= iPlayerNumber && selPlayers[ifixedPlayer] == unset)
				{
					tmpPlayers[iPlayer] = ifixedPlayer;
					if (setLinksForOnePlayer(tmpPlayers, iPlayer, 1))
						return ifixedPlayer;
					tmpPlayers[iPlayer] = unset;
				}
				return m_numPlayers;
			}
		}
		else if (iDay == 3)
		{
			if (iPlayer == 2 && iPlayerNumber <= result(iDay - 1)[iPlayer])
				iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
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
			else if (iPlayer == 1 && !m_bCheckLinkV) // iDay != 0
			{
				if (iPlayerNumber <= result(iDay - 1)[iPlayer])
					iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
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
			return firstNotSel;
	}
	else
	{
		char i0 = iPlayer > 0 ? tmpPlayers[iPlayer - 1] : 0;
		char i1 = iPlayer > 1 ? tmpPlayers[iPlayer - 2] : 0;
		char* l0 = links(i0);
		char* l1 = links(i1);
		char day = iDay;
		for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
		{
			if (selPlayers[iPlayerNumber] != unset)
				continue;
			if (m_groupSize == 3) {
				/**  3 4 5 7 <= x1 < x2 < ... < xd **/
				if (iPlayer == 1 && iDay > 0)
				{
					if (iDay >= 4)
					{
						if (iPlayerNumber < 7)
							iPlayerNumber = 7;
						if (iPlayerNumber <= result(iDay - 1)[1])
							iPlayerNumber = result(iDay - 1)[1];
					}
					else
					{
						if (iPlayerNumber <= iDay + 2)
							iPlayerNumber = iDay + 2;
						else
							return m_numPlayers;
					}
				}
				/** p7 > p4 in the Second day
				4 <= z1 <= 11, 5 <= z2 <= 14,
				and if z1 <= 8, then z2 <= 11 **/
				if (iDay == 1)
				{
					if (iPlayer == 4)
					{
						if (iPlayerNumber < 4)
							iPlayerNumber = 4;
						if (iPlayerNumber > 11)
							return m_numPlayers;
					}
					else if (iPlayer == 7)
					{
						if (iPlayerNumber <= tmpPlayers[4])
							iPlayerNumber = tmpPlayers[4] + 1;
						if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
							return m_numPlayers;
						if (iPlayerNumber > 14)
							return m_numPlayers;
					}
					if (selPlayers[iPlayerNumber] != unset)
						continue;
				}
				if (l0[iPlayerNumber] != unset)
					continue;
				char* li = links(iPlayerNumber);
				if (m0 == 2)
				{
					if (l1[iPlayerNumber] != unset)
						continue;
					l1[iPlayerNumber] = li[i1] = day;
				}
				l0[iPlayerNumber] = li[i0] = day;
				if (m0 == 2)
				{
					bool bV = true, bH = true;
					if (iPlayer > m_numPlayers / 2 && m_bCheckLinkH)
					{
						int nUnselected = m_numPlayers - iPlayer - 1;
						if (nUnselected > 0 && iDay > 2)
						{
							selPlayers[iPlayerNumber] = iPlayer;
							tmpPlayers[iPlayer] = iPlayerNumber;
							getUnselected(m_h, nUnselected);
							bH = m_pCheckLink->checkLinksH(links(), m_h, nUnselected, unset, unset, m_ho);
						}
					}
					if (m_bCheckLinkV && bH == true)
					{
						bV = m_pCheckLink->checkLinks(links(), iDay);
					}
					if (!bH)
						bH = bH;
					if (!bH || !bV)
					{
						selPlayers[iPlayerNumber] = unset;
						tmpPlayers[iPlayer] = unset;
						li[i0] = li[i1] = unset;
						l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
						continue;
					}
				}
				break;
			}
			else {
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