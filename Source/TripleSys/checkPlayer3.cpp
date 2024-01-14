#include <iostream>
#include "TripleSys.h"

int alldata::checkPlayer3(int iPlayerNumber, int lastPlayer)
{
	int m0 = iPlayer % GroupSize;
	if (m0 == 0)
		abort();
	int m1 = m0 == 0 ? GroupSize : 1;

	char i0 = iPlayer > 0 ? tmpPlayers[iPlayer - 1] : 0;
	char i1 = iPlayer > 1 ? tmpPlayers[iPlayer - 2] : 0;
	char* l0 = links(i0);
	char* l1 = links(i1);
	char day = iDay;

	if (l0[iPlayerNumber] != unset)
		return iPlayerNumber + 1;
#if UseSS == 0
	if (iPlayer == 1)
	{
		// AI statement #15
		for (int i = 1; i < iPlayerNumber; i++)
		{
			if (l0[i] == unset)
				return m_numPlayers;
		}
	}
#endif
	char* li = links(iPlayerNumber);
	if (m0 == 2)
	{
		if (l1[iPlayerNumber] != unset)
			return iPlayerNumber + 1;
		l1[iPlayerNumber] = li[i1] = day;
	}
	l0[iPlayerNumber] = li[i0] = day;
	if (m0 == 2)
	{
		bool bV = true, bH = true, bHdone = false;
#if UseSS == 0
		if (iPlayer > m_numPlayers / 2 && m_bCheckLinkH)
		{
			int nUnselected = m_numPlayers - iPlayer - 1;
			if (nUnselected > 0 && iDay > 2)
			{
				selPlayers[iPlayerNumber] = iPlayer;
				tmpPlayers[iPlayer] = iPlayerNumber;
				getUnselected(m_h, nUnselected);
				bHdone = bH = m_pCheckLink->checkLinksH(links(), NULL, NULL, 0, m_h, nUnselected, nUnselected, unset, unset, m_ho);
			}
		}
#endif
		if (m_bCheckLinkV && bH == true)
		{
			bV = m_pCheckLink->checkLinks(links(), iDay);
		}
		if (!bH || !bV)
		{
			selPlayers[iPlayerNumber] = unset;
			tmpPlayers[iPlayer] = unset;
			li[i0] = li[i1] = unset;
			l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
			return iPlayerNumber + 1;
		}
		if (bHdone)
		{
			for (int i = iPlayer + 1, j = 0; i < m_numPlayers; i++, j++)
			{
				char v = m_ho[j];
				if (!setLinksForOnePlayer(tmpPlayers, i, v))
					abort();
				indexPlayer[i] = v;
				if ((i % 3) == 2)
				{
					if (m_bCheckLinkV)
					{
						if (!m_pCheckLink->checkLinks(links(), iDay))
						{
							if (!unsetLinksForOnePlayer(tmpPlayers, i))
								abort();
							tmpPlayers[i] = unset;
							iPlayer = i;
							return -1;
						}
					}
				}

				selPlayers[v] = i;
			}
			iPlayer = m_numPlayers;
			return -1;
		}
	}
	return iPlayerNumber;
}