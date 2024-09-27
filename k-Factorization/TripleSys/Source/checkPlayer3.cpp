#include "TripleSys.h"

CC int alldata::checkPlayer3(int iPlayerNumber, int lastPlayer)
{
	const auto m0 = iPlayer % m_groupSize;
	ASSERT(!m0);

	auto i0 = iPlayer > 0 ? tmpPlayers[iPlayer - 1] : 0;
	auto i1 = iPlayer > 1 ? tmpPlayers[iPlayer - 2] : 0;
	auto* l0 = links(i0);
	auto* l1 = links(i1);
	char day = iDay;

	if (l0[iPlayerNumber] != unset)
		return iPlayerNumber + 1;
	if (iPlayer == 1)
	{
		// AI statement #15
		for (int i = 1; i < iPlayerNumber; i++)
		{
			if (l0[i] == unset)
				return m_numPlayers;
		}
	}

	auto* li = links(iPlayerNumber);
	if (m0 == 2)
	{
		if (l1[iPlayerNumber] != unset)
			return iPlayerNumber + 1;
		l1[iPlayerNumber] = li[i1] = day;
	}
	l0[iPlayerNumber] = li[i0] = day;
	if (m0 == 2)
	{
		bool bH = true, bHdone = false;
		int ipToUseLinksH = m_numPlayers - m_numPlayers / 3;
		//if (!m_p1f && iPlayer > m_numPlayers / 2 && iDay > 2 && m_pCheckLinksH)
		if (iPlayer > ipToUseLinksH && iDay > 2 && m_pCheckLinksH)
		{
			const int nUnselected = m_numPlayers - iPlayer - 1;
			if (nUnselected > 0)
			{
				selPlayers[iPlayerNumber] = iPlayer;
				tmpPlayers[iPlayer] = iPlayerNumber;
				if (iPlayer + 1 < m_numPlayers && iPlayerNumber != indexPlayer[iPlayer])
					indexPlayer[iPlayer + 1] = m_indexPlayerMin[iPlayer + 1];
				indexPlayer[iPlayer] = iPlayerNumber;
				getUnselected(m_h, nUnselected);
				ASSERT(nUnselected != m_numPlayers - iPlayer - 1);
				memcpy(m_ho, tmpPlayers, iPlayer + 1);
				int pind = m_playerIndex = m_numPlayers * iDay + iPlayer;
				bHdone = bH = (this->*m_pCheckLinksH)(m_h, nUnselected, nUnselected, -1, -1, m_ho + iPlayer + 1);

				if (!bH)
				{
					tmpPlayers[iPlayer] = selPlayers[iPlayerNumber] = unset;
					li[i0] = li[i1] = unset;
					l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
					if (pind <= m_playerIndex)
						return iPlayerNumber + 1;
					while (iDay * m_numPlayers + iPlayer > m_playerIndex) // iPlayer changed by getPrevPlayer()
						getPrevPlayer();
					return -1;
				}
			}
		}
		if (m_bCheckLinkV && !checkLinks(links(), iDay))
		{
			selPlayers[iPlayerNumber] = unset;
			tmpPlayers[iPlayer] = unset;
			li[i0] = li[i1] = unset;
			l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
			return iPlayerNumber + 1;
		}
		if (bHdone)
		{
			for (int i = iPlayer + 1; i < m_numPlayers; i++)
			{
				const auto v = m_ho[i];
				auto linksOK = setLinksForOnePlayer(iDay, links(), tmpPlayers, i, v);
				ASSERT(!linksOK);
				indexPlayer[i] = v;
				if ((i % 3) == 2)
				{
					if (m_bCheckLinkV)
					{
						if (!checkLinks(links(), iDay))
						{
							linksOK = unsetLinksForOnePlayer(tmpPlayers, i);
							ASSERT(!linksOK);
							tmpPlayers[iPlayer = i] = unset;
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
