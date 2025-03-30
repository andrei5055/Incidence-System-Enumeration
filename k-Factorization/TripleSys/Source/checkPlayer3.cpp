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
	if (iPlayer == 1 && (m_useRowsPrecalculation != eCalculateRows || iDay != 3))
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
		l0[iPlayerNumber] = li[i0] = day;
		if (m_bCheckLinkV && m_useRowsPrecalculation != eCalculateRows /* && iPlayerNumber == 11 */ && !checkLinks(links(), iDay))
		{
			selPlayers[iPlayerNumber] = unset;
			tmpPlayers[iPlayer] = unset;
			li[i0] = li[i1] = unset;
			l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
			return iPlayerNumber + 1;
		}
	}
	else
		l0[iPlayerNumber] = li[i0] = day;
	return iPlayerNumber;
}
