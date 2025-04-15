#include "TripleSys.h"

CC int alldata::checkPlayer1(int iPlayerNumber)
{
	int m0 = iPlayer % m_groupSize;
	int m1 = m0 == 0 ? m_groupSize : 1;
	if (iDay <= 0)
	{
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
		return iPlayer;
	}
	// iDay > 0
	// for all group sizes:
	auto prevPlayer = iPlayer ? tmpPlayers[iPlayer - 1] : 0;
	if (iPlayer >= m1)
	{
		// AI #1 and #2 (for all groups)
		if (iPlayerNumber <= tmpPlayers[iPlayer - m1])
			iPlayerNumber = tmpPlayers[iPlayer - m1] + 1;
	}
	if (m0 == 0)
	{
		for (; m_firstNotSel < m_numPlayers; m_firstNotSel++)
		{
			if (selPlayers[m_firstNotSel] == unset)
				break;
		}
		if (iPlayerNumber > m_firstNotSel)
			return m_numPlayers;
		return m_firstNotSel;
	}
	else
	{
		iPlayerNumber = MAX2(prevPlayer + 1, iPlayerNumber);
		// do we need a separate cycle for m0 == 2?
		tchar* lnk = links(prevPlayer);
		const auto maxp = m_indexPlayerMax[iPlayer];
		for (; iPlayerNumber <= maxp; iPlayerNumber++)
			if (lnk[iPlayerNumber] == unset && selPlayers[iPlayerNumber] == unset)
				goto PlayerOk1;
		return m_numPlayers;
	}
PlayerOk1:
	switch (m_groupSize) {
	case 2:
		if (iPlayer < 3)// for indices=0,1,2 and group size=2 m_indexPlayerMin equal m_indexPlayerMax
			return iPlayerNumber;
		if (m0 && m_checkForUnexpectedCycle)
			iPlayerNumber = checkForUnexpectedCycle(iPlayerNumber, iPlayer, m_numPlayers, links(), tmpPlayers);
		break;
	case 3:
		if (iDay == 1)
		{
#if 0
			//if player[1, 4] == 4 ==> player[1, 5] <= player[0, { 5 }]
			if (iPlayer > 5 && tmpPlayers[4] == 4 && iPlayerNumber == 5 && iPlayer < tmpPlayers[5])
				return 6;
#endif
			if (!param(t_bipartiteGraph)) {
				// AI #4 (part)
				switch (iPlayer)
				{
				case 2: /*~1sec*/ return (iPlayerNumber <= 6) ? 6 : m_numPlayers;
				case 4: /* AI #17 */ return (iPlayerNumber <= 4) ? 4 : (m_numPlayers > 9 && iPlayerNumber <= 9) ? 9 : m_numPlayers;
				case 5: /* AI #14 */
					if (prevPlayer == 4)
						return (iPlayerNumber <= 7) ? 7 : (iPlayerNumber <= 9) ? 9 : m_numPlayers;
					break;
				case 7:
					// AI #7 
					if (iPlayerNumber <= tmpPlayers[4])
						iPlayerNumber = tmpPlayers[4] + 1;
					if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
						return m_numPlayers;
					// AI #19
					if (tmpPlayers[4] != 4)
					{
						if (iPlayerNumber <= 12)
							iPlayerNumber = 12;
					}
					break;
				case 8:
					if (prevPlayer == 5 && tmpPlayers[4] == 4 && iPlayerNumber <= tmpPlayers[5])
					{
						// AI #10
						iPlayerNumber = tmpPlayers[5] + 1;
					}
					break;
				case 9:
					// AI #8,9a 
					if (tmpPlayers[4] == 4 && tmpPlayers[7] != 5)
					{
						if (iPlayerNumber <= 5)
							return 5; // not happened
						return m_numPlayers; // not happened
					}
					break;
				}
			}
		}
	}
	if (!param(t_bipartiteGraph) && iDay == 1 && (iPlayerNumber % m_groupSize) && selPlayers[iPlayerNumber - 1] == unset)
		return iPlayerNumber + 1;
	return iPlayerNumber;
}