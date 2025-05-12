#include "TripleSys.h"

CC int alldata::getNextPlayer()
{
	int iRet;
	int m0 = m_groupSizeRemainder[iPlayer];
	int iPrevPlayer = m0 ? tmpPlayers[iPlayer - 1] : iPlayer ? tmpPlayers[iPlayer - m_groupSize] : -1;
	const auto cbmpGraph = !completeGraph();
	int iPlayerNumber = MAX2(indexPlayer[iPlayer], iPrevPlayer + 1);
	int iPlayerMax = m_indexPlayerMax[iPlayer];
	for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
	{
	checkPlayerNumber:
		if (cbmpGraph) {
			checkPlayerNumber1:
			int m1 = m_groupSizeRemainder[iPlayerNumber];
			for (int k = 1; k <= m0; k++) {
				if (m_groupSizeRemainder[tmpPlayers[iPlayer - k]] == m1) {
					iPlayerNumber++;
					goto checkPlayerNumber1;
				}
			}
		}
		
		if (iPlayerNumber > iPlayerMax)
			return m_numPlayers;

		if (selPlayers[iPlayerNumber] != unset)
			continue;

		if ((iRet = checkPlayer1(iPlayerNumber)) > iPlayerMax)
				return m_numPlayers;

		if (iPlayerNumber != iRet)
		{
			iPlayerNumber = iRet;
			goto checkPlayerNumber;
		}

		if (m0 != 0)
		{
			if (m_groupSize == 2)
			{
				tchar* lnk = links(iPrevPlayer);
				if (lnk[iPlayerNumber] != unset)
					continue;
				tmpPlayers[iPlayer] = iPlayerNumber;
				lnk[iPlayerNumber] = iDay;
				links(iPlayerNumber)[iPrevPlayer] = iDay;
				break;
			}
			else if (m_groupSize == 3)
			{
				if ((iRet = checkPlayer3(iPlayerNumber, m_numPlayers)) > iPlayerMax)
					return m_numPlayers;
				if (iPlayerNumber == iRet)
					break;
				iPlayerNumber = iRet;
				goto checkPlayerNumber;
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
