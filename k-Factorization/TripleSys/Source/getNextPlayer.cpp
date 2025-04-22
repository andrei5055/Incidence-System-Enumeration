#include "TripleSys.h"

CC int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int iRet;
	int m0 = iPlayer % m_groupSize;
	const auto cmpGraph = param(t_CMP_Graph);
	for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
	{
	checkPlayerNumber:
		if (cmpGraph) {
			checkPlayerNumber1:
			int m1 = iPlayerNumber % m_groupSize;
			for (int k = 1; k <= m0; k++) {
				if ((tmpPlayers[iPlayer - k] % m_groupSize) == m1) {
					iPlayerNumber++;
					goto checkPlayerNumber1;
				}
			}
		}
		
		if (iPlayerNumber > m_indexPlayerMax[iPlayer])
			return m_numPlayers;

		if (selPlayers[iPlayerNumber] != unset)
			continue;

		if ((iRet = checkPlayer1(iPlayerNumber)) >= m_numPlayers)
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
				tchar* lnk = links(tmpPlayers[iPlayer - 1]);
				if (lnk[iPlayerNumber] != unset)
					continue;
				tmpPlayers[iPlayer] = iPlayerNumber;
				lnk[iPlayerNumber] = iDay;
				links(iPlayerNumber)[tmpPlayers[iPlayer - 1]] = iDay;
				break;
			}
			else if (m_groupSize == 3)
			{
				if ((iRet = checkPlayer3(iPlayerNumber, m_numPlayers)) >= m_numPlayers)
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
