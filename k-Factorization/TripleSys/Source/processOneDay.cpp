#include "TripleSys.h"
CC bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	//if (iPlayer == 0)
	//	updateIndexPlayerMinMax();
	while (iPlayer < m_numPlayers)
	{
		if (iPlayer < 0)
			return false;
		if (iPlayer == 1)
			iPlayer = iPlayer;
		const auto iPlayerNumber = getNextPlayer();
		if (iPlayerNumber >= m_numPlayers)
		{
			getPrevPlayer();
			continue;
		}
		else if (iPlayerNumber >= 0)
		{
			if (iPlayer + 1 < m_numPlayers)// && iPlayerNumber != indexPlayer[iPlayer])
				indexPlayer[iPlayer + 1] = m_indexPlayerMin[iPlayer + 1];
			tmpPlayers[iPlayer] = indexPlayer[iPlayer] = iPlayerNumber;
			selPlayers[iPlayerNumber] = iPlayer++; // check values of selPlayers only for equal or not to unset (-1)
		}
		else
		{
			iPlayer = m_numPlayers;
			break;
		}
		if (iPlayer == m_numPlayers && iDay > 0)
		{
			u1fSetTableRow(neighbors(iDay), tmpPlayers);
			if (m_use2RowsCanonization || param(t_u1f))
			{
				memcpy(result(iDay), tmpPlayers, m_numPlayers);
				m_playerIndex = m_numPlayers * (iDay + 1) - m_groupSize - 1;
				if (/**m_groupSize <= 3 && **/ !matrixStat(neighbors(), iDay + 1))
				{
					while (iDay * m_numPlayers + iPlayer > m_playerIndex)
						getPrevPlayer();
				}
				m_playerIndex = 0;
			}
		}
	}
	return true;
}
