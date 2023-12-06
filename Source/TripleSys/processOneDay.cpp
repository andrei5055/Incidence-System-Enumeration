#include <iostream>
#include "TripleSys.h"
bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	while (iPlayer < m_numPlayers)
	{
		//if (iDay == 3)
		//	printf("%d p=%d\n", iDay, iPlayer);

		//if (iDay == 3 && iPlayer == 1)
		//	iDay = iDay;
		int iPlayerNumber = 0;
		if (iPlayer < 0)
			return false;

		if (UseLastSix && iPlayer == m_numPlayers - 6)
		{
			iPlayerNumber = processLastSix();
			if (iPlayer == m_numPlayers)
				break;
		}
		else
			iPlayerNumber = getNextPlayer();

		if (iPlayerNumber >= m_numPlayers)
		{
			//printf("d=%d p=%d go back\n", iDay, iPlayer);
			getPrevPlayer();
			continue;
		}
		//if (iDay == 3)
		//	printf("%d p=%d v=%d\n", iDay, iPlayer, iPlayerNumber);
		if (iPlayer + 1 < m_numPlayers && iPlayerNumber != indexPlayer[iPlayer])
			indexPlayer[iPlayer + 1] = 0;
		indexPlayer[iPlayer] = iPlayerNumber;
		tmpPlayers[iPlayer] = iPlayerNumber;
		selPlayers[iPlayerNumber] = iPlayer; // check values of selPlayers only for equal or not to unset (-1)
		iPlayer++;
	}
	//if (iDay == 3)
	//	printf("return d=%d p=%d\n", iDay, iPlayer);
	return true;
}
