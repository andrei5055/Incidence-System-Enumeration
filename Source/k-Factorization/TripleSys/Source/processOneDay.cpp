#include <iostream>
#include "TripleSys.h"
bool alldata::processOneDay()
{   // returns: false - go to prev day, true - day processed
	while (iPlayer < m_numPlayers)
	{
#if 0
			static char c[] = { 0,   3,   6,    1,   4,   7,    2,   5 };
			if (iDay == 1 && iPlayer == 7)
			{
				if (memcmp(tmpPlayers, c, 7) == 0)
				{
					iDay = iDay;
				}
				else
				{
					iDay = iDay;
				}
			}
		//	printf("%d p=%d\n", iDay, iPlayer);
#endif

		if (iPlayer < 0)
			return false;

		int iPlayerNumber = getNextPlayer();
		
		if (iPlayerNumber >= m_numPlayers)
		{
			//printf("d=%d p=%d go back\n", iDay, iPlayer);
			getPrevPlayer();
			continue;
		}
		else if (iPlayerNumber >= 0)
		{
			//if (iDay == 3)
			//	printf("%d p=%d v=%d\n", iDay, iPlayer, iPlayerNumber);
			if (iPlayer + 1 < m_numPlayers && iPlayerNumber != indexPlayer[iPlayer])
				indexPlayer[iPlayer + 1] = 0;
			indexPlayer[iPlayer] = iPlayerNumber;
			tmpPlayers[iPlayer] = iPlayerNumber;
			selPlayers[iPlayerNumber] = iPlayer; // check values of selPlayers only for equal or not to unset (-1)
			iPlayer++;
		}
		else
		{
			iPlayer = m_numPlayers;
			return true;
		}
	}
	//if (iDay == 3)
	//	printf("return d=%d p=%d\n", iDay, iPlayer);
	return true;
}
