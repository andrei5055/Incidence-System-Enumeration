#include <iostream>
#include "TripleSys.h"

int alldata::checkPlayer1(int iPlayerNumber)
{
	int m0 = iPlayer % GroupSize;
	int m1 = m0 == 0 ? GroupSize : 1;
	int ifixedPlayer = -1;
	if (m0 == 0)
	{
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
	}
	if (iDay <= 0)
	{
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
		return iPlayer;
	}
	// iDay > 0
	if (iPlayer == 0)
	{
		// AI statement #4 (part 1)
		if (iPlayerNumber != 0)
			return m_numPlayers;
		return 0;
	}
	if (iPlayer >= m1)
	{
		// AI statement #1 and #2 (for all groups)
		if (iPlayerNumber <= tmpPlayers[iPlayer - m1])
			iPlayerNumber = tmpPlayers[iPlayer - m1] + 1;
	}
	if (iPlayer == 1)
	{
		// AI statement #6
		if (iPlayerNumber <= result(iDay - 1)[1])
			iPlayerNumber = result(iDay - 1)[1] + 1;
		if (GroupSize == 2)
		{
			return (iPlayerNumber > iDay + 1) ? m_numPlayers : iDay + 1;
		}
	}
	if (GroupSize == 2)
	{
		if (iPlayer == 2)
			return (iPlayerNumber > 1) ? m_numPlayers : 1;
		if (iDay == 1 && m_numPlayers == 12)
		{
			if (iPlayer == m_numPlayers - 2)
				return (iPlayerNumber > 9) ? m_numPlayers : 9;
			else if (iPlayer == m_numPlayers - 1)
				return (iPlayerNumber > 11) ? m_numPlayers : 11;
		}
	}
	else if (GroupSize == 3)
	{
#if UseSS != 2
		if (iPlayer < 7 && m0 == 0)
		{
			// AI statement #4 (part 2)
			ifixedPlayer = iPlayer / 3;
			if (ifixedPlayer >= iPlayerNumber)
				return ifixedPlayer;
			return m_numPlayers;
		}
//#if UseSS == 0
		{
			if (iPlayer == 1)
			{
				// AI statement #6
				if (iDay >= 4)
				{
					if (iPlayerNumber < 7)
						iPlayerNumber = 7;  // not happend
				}
				else
				{
					if (iPlayerNumber <= iDay + 2)
						iPlayerNumber = iDay + 2;
					else
						return m_numPlayers;
				}
				//new
				if (iPlayerNumber > m_numPlayers - m_numDays + iDay - 3) //for n=15 in last day player#1 cant be max 11, for 21 - 17
					return m_numPlayers;
			}
			else if (iDay == 1)
			{
#if 0
				{
					static char exp[] = { 0,   3,   6,    1,   4,   9,   2,   7,  12,    5,  10,  13,    8,  11,  14 };
					if (memcmp(exp, tmpPlayers, 14) == 0)
						iPlayer = iPlayer;
				}
#endif
#if 1
				//if player[1, 4] == 4 ==> player[1, 5] <= player[0, { 5 }]
				if (iPlayer > 5 && tmpPlayers[4] == 4 && iPlayerNumber == 5 && iPlayer < tmpPlayers[5])
					return 6;
#endif
				// AI statement #4 (part)
				switch (iPlayer) 
				{
					case 2:
					{
						ifixedPlayer = 6;
						if (ifixedPlayer >= iPlayerNumber)
						{
							return ifixedPlayer;
						}
						return m_numPlayers;
					}
#if 0
					case 4:
					{
						if (iPlayerNumber < 4)
							iPlayerNumber = 4;
						if (iPlayerNumber > 11)
							return m_numPlayers;
						break;
					}
#else
					case 4:
					{
						// ai 17
						if (iPlayerNumber <= 4)
							return 4;
						if (m_numPlayers > 9 && iPlayerNumber <= 9)
							return 9;
						return m_numPlayers;
					}
#endif
					case 5:
					{
#if 1
						// AI statement #14
						if (tmpPlayers[4] == 4)
						{
							if (iPlayerNumber <= 7)
								return 7;
							else if (iPlayerNumber <= 9)
								return 9; // good for n=9?
							else
								return m_numPlayers;
						}
#endif
						break;
					}
					case 7:
					{
						// AI statement #7 
						if (iPlayerNumber <= tmpPlayers[4])
							iPlayerNumber = tmpPlayers[4] + 1;
						if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
							return m_numPlayers;
						if (iPlayerNumber > 14)
							return m_numPlayers; // not happed
						// ai 19
						if (tmpPlayers[4] != 4)
						{
							if (iPlayerNumber <= 12)
								iPlayerNumber = 12;
						}
						break;
					}
					case 8:
					{
						if (tmpPlayers[7] == 5 && tmpPlayers[4] == 4 && iPlayerNumber <= tmpPlayers[5])
						{
							// AI statement #10
							iPlayerNumber = tmpPlayers[5] + 1;
						}
						/** covered in case 7
						if (iPlayerNumber <= tmpPlayers[4])
							iPlayerNumber = tmpPlayers[4] + 1; // not happend **/
						//if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
						//	return m_numPlayers;
						if (iPlayerNumber > 14)
							return m_numPlayers;
						break;
					}
					case 9:
					{
						// AI statement #8 
						if (tmpPlayers[4] != 4)
						{
							if (iPlayerNumber <= 4)
								return 4;
							else
								return m_numPlayers;
						}

						// AI statement #9 part a
						else if (tmpPlayers[7] != 5)
						{
							if (iPlayerNumber <= 5)
								return 5; // not happend
							return m_numPlayers; // not happend
						}
						break;
					}
					case 12:
					{
						if (tmpPlayers[9] == 4)
						{
							// AI statement #9 part b
							if (iPlayerNumber <= 5)
								return 5;
							return m_numPlayers;
						}
						break;
					}
					case 13:
					{
						if (tmpPlayers[12] == 5 && tmpPlayers[9] == 4 && iPlayerNumber <= tmpPlayers[10])
						{
							// AI statement #11 ???
							iPlayerNumber = tmpPlayers[10] + 1; // not happend
						}
						break;
					}
				}
				/**/
				// AI statement #18 (for group size == 3)
				if ((iPlayerNumber % 3) == 2)
				{
					if (selPlayers[iPlayerNumber - 2] == unset || selPlayers[iPlayerNumber - 1] == unset ||
						selPlayers[iPlayerNumber - 2] >= selPlayers[iPlayerNumber - 1])
						return iPlayerNumber + 1;
				}
				/**/
			}
		}
#endif
	}

	if (m0 == 0)
	{
		int firstNotSel = 0;
		for (int i = 0; i < m_numPlayers; i++)
		{
			if (selPlayers[i] == unset)
			{
				firstNotSel = i;
				break;
			}
		}
		// new 3
		if (iPlayerNumber < iPlayer / m_groupSize)
			iPlayerNumber = iPlayer / m_groupSize; // not happen

		if (iPlayerNumber > firstNotSel)
			return m_numPlayers; // happen for 21, 10
		else
			iPlayerNumber = firstNotSel; // happend for 21, 10
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers; // not happen if at the end
		//new 2
		if (GroupSize == 3)
		{
			if (iPlayerNumber > m_numPlayers - 6 - (m_numPlayers - iPlayer) / 3)
				return m_numPlayers; // not happen if at the end
		}
	}
	return iPlayerNumber;
}