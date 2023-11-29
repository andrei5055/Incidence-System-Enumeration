
#include <iostream>
#include "TripleSys.h"


void alldata::initCurrentDay()
{
	iPlayer = 0;
	memset(indexPlayer, 0, m_numPlayers);
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);
	//checkbmask(selPlayers, bmask);
}


bool alldata::setLinksAndDevCounts(char* p, int ip, char iset)
{
	const int i = ip % GroupSize;
	if (i == 0)
		return true;
	//	if (iset != 1 && p[ip] < 0)
	//		return true;
	char bset = iset == 1 ? iDay : unset;
	char i1 = p[ip];
	//	if (i1 < 0 || i1 >= m_numPlayers)
	//		i1 = i1;
	auto* linkPtr = links(i1);
	if (bset != unset)
	{
		for (int j = 1; j <= i; j++)
		{
			char i2 = p[ip - j];
			//		if (i2 < 0 || i2 >= m_numPlayers)
			//			i2 = i2;
			if (linkPtr[i2] != unset)
				return false;

		}
	}
	for (int j = 1; j <= i; j++)
	{
		const char i2 = p[ip - j];
		linkPtr[i2] = *(links(i2) + i1) = bset;
	}
	return true;
}

void alldata::initPrevDayProcess()
{
	if (m_pCheckLink)
	{
		if (iDay < 0)
			return;
	}
	else
	{
		if (iDay < 1)  // keep first line (first day result)
		{
			iDay = -1;
			return;
		}
	}

	auto* const pRes = result(iDay);
	memcpy(indexPlayer, pRes, m_numPlayers);
	memcpy(tmpPlayers, pRes, m_numPlayers);

	if (UseLastSix)
	{
		iPlayer = m_numPlayers - GroupSize * 2;
		indexPlayer[iPlayer] = index6[iDay]; // only one of six indices used, last 5 values of array indexPlayers not used
	}
	else
		iPlayer = m_numPlayers - GroupSize - 1;

	if (iPlayer < 0)
	{
		iDay = -1;
		return;
	}

	int ind = indexPlayer[iPlayer];
	if (ind < 0)
		abort();
	else
	{
		memset(selPlayers, 1, m_numPlayers);

		for (int j = iPlayer; j < m_numPlayers; j++)
		{
			if (*(pRes + j) < 0)
				abort();
			if (!setLinksAndDevCounts(pRes, j, unset))
				abort();
			int k = tmpPlayers[j];
			tmpPlayers[j] = selPlayers[k] = unset;
			indexPlayer[j] = 0;
		}
	}
	if (UseLastSix && iPlayer == m_numPlayers - 6)
		index6[iDay] = ind + 1;
	indexPlayer[iPlayer] = ind + 1;
}

bool alldata::initStartValues(const char* ivcb, bool printStartValues)
{
	char* iv = m_co; // We can use existing array m_co
	int v;
	int ind = 0;
	int id = iDay = 0;
	memset(iv, unset, m_nLenResults);

	for (ind = 0; ; ind++)
	{
		while (*ivcb == ' ' || *ivcb == '\n')
		{
			if (*ivcb++ == '\n')
			{
				if (id >= m_numDays - 1)
					goto doneInit;
				ind = 0;
				id++;
			}
		}

		if (sscanf_s(ivcb, "%d", &v) != 1)
			break;
		if (ind >= m_numPlayers)
		{
			ind = 0;
			if (++id >= m_numDays)
				break;
		}
		*(iv + id * m_numPlayers + ind) = (char)v;
		while ((*ivcb >= '0' && *ivcb <= '9') || *ivcb == '-')
			ivcb++;
	}
doneInit:
	if (ind <= 0 && id <= 0)
		return false;

	char* iv_id = iv;
	auto* res = result(0);
	for (int i = 0; i < id; i++, iv_id += m_numPlayers, res += m_numPlayers)
	{
		iDay = i;
		for (int j = 0; j < m_numPlayers; j++)
		{
			const auto ivId = iv_id[j];
			if (ivId == unset)
			{
				printf("Init: value for day %d position %d not devined\n", i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}

			*(res + j) = ivId;
			if (!setLinksAndDevCounts(res, j, 1))
			{
				printf("Init: value of %d (for day %d position %d) already devined in links table\n",
					ivId, i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}
		}
		index6[i] = getLastSixIndex(res);
	}

	iDay = id;
	if (printStartValues) {
		printTable("Result", result(), numDays(), numPlayers());
		printTable("Links", links(0), numPlayers(), numPlayers());
	}
	return true;
}

void alldata::getPrevPlayer()
{
	if (iPlayer >= m_numPlayers)
		abort();

	indexPlayer[iPlayer] = 0;
	while (--iPlayer >= 0)
	{
		int iPlayerNumber;
		if (!setLinksAndDevCounts(tmpPlayers, iPlayer, unset))
			abort();
		if (tmpPlayers[iPlayer] < 0 || tmpPlayers[iPlayer] >= m_numPlayers)
			abort();
		iPlayerNumber = tmpPlayers[iPlayer];
		//checkbmask(selPlayers, bmask);
		selPlayers[iPlayerNumber] = unset;
		//checkbmask(selPlayers, bmask);
		tmpPlayers[iPlayer] = unset;

		if (UseLastSix && iPlayer >= m_numPlayers - 6)
			abort();

		if (iPlayer >= m_numPlayers - GroupSize || iPlayerNumber >= m_numPlayers)
		{
			indexPlayer[iPlayer] = 0;
			continue;
		}

		iPlayerNumber++;
		indexPlayer[iPlayer] = iPlayerNumber;
		break;
	}
}

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int m0 = iPlayer % GroupSize;
	int m1 = m0 == 0 ? GroupSize : 1;
	int ifixedPlayer = -1;

	if (!m_pCheckLink)
	{
		if (iDay == 0)
		{
			if (iPlayerNumber > iPlayer)
				return m_numPlayers;
			iPlayerNumber = iPlayer;
		}
		else if (GroupSize == 3)
		{
			if (iDay == 1 && iPlayer > 0 && iPlayer < 3)
			{
				ifixedPlayer = iPlayer * 3;
				if (ifixedPlayer >= iPlayerNumber)
				{
					tmpPlayers[iPlayer] = ifixedPlayer;
					if (setLinksAndDevCounts(tmpPlayers, iPlayer, 1))
						return ifixedPlayer;
					tmpPlayers[iPlayer] = unset;
				}
				return m_numPlayers;
			}
			if (iPlayer < 7 && m0 == 0 && iDay > 0)
			{
				ifixedPlayer = iPlayer / 3;
				if (ifixedPlayer >= iPlayerNumber)
					return ifixedPlayer;
				return m_numPlayers;
			}
			else if (iDay == 3 && iPlayer == 2 && iPlayerNumber <= result(iDay - 1)[iPlayer])
				iPlayerNumber = result(iDay - 1)[iPlayer] + 1;

			/** attemt to sort days by the second columnt **/
			else if (iPlayer == 1 && iDay > 0)
			{
				if (iPlayerNumber <= result(iDay - 1)[iPlayer])
					iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
			}
			/**/
		}
	}

	if (iPlayer >= m1)
	{
		// the following 2 lines makes signinficant change in speed
		if (iPlayerNumber <= tmpPlayers[iPlayer - m1])
			iPlayerNumber = tmpPlayers[iPlayer - m1] + 1;


		/** p7 > p4 in the Second day
		4 <= z1 <= 11, 5 <= z2 <= 14,
		and if z1 <= 8, then z2 <= 11 **/
		if (GroupSize == 3 && iDay == 1)
		{
			if (iPlayer == 4)
			{
				if (iPlayerNumber < 4)
					iPlayerNumber = 4;
				if (iPlayerNumber > 11)
					return m_numPlayers;
			}
			else if (iPlayer == 7)
			{
				if (iPlayerNumber <= tmpPlayers[4])
					iPlayerNumber = tmpPlayers[4] + 1;
				if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
					return m_numPlayers;
				if (iPlayerNumber > 14)
					return m_numPlayers;
			}
		}
	}

	if (m0 == 0)
	{
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
		int firstNotSel = 0;
		for (int i = 0; i < m_numPlayers; i++)
		{
			if (selPlayers[i] == unset)
			{
				firstNotSel = i;
				break;
			}
		}
		if (iPlayerNumber > firstNotSel)
			return m_numPlayers;
		else
			return firstNotSel;
	}
	else
	{
		char i0 = iPlayer > 0 ? tmpPlayers[iPlayer - 1] : 0;
		char i1 = iPlayer > 1 ? tmpPlayers[iPlayer - 2] : 0;
		char* l0 = links(i0);
		char* l1 = links(i1);
		char day = iDay;
		for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
		{
			if (selPlayers[iPlayerNumber] != unset)
				continue;
			if (m_numPlayers) { // GroupSize == 3
				if (l0[iPlayerNumber] != unset)
					continue;
				char* li = links(iPlayerNumber);
				if (m0 == 2)
				{
					if (l1[iPlayerNumber] != unset)
						continue;
					l1[iPlayerNumber] = li[i1] = day;
				}
				l0[iPlayerNumber] = li[i0] = day;

				/**/
				if (m0 == 2 && m_pCheckLink)
				{
					if (!m_pCheckLink->checkLinks(links(), iDay))
					{
						li[i0] = li[i1] = unset;
						l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
						continue;
					}
				}
				/**/
				break;
			}
			else {
				tmpPlayers[iPlayer] = iPlayerNumber;
				if (setLinksAndDevCounts(tmpPlayers, iPlayer, 1))
				{
					break;
				}
				tmpPlayers[iPlayer] = unset;
			}
		}
	}

	return iPlayerNumber;
}

