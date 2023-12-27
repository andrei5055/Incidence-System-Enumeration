
#include <iostream>
#include "TripleSys.h"


bool alldata::initCurrentDay()
{
	iPlayer = 0;
	memset(indexPlayer, 0, m_numPlayers);
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);

#if UseSS == 0

	if (m_bCheckLinkH)
	{
		int np = UseLastSix ? m_numPlayers - 6 : m_numPlayers;
		//int np = m_numPlayers;
		for (int i = 0; i < m_numPlayers; i++)
			m_h[i] = i;
		if (iDay == 0)
		{
			memcpy(m_ho, m_h, np);
		}
		else
		{/**
			double cnt = 0;
			m_pCheckLink->checkLinksH(links(), m_h, m_numPlayers, np, unset, result(iDay-1)[1], m_ho, &cnt);
			printf("d = % d n = % .0f\n", iDay, cnt);
			**/
			if (!m_pCheckLink->checkLinksH(links(), m_h, m_numPlayers, np, unset, result(iDay - 1)[1], m_ho))
			{
				bPrevResult = true;
				//printf("day=%d\n", iDay);
				//printTable("no result", result(iDay), 1, np);
				return false;
			}
		}

		memcpy(result(iDay), m_ho, np);
		auto* const pRes = result(iDay);
		memcpy(indexPlayer, pRes, m_numPlayers);
		memcpy(tmpPlayers, pRes, m_numPlayers);

		if (UseLastSix)
		{
			iPlayer = m_numPlayers - GroupSize * 2;
			indexPlayer[iPlayer] = index6[iDay] = 0; // only one of six indices used, last 5 values of array indexPlayers not used
		}
		else if (m_numPlayers > GroupSize)
			iPlayer = m_numPlayers - GroupSize - 1;
		else
			iPlayer = m_numPlayers - 1;
   
		memset(selPlayers, unset, m_numPlayers);

		for (int j = 0; j < m_numPlayers; j++)
		{
			int k = tmpPlayers[j];
			if (j < iPlayer)
			{
				if (!setLinksForOnePlayer(pRes, j, 1))
				{
					if (iDay == 0)
					{
						bPrevResult = true;
						return false;
					}
					abort();
				}
				selPlayers[k] = j;
			}
			else
			{
				tmpPlayers[j] = selPlayers[k] = unset;
				indexPlayer[j] = 0;
			}
		}
	}
#endif // UseSS
	return true;
}


bool alldata::setLinksForOnePlayer(const char* p, int ip, char iset) const
{
	const int i = ip % GroupSize;
	if (i == 0)
		return true;
	char bset = iset == 1 ? iDay : unset;
	char i1 = p[ip];
	auto* linkPtr = links(i1);
	if (bset != unset)
	{
		for (int j = 1; j <= i; j++)
		{
			char i2 = p[ip - j];
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

bool alldata::initPrevDay()
{
	if (iDay >= 0 && iDay < m_numDays)
		memset(result(iDay), 0, m_numPlayers);
	iDay--;
	bPrevResult = false;

	if (m_bCheckLinkV)
	{
		if (iDay < 1) //0)
			return false;
	}
	else
	{
		if (iDay < 1)  // keep first line (first day result)
		{
			iDay = -1;
			return false;
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
	else if (m_numPlayers > GroupSize)
		iPlayer = m_numPlayers - GroupSize - 1;
	else
		iPlayer = m_numPlayers - 1;

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
			if (!setLinksForOnePlayer(pRes, j, unset))
				abort();
			int k = tmpPlayers[j];
			tmpPlayers[j] = selPlayers[k] = unset;
			indexPlayer[j] = 0;
		}
	}
	if (UseLastSix)
		index6[iDay] = ind + 1;
	indexPlayer[iPlayer] = ind + 1;
	return true;
}

void alldata::getPrevPlayer()
{
	if (iPlayer >= m_numPlayers)
		abort();

	if (UseLastSix && iPlayer > m_numPlayers - 6)
		abort();

	indexPlayer[iPlayer] = 0;
	while (--iPlayer >= 0)
	{
		int iPlayerNumber;
		if (!setLinksForOnePlayer(tmpPlayers, iPlayer, unset))
			abort();
		if (tmpPlayers[iPlayer] < 0 || tmpPlayers[iPlayer] >= m_numPlayers)
			abort();
		iPlayerNumber = tmpPlayers[iPlayer];
		//checkbmask(selPlayers, bmask);
		selPlayers[iPlayerNumber] = unset;
		//checkbmask(selPlayers, bmask);
		tmpPlayers[iPlayer] = unset;

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


