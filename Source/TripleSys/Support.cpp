#include <iostream>
#include "TripleSys.h"

bool alldata::initCurrentDay()
{
	iPlayer = 0;
	memset(indexPlayer, 0, m_numPlayers);
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);

#if UseSS == 0

	if (m_bCheckLinkH && iDay > 1)
	{
		const int np = m_numPlayers;
		for (int i = 0; i < m_numPlayers; i++)
			m_h[i] = i;
		if (iDay == 0)
		{
			memcpy(m_ho, m_h, np);
		}
		else
		{
			if (!m_pCheckLink->checkLinksH(links(), m_h, m_numPlayers, np, unset, result(iDay - 1)[1], m_ho))
			{
				bPrevResult = true;
				return false;
			}
			//printf("day=%d ", iDay);
			//printTable("m_ho", m_ho, 1, m_numPlayers);
		}
		memcpy(tmpPlayers, m_ho, np);
		memcpy(indexPlayer, m_ho, np);
		iPlayer = m_numPlayers;

		for (int j = 0; j < m_numPlayers; j++)
		{
			char k = tmpPlayers[j];
			if (!setLinksForOnePlayer(iDay, m_numPlayers, links(), tmpPlayers, j, k))
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
	}
#endif // UseSS
	return true;
}
// setLinksForOnePlayer need to be available outside of alldata
bool setLinksForOnePlayer(int id, int np, char* lnk, char* p, int ip, char v)
{
	const int i = ip % GroupSize;
	if (i != 0)
	{
		char bset = id;
		auto* linkPtr = lnk + v * np;
		for (int j = 1; j <= i; j++)
		{
			char i2 = p[ip - j];
			if (linkPtr[i2] != unset)
				return false;
		}
		for (int j = 1; j <= i; j++)
		{
			const char i2 = p[ip - j];
			linkPtr[i2] = *(lnk + i2 * np + v) = bset;
		}
	}
	p[ip] = v;
	return true;
}

bool alldata::unsetLinksForOnePlayer(char* p, int ip) const
{
	const int i = ip % GroupSize;
	if (i == 0)
		return true;
	char i1 = p[ip];
	auto* linkPtr = links(i1);
	for (int j = 1; j <= i; j++)
	{
		const char i2 = p[ip - j];
		linkPtr[i2] = *(links(i2) + i1) = unset;
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

	if (m_numPlayers > GroupSize)
		iPlayer = m_numPlayers - GroupSize - 1;
	else
		iPlayer = m_numPlayers - 1;

	int ind = indexPlayer[iPlayer];
	if (ind < 0)
		abort();
	else
	{
		for (int j = 0; j < m_numPlayers; j++)
		{
			if (*(pRes + j) < 0)
				abort();
			int k = tmpPlayers[j];
			if (j < iPlayer)
			{
				selPlayers[k] = j;
			}
			else
			{
				if (!unsetLinksForOnePlayer(pRes, j))
					abort();
				tmpPlayers[j] = selPlayers[k] = unset;
				indexPlayer[j] = 0;
			}
		}
	}
	indexPlayer[iPlayer] = ind + 1;
	return true;
}

void alldata::getPrevPlayer()
{
	if (iPlayer > m_numPlayers)
		abort();
	if (iPlayer < m_numPlayers)
	    indexPlayer[iPlayer] = 0;
	while (--iPlayer >= 0)
	{
		int iPlayerNumber = tmpPlayers[iPlayer];
		if (!unsetLinksForOnePlayer(tmpPlayers, iPlayer))
			abort();
		tmpPlayers[iPlayer] = selPlayers[iPlayerNumber] = unset;

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


