
#include <iostream>
#include "TripleSys.h"


void alldata::initCurrentDay()
{
	iPlayer = 0;
	memset(indexPlayer, 0, m_numPlayers);
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);

	return;

	if (m_bCheckLinkH)
	{
		static int mx = -1;
		double c = 0.0;
		//if (mx != iDay)
		{
			mx = iDay;
			for (int i = 0; i < m_numPlayers; i++)
				m_h[i] = i;
			if (iDay == 0)
			{
				memcpy(m_ho, m_h, m_numPlayers);
			}
			else if (!m_pCheckLink->checkLinksH(links(), m_h, m_numPlayers, unset, unset, m_ho))
			{
				bPrevResult = true;
				initPrevDay();
				//printf("day=%d\n", iDay);
				//printTable("no result", result(iDay), 1, m_numPlayers);
				return;
			}
		}
	}
}


bool alldata::setLinksForOnePlayer(char* p, int ip, char iset)
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
		if (iDay < 0)
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
	else
		iPlayer = m_numPlayers - GroupSize - 1;

	if (iPlayer < 0)
	{
		iDay = -1;
		return false;
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
			if (!setLinksForOnePlayer(pRes, j, unset))
				abort();
			int k = tmpPlayers[j];
			tmpPlayers[j] = selPlayers[k] = unset;
			indexPlayer[j] = 0;
		}
	}
	if (UseLastSix && iPlayer == m_numPlayers - 6)
		index6[iDay] = ind + 1;
	indexPlayer[iPlayer] = ind + 1;
	return true;
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
				printf("Init: value for day %d position %d not defined\n", i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}

			*(res + j) = ivId;
			if (!setLinksForOnePlayer(res, j, 1))
			{
				printf("Init: value of %d (for day %d position %d) already defined in links table\n",
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
		if (!setLinksForOnePlayer(tmpPlayers, iPlayer, unset))
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


