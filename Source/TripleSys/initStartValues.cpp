#include <iostream>
#include "TripleSys.h"

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
				if (id >= m_numDays)
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
	maxDays = iDay;
	memcpy(maxResult, result(0), m_nLenResults);
	if (printStartValues) {
		printTable("Result", result(), numDays(), numPlayers());
		printTableColor("Links", links(0), numPlayers(), numPlayers());
	}
	return true;
}

