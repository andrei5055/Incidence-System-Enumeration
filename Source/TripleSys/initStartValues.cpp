#include <iostream>
#include "TripleSys.h"
void alldata::linksFromMatrix(char* iv, int id)
{
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

			if (!setLinksForOnePlayer(res, j, ivId))
			{
				printf("Init: value of %d (for day %d position %d) already defined in links table\n", ivId, i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}
		}
	}
}

bool alldata::initStartValues(const char* ivcb, bool printStartValues)
{
	char* iv = m_co; // We can use existing array m_co
	int v;
	int ind = 0, lastInd = 0;
	int id = iDay = 0;
	memset(iv, unset, m_nLenResults);

	for (ind = 0; ; ind++)
	{
		while (*ivcb == ' ' || *ivcb == '\n')
		{
			if (*ivcb++ == '\n')
			{
				if (id + 1 >= m_numDays)
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
			if (id + 1 >= m_numDays)
			{
				printf("Init: too many input values, data ignored after pos=%d, day=%d\n", m_numPlayers - 1, id);
				goto doneInit;
			}
			id++;
		}
		*(iv + id * m_numPlayers + ind) = (char)v;
		lastInd = ind;
		while ((*ivcb >= '0' && *ivcb <= '9') || *ivcb == '-')
			ivcb++;
	}
doneInit:
	if (ind <= 0 && id <= 0)
		return false;

	if (lastInd == m_numPlayers - 1)
		id++;
	else
	{
		printf("Init: values for day %d positions %d-%d not defined\n", id, lastInd + 1, m_numPlayers - 1);
		printTable("Initial result", result(0), m_numDays, m_numPlayers);
		exit(0);
	}
	linksFromMatrix(iv, id);
	iDay = id;
	maxDays = iDay;
	memcpy(maxResult, result(0), m_nLenResults);
	if (printStartValues) {
		printTable("Start Result", result(), numDays(), numPlayers());
		printTableColor("Start Links", links(0), numPlayers(), numPlayers());
	}
	return true;
}

