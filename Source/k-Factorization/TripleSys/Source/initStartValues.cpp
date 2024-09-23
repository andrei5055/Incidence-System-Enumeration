#include <iostream>
#include "TripleSys.h"
void linksFromMatrix(char* lnk, char* iv, int nr, int np)
{
	char* iv_id = iv;
	memset(lnk, unset, np * np);
	for (int i = 0; i < nr; i++, iv_id += np)
	{
		for (int j = 0; j < np; j++)
		{
			const auto ivId = iv_id[j];
			if (ivId == unset)
			{
				printf("Init: value for day %d position %d not defined\n", i, j);
				printTable("Initial result", iv, nr, np);
				exit(0);
			}

			if (!setLinksForOnePlayer(i, np, lnk, iv_id, j, ivId))
			{
				printf("Init: value of %d (for day %d position %d) already defined in links table\n", ivId, i, j);
				printTable("Initial result", iv, nr, np);
				exit(0);
			}
		}
	}
}

bool alldata::initStartValues(const char* ivcb, bool printStartValues)
{
	char* iv = result();
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
	linksFromMatrix(links(), result(), id, m_numPlayers);
	iDay = id;
	maxDays = iDay;
	if (printStartValues) {
		printTable("Start Result", result(), numDays(), numPlayers());
		printTableColor("Start Links", links(0), numPlayers(), numPlayers());
	}
	return true;
}

