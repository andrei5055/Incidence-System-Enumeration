#include "TripleSys.h"
CC void SizeParam::linksFromMatrix(tchar* lnk, ctchar* iv, int nr) const
{ // ANDREI ??? ctchar* iv
	auto* iv_id = (tchar*)iv;
	const auto np = m_numPlayers;
	memset(lnk, unset, np * np);
	for (int i = 0; i < nr; i++, iv_id += np)
	{
		for (int j = 0; j < np; j++)
		{
			const auto ivId = iv_id[j];
			ASSERT_(ivId == unset,
				printfRed("*** Init: value for day %d position %d not defined\n", i, j);
				printTable("Initial result", iv, nr, np, m_groupSize);
				myExit(1)
			)

			const auto linksOK = setLinksForOnePlayer(i, lnk, iv_id, j, ivId);
			if (!linksOK) {
				for (int k = 0; k < j % m_groupSize; k++) {
					const auto iv0 = iv_id[j - k - 1];
					const auto lDay = *(lnk + ivId * np + iv0);
					if (lDay != unset) {
#if !USE_CUDA
						printfRed("*** Init: pair (%d,%d) in day %d already defined in day %d\n", iv0, ivId, i, lDay);
						printTableColor("Initial result", iv, nr, np, m_groupSize);
						myExit(1);
#endif
					}
				}
			}
		}
	}
}

bool alldata::initStartValues(const char* ivcb, bool printStartValues)
{
	auto* iv = result();
	int v;
	int ind = 0, lastInd = 0;
	int id = iDay = 0;
	memset(iv, unset, m_nLenResults);
	if (memcmp(ivcb, " 0 1 ", 4) == 0) {
		tchar cv = 0, ip0, ip1;
		for (id = 0; id < m_numDays; id++) {
			for (ind = 0; ind < m_numPlayers; ind++) {
				ip0 = *ivcb;
				if (ip0 != ' ' && (ip0 < '1' || ip0 > '9'))
					goto doneInit;
				ip1 = *(ivcb + 1);
				if (ip1 != ' ' && (ip1 < '0' || ip1 > '9'))
					goto doneInit;
				cv = (ip0 == ' ' ? 0 : ip0 - '0') * 10 + (ip1 - '0');
				if (cv >= m_numPlayers)
					goto doneInit;
				result(id)[ind] = cv;
				lastInd = ind;
				ivcb += 2;
				if (*ivcb == 0 || *(ivcb + 1) == 0)
					goto doneInit;
				if (!((ind + 1) % m_groupSize))
					ivcb++;
			}
		}
		goto doneInit;
	}

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
		printfRed("*** Init: values for day %d positions %d-%d not defined\n", id, lastInd + 1, m_numPlayers - 1);
		printTable("Initial result", result(0), m_numDays, m_numPlayers, m_groupSize);
		myExit(1);
	}
	linksFromMatrix(links(), result(), id);
	maxDays = iDay = id;
	if (printStartValues) {
		printTable("Start Matrix from data.h", result(), iDay, numPlayers(), m_groupSize);
		//printTableColor("Start Links", links(0), numPlayers(), numPlayers(), m_groupSize);
	}
	return true;
}

