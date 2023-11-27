#include "TripleSys.h"
#include <iostream>
#include <chrono>

bool oktoset1(char* c, char av)
{
	int cnt = 0;
	for (int i = 0; i < nPlayers; i++)
	{
		if (c[i] == av)
		{
			cnt++;
			if (cnt > 1)
				return false;
		}
	}
	return true;
}
bool oktoset2(char* c, char av)
{
	for (int i = 0; i < nPlayers; i++)
	{
		if (c[i] == av)
			return false;
	}
	return true;
}

bool CChecklLink::checkLinks27(char* c, int id)
{
#define cntm 10000000
	auto *faultsID = faults + id * m_numPlayers;
	std::chrono::steady_clock::time_point start0, start;
	const auto len1 = m_numPlayers * m_numPlayers;
	char *cb = new char[len1];
	bool ret = true;
	int idfin = 10;
	cnt++;
	icnt++;

	/**
	char ca[nDays][nPlayers];
	char cc[nDays][nPlayers];
	char cl[nPlayers][nPlayers];
	memset(&ca, unset, sizeof(ca));
	memset(&cc, unset, sizeof(cc));
	memset(&cl, unset, sizeof(cl));
	for (int j = 0; j < nPlayers; j++)
	{
		for (int i = 0; i < nDays; i++)
		{
			int k, k0;
			for (k0 = 0; k0 < nPlayers; k0++)
			{
				k = k0; //k = ((j % 2) == 0) ? k0 : nPlayers - k0 - 1;
				if (i == 3 && j == 4 && k == 5)
					i = i;
				if (ca[i][k] != unset)
					continue;

				if ((j % 3) != 0)
				{
					if (cl[cc[i][j - 1]][k] != unset)
						continue;
					if ((j % 3) == 1)
					{
						//if (!oktoset2(cl[j], i))
						//	continue;
					}
					if ((j % 3) == 2)
					{
						//if (!oktoset1(cl[j], i))
							//continue;
						if (cl[cc[i][j - 2]][k] != unset)
							continue;
						cl[k][cc[i][j - 2]] = i;
						cl[cc[i][j - 2]][k] = i;
					}
					cl[k][cc[i][j - 1]] = i;
					cl[cc[i][j - 1]][k] = i;
				}
				ca[i][k] = j;
				cc[i][j] = k;
				break;
			}
		}
	}
	printTable("Links", cl[0], nPlayers, nPlayers);
	char co[nDays][nPlayers];
	convertLinksToResult(cl, co);
	printTable("result from link", co[0], nDays, nPlayers);
	printTable("result", cc[0], nDays, nPlayers);
	exit(0);
	**/

	//if (id > nDays - 4)
	//	return true;
	memcpy(cb, c, len1);

	if ((icnt % cntm) == 0)
	{
		if (cnt != 0.0)
		{
			if (PrintLinksStat)
			{
				printf("\nTotal calls to 'checkLinks'=%.0fM, 'False' returns=%.2f%%, 'OK' returns=%.2f%%\n",
					cnt / 1000000, cntErr * 100.0 / cnt, cntOk * 100 / cnt);
				/**
				printTable("'False' returns per player per day (% of total calls)", counts[0],
					nDays, nPlayers, 0, 0, false, 100.0 / cnt);
				if (PrintLinksStatTime)
				{
					printf("\nTotal time for 'checkLinks'=%.0fsec, 'False' returns time=%.2f%%, 'OK' returns time=%.2f%% \n",
						tmtotal / 1000000, tmtotalFalse * 100.0 / tmtotal, tmtotalOk * 100.0 / tmtotal);
					printTable("'False checkLinks' times per player per day (% of total time)",
						tmfalse[0], nDays, nPlayers, 0, 0, false, 100.0 / tmtotal);
					printTable("'OK' checkLinks' times per player per day (% of total time)",
						tmok[0], nDays, nPlayers, 0, 0, false, 100.0 / tmtotal);
				}
				**/
				printTableColor("'LinksCheck': 1-Fault, 2-OK, 3-Mix",
					faults, m_numDays, m_numPlayers, 0, 0, false);
			}
		}
		icnt = 0;
		/**
		cnt = 0;
		cntErr = 0;
		tmtotal = 0;
		tmtotalFalse = 0;
		tmtotalOk = 0;
		memset(faults[0], 0, sizeof(faults));
		memset(counts[0], 0, sizeof(counts));
		memset(tmfalse[0], 0, sizeof(tmfalse));
		memset(tmfalse[0], 0, sizeof(tmok));
		**/
	}

	int i = 0;
	if (PrintLinksStatTime)
		start0 = std::chrono::high_resolution_clock::now();
	int ie = id < 3 ? nPlayers / 3 : nPlayers;
	ie = 9;
	for (int i0 = 0; i0 < ie; i0++)
	{
		i = (i0 + 5) % m_numPlayers;
		i = i0;
		//int i = ((i0 % 2) == 0) ? i0 / 2 : numPlayers - 1 - i0 / 2;
		if (PrintLinksStatTime)
		{
			start = std::chrono::high_resolution_clock::now();
		}

		char v[nPlayers], vo[nPlayers];
		char* lnk = c + i * m_numPlayers;
		int nv = 0, ns = 0;

		for (int j = 0; j < m_numPlayers; j++)
		{
			if (lnk[j] == unset && i != j)
			{
				v[nv] = j;
				nv++;
			}
		}
		if (nv == 0)
			continue;
			continue;
		if (nv >= (m_numDays - 1) * 2) // 6 times faster
			continue;
		if ((nv % 2) != 0)
			abort();
		/**
		if (id > 4 && nv / 2 < nDays - id - 1)
			goto fltPlayer;
			**/
		//memset(vo, unset, sizeof(vo));
		if (nv < 0)
			goto okplayer;

		if (checkLinksV(c, v, nv, -1, vo))
		{
			int idd = 0;
			for (int n = 0; n < nv; n += 2)
			{
				char a = vo[n];
				char b = vo[n + 1];

				if (a == unset || b == unset)
					abort();
				auto* ca = c + a * m_numPlayers;
				if (lnk[a] != unset || lnk[b] != unset || ca[b] != unset)
					abort();
				idd = -2;
				/**/
				auto* cb = c + b * m_numPlayers;
				if (id > 7)
				{
					int k;
					for (k = id; k < m_numDays; k++)
					{
						idd = k;
						if (oktoset2(lnk, idd) && oktoset2(ca, idd) && oktoset2(cb, idd))
							break;

					}
					if (k >= m_numDays)
						goto fltPlayer;
				}
				/**/
				 

				lnk[a] = ca[i] = idd;
				lnk[b] = cb[i] = idd;
				
				if (1)
					cb[a] = ca[b] = idd;
			}
			goto okplayer;
		}
	fltPlayer:
		faultsID[i] |= 1;
		++*(counts + id * m_numPlayers + i);
		cntErr++;
		if (PrintLinksStatTime)
		{
			auto elapsed = std::chrono::high_resolution_clock::now() - start;
			long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
			*(tmfalse + id * m_numPlayers + i) += (double)microseconds;
			tmtotalFalse += (double)microseconds;
		}
		ret = false;
		break;
	okplayer:

		faultsID[i] |= 2;

		if (PrintLinksStatTime)
		{
			auto elapsed = std::chrono::high_resolution_clock::now() - start;
			long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
			*(tmok + id * m_numPlayers + i) += (double)microseconds;
			tmtotalOk += (double)microseconds;
		}
	}
	if (ret)
	{
		cntOk++;
		/**/
		if (id > 7)
		{
			printTable("Links", c, m_numPlayers, m_numPlayers);
			convertLinksToResult(c);
			printTable("result", m_co, m_numDays, m_numPlayers);
			i = i;
			exit(0);
		}
		/**/
	}
	if (ret)
		ret = ret;
	memcpy(c, cb, len1);

	if (PrintLinksStatTime)
	{
		auto elapsed0 = std::chrono::high_resolution_clock::now() - start0;
		long long microseconds0 = std::chrono::duration_cast<std::chrono::microseconds>(elapsed0).count();
		tmtotal += (double)microseconds0;
	}

	delete[] cb;
	return ret;
}