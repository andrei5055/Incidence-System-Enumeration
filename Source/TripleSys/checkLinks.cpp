#include "TripleSys.h"
#include <iostream>
#include <chrono>

void alldata::setCheckLinks() {

}

CChecklLink::CChecklLink(int numDays, int numPlayers) :
		SizeParam(numDays, numPlayers) {
	const auto len = numDays * numPlayers;
	initArray(&counts, len);
	initArray(&faults, len);
	initArray(&tmfalse, len);
	initArray(&tmok, len);
	initArray(&cb, numPlayers * numPlayers);
#if PrintNVminmax
	initArray(&nvmn, len, (char)99);
	initArray(&nvmx, len, (char)(-1));
#endif
}

CChecklLink::~CChecklLink() {
	delete[] counts;
	delete[] faults;
	delete[] tmfalse;
	delete[] tmok;
	delete[] cb;
#if PrintNVminmax
	delete[] nvmn;
	delete[] nvmx;
#endif
}

#if PrintNVminmax
void CChecklLink::setNV_MinMax(int id, int idx, char nv) {
	idx += id * m_numPlayers;
	if (nvmn[idx] > nv)
		nvmn[idx] = nv;
	if (nvmx[idx] < nv)
		nvmx[idx] = nv;
}
#endif


bool CChecklLink::checkLinks(char *c, int id, bool printLinksStatTime)
{
#define cntm 10000000
	std::chrono::steady_clock::time_point start0, start;
	bool ret = true;

//	if (m_numPlayers == 27)
//		return checkLinks27(c, id);
	const auto len = m_numPlayers * m_numPlayers;
	memcpy(cb, c, len);

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
					m_numDays, m_numPlayers, 0, 0, false, 100.0 / cnt);
				if (printLinksStatTime)
				{
					printf("\nTotal time for 'checkLinks'=%.0fsec, 'False' returns time=%.2f%%, 'OK' returns time=%.2f%% \n",
						tmtotal / 1000000, tmtotalFalse * 100.0 / tmtotal, tmtotalOk * 100.0 / tmtotal);
					printTable("'False checkLinks' times per player per day (% of total time)",
						tmfalse[0], m_numDays, m_numPlayers, 0, 0, false, 100.0 / tmtotal);
					printTable("'OK' checkLinks' times per player per day (% of total time)",
						tmok[0], m_numDays, m_numPlayers, 0, 0, false, 100.0 / tmtotal);
				}
				**/
				printTableColor("'LinksCheck': 1-Fault, 2-OK, 3-Mix",
					faults, m_numDays, m_numPlayers, 0, 0, false);
			}

#if PrintNVminmax
			printTable("nv min", nvmn, m_numDays, m_numPlayers);
			printTable("nv max", nvmx, m_numDays, m_numPlayers);
			icnt = 0;
#endif
		}
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

	cnt++;
	icnt++;

	if (printLinksStatTime)
		start0 = std::chrono::high_resolution_clock::now();

	char vBuffer[100];
	char *v = m_numPlayers <= sizeof(vBuffer) / 2 ? vBuffer : new char[2 * m_numPlayers];
	char *vo = v + m_numPlayers;
	const auto idx = id * m_numPlayers;
	auto *faults_id = faults + idx;
	auto* counts_id = counts + idx;
	for (int i0 = 0; i0 < m_numPlayers; i0++)
	{
		//i = i0;
		int i = (i0 + 5) % m_numPlayers;
		auto *ci = c + i * m_numPlayers;
		if (printLinksStatTime)
			start = std::chrono::high_resolution_clock::now();

		char* lnk = c + i * m_numPlayers;
		int nv = 0, ns = 0;

		for (int j = 0; j < m_numPlayers; j++)
		{
			if (lnk[j] == unset && i != j)
				v[nv++] = j;
		}
		if (nv == 0)
			continue;

		if (nv >= (m_numDays - 1) * 2) // this check make it faster
			continue;

		//if (i0 < 3 && nv / 2 < m_numDays - id - 1)
		//	goto fltPlayer;
		///continue;
		if ((nv % 2) != 0)
		{
			//continue;
			abort();
		}
		//memset(vo, unset, sizeof(vo));
		if (nv < 0)
			goto okplayer;

		if (checkLinksV(c, v, nv, -1, vo))
		{
			int idd = 0;
			for (int n = 0; n < nv; n += 2)
			{
				const auto a = vo[n];
				const auto b = vo[n + 1];

				if (a == unset || b == unset)
					abort();
				auto* ca = c + a * m_numPlayers;
				if (ci[a] != unset || ci[b] != unset || ca[b] != unset)
					abort();
				idd = -2;//m_numDays - nv / 2 + n;

				auto* cb = c + b * m_numPlayers;
				ci[a] = ca[i] = idd;
				ci[b] = cb[i] = idd;
				/**/
				cb[a] = ca[b] = idd;
				/**/
			}
			goto okplayer;
		}
	//fltPlayer:
		setNV_MinMax(id, i, nv);
		faults_id[i] |= 1;
		counts_id[i]++;
		cntErr++;
		if (printLinksStatTime)
		{
			auto elapsed = std::chrono::high_resolution_clock::now() - start;
			long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
			*(tmfalse + idx + i) += (double)microseconds;
			tmtotalFalse += (double)microseconds;
		}
		ret = false;
		break;
	okplayer:
		setNV_MinMax(id, i, nv);
		faults_id[i] |= 2;

		if (printLinksStatTime)
		{
			auto elapsed = std::chrono::high_resolution_clock::now() - start;
			long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
			*(tmok + idx + i) += (double)microseconds;
			tmtotalOk += (double)microseconds;
		}
	}
	if (ret)
	{
		cntOk++;
		/**/
		if (id > 33)
		{
			printTable("Links", c, m_numPlayers, m_numPlayers);
			convertLinksToResult(c);
			printTable("result", m_co, m_numDays, m_numPlayers);
		}
		/**/
	}

	memcpy(c, cb, len);

	if (printLinksStatTime)
	{
		auto elapsed0 = std::chrono::high_resolution_clock::now() - start0;
		long long microseconds0 = std::chrono::duration_cast<std::chrono::microseconds>(elapsed0).count();
		tmtotal += (double)microseconds0;
	}

	if (v != vBuffer)
		delete[] v;

	return ret;
}