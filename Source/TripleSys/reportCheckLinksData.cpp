#include "TripleSys.h"
#include <iostream>
#include <chrono>
void CChecklLink::reportCheckLinksData()
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
			printTableColor("CheckLinks: 1-Fault, 2-OK, 3-Mix",
				faults, m_numDays, m_numPlayers, 0, 0, false);
			//printTableColor("CheckLinks Links",
			//	m_pLinksCopy, m_numPlayers, m_numPlayers, 0, 0, false);
		}

#if PrintNVminmax
		printTable("nv min", nvmn, m_numDays, m_numPlayers);
		printTable("nv max", nvmx, m_numDays, m_numPlayers);
		icnt = 0;
#endif
	}
	/**/
	cnt = 0;
	cntOk = 0;
	cntErr = 0;
	tmtotal = 0;
	tmtotalFalse = 0;
	tmtotalOk = 0;
	const auto len = m_numDays * m_numPlayers;
	const auto dlen = len * sizeof(double);
	memset(faults, 0, len);
	memset(counts, 0, dlen);
	memset(tmfalse, 0, dlen);
	memset(tmfalse, 0, dlen);
	/**/
}