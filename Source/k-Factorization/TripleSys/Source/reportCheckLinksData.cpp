#include "TripleSys.h"
#include <iostream>
void CChecklLink::reportCheckLinksData()
{
	if (m_cnt != 0)
	{
		if (PrintLinksStat)
		{
			printf("\nTotal calls to 'checkLinks'=%.0fM, 'False' returns=%.2f%%, 'OK' returns=%.2f%%\n",
				m_cnt / 1000000.0, m_cntErr * 100.0 / m_cnt, m_cntOk * 100.0 / m_cnt);
			/**
			printTable("'False' returns per player per day (% of total calls)", m_counts[0],
				m_numDays, m_numPlayers, 0, 0, false, 100.0 / m_cnt);
			if (printLinksStatTime)
			{
				printf("\nTotal time for 'checkLinks'=%.0fsec, 'False' returns time=%.2f%%, 'OK' returns time=%.2f%% \n",
					m_tmtotal / 1000000, m_tmtotalFalse * 100.0 / m_tmtotal, m_tmtotalOk * 100.0 / m_tmtotal);
				printTable("'False checkLinks' times per player per day (% of total time)",
					m_tmfalse[0], m_numDays, m_numPlayers, 0, 0, false, 100.0 / m_tmtotal);
				printTable("'OK' checkLinks' times per player per day (% of total time)",
					m_tmok[0], m_numDays, m_numPlayers, 0, 0, false, 100.0 / tmtotal);
			}
			**/
			printTableColor("CheckLinks: 1-Fault, 2-OK, 3-Mix",
				m_faults, m_numDays, m_numPlayers, 0, 0, false);
			//printTableColor("CheckLinks Links",
			//	m_pLinksCopy, m_numPlayers, m_numPlayers, 0, 0, false);
		}

#if PrintNVminmax
		printTable("nv min", m_nvmn, m_numDays, m_numPlayers);
		printTable("nv max", m_nvmx, m_numDays, m_numPlayers);
		icnt = 0;
#endif
	}
	/**/
	m_cnt = 0;
	m_cntOk = 0;
	m_cntErr = 0;
	m_tmtotal = 0;
	m_tmtotalFalse = 0;
	m_tmtotalOk = 0;
	const auto len = m_numDays * m_numPlayers;
	const auto dlen = len * sizeof(sLongLong);
	memset(m_faults, 0, len);
	memset(m_counts, 0, dlen);
	memset(m_tmfalse, 0, dlen);
	memset(m_tmfalse, 0, dlen);
	/**/
}