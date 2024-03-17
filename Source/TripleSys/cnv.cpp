#include <iostream>
#include "TripleSys.h"
int factorial(int n) {
	return n > 2 ? n * factorial(n - 1) : n;
}
bool checkSet(char* tr, int nt)
{
	int n = 1<<tr[0];
	for (int i = 1; i < nt; i++)
	{
		int m = 1<<tr[i];
		if ((n & m) != 0)
			return false;
		n |= m;
	}
	return true;
}

bool cnvInit1(char* t, int nt, int nv, int mv, int gs, bool bCheck)
{
	char* tt = t + nv; 
	for (int i = 0; i < nv; i++)
	{
		t[i] = bCheck ? i : 0;
	}
	for (int i = 1; i < nt; i++, tt+=nv)
	{
		memcpy(tt, tt - nv, nv);
		while (1)
		{
			for (int j = 1; j <= nv; j++)
			{
				tt[nv - j]++;
				if (tt[nv - j] < mv)
					break;
				tt[nv - j] = 0;
			}
			if (!bCheck || checkSet(tt, nv))
				break;
		}
	}
	for (int i = 0; i < nt * nv; i++)
	{
		t[i] *= gs;
	}
	return true;
}
void alldata::cnvInit()
{
	cnvInit1(m_allTr, m_nallTr, m_nGroups, m_nGroups, m_groupSize, true);
	cnvInit1(m_allTg, m_nallTg, m_nGroups, m_groupSizeFactorial, m_groupSize, false);
	for (int i = 0; i < m_groupSize; i++)
	{
		m_groups[i] = i;
	}
	for (int i = 1; i < m_groupSizeFactorial; i++)
	{
		int ig = i * m_groupSize;
		memcpy(m_groups + ig, m_groups + ig - m_groupSize, m_groupSize);
		while (1)
		{
			for (int j = m_groupSize - 1; j >= 0; j--)
			{
				if (m_groups[ig + j] < m_groupSize - 1)
				{
					m_groups[ig + j]++;
					break;
				}
				m_groups[ig + j] = 0;
			}
			if (checkSet(m_groups + ig, m_groupSize))
				break;
		}
	}
}
bool alldata::cnvCheckKm1(char* tr)
{
	bool ret = true;
	char ttr[27];
	int npm = iDay * m_numPlayers;
	char* mo = m_Km;
	char* res = result();
	char* mo2 = mo + m_numPlayers;
	char* mo3 = mo + npm;
	char* res2 = res + m_numPlayers;
	int npm2 = npm - m_numPlayers;
	for (int n = 0; n < iDay; n++)
	{
		char* resn = result(n);
		for (int i = 0; i < m_numPlayers; i++)
		{
			ttr[resn[i]] = tr[i];
		}
#if 0
		if (n == 1)
		{
			printTable("Tr source", tr, 1, m_numPlayers);
			printTable("Tr actual", ttr, 1, m_numPlayers);
			printTable("Original", res, iDay, m_numPlayers);
		}
#endif
		if (m_groupSize == 2)
		{
			kmTranslate(mo3, res, ttr, iDay, m_numPlayers);
			kmFullSort2(mo, mo3, iDay, m_numPlayers);
		}
		else
		{
			kmTranslate(mo, res, ttr, iDay, m_numPlayers);
			kmFullSort(mo, iDay, m_numPlayers, m_groupSize);
		}
#if 0
		//if (n == 1)
		{
			static int icnt = 0,icnt2 = 0;
			static char ttrd[11][12] ;
			if (icnt == 0)
				memset(ttrd[0], -1, sizeof(ttrd));
			if (ttrd[n][0] == unset)
				memcpy(ttrd[n], ttr, m_numPlayers);
			else
			{
				char ttrn[24];
				int j = 0;
				for (int i = 0;i < m_numPlayers;i++)
				{
					if (ttrd[n][i] != ttr[i])
					{
						if (j < 3)
						{
							ttrn[j] = ttrd[n][i];
							ttrn[j+1] = ttr[i];
						}
						j += 2;
					}
				}
				if (j < 25)
				{
					icnt2++;
					if ((ttrn[0] != ttrn[3]) != 0 || ttrn[1] != ttrn[2])
						j = j;
					if (ttrn[0] != ttrn[1] + 1 && ttrn[0] != ttrn[1] - 1)
						j = j;
					if ((min(ttrn[0], ttrn[1]) & 1) != 0)
						j = j;
				}
			}
			memcpy(ttrd[n], ttr, m_numPlayers);
			icnt++;
			if ((icnt % 100000) == 0)
				printf("icnt=%d, %d\n", icnt, icnt2);
			//printTable("Tr actual", ttr, 1, m_numPlayers);
			//printTable("Translated", mo, iDay, m_numPlayers);
		}
#endif
		if (memcmp(mo2, res2, npm2) < 0)
		{
#if 0
			printf("Calculated Matrix %.0f can be improved. See below (nKm=%d row=%d)\n", nLoops, m_finalKMindex, n);
			printTable("Tr source", tr, 1, m_numPlayers);
			printTable("Tr actual", ttr, 1, m_numPlayers);
			printTable("Original", res, iDay, m_numPlayers);
			printTable("Translated", mo, iDay, m_numPlayers);
#endif
			//printf(" d%d", n);
			ret = false;
			break;
		}	
	}
	return ret;
}
bool alldata::cnvCheckKm(char* tr, char* tg, int gfs)
{
	int ii = 0;
	for (int i = 0; i < m_nGroups; i++, ii+=m_groupSize)
	{
		int itr = tr[i];
		int itg = tg[i];
		for (int j = 0; j < m_groupSize; j++)
		{
			m_trmk[ii + j] = itr + m_groups[itg + j];
		} 
	}
#if 0
	printTable("Tr", tr, 1, m_nGroups);
	printTable("Tg", tg, 1, m_nGroups);
#endif
	bool ret = cnvCheckKm1(m_trmk);
#if 0
	if (!ret)
	{
		printTable("Tr", tr, 1, m_nGroups);
		printTable("Tg", tg, 1, m_nGroups);
	}
#endif
	return ret;
}
bool alldata::cnvCheckTg(char* tr, char* tg, int ntg, int gsf)
{
	int nm = ntg > 1000 ? 1000 : ntg - 1;
	char* ttg = tg;
	for (int i = 0; i < ntg; i++, ttg += m_nGroups)
	{
		if (!cnvCheckKm(tr, ttg, gsf))
		{
#if 1  // a little bit faster with 1
			char t[16];
			if (i > m_bestTg && m_bestTg < nm)
			{
				memcpy(t, ttg, m_nGroups);
				memcpy(ttg, tg + m_bestTg * m_nGroups, m_nGroups);
				memcpy(tg + m_bestTg * m_nGroups, t, m_nGroups);
				m_bestTg++;
				//printf("i=%d bestTg=%d\n", i, m_bestTg);
			}
#endif
			//printf(" g%d", i);
			return false;
		}
	}
	//printf(" BestTg=%d\n", m_bestTg);
	return true;
}
bool alldata::cnvCheck()
{
	char t[16];
	int nm = m_nallTr > 1000 ? 1000 : m_nallTr - 1;
	char* ttr = m_allTr;
	if (m_nGroups > sizeof(t))
		abort();
	for (int i = 0; i < m_nallTr; i++, ttr += m_nGroups)
	{
		if (!cnvCheckTg(ttr, m_allTg, m_nallTg, m_groupSizeFactorial))
		{
#if 1 // a little bit faster with 1
			if (i > m_bestTr && m_bestTr < nm)
			{
				memcpy(t, ttr , m_nGroups);
				memcpy(ttr, m_allTr + m_bestTr * m_nGroups, m_nGroups);
				memcpy(m_allTr + m_bestTr * m_nGroups, t, m_nGroups);
				m_bestTr++;
				//printf("i=%d best=%d\n", i, m_bestTr);
			}
#endif
			//printf(" t%d ", i);
			return false;
		}
	}
	//printf(" BestTr=%d\n", m_bestTr);
	return true;
}