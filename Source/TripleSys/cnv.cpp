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

bool cnvInit1(char* t, int nt, int nv, int mv, bool bCheck)
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
	return true;
}
void alldata::cnvInit()
{
	cnvInit1(m_allTr, m_nallTr, m_nGroups, m_nGroups, true);
	cnvInit1(m_allTg, m_nallTg, m_nGroups, m_groupSizeFactorial, false);
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
		kmTranslate(mo3, res, ttr, iDay, m_numPlayers);
		kmFullSort(mo, mo3, iDay, m_numPlayers, m_groupSize);
#if 0
		if (n == 1)
		{
			printTable("Tr actual", ttr, 1, m_numPlayers);
			printTable("Translated", mo, iDay, m_numPlayers);
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
	for (int i = 0; i < m_nGroups; i++)
	{
		int ii = i * m_groupSize;
		int itr = tr[i] * m_groupSize;
		int itg = tg[i] * m_groupSize;
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
	char* ttg = tg;
	for (int i = 0; i < ntg; i++, ttg += m_nGroups)
	{
		if (!cnvCheckKm(tr, ttg, gsf))
		{
			//printf(" g%d", i);
			return false;
		}
	}
	return true;
}
bool alldata::cnvCheck()
{
	char* ttr = m_allTr;
	for (int i = 0; i < m_nallTr; i++, ttr += m_nGroups)
	{
		if (!cnvCheckTg(ttr, m_allTg, m_nallTg, m_groupSizeFactorial))
		{
			//printf(" t%d ", i);
			return false;
		}
	}
	return true;
}