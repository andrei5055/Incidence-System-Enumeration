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
	if (gs != 1)
	{
		for (int i = 0; i < nt * nv; i++)
		{
			t[i] *= gs;
		}
	}
	return true;
}
void alldata::cnvInit()
{
	//cnvInit1(m_allTr, m_nallTr, m_nGroups, m_nGroups, m_groupSize, true);
	//cnvInit1(m_allTg, m_nallTg, m_nGroups, m_groupSizeFactorial, 1, false);
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
int alldata::cnvCheckKm1(char* tr, int nrows, unsigned char* pOrbits) const
{
	int ret = 1;
	char ttr[27];
	const char* res = result();
	const char* resSecondRow = res + m_numPlayers;
	const int npmMinus1Row = (nrows - 1) * m_numPlayers;
	for (int day = 0; day < m_NumDaysToTransform; day++)
	{
		const auto n = m_DayIdx[day];
		const char* resn = result(n);
		for (int i = 0; i < m_numPlayers; i++)
		{
			ttr[resn[i]] = tr[i];
		}
#if 0
		if (n == 1)
		{
			printTable("Tr source", tr, 1, m_numPlayers);
			printTable("Tr actual", ttr, 1, m_numPlayers);
			printTable("Original", res, nrows, m_numPlayers);
		}
#endif
		if (m_groupSize == 2 && nrows == m_numDays)
		{
			kmTranslate(m_Km2, res, ttr, nrows, m_numPlayers);
			kmFullSort2(m_Km, m_Km2, nrows, m_numPlayers);
#if 0
			if (n == 1)
			{
				printf("n=%d ", n);
				printTable("Tr actual", ttr, 1, m_numPlayers);
				printTable("Translated", m_Km, nrows, m_numPlayers);
		}
#endif
		}
		else
		{
			kmTranslate(m_Km, res, ttr, nrows, m_numPlayers);
			kmFullSort(m_Km, nrows, m_numPlayers, m_groupSize);
		}
		const int icmp = memcmp(m_KmSecondRow, resSecondRow, npmMinus1Row);
#if 0
		static int a, b;
		a++;
		if (icmp == 0)
			b++;
		if ((a % 1000000) == 0)
			printf("%d %d\n", a, b);
#endif
#if CHECK_WITH_GROUP
		if (icmp > 0)
			continue;

		if (!icmp) {
#if 0
			if (pOrbits) {
				if (!ret) {
					;
				}
				else
					memcpy(pOrbits, ttr, m_numPlayers);
			}
#endif
			ret = 0;
			if (day) 
				m_DayIdx[day--] = m_DayIdx[--m_NumDaysToTransform];

            continue;
		}
#else
		if (icmp >= 0)
			continue;
#endif
#if PRINT_TRANSFORMED
		//extern bool flg;
		if (/*flg && */icmp < 0)
			printTransformed(nrows, m_numPlayers, tr, ttr, res, m_Km, n, nLoops, m_finalKMindex);
#endif
		//printf(" d%d", n);
		return -1;
	}
	return ret;
}
bool alldata::cnvCheckKm(char* tr, char* tg)
{
	int ii = 0;
	for (int i = 0; i < m_nGroups; i++, ii+=m_groupSize)
	{
		int itr = tr[i];
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
	bool ret = cnvCheckKm1(m_trmk, iDay) >= 0;
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
	bool ret = true;
	for (int i = 0; i < ntg; i++, ttg += m_nGroups)
	{
		if (!cnvCheckKm(tr, ttg))
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
			ret = false;
			break;
		}
	}
	//printf(" BestTg=%d\n", m_bestTg);
	return ret;
}
bool alldata::cnvCheckTgNew(char* tr, int gsf)
{
	char tg[16];
	char* ttg = tg;
	bool ret = true;
	register char ng = m_nGroups;
	register char gs = gsf;
	register int i, cnt = 0;
	memset(tg, 0, ng);
	while(1)
	{
		if (!cnvCheckKm(tr, tg))
		{
			//printf(" g%d", i);
			ret = false;
			break;
		}
		i = ng;
		while (--i >= 0)
		{
			tg[i] += 1;
			if (tg[i] < gs)
			{
				break;
			}
			tg[i] = 0;
		}
		cnt++;
		if (i < 0)
			break;
	}
	//if (ret == true)printf(" BestTg=%d\n", cnt);
	return ret;
}
bool alldata::cnvCheck()
{
	char t[16];
	int nm = m_nallTr > 1000 ? 1000 : m_nallTr - 1;
	char* ttr = m_allTr;
	bool ret = true;
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
			ret = false;
			break;
		}
	}
	//printf(" BestTr=%d\n", m_bestTr);
	return ret;
}
bool alldata::cnvCheckNew()
{
	// Head Permutations Using a Linear Array Without Recursion by Phillip Paul Fuchs
	char a[16], p[16];
	register char i, j, tmp; // Upper Index i; Lower Index j
	int itr = 0;
	bool ret = true;
#if USE_EQUAL
	initDayIdx(iDay);
#endif
	if (m_nGroups + 1 > sizeof(p))
		abort();
	for (i = 0; i < m_nGroups; i++)   // initialize arrays; a[N] can be any type
	{
		a[i] = i * m_groupSize;   // a[i] value is not revealed and can be arbitrary
		p[i] = i;
	}
	if (!cnvCheckTgNew(a, m_groupSizeFactorial))
	{
		//printf(" t%d ", itr);
		return false;
	}
	itr++;
	p[m_nGroups] = m_nGroups; // p[N] > 0 controls iteration and the index boundary for i
	i = 1;   // setup first swap points to be 1 and 0 respectively (i & j)
	while (i < m_nGroups)
	{
		p[i]--;             // decrease index "weight" for i by one
		j = i % 2 * p[i];   // IF i is odd then j = p[i] otherwise j = 0
		tmp = a[j];         // swap(a[j], a[i])
		a[j] = a[i];
		a[i] = tmp;
		if (!cnvCheckTgNew(a, m_groupSizeFactorial))
		{
			//printf(" t%d ", itr);
			ret = false;
			break;
		}
		itr++;
		i = 1;              // reset index i to 1 (assumed)
		while (!p[i])       // while (p[i] == 0)
		{
			p[i] = i;        // reset p[i] zero value
			i++;             // set new index value for i (increase by one)
		} // while(!p[i])
	} // while(i < m_nGroups)
	//if (ret == true)printf(" itr=%d\n", itr);
	return ret;
}
