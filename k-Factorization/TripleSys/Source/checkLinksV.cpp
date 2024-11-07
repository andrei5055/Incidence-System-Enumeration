#include "TripleSys.h"
CC bool CChecklLink::checkLinksV2(ctchar* lnks, int nr) const
{
	return true;
	int nrr = m_numPlayers - nr - 1;
	if (nrr <= 0 || nr < 13 || m_groupSize != 2)
		return true;
	bool ok = true;
	tchar p[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	tchar sel[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	int nOdd = 0, nEven = 0;
	memset(sel, unset, nrr * m_numPlayers);
	memset(p, unset, nrr * m_numPlayers);
	// check p1fTable cells in column i and "corresponding" cells (p1fTable(k)[p1fTable(k)[i]] = i)
	//tchar* pi = p1fTable(nr + 1);
	//tchar* pk = p1fTable(nr + 1);
	// set columns from 1 to (m_numPlayers-1) and corresponding cells
	for (int i = 0; i < m_numPlayers; i++)
	{
		tchar* pp = p;
		//tchar c[MAX_PLAYER_NUMBER];
		int nc = 0;
		// get unset values for column i;
		ctchar* lnk = lnks + i * m_numPlayers;
		for (int k = m_numPlayers - 1; k >= 0; k--)
		{
			if (lnk[k] == unset && k != i)// && k > i)
			{
				if ((i - k) & 1)
					nOdd++;
				else
					nEven++;
#if 0
			}
			continue;
			{
#endif
				int m = 0;
				nc = (nrr - 1) * m_numPlayers;
				while ((pp[nc + i] != unset && pp[nc + i] != k) && m < nrr)
				//while (((pp[nc + i] != unset && pp[nc + i] != k) || sel[nc + k] != unset) && m < nrr)
				{
					m++;
					nc -= m_numPlayers;
				}
				if (m >= nrr)
				{
					printf("i=%d k=%d\n", i, k);
					printTable("neighbors", neighbors(), nr, 16, 2);
					printTable("rest", p, nrr, 16, 2);
					ok = false;
				}
				//ASSERT(m >= nrr);
				pp[nc + i] = k;
				//if (pp[nc + k] == unset)
				//	pp[nc + k] = i;
				sel[nc + k] = i;
			}
		}
	}
	printTable("neighbors", neighbors(), nr, 16, 2);
	printTable("rest", p, nrr, 16, 2);
	printf("nOdd=%d nEven=%d nrr=%d\n", nOdd/2, nEven/2, nrr);
	printTable("links", lnks, 16, 16, 2);
	return ((nOdd/2) & 1) == 0;
	if (!ok)
		return false;
	for (int i = 0; i < nrr; i++)
	{
		tchar* pp = p + i * m_numPlayers;
		for (int k = 0; k < m_numPlayers; k++)
		{
			tchar pk = pp[k];
			ASSERT(pk == unset);
			int j = 0;
			for (; j < nrr; j++)
			{
				if (p[pk + j * m_numPlayers] == k)
					break;
			}
			if (j == nrr)
			{
				//printTable("neighbors", neighbors(), nr, 16, 2);
				//printTable("rest", p, nrr, 16, 2);
				return false;
			}
		}
	}
	printTable("neighbors", neighbors(), nr, 16, 2);
	printTable("rest", p, nrr, 16, 2);
	return true;
}
CC bool CChecklLink::checkLinksV(ctchar* c, ctchar* v, int nv, int ind, tchar* vo) const
{
	if (nv <= 0)
		return true;
	if (nv == 2 && ind != -1)
	{
		ASSERT(vo + 1 - m_vo >= m_numPlayers,
		    printf("vo + 1 - m_vo=%d\n", (int)(vo + 1 - m_vo));
		)
		switch (ind) {
		case 0:
			if (c[v[1] * m_numPlayers + v[2]] != unset)
				return false;
			*vo = v[1];
			*(vo + 1) = v[2];
			return true;
		case 1:
			if (c[v[0] * m_numPlayers + v[2]] != unset)
				return false;
			*vo = v[0];
			*(vo + 1) = v[2];
			return true;
		case 2:
			if (c[v[0] * m_numPlayers + v[1]] != unset)
				return false;
			*vo = v[0];
			*(vo + 1) = v[1];
			return true;
		}
		return false;
	}
	tchar t[MAX_PLAYER_NUMBER];
	if (ind == -1)
		memcpy(t, v, nv);
	else if (ind == 0)
		memcpy(t, v + 1, nv);
	else
	{
		memcpy(t, v, ind);
		if (nv > ind)
			memcpy(t + ind, v + ind + 1, nv - ind);
	}
	const auto* ct0 = c + t[0] * m_numPlayers;
	for (int i = 1; i < nv; i++)
	{
		if (ct0[t[i]] == unset && checkLinksV(c, t + 1, nv - 2, i - 1, vo + 2))
		{
			ASSERT(vo + 1 - m_vo >= m_numPlayers);
			*vo = t[0];
			*(vo + 1) = t[i];
			return true;
		}
	}

	return false;
}
