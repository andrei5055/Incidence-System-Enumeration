#include "TripleSys.h"
#include <iostream>

bool CChecklLink::checkLinksH(const char* c, char *st, char *s, int ipos, const char* v, int nv, int nvo, int ind1, int ind2, char* vo, double* counter)
{
	if (nvo <= 0)
	{/**
		if (counter != NULL)
		{
			if ((int)(*counter / 10000000.0) * 10000000.0 == *counter)
				printf("c=%.0fM\n", *counter / 1000000.);
		}**/
		return true;
	}
	char t[100]; // more time with malloc
	char sel[100];
	//char* t = (char *)malloc(nv);

	if (s != NULL)
	{
		if (st == NULL)
			memcpy(sel, s, m_numPlayers);
		else
			memcpy(sel, st, m_numPlayers);
	}

	if (ind1 == unset)
	{
		if (counter != NULL)
			*counter = 0;
		memcpy(t, v, nv);
	}
	else 
	{
		memcpy(t, v, ind1);
		memcpy(t + ind1, v + ind1 + 1, ind2 - ind1 - 1);
		if (nv >= ind2)
		    memcpy(t + ind2 - 1, v + ind2 + 1, nv - ind2 + 1);
	}

	char t0 = t[0];
	const char* ct0 = c + t0 * m_numPlayers;
	char is = 0;
	char ie = m_numPlayers;
	if (ind1 == unset && ind1 != unset)
	{
		is = ind2;
		ie = m_numPlayers - (nv - 3) / 3;
	}
	for (int i = 1; i < nv - 1; i++)
	{
		//if (t[0] == 0 && i >= 7  && i <= 10 && m_numPlayers == 21)
		//	i = 11;
		char ti = t[i];
		if (ti <= is)
			continue;
		if (ti >= ie)
			break;
		if (ct0[ti] == unset)
		{
			const char* cti = c + ti * m_numPlayers;
			for (int j = i + 1; j < nv; j++)
			{
				char tj = t[j];
				if ((s == NULL || (tj % 3) != 2 || (sel[tj - 2] != unset && sel[tj - 1] != unset && sel[tj - 2] < sel[tj - 1])) &&
					ct0[tj] == unset &&
					cti[tj] == unset)
				{
					if (sel != NULL)
					{
						sel[t0] = ipos;
						sel[ti] = ipos + 1;
						sel[tj] = ipos + 2;
					}
					if (checkLinksH(c, sel, s, ipos + 3, t + 1, nv - 3, nvo - 3, i - 1, j - 1, vo + 3, counter))
					{
						*vo = t0;
						*(vo + 1) = ti;
						*(vo + 2) = tj;
						if (counter != NULL)
							(*counter)++;
						else
						{
							if (s != 0)
							{
								s[t0] = ipos;
								s[ti] = ipos + 1;
								s[tj] = ipos + 2;
							}
							//free(t);
							return true;
						}
					}
					if (sel != NULL)
					{
						sel[t0] = unset;
						sel[ti] = unset;
						sel[tj] = unset;
					}
				}
			}
		}
	}
	//free(t);
	return false;
}