#include "TripleSys.h"
#include <iostream>

bool CChecklLink::checkLinksH2(const char* lnks, const char* v, int nv, int nvo, int ind1, int ind2, char* vo)
{
	// v - array with players
	// nv - number of pairs in v
	// nvo - number of pairs for result
	// ind1 index of the second player in array nv (selected in the prev iteration), or unset
	// ind2 equival to ind1, or (if ind1==unset) value of player in the second pos in the prev day
	if (nvo <= 0)
	{
		return true;
	}
	char t[32]; // more time with malloc

	if (ind1 == unset)
	{
		memcpy(t, v, nv);
	}
	else
	{
		memcpy(t, v, ind1);
		if (nv >= ind1)
			memcpy(t + ind1, v + ind1 + 1, nv - ind1 + 1);
	}

	char t0 = t[0];
	const char* ct0 = lnks + t0 * m_numPlayers;
	char is = 0;
	char ie = m_numPlayers;
	if (ind1 == unset && ind2 != unset)
	{
		is = ind2;
		ie = ind2 + 2;
	}
	for (int i = 1; i < nv; i++)
	{
		char ti = t[i];
		if (ti <= is)
			continue;
		if (ti >= ie)
			break;
		if (ct0[ti] == unset)
		{
			if (checkLinksH2(lnks, t + 1, nv - 2, nvo - 2, i - 1, i - 1, vo + 2))
			{
				*vo = t0;
				*(vo + 1) = ti;
				return true;
			}
		}
	}
	return false;
}
bool CChecklLink::checkLinksH(const char* lnks, const char* v, int nv, int nvo, int ind1, int ind2, char* vo)
{
	if (nvo <= 0)
		return true;
	if (m_groupSize == 2)
		return checkLinksH2(lnks, v, nv, nvo, ind1, ind2, vo);
	char t[32];

	if (ind1 == unset)
	{
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
	const char* ct0 = lnks + t0 * m_numPlayers;
	char is = 0;
	char ie = m_numPlayers;
	if (ind1 == unset && ind2 != unset)
	{
		is = ind2;
		ie = m_numPlayers - (nv - 3) / 3 + 1;
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
			const char* cti = lnks + ti * m_numPlayers;
			for (int j = i + 1; j < nv; j++)
			{
				char tj = t[j];
				if (ct0[tj] == unset && cti[tj] == unset)
				{
					if (checkLinksH(lnks, t + 1, nv - 3, nvo - 3, i - 1, j - 1, vo + 3))
					{
						*vo = t0;
						*(vo + 1) = ti;
						*(vo + 2) = tj;
						return true;
					}
				}
			}
		}
	}
	return false;
}