#include "TripleSys.h"
#include <iostream>

bool CChecklLink::checkLinksTR(const char* v, int nvAll, int nv, int ind)
{
	if (nvAll <= 0)
		return true;
	if (nv <= 0)
	{
		if (ind != unset)
			v++;
		nv = *v;
		v++;
		nvAll--;
		ind = unset;
	}
	char t[256];
	if (ind == unset)
		memcpy(t, v, nvAll);
	else if (ind == 0)
		memcpy(t, v + 1, nvAll);
	else
	{
		memcpy(t, v, ind);
		if (nvAll > ind)
			memcpy(t + ind, v + ind + 1, nvAll - ind);
	}
	char* lnkt0 = m_pLinksCopy + t[0] * m_numPlayers;
	for (int i = 1; i < nv; i++)
	{
		char* lnk = lnkt0 + t[i];
		if (*lnk == unset)
		{
			*lnk = 1;
			if (checkLinksTR(t + 1, nvAll - 2, nv - 2, i - 1))
				return true;
			*lnk = unset;
		}
	}
	return false;
}
bool CChecklLink::checkLinksV(const char* c, const char* v, int nv, int ind, char* vo)
{
	if (nv <= 0)
		return true;
	if (nv == 2 && ind != unset)
	{
		if (vo + 1 - m_vo >= m_numPlayers)
		{
			printf("vo + 1 - m_vo=%d\n", (int)(vo + 1 - m_vo));
			abort();
		}
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
	char t[100]; // malloc required additional ~50% cpu
	//char* t = (char*)malloc(nv);
	if (ind == unset)
		memcpy(t, v, nv);
	else if (ind == 0)
		memcpy(t, v + 1, nv);
	else
	{
		memcpy(t, v, ind);
		if (nv > ind)
			memcpy(t + ind, v + ind + 1, nv - ind);
	}
	const char* ct0 = c + t[0] * m_numPlayers;
	for (int i = 1; i < nv; i++)
	{
		if (ct0[t[i]] == unset && checkLinksV(c, t + 1, nv - 2, i - 1, vo + 2))
		{
			if (vo + 1 - m_vo >= m_numPlayers)
				abort();
			*vo = t[0];
			*(vo + 1) = t[i];
			//free(t);
			return true;
		}
	}
	//free(t);
	return false;
}