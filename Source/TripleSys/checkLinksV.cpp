#include "TripleSys.h"
#include <iostream>

bool CChecklLink::checkLinksV(const char *c, const char *v, int nv, int ind, char *vo)
{
	char t[nPlayers];
	if (nv <= 0)
		return true;
	if (ind == unset)  
		memcpy(t, v, nv);
	else if (nv == 2)
	{
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
	else if (ind == 0)
		memcpy(t, v + 1, nv);
	else
	{
		memcpy(t, v, ind);
		if (nv > ind)
		    memcpy(t + ind, v + ind + 1, nv - ind);
	}
	nv -= 2;
	const char *ct0 = c + t[0] * m_numPlayers;
	for (int i = 0; i <= nv; i++)
	{
		if (ct0[t[i+1]] == unset && checkLinksV(c, t + 1, nv, i, vo + 2))
		{
			*vo = t[0];
			*(vo + 1) = t[i+1];
			return true;
		}
	}
	return false;
}