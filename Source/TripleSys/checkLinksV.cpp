#include "TripleSys.h"
#include <iostream>

bool CChecklLink::checkLinksV(const char *c, const char *v, int nv, int ind, char *vo)
{
	char t[nPlayers];
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

	nv -= 2;
	const char *ct0 = c + t[0] * m_numPlayers;
	for (int i = 0; i <= nv; i++)
	{
		if (ct0[t[i+1]] == unset && (nv <= 0 || (nv > 1) && checkLinksV(c, t + 1, nv, i, vo + 2)))
		{
			*vo = t[0];
			*(vo + 1) = t[i+1];
			return true;
		}
	}
	return false;
}