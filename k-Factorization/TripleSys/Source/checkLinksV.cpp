#include "TripleSys.h"

CC bool CChecklLink::checkLinksV(const tchar* c, const tchar* v, int nv, int ind, tchar* vo) const
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
