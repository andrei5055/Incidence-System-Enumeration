#include "TripleSys.h"
CC bool CChecklLink::checkLinksV(ctchar* links, ctchar* v, int nv, int ind, tchar* vo) const
{
	if (nv <= 0) {
		//printTableColor("v", m_vo, 1, vo + 1 - m_vo, 0);
		//return false;
		return true;
	}
	const auto cbmpGraph = !completeGraph();
	if (nv == 2 && ind != -1)
	{
		ASSERT_IF(vo + 1 - m_vo >= m_numPlayers,
			printf("vo + 1 - m_vo=%d\n", (int)(vo + 1 - m_vo));
			)
			switch (ind) {
			case 0:
				if (links[v[1] * m_numPlayers + v[2]] != unset)
					return false;
				if (cbmpGraph) {
					if (m_remainder3[v[1]] == m_remainder3[v[2]])
						return false;
				}
				*vo = v[1];
				*(vo + 1) = v[2];
				break;
			case 1:
				if (links[v[0] * m_numPlayers + v[2]] != unset)
					return false;
				if (cbmpGraph) {
					if (m_remainder3[v[0]] == m_remainder3[v[2]])
						return false;
				}
				*vo = v[0];
				*(vo + 1) = v[2];
				break;
			case 2:
				if (links[v[0] * m_numPlayers + v[1]] != unset)
					return false;
				if (cbmpGraph) {
					if (m_remainder3[v[0]] == m_remainder3[v[1]])
						return false;
				}
				*vo = v[0];
				*(vo + 1) = v[1];
				break;
			default:
				ASSERT_IF(1);
				EXIT_(104);
				return false;
			}
		//printTableColor("v", m_vo, 1, vo + 2 - m_vo, 0);
		//return false;
		return true;
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
	const auto* ct0 = links + t[0] * m_numPlayers;
	if (cbmpGraph) {
		const auto t0r3 = m_remainder3[t[0]];
		for (int i = 1; i < nv; i++)
		{
			if (t0r3 == m_remainder3[t[i]])
				continue;
			if (ct0[t[i]] == unset) {
				ASSERT_IF(vo + 1 - m_vo >= m_numPlayers);
				//*vo = t[0];
				//*(vo + 1) = t[i];
				if (checkLinksV(links, t + 1, nv - 2, i - 1, vo + 2)) {
					*vo = t[0];
					*(vo + 1) = t[i];
					return true;
				}
			}
		}
	}
	else {
		for (int i = 1; i < nv; i++)
		{
			if (ct0[t[i]] == unset && checkLinksV(links, t + 1, nv - 2, i - 1, vo + 2))
			{
				ASSERT_IF(vo + 1 - m_vo >= m_numPlayers);
				*vo = t[0];
				*(vo + 1) = t[i];
				return true;
			}
		}
	}
	return false;
}
