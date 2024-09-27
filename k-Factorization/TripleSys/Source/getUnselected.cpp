#include "TripleSys.h"

CC void alldata::getUnselected(tchar* v, int nv) const
{
	int j = 0;
	for (int i = 0; i < m_numPlayers; i++)
	{
		if (selPlayers[i] == unset)
		{
			v[j] = i;
			if (++j == nv)
				return;
		}
	}
	ASSERT(1);
}

