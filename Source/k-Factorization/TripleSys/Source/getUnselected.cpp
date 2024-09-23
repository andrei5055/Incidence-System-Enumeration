#include <iostream>
#include "TripleSys.h"

void alldata::getUnselected(char* v, int nv)
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
	abort();
}

