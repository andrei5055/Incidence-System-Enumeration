#include <iostream>
#include "TripleSys.h"

void SizeParam::convertLinksToResult(const char *ci)
{
	const auto len = m_numPlayers * m_numPlayers;
	char* co = m_co;
	char *c = new char[len];
	memcpy(c, ci, len);
	char* a = m_pBuf;
	for (int i = 0; i < m_numDays; i++, co += m_numPlayers)
	{
		int m = 0;

		memset(co, unset, m_numPlayers);
		auto* cj = c;
		const auto gs_1 = m_groupSize - 1;
		for (int j = 0; j < m_numPlayers; j++, cj += m_numPlayers)
		{
			int cnt = 0;

			for (int k0 = 0; k0 < m_numPlayers; k0++)
			{
				if (cj[k0] == i)
				{
					if (cnt >= m_groupSize)
					{
						printf("too many players on day %d (%d, %d)\n", i, j, k0);
						break;
					}
					a[cnt++] = k0;
				}
			}
			if (cnt == gs_1)
			{
				if (m >= m_numPlayers - cnt)
				{
					printf("error on day %d (%d)\n", i, m);
					break;
				}

				co[m] = j;
				memcpy(co + m + 1, a, gs_1);
				for (int i0 = 0; i0 < gs_1; i0++)
				{
					const auto idx0 = co[m + i0];
					if (idx0 == unset)
						continue;

					for (int i1 = i0 + 1; i1 < m_groupSize; i1++)
					{
						const auto idx1 = co[m + i1];
						if (idx1 != unset) {
							c[idx0 * m_numPlayers + idx1] = unset;
							c[idx1 * m_numPlayers + idx0] = unset;
						}
					}
				}
				m += m_groupSize;
			}
		}
	}
	delete[] c;
}