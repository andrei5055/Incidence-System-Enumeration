#include "TripleSys.h"

//no recursion in checkLinsH2
#define Set1LoopH2(ct1, v0, v1, t0, t1, i1, nv) \
	const auto v0 = t0[0]; \
	const auto* ct1 = links(v0); \
	tchar* t1 = t0 + m_numPlayers; \
    memcpy(t1, t0 + 2, nv - 2);	\
    for (tchar i1 = 1; i1 < nv; i1++)  { \
		const auto v1 = t0[i1]; \
		if (ct1[v1] != unset) continue; \
		if (i1 > 1) memcpy(t1, t0 + 1, i1 - 1);
#define SetFirstLoopH2(nv) \
	Set1LoopH2(ct2, v2, v3, t1, t2, i2, nv - 2)
#define SetFirst2LoopsH2(nv) \
	SetFirstLoopH2(nv); \
	Set1LoopH2(ct3, v4, v5, t2, t3, i3, nv - 4)
#define SetFirst3LoopsH2(nv) \
	SetFirst2LoopsH2(nv); \
	Set1LoopH2(ct4, v6, v7, t3, t4, i4, nv - 6)
#define Set2ValuesH2(n, v0, v1) vo[n] = v0; vo[(n+1)] = v1
#define Set2ResultsH2() vo[0] = 0; vo[1] = nextDay
#define Set4ResultsH2() Set2ResultsH2(); \
						Set2ValuesH2(2, v2, v3)
#define Set6ResultsH2() Set4ResultsH2(); \
						Set2ValuesH2(4, v4, v5)
#define Set8ResultsH2() Set6ResultsH2(); \
						Set2ValuesH2(6, v6, v7)

#define Set2LastValues(n, t4, action) \
	if (*(links(t4[0]) + t4[1]) != unset) \
		action; \
	vo[n] = t4[0]; vo[n+1] = t4[1]
#define Set2LastValuesH2(n, t4)  Set2LastValues(n, t4, continue)
#define Set2LastValuesRH2(n, t4) Set2LastValues(n, t4, return false)

#define P1FCheck() \
	p1fSetTableRow(p1ftable(nextDay - 1), vo); \
	if (!m_p1f || p1fCheck(nextDay, vo) < 0) \
		return true

CC bool alldata::checkLinksH2(const tchar* v, int nv, int nvo, int ind1, int iday, tchar* vo)
{
	ASSERT(nv < 4 || nv > 16 || (nv & 1));
	ASSERT(v[0]);

	tchar* t1 = m_tx;
	memcpy(t1, v + 1, iday);
	if (nv - iday - 2 > 0)
		memcpy(t1 + iday, v + iday + 2, nv - iday - 2);
	const auto nextDay = iday + 1;
	switch (nv) {
		case 4:
		{
			Set2LastValuesRH2(2, t1);
			Set2ResultsH2();
			P1FCheck();
			break;
		}
		case 6:
		{
			SetFirstLoopH2(6);
			Set2LastValuesH2(4, t2);
			Set4ResultsH2();
			P1FCheck();
			}
			break;
		}
		case 8:
		{
			SetFirst2LoopsH2(8);
			Set2LastValuesH2(6, t3);
			Set6ResultsH2();
			P1FCheck();
			}}
			break;
		}
		case 10:
		{
			SetFirst3LoopsH2(10);
			Set2LastValuesH2(8, t4);
			Set8ResultsH2();
			P1FCheck();
			}}}
			break;
		}
		case 12:
		{
			SetFirst3LoopsH2(12);
			Set1LoopH2(ct5, v8, v9, t4, t5, i5, 4);
			Set2LastValuesH2(10, t5);
			Set8ResultsH2();
			Set2ValuesH2(8, v8, v9);
			P1FCheck();
			}}}}
			break;
		}
		case 14:
		{
			SetFirst3LoopsH2(14);
			Set1LoopH2(ct5, v8, v9, t4, t5, i5, 6);
			Set1LoopH2(ct6, va, vb, t5, t6, i6, 4);
			Set2LastValuesH2(12, t6);
			Set8ResultsH2();
			Set2ValuesH2(8, v8, v9);
			Set2ValuesH2(10, va, vb);
			P1FCheck();
			}}}}}
			break;
		}
		case 16:
		{
			SetFirst3LoopsH2(16);
			Set1LoopH2(ct5, v8, v9, t4, t5, i5, 8);
			Set1LoopH2(ct6, va, vb, t5, t6, i6, 6);
			Set1LoopH2(ct7, vc, vd, t6, t7, i7, 4);
			Set2LastValuesH2(14, t7);
			Set8ResultsH2();
			Set2ValuesH2(8, v8, v9);
			Set2ValuesH2(10, va, vb);
			Set2ValuesH2(12, vc, vd);
			P1FCheck();
			}}}}}}
			break;
		}
	}

	return false;
}
CC bool alldata::checkLinksH(const tchar* v, int nv, int nvo, int ind1, int ind2, tchar* vo)
{
	if (nvo <= 0)
	{
		p1fSetTableRow(p1ftable(iDay), vo -= m_numPlayers);
		if (m_p1f)
		{
			memcpy(result(iDay), vo, m_numPlayers);
			memcpy(tmpPlayers, vo, m_numPlayers);
			m_playerIndex = m_numPlayers * iDay + iPlayer - 1;
			const bool ret = matrixStat(p1ftable(), iDay + 1);
			if (!ret)
				memset(result(iDay), 0, m_numPlayers);
			return ret;
		}
		return true;
	}

	tchar t[MAX_PLAYER_NUMBER];
	if (ind1 == -1)
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

	const auto t0 = t[0];
	const auto* ct0 = links(t0);
	nvo -= 3;
	const auto nvMinus3 = nv - 3;
	tchar is = 0;
	auto ie = m_numPlayers;
	if (ind1 == -1 && ind2 != -1)
	{
		is = ind2;
		ie -= nvMinus3 / 3 - 1;
	}

	for (int i = 0; i <= nvMinus3; i++)
	{
		const auto ti = t[i+1];
		if (ti <= is)
			continue;
		if (ti >= ie)
			break;
		if (ct0[ti] != unset)
			continue;

		const auto* cti = links(ti);
		for (int j = i + 2; j < nv; j++)
		{
			const auto tj = t[j];
			if (ct0[tj] == unset && cti[tj] == unset)
			{
				* vo = t0;
				*(vo + 1) = ti;
				*(vo + 2) = tj;
				if (checkLinksH(t + 1, nvMinus3, nvo, i, j - 1, vo + 3))
				{
					return true;
				}
			}
		}
	}

	return false;
}
