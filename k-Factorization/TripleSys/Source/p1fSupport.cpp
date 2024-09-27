#include "TripleSys.h"

CC bool alldata::create2P1FTr(tchar* tr, tchar kStart, ctchar* pf0, ctchar* pf1, ctchar* pfi, ctchar* pfj) const
{
	tchar v1 = kStart, v2;
	tchar itr = 0;
	for (int n = 0; n < m_numPlayers; n += 2)
	{
		tr[v1] = itr;
		tr[v2 = pfi[v1]] = itr = pf0[itr];
		itr = pf1[itr];
		v1 = pfj[v2];
	}
	return true;
}
CC bool alldata::create3P1FTr1(tchar* tr, tchar k0Start, tchar k1Start, ctchar* v0, ctchar* v1,
	ctchar* t0, ctchar* t1, ctchar* res1, tchar ir0, tchar ir1, int idir, int iPrint) const
{
	tchar a0, b0, c0, a1, b1, c1;
	tchar* ti0 = p1ftable(ir0);
	tchar* ti1 = p1ftable(ir1);
	ctchar* res0 = result();
	auto* resi0 = result(ir0);
	auto* resi1 = result(ir1);

	tchar tv1[MAX_GROUP_NUMBER];
	tchar tvi1[MAX_GROUP_NUMBER];
	m_numCycles = 1;
	for (int i = 0; i < m_nGroups; i++)
	{
		tvi1[ti1[v1[i]] / 3] = tv1[t1[v0[i]] / 3] = i;
	}

	tchar k0 = k0Start;
	tchar k1 = k1Start;
	//SWAP(k0, k1);
	b0 = v0[k0];
	b1 = v1[k1];
	tchar posb = t0[b0];
	int mask = 1;
	if (idir & mask) // select a0 < c0 (clockwise), or a0 > c0 (contrclockwise)
	{
		// a0 < b0
		switch (posb % 3) {
		case 0: a0 = res0[posb + 1]; break;
		case 1: a0 = res0[posb - 1]; break;
		case 2: a0 = res0[posb - 2]; break;
		}
	}
	else
	{
		switch (posb % 3) {
		case 0: a0 = res0[posb + 2]; break;
		case 1: a0 = res0[posb + 1]; break;
		case 2: a0 = res0[posb - 1]; break;
		}
	}
	posb = ti0[b1];
	switch (posb % 3) {
	case 0: a1 = resi0[posb + 2]; break;
	case 1: a1 = resi0[posb + 1]; break;
	case 2: a1 = resi0[posb - 1]; break;
	}
	tchar v0t[MAX_GROUP_NUMBER];
	memcpy(v0t, v0, m_nGroups);
	tchar v1t[MAX_GROUP_NUMBER];
	memcpy(v1t, v1, m_nGroups);
	tchar v00[MAX_GROUP_NUMBER];
	memcpy(v00, v0, m_nGroups);
	tchar v01[MAX_GROUP_NUMBER];
	memcpy(v01, v1, m_nGroups);
	memset(tr, unset, m_numPlayers);
	int i;
	for (i = 0; i < m_nGroups; i++)
	{
		b0 = v0t[k0];
		b1 = v1t[k1];

		c0 = res0[k0 * 9 + 3 - t0[a0] - t0[b0]];
		c1 = resi0[k1 * 9 + 3 - ti0[a1] - ti0[b1]];

		if ((tr[a1] != unset && tr[a1] != a0) || (tr[b1] != unset && tr[b1] != b0) || (tr[c1] != unset && tr[c1] != c0))
			break;

		if (tr[a1] == a0 && tr[b1] == b0 && tr[c1] == c0)
		//if (tr[a1] == a0 && tr[b1] == b0 && tr[c1] == c0)
		//if (tr[a1] == a0 && tr[c1] == c0)
		{
			for (k0 = 0; k0 < m_nGroups; k0++)
			{
				if (v00[k0] != unset)
					break;
			}
			if (k0 == m_nGroups)
				break;
			for (k1 = 0; k1 < m_nGroups; k1++)
			{
				if (v01[k1] != unset)
					break;
			}
			if (k1 == m_nGroups)
				break;
			b0 = v0t[k0];
			b1 = v1t[k1];
			posb = t0[b0];
			mask = mask << 1;
			if (idir & mask) // select a0 < c0 (clockwise), or a0 > c0 (contrclockwise)
			{
				// a0 < b0
				switch (posb % 3) {
				case 0: a0 = res0[posb + 1]; break;
				case 1: a0 = res0[posb - 1]; break;
				case 2: a0 = res0[posb - 2]; break;
				}
			}
			else
			{
				switch (posb % 3) {
				case 0: a0 = res0[posb + 2]; break;
				case 1: a0 = res0[posb + 1]; break;
				case 2: a0 = res0[posb - 1]; break;
				}
			}
			posb = ti0[b1];
			switch (posb % 3) {
			case 0: a1 = resi0[posb + 2]; break;
			case 1: a1 = resi0[posb + 1]; break;
			case 2: a1 = resi0[posb - 1]; break;
			}
			c0 = res0[k0 * 9 + 3 - t0[a0] - t0[b0]];
			c1 = resi0[k1 * 9 + 3 - ti0[a1] - ti0[b1]];
			m_numCycles++;
			i--;
			if (iPrint)
				printf(" *");
			continue;
		}
		v00[k0] = v01[k1] = unset;
		if (tr[b1] != unset && tr[b1] != b0)
			break;
#if 0 // not happend
		if (tr[a1] != unset && tr[a1] != a0)
			break;
		if (tr[c1] != unset && tr[c1] != c0)
			break;
#endif
		ASSERT(a1 >= m_numPlayers || b1 >= m_numPlayers || c1 >= m_numPlayers);
		switch (iPrint){
		case 1: printf("  %2d %2d %2d", a1, b1, c1); break;
		case 2: printf("  %2d %2d %2d", a0, b0, c0); break;
		case 3: printf("  %2d    %2d", a1, c1); break;
		case 4: printf("  %2d    %2d", a0, c0); break;
		}
		tr[a1] = a0;
		tr[b1] = b0;
		tr[c1] = c0;
		tchar da0 = t1[c0] / 3;
		tchar da1 = ti1[c1] / 3;
		a0 = c0;
		a1 = c1;
		k0 = tv1[da0];
		k1 = tvi1[da1];
		b0 = v0t[k0];
		b1 = v1t[k1];
		tchar indc0 = da0 * 9 + 3 - t1[a0] - t1[b0];
		tchar indc1 = da1 * 9 + 3 - ti1[a1] - ti1[b1];
		ASSERT(indc0 >= m_numPlayers || indc1 >= m_numPlayers);
		c0 = res1[indc0];
		c1 = resi1[indc1];
		switch (iPrint) {
		case 1: printf("  %2d %2d %2d", a1, b1, c1); break;
		case 2: printf("  %2d %2d %2d", a0, b0, c0); break;
		case 3: printf("     %2d   ", b1); break;
		case 4: printf("     %2d   ", b0); break;
		}
		ASSERT(b1 >= m_numPlayers);
		if (tr[b1] != unset && tr[b1] != b0)
		{
			break;
		}
		tr[b1] = b0;
		a0 = c0;
		a1 = c1;
		k0 = t0[a0] / 3;
		k1 = ti0[a1] / 3;
	}
	if (iPrint)
		printf("\n");
	if (i == m_nGroups)
	{
		int s = 0;
		for (i = 0; i < m_numPlayers; i++)
		{
			if (tr[i] == unset)
			{
#if 0
				printf("Ks=%d i=%d j=%d clock=%d\n", kStart, ir0, ir1, idir);
				printTable("src", result(), iDay, m_numPlayers, 3);
				printTable("v0", v0, 1, 5, 1);
				printTable("v1", v1, 1, 5, 1);
				printTable("tr", tr, 1, m_numPlayers, 3);
#endif
				goto falseret;
			}
			s += 1 << tr[i];
		}
		if (s != (1 << m_numPlayers) - 1)
			goto falseret;
		return true;
	}
	falseret: m_numCycles = 0;
	return false;
}

CC bool alldata::getCyclesAndPath3(TrCycles* trc, ctchar* v, ctchar* t0, ctchar* t1, ctchar* res0, ctchar* res1)
{
	memset(trc, 0, sizeof(TrCycles));
	tchar tt1[MAX_PLAYER_NUMBER], tt2[MAX_PLAYER_NUMBER];
	tchar tt3[MAX_PLAYER_NUMBER], tt4[MAX_PLAYER_NUMBER];

	getTT14ForG3(tt1, tt2, tt3, tt4, v, t0, t1, res0, res1, m_nGroups);
	return getCyclesAndPath(trc, 1, tt1, tt2, tt3, tt4);
}

CC bool alldata::create3P1FTr(tchar* tr, tchar k0Start, tchar k1Start, ctchar* v0, ctchar* v1,
	ctchar* t0, ctchar* t1, ctchar* res1, tchar ir0, tchar ir1, int idir, int iPrint)
{
	TrCycles trCycles01;
	TrCycles trCycles;
	m_numCycles = 0;
	if (!getCyclesAndPath3(&trCycles01, v0, t0, t1, result(), res1))
		return false;
	if (!getCyclesAndPath3(&trCycles, v1, p1ftable(ir0), p1ftable(ir1), result(ir0), result(ir1)))
		return false;
	if (MEMCMP(trCycles01.length, trCycles.length, MAX_CYCLES_PER_SET) != 0)
		return false;
	int ntr = 0;
	memset(tr, unset, m_numPlayers);
	//printTable("t1", trCycles.fullPath, 1, m_numPlayers * 2, 3);
	//printTable("t2", m_TrCycles.fullPath, 1, m_numPlayers * 2, 3);
	for (int i = 0; i < MAX_CYCLES_PER_SET; i++)
	{
		if (!trCycles.length[i])
			break;
		m_numCycles++;
		int cycleLength = trCycles.length[i] / 3;
		int ioff = k0Start % cycleLength;
		ioff = ioff * 6;
		k0Start /= cycleLength;
		int ip = (idir & (1 << i)) ? ioff + 2: ioff;
		int istep = (idir & (1 << i)) ? -1 : 1;
		cycleLength *= 6;
		for (int j = 0; j < cycleLength; j += 1, ip += istep)
		{
			ip = ip >= cycleLength ? ip - cycleLength : ip < 0 ? ip + cycleLength : ip;
			int i0 = trCycles.fullPath[j + trCycles.start[i] * 2];
			int i1 = trCycles01.fullPath[ip + trCycles01.start[i] * 2];
			if (tr[i0] != i1)
			{
				if (tr[i0] != unset)
					break;
				tr[i0] = i1;
				ntr++;
			}
		}
	}
	if (ntr != m_numPlayers)
	{
		m_numCycles = 0;
		return false;
	}
	return true;
}
