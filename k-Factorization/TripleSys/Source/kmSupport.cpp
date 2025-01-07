#include "TripleSys.h"
#include "p1fCheck.h"

#define SetUpToTwoGroups2(i, m, j, ta, tb) \
		if (ta[i]) \
		{m[j] = ta[i]; m[j + 1] = tb[i]; j += 2;} \
		if (ta[i + 1]) \
		{m[j] = ta[i + 1]; m[j + 1] = tb[i + 1]; j += 2;}
#define SetUpToThreeGroups3(i, m, j, ta, tb, tc) \
		if (ta[i]) \
		{m[j] = ta[i]; m[j + 1] = tb[i]; m[j + 2] = tc[i]; j += 3;} \
		if (ta[i + 1]) \
		{m[j] = ta[i + 1]; m[j + 1] = tb[i + 1]; m[j + 2] = tc[i + 1]; j += 3;} \
		if (ta[i + 2]) \
		{m[j] = ta[i + 2]; m[j + 1] = tb[i + 2]; m[j + 2] = tc[i + 2]; j += 3;}

#define sort3bytesInMemory(p)	sort3bytes(*p, *(p + 1), *(p + 2))

#define SetTaTb(i, m, ta, tb, a, b) \
		if ((a=tr[m[i]]) > (b=tr[m[i+1]])) tb[ta[b] = b] = a; else tb[ta[a] = a] = b;

#define SetTaTbTc(i, m, ta, tb, tc, a, b, c) a=tr[m[i]]; b=tr[m[i+1]]; c=tr[m[i+2]]; \
		sort3bytes(a, b, c); \
		ta[a] = a; tb[a] = b; tc[a] = c

CC void kmSortGroupsByFirstValue(tchar* mo, ctchar* mi, tchar nr, tchar nc, tchar np)
{
	ctchar* pmi[MAX_PLAYER_NUMBER];
	memset(pmi, 0, np * sizeof(pmi[0]));
	auto* mic = mi;

	for (tchar ir = 0; ir < nr; ir++)
	{
		pmi[*mic] = mic;
		mic += nc;
	}
	const auto iMax = np - nc + 1;
	switch (nc)
	{
		case 2: {
			auto* mos = (short int*)mo;
			for (tchar i = 0; i < iMax; i++)
			{
				if (auto *mis = (short int*)pmi[i])
				{
					*mos++ = *mis;
				}
			}
			break;
		}
		case 3: {
			auto* moc = mo;
			for (tchar i = 0; i < iMax; i++)
			{
				if (auto* mii = (int*)(pmi[i]))
				{
					*(int *)moc = *mii;
					moc += 3;
				}
			}
			break;
		}
		case 4: {
			auto* moi = (int*)mo;
			for (tchar i = 0; i < iMax; i++)
			{
				if (auto* mii = (int*)(pmi[i]))
				{
					*moi++ = *mii;
				}
			}
			break;
		}
		default: {
			auto* moc = mo;
			for (tchar i = 0; i < iMax; i++)
			{
				if (mic = pmi[i])
				{
					memcpy(moc, mic, nc);
					moc += nc;
				}
			}
		}
	}
}
CC void kmSortRowsBy2ndValue(tchar* mo, ctchar* mi, char nr, char nc, tchar* tm)
{
	ctchar* pmi[MAX_PLAYER_NUMBER];
	memset(pmi, 0, nc * sizeof(pmi[0]));
	auto* mic = mi;

	for (char ir = 0; ir < nr; ir++)
	{
		const auto iv = *(mic + 1);
		ASSERT(iv >= nc || iv < 0);// || pmi[iv] != NULL);
		tm[iv] = ir;
		pmi[iv] = mic;
		mic += nc;
	}
	auto* moc = mo;
	for (char i = 0, j = 0; i < nc; i++)
	{
		mic = pmi[i];
		if (mic)
		{
			tm[j++] = tm[i];
			memcpy(moc, mic, nc);
			moc += nc;
		}
	}
}

CC void alldata::kmSortGroups2(tchar* mi, int nr) const
{
	for (int i = nr * m_nGroups; i--; mi += 2)
	{
		const auto t = *mi;
		if (t > mi[1])
		{
			*mi = mi[1]; mi[1] = t;
		}
	}
}
CC void alldata::kmSortGroups(tchar* mi, int nr) const
{
	for (int i = nr * m_nGroups; i--; mi += m_groupSize)
	{
		for (int j = 1; j < m_groupSize; j++)
		{
			bool ex = true;
			for (int k = 0; k < m_groupSize - j; k++)
			{
				if (mi[k] > mi[k + 1])
				{
					ex = false;
					SWAP(mi[k], mi[k + 1]);
				}
			}
			if (ex)
				break;
		}
	}
}
CC void alldata::kmSortGroups3(tchar* mi, int nr) const
{
	for (int i = nr * m_nGroups; i--; mi += 3)
	{
		sort3bytesInMemory(mi);
	}
}
CC void kmTranslate(tchar* mo, ctchar* mi, ctchar* tr, int len)
{
	for (int i = len; i--;)
	{
		mo[i] = tr[mi[i]];
	}
}
CC int alldata::kmProcessMatrix(ctchar* mi, ctchar* tr, int nr, tchar ind) const
{
	tchar* mo = m_Km;
	const auto nc = m_numPlayers;
	const auto len = nc * nr;
	if (tr)
		kmTranslate(mo, mi, tr, len);
	else 
		memcpy(mo, mi, len);

	(this->*m_pSortGroups)(mo, nr);

	auto* cii = mo;
	auto* coi = m_Ktmp;
	const auto numGroups = nc / m_groupSize;
	for (int i = 0; i < nr; i++, coi += nc, cii += nc)
	{
		kmSortGroupsByFirstValue(coi, cii, numGroups, m_groupSize, nc);
	}
	if (param(t_nestedGroups) > 1)
	{
		if (MEMCMP(mi, m_Ktmp, len) == 0)
			return 0;
	}
	// result of the loop above is in m_Ktmp, sort and send it to mo
	tchar tm[MAX_PLAYER_NUMBER];
	kmSortRowsBy2ndValue(mo, m_Ktmp, nr, nc, tm);
	auto dayMax = tm[0];
	auto miFrom = mi;
	coi = mo;
	int nrr = param(t_useRowsPrecalculation);
	bool bPrecalcRow = m_useRowsPrecalculation == eCalculateRows && nr > nrr && *(mi + nc * nrr + 1) != nrr + 1;
	if (bPrecalcRow)
		nr = nrr;
	for (int i = 0; i < nr; i++, coi += nc, mi += nc)
	{
		ASSERT(tm[i] >= nc);
		switch (MEMCMP(coi, mi, nc))
	    {
		case -1: setPlayerIndex(tr, dayMax, tm[i], coi, mi, miFrom + nc * tm[i], nc); return -1;
		case 0: if (dayMax < tm[i]) { dayMax = tm[i]; } break;
		case 1: return coi[1] == mi[1] ? 1 : 2;
		}
	}
	return (bPrecalcRow || param(t_nestedGroups) > 1) ? 1 : 0;
}
CC void alldata::setPlayerIndexByPos(ctchar* tr, ctchar* co, ctchar* ciFrom, int iDayMax, int iDayCurrent, int nc, int ip) const
{
	if (!tr)
		return;
	if (iDayMax <= iDayCurrent)
	{
		tchar ttr[MAX_PLAYER_NUMBER];
		tchar tpr[MAX_PLAYER_NUMBER];
		int i, j = 0;
		for (i = 0; i < nc; i++)
		{
			tpr[ciFrom[i]] = ttr[tr[i]] = i;
		}
		for (i = 0; i <= ip; i++)
		{
			if (j < tpr[ttr[co[i]]])
				j = tpr[ttr[co[i]]];
		}
		int playerIndex = iDayCurrent * m_numPlayers + j;
		if (m_playerIndex > playerIndex)
			m_playerIndex = playerIndex;
#define ProfilePlayerIndex 0 && !USE_CUDA
#if ProfilePlayerIndex 
		int ipMax = (iDayCurrent + 1) * m_numPlayers - m_groupSize - 1;
		int n = ipMax - m_playerIndex;
		Stat("0", 6, n == 0);
		Stat("1-4", 7, n > 0 && n < 5);
		Stat("5-10", 8, n > 4 && n < 11);
		Stat("11-15", 9, n > 10 && n < 16);
		Stat(">15", 10, n > 15);
#endif
	}
}
CC void alldata::setPlayerIndex(ctchar* tr, int iDayMax, int iDayCurrent, ctchar* co, ctchar* ci, ctchar* ciFrom, int nc) const
{
	if (iDayMax <= iDayCurrent)
	{
		int i = 0;
		for (; i < nc - m_groupSize - 1; i++)
		{
			if (co[i] < ci[i])
				break;
		}
		setPlayerIndexByPos(tr, co, ciFrom, iDayMax, iDayCurrent, nc, i);
	}
}
CC int alldata::kmProcessOneNot1stRow2(ctchar* mi, int mind, tchar* tb, ctchar* tr, int nr, int irow) const
{
	const auto nc = m_numPlayers;
	auto* mii = mi + mind * nc;
	auto* ta = m_tx;
	memset(ta, 0, nc);
	tchar va, vb;
	switch (nc)
	{
	case 18:
		SetTaTb(16, mii, ta, tb, va, vb)
	case 16:
		SetTaTb(14, mii, ta, tb, va, vb)
	case 14:
		SetTaTb(12, mii, ta, tb, va, vb)
	case 12:
		SetTaTb(10, mii, ta, tb, va, vb)
	case 10:
		SetTaTb(8, mii, ta, tb, va, vb)
		SetTaTb(6, mii, ta, tb, va, vb)
		SetTaTb(4, mii, ta, tb, va, vb)
		SetTaTb(2, mii, ta, tb, va, vb)
		SetTaTb(0, mii, ta, tb, va, vb)
		break;
	default:
		for (int j = 0; j < nc; j += 2)
		{
			SetTaTb(j, mii, ta, tb, va, vb)
		}
	}
	const auto row2ndValue = tb[0];

	const int ncc = (row2ndValue - 1) * nc;
	if (row2ndValue == irow)
	{
		if (tb[1] > mi[ncc + 3])
		{
			return 1;
		}
	}
	auto* mon = m_Ktmp + ncc; // needed for macros below
	tchar m = 4;
	int j;
	if (nc >= 4)
	{
		// this is not first row, first 4 output values always in ta[0,1], tb[0,1]:
		mon[0] = ta[0];
		mon[1] = tb[0];
		mon[2] = ta[1];
		mon[3] = tb[1];
	}
	switch (nc)
	{
	case 10:
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		break;
	case 12:
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		SetUpToTwoGroups2(10, mon, m, ta, tb)
		break;
	case 14:
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		SetUpToTwoGroups2(10, mon, m, ta, tb)
		SetUpToTwoGroups2(12, mon, m, ta, tb)
		break;
	case 16:
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		SetUpToTwoGroups2(10, mon, m, ta, tb)
		SetUpToTwoGroups2(12, mon, m, ta, tb)
		SetUpToTwoGroups2(14, mon, m, ta, tb)
		break;
	default:
		j = nc >= 4 ? 2 : 0;
		m = 2 * j;
		for (; j < nc; j += 2)
		{
			SetUpToTwoGroups2(j, mon, m, ta, tb)
		}
	}
	if (row2ndValue == irow)
		return MEMCMP(mon, mi + ncc, nc);
	return 0;
}
CC int alldata::kmProcessOneNot1stRow3(tchar* mo, ctchar* mi, int mind, tchar* tb, tchar* tc, ctchar* tr, int nr) const
{
	int nc = m_numPlayers;
	auto* mii = mi + mind * nc;
	auto* ta = m_tx;
	memset(ta, 0, nc);
	tchar va, vb, vc;
	for (int j = 0; j < nc; j += 3)
	{
		SetTaTbTc(j, mii, ta, tb, tc, va, vb, vc);
	}

	mi += nc;
	const auto row2ndValue = tb[0];
	if (row2ndValue == 3)
	{
		const char diff = (tc[0] == mi[2]) ? tb[1] - mi[4] : tc[0] - mi[2];
		if (diff > 0)
			return 1;
	}
	const auto ncc = (row2ndValue - 1) * nc;
	auto* mon = mo + ncc; // needed for macros below
	tchar m;
	int j;
	if (nc >= 9)
	{
		m = 9;
		j = 3;
		// this is not first row, first 3 groups always present in ta, tb, tc [0,1,2]
		mon[0] = ta[0]; mon[1] = tb[0]; mon[2] = tc[0];
		mon[3] = ta[1]; mon[4] = tb[1]; mon[5] = tc[1];
		mon[6] = ta[2]; mon[7] = tb[2]; mon[8] = tc[2];
	}
	else
	{
		j = m = 0;
	}
	nc -= 2;
	for (; j < nc; j += 3)
	{
		SetUpToThreeGroups3(j, mon, m, ta, tb, tc)
	}

	return row2ndValue == 3 ? MEMCMP(mon, mi, nc) : 0;
}
CC int alldata::kmProcessMatrix2p1f(tchar* tr, int nr, int ind0, int ind1)
{
	int iRet;
	tchar tb[MAX_PLAYER_NUMBER], tm[MAX_PLAYER_NUMBER];
	bool bProc3 = false;
	const auto nc = m_numPlayers;
	const auto mi = result();
	tb[0] = 1; // to make compiler happy
	memset(tm, unset, nc);
	auto rowMax = (tchar)(MAX2(ind0, ind1));
	char row2ndValue = 0;

	int nrr = param(t_useRowsPrecalculation);
	bool bPrecalcRow = m_useRowsPrecalculation == eCalculateRows && nr > nrr && *(mi + nc * nrr + 1) != nrr + 1;
	for (tchar i = 0; i < nr; i++)
	{
		if (i == ind0 || i == ind1)
			continue;
		iRet = kmProcessOneNot1stRow2(mi, i, tb, tr, nr, 3);
		row2ndValue = tb[0];
		if (row2ndValue == 3 && (!bPrecalcRow || nrr > 2))
		{
			if (iRet)
			{
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, m_Ktmp + nc * 2, mi + nc * 2, mi + nc * i, nc);
					return  -1;
				}
	//			printf("%d %d %d:", i, ind0, ind1);
				return iRet;
			}
			rowMax = MAX2(rowMax, i);
			bProc3 = true;
		}
		else if (row2ndValue <= nr && !bPrecalcRow) // cant use quick check if row2ndValue > nr, or precalculation of rows nrr and up
		{
			if (bProc3 && row2ndValue == 4)
			{
				iRet = MEMCMP(m_Ktmp + nc * 3, mi + nc * 3, nc);
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, m_Ktmp + nc * 3, mi + nc * 3, mi + nc * i, nc);
					return -1;
				}
				else if (iRet == 1)
					return 1;
			}
		}
		tm[row2ndValue] = i;
	}
	if (!bProc3 && (!bPrecalcRow || nrr > 2))
	{
		return 2;
	}

	if (bPrecalcRow)
		nr = nrr;
	
	for (int i = 4; i <= nr; i++)
	{
		if (tm[i] == unset) // all values of tm are >= 0; unset indicates that row is missing
			return 2;
		const auto shift = nc * (i - 1);
		iRet = MEMCMP(m_Ktmp + shift, mi + shift, nc);
		switch (iRet)
		{
		case -1: setPlayerIndex(tr, rowMax, tm[i], m_Ktmp + shift, mi + shift, mi + nc * tm[i], nc); return -1;
		case 0: rowMax = MAX2(tm[i], rowMax); continue;
		case 1: return 1;
		}
	}
	return bPrecalcRow ? 1 : 0;
}
CC int alldata::kmProcessMatrix2(ctchar* mi, ctchar* tr, int nr, tchar ind) const
{
	int iRet;
	tchar tb[MAX_PLAYER_NUMBER], tm[MAX_PLAYER_NUMBER];
	bool bProc2 = false;
	const auto nc = m_numPlayers;
	tb[0] = 1; // to make compiler happy
	memset(tm, unset, nc);
	auto rowMax = tm[1] = ind;
	auto r2ind = m_Km2ndRowInd[ind];
	if (r2ind >= nr)
		r2ind = 1;
	if (r2ind == ind)
		r2ind = (ind + 1) % nr; // function kmProcessOneNot1stRow2 cant be used for 1st row

	iRet = kmProcessOneNot1stRow2(mi, r2ind, tb, tr, nr);
	auto row2ndValue = tb[0];
	if (row2ndValue == 2)
	{
		if (iRet)
		{
			if (iRet == -1)
			{
				setPlayerIndex(tr, rowMax, r2ind, m_Ktmp + nc, mi + nc, mi + nc * r2ind, nc);
				return -1;
			}
			return iRet;
		}
		rowMax = MAX2(rowMax, r2ind);
		bProc2 = true;
	}
	ASSERT(row2ndValue >= nc);
	tm[row2ndValue] = r2ind;

	for (tchar i = 0; i < nr; i++)
	{
		if (i == r2ind || i == ind)
			continue;
		iRet = kmProcessOneNot1stRow2(mi, i, tb, tr, nr);
		if ((row2ndValue = tb[0]) == 2)
		{
			m_Km2ndRowInd[ind] = i;
			if (iRet)
			{
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, m_Ktmp + nc, mi + nc, mi + nc * i, nc);
					return -1;
				}
				return iRet;
			}
			rowMax = MAX2(rowMax, i);
			bProc2 = true;
		}
		else if (row2ndValue <= nr) // cant return if row2ndValue > nr
		{
			if (bProc2 && row2ndValue == 3)
			{
				iRet = MEMCMP(m_Ktmp + nc * 2, mi + nc * 2, nc);
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, m_Ktmp + nc * 2, mi + nc * 2, mi + nc * i, nc);
					return -1;
				}
				else if (iRet == 1)
					return 1;
			}
		}
		tm[row2ndValue] = i;
	}
	if (!bProc2)
	{
		return 2;
	}
	for (int i = 3; i <= nr; i++)
	{
		if (tm[i] == unset) // all values of tm are >= 0; unset indicates that row is missing
			return 2;
		const auto shift = nc * (i - 1);
		switch (MEMCMP(m_Ktmp + shift, mi + shift, nc))
		{
		case -1: setPlayerIndex(tr, rowMax, tm[i], m_Ktmp + shift, mi + shift, mi + nc * tm[i], nc); return -1;
		case 0: rowMax = MAX2(tm[i], rowMax); continue;
		case 1: return 1;
		}
	}
	return 0;
}
CC int alldata::kmProcessMatrix3(ctchar* mi, ctchar* tr, int nr, tchar ind) const
{
	int iRet;
	tchar tb[MAX_PLAYER_NUMBER], tc[MAX_PLAYER_NUMBER], tm[MAX_PLAYER_NUMBER];
	bool bProc2 = false;
	tb[0] = 1; // to make compiler happy
	const auto nc = m_numPlayers;
	memset(tm, unset, nc);
	auto rowMax = tm[0] = ind; // indices of input rows in result matrices
	auto r2ind = m_Km2ndRowInd[ind];
	if (r2ind >= nr)
		r2ind = 1;
	if (r2ind == ind)
		r2ind = (ind + 1) % nr; // function kmProcessOneNot1stRow3 cant be used for 1st row
	tchar* mo = m_Ktmp;
	iRet = kmProcessOneNot1stRow3(mo, mi, r2ind, tb, tc, tr, nr);
	auto row2ndValue = tb[0];
#if 0
	Stat("kmProcessMatrix3", 0, true);
	Stat("2nd row first", 1, row2ndValue == 3);
	Stat("fast returns", 2, row2ndValue == 3 && iRet);
#endif
	if (row2ndValue == 3) // 2nd row
	{
		if (iRet)
		{
			if (iRet == -1)
			{
				setPlayerIndex(tr, rowMax, r2ind, mo + nc * 2, mi + nc, mi + nc * r2ind , nc);
				return -1;
			}
			return iRet;
		}
		rowMax = MAX2(rowMax, r2ind);
		bProc2 = true;
	}
	ASSERT(row2ndValue < 2);
	tm[row2ndValue - 2] = r2ind;

	for (tchar i = 0; i < nr; i++)
	{
		if (i == r2ind || i == ind)
			continue;
		iRet = kmProcessOneNot1stRow3(mo, mi, i, tb, tc, tr, nr);
		row2ndValue = tb[0];
		ASSERT(row2ndValue >= nc);

		if (row2ndValue == 3) // 2nd row
		{
			m_Km2ndRowInd[ind] = i; // update index for next time
			if (iRet)
			{
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, mo + nc * 2, mi + nc, mi + nc * i, nc);
					return -1;
				}
				return iRet;
			}
			rowMax = MAX2(rowMax, i);
			bProc2 = true;
		}
		else if (row2ndValue <= nr) // cant return if row2ndValue > nr
		{
			if (bProc2 && row2ndValue == 4) // 3rd row and 2nd row was processed cmp was equal 0
			{
				iRet = MEMCMP(mo + nc * 3, mi + nc * 2, nc);
				if (iRet == -1)
				{
					setPlayerIndex(tr, rowMax, i, mo + nc * 3, mi + nc * 2, mi + nc * i, nc);
					return -1;
				}
				else if (iRet == 1)
					return 1;
				rowMax = MAX2(rowMax, i);
			}
		}
		ASSERT(row2ndValue < 3)
		tm[row2ndValue - 2] = i;
	}
	if (!bProc2)
	{
		return 2;
	}
	rowMax = MAX2(tm[0], tm[1]); // tm[0] - first row index, tm[1] - second row index
	auto* moi = mo + nc * 3;
	auto* mii = mi + nc * 2;
	// following loop checks all rows starting from row 3 
	for (int i = 2, j = 2; i < nc - 3 && j < nr; i++, moi += nc)
	{
		if (tm[i] == unset) // all values of tm are >= 0; unset indicates that row is missing
			continue;
		if (*(moi + 1) != *(mii + 1)) //????
			return 2;

		switch (MEMCMP(moi, mii, nc))
		{
		case -1: setPlayerIndex(tr, rowMax, tm[i], moi, mii, mi + nc * tm[i], nc); return -1;
		case 0: rowMax = MAX2(tm[i], rowMax); break;
		case 1: return 1;
		}
		j++;
		mii += nc;
	}
	return 0;
}
