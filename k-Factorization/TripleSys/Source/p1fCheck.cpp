#include "TripleSys.h"
#include "p1fCheck.h"
inline void p1ftmp14(tchar* ro, const tchar* ri, int nc)
{
	P1FT(0);
	P1FT(2);
	P1FT(4);
	P1FT(6);
	P1FT(8);
	P1FT(10);
	P1FT(12);
}
inline void p1ftmp16(tchar* ro, const tchar* ri, int nc)
{
	P1FT(0);
	P1FT(2);
	P1FT(4);
	P1FT(6);
	P1FT(8);
	P1FT(10);
	P1FT(12);
	P1FT(14);
}
CC void CChecklLink::p1fSetTableRow(tchar* ro, const tchar* ri) const
{
	if (m_groupSize == 3)
	{
		for (int i = 0; i < m_numPlayers; i += 3)
		{
			P1FT3(i);
		}
		if (ro[7] == ro[6])
		{
			printTable("ri", ri, 1, 15, 3);
			printTable("ro", ro, 1, 15, 3);
		}
	}
	else 
	if (m_groupSize == 2)
	{
#if !USE_CUDA
		switch (m_numPlayers) {
			case 14: p1ftmp14(ro, ri, m_numPlayers); break;
			case 16: p1ftmp16(ro, ri, m_numPlayers); break;
			default:
#endif
				for (int i = 0; i < m_numPlayers; i += 2)
				{
					P1FT(i);
				}
#if !USE_CUDA
				break;
		}
#endif
	}
}
CC int u1fCheck(ctchar* u1f, int nu1f, const tchar* rowir, int ir, int nc, const tchar* t1)
{
	// returns rejected row number (0:nr-1), or -1(no rejections)
	tchar t2[MAX_PLAYER_NUMBER];
	tchar k = 0, k1;
	memcpy(t2, rowir, nc);
	int i = ir; // needed for macros below
	for (int iu1f = 0; iu1f < nu1f; iu1f++)
	{
		if (!u1f[iu1f])
			break;
		if (iu1f)
		{
			for (k1 = 0; k1 < nc; k1++)
				if (t2[k1] != unset)
					break;
			if (k1 >= nc)
				return i;
			k = t2[k1];
		}
		for (int j = 0; j < u1f[iu1f] / 2; j++)
		{
			P1F();
		}
		P1FN();
	}
	return -1;
}
CC void alldata::sortCycles(tchar* length, tchar* start, int ncycles)
{
	for (int j = 1; j < ncycles; j++) {
		for (int i = 0; i < ncycles - j; i++) {
			if (length[i] > length[i + 1])
			{
				SWAP(length[i], length[i + 1]);
				SWAP(start[i], start[i + 1]);
			}
		}
	}
}
CC int alldata::collectCyclesAndPath(TrCycles* trc)
{
	int ncr = MAX_CYCLE_SETS;
	int iLength = sizeof(TrCycles);
	const int ncc = MAX_CYCLES_PER_SET;
	int j;
	for (j = 0; j < ncr; j++)
	{
		if (m_TrCyclesAll[j].counter == 0)
		{
			memcpy(&m_TrCyclesAll[j], trc, iLength);
			m_TrCyclesAll[j].counter = 1;
			break;
		}

		// add cycle to array of all cycles
		switch (MEMCMP(trc->length, m_TrCyclesAll[j].length, ncc)) {
		case -1:
			for (int i = ncr - 2; i >= j; i--) 
				memcpy(&m_TrCyclesAll[i + 1], &m_TrCyclesAll[i], iLength);
			memcpy(&m_TrCyclesAll[j], trc, iLength);
			m_TrCyclesAll[j].counter = 1;
			break;
		case 0: m_TrCyclesAll[j].counter++; break;
		case 1: continue;
		}
		break;
	}

	return j; // trc position in collection
}
CC void alldata::adjustPlayerPosition(tchar* path, tchar length, tchar nrows)
{
	if (length / 2 == m_numPlayers)
		return; // no adjustment can be made if full cycle
	tchar sp[MAX_PLAYER_NUMBER], *res = result(nrows - 1);
	for (int i = 0; i < m_numPlayers; i++)
		sp[res[i]] = i;

	tchar maxPlayer = 0;
	for (int i = 0; i < length; i++)
	{
		tchar ip = path[i];
		tchar ipPos = sp[ip];
		if (maxPlayer < ipPos)
		{
			//printf(" %d.%d", maxPlayer, ipPos);
			maxPlayer = ipPos;
		}
	}
	int diff = m_playerIndex - (nrows - 1) * m_numPlayers - maxPlayer;
#if 0
	StatAdd("CErr", 1, true);
	StatAdd(">0", 2, diff > 0);
	StatAdd(">3", 3, diff > 3);
	StatAdd(">6", 4, diff > 6);
	StatAdd(">12", 5, diff > 12);
#endif
	if (diff > 0)
		m_playerIndex -= diff;
}
CC int alldata::getCyclesAndPath(TrCycles* trc, int ncr, ctchar* tt1, ctchar* tt2, ctchar* tt3, ctchar* tt4)
{
	// calculate cycle(s) between two rows for group size 2 or 3.
	// return number of cycles calculated
	// tt1, tt2 contains values of "neighbor of player" (for not common elements), or value of common element
	// tt3, tt4 contains values of same group common element for each player

	const auto nc = m_numPlayers;
	const int ncc = MAX_CYCLES_PER_SET;
	tchar tt2tmp[MAX_PLAYER_NUMBER];
	int ncycles = 0;
	memset(trc, 0, sizeof(TrCycles));
	memcpy(tt2tmp, tt2, nc);
	tchar k = 0, k0, k1 = 0, ip = 0;
	for (k0 = 0; k0 < nc && ncycles < ncc; k0++)
	{
		if (tt2tmp[k0] != unset && k0 != tt2tmp[k0]) // not used before and not common value
		{
			trc->start[ncycles] = ip / 2;
			k1 = tt2tmp[k = k0];
			int i = 0;
			for (; i < nc; i += 2)
			{
				if (k == unset || k1 == unset)
					break;
				trc->fullPath[ip] = k;
				if (m_groupSize == 3)
					trc->fullPath[++ip] = tt3[k];
				trc->fullPath[ip + 1] = tt1[k];
				ip += 2;
				tt2tmp[k] = tt2tmp[k1] = unset;
				k = tt2tmp[k1 = tt1[k]];

				trc->fullPath[ip] = k1;
				if (m_groupSize == 3)
					trc->fullPath[++ip] = tt4[k1];
				trc->fullPath[ip + 1] = tt2[k1];
				ip += 2;
			}
			tchar length = (m_groupSize == 3) ? i + i / 2 : i;
			if (ncr == 1 && !cycleLengthOk(length))
			{
				// in check mode (ncr == 1) report only one (first) cycle with incorrect length
				trc->start[0] = trc->start[ncycles];
				trc->length[0] = length;
				trc->ncycles = 1;
				return 1;
			}
			trc->length[ncycles++] = length;
		}
	}
	ASSERT(ncycles == 0);
	trc->ncycles = ncycles;
	sortCycles(trc->length, trc->start, ncycles);
	if (ncr > 1)
		collectCyclesAndPath(trc);
	return ncycles;
}
CC void getTT14ForG3(tchar* tt1, tchar* tt2, tchar* tt3, tchar* tt4, const tchar* v, const tchar* t1, const tchar* t2, const tchar* res1, const tchar* res2, int gn)
{
	for (int i = 0; i < gn; i++)
	{
		const auto iv = v[i];
		auto ig = t1[iv];
		auto res = res1 + ig;
		//tt1[*res] = unset;
		tt1[iv] = iv;
		tt3[iv] = iv;
		switch (ig % 3) {
		case 0: tt1[tt1[*(res + 1)] = *(res + 2)] = *(res + 1); tt3[*(res + 1)] = tt3[*(res + 2)] = iv; break;
		case 1: tt1[tt1[*(res - 1)] = *(res + 1)] = *(res - 1); tt3[*(res - 1)] = tt3[*(res + 1)] = iv; break;
		case 2: tt1[tt1[*(res - 2)] = *(res - 1)] = *(res - 2); tt3[*(res - 2)] = tt3[*(res - 1)] = iv; break;
		}
		//printf("i=%d iv=%d ig=%d %d %d %d\n", i, iv, ig, tt1[res1[ig / 3 * 3]], tt1[res1[ig / 3 * 3 + 1]], tt1[res1[ig / 3 * 3 + 2]]);
		ig = t2[iv];
		res = res2 + ig;
		//tt2[*res] = unset;
		tt2[iv] = iv;
		tt4[iv] = iv;
		switch (ig % 3) {
		case 0: tt2[tt2[*(res + 1)] = *(res + 2)] = *(res + 1); tt4[*(res + 1)] = tt4[*(res + 2)] = iv; break;
		case 1: tt2[tt2[*(res - 1)] = *(res + 1)] = *(res - 1); tt4[*(res - 1)] = tt4[*(res + 1)] = iv; break;
		case 2: tt2[tt2[*(res - 2)] = *(res - 1)] = *(res - 2); tt4[*(res - 2)] = tt4[*(res - 1)] = iv; break;
		}
		//printf("i=%d iv=%d ig=%d %d %d %d\n", i, iv, ig, tt2[res2[ig / 3 * 3]], 
		// tt2[res2[ig / 3 * 3 + 1]], tt2[res2[ig / 3 * 3 + 2]]);
	}
}
CC int alldata::p3Cycles(int ncr, ctchar* t1, ctchar* t2, ctchar* v, ctchar* res1, ctchar* res2)
{
	ASSERT(m_groupSize != 3);
	tchar tt1[MAX_PLAYER_NUMBER], tt2[MAX_PLAYER_NUMBER];
	tchar tt3[MAX_PLAYER_NUMBER], tt4[MAX_PLAYER_NUMBER];

	getTT14ForG3(tt1, tt2, tt3, tt4, v, t1, t2, res1, res2, m_nGroups);
	int iret;
	iret = getCyclesAndPath(&m_TrCycles, ncr, tt1, tt2, tt3, tt4);
#if 0
	printf("\niret=%b", iret);
	printTable("ri", res1, 1, m_numPlayers, m_groupSize);
	printTable("rm", res2, 1, m_numPlayers, m_groupSize);
	printTable("v", v, 1, m_nGroups, 0);
	printTable("fullPath", m_TrCycles.fullPath, 1, m_numPlayers * 2, m_groupSize);
	printTable("length", m_TrCycles.length, 1, MAX_CYCLES_PER_SET, 0);
	printTable("start", m_TrCycles.start, 1, MAX_CYCLES_PER_SET, 0);
#endif
	return iret;
}
CC int alldata::u1fGetCycleLength(int ncr, const tchar* t1, const tchar* t2, const tchar* res1, const tchar* res2, int ind)
{
	// calculate cycle(s) length for rows res1, res2.
	// t1, t2 - precalculated arrays with 
	// for group size = 2: neighbor for each player (t1[7] - neighbor of player 7 in row res1)
	// for group size = 3: position of player (t1[7] - position of player 7 in row res2)
	int ncycles = 0;
	switch (m_groupSize) {
		case 2:
			return getCyclesAndPath(&m_TrCycles, ncr, t1, t2);
		case 3: {
			tchar us[MAX_PLAYER_NUMBER];
			tchar v[MAX_GROUP_NUMBER];
			tchar t2d3[MAX_PLAYER_NUMBER];
			for (int it2 = m_numPlayers; it2--;)
				t2d3[it2] = t2[it2] / 3; // used in macros below
			//printTable("t1", t1, 1, m_numPlayers, 3);
			//printTable("t2", t2, 1, m_numPlayers, 3);
			memset(us, 0, m_numPlayers);
			switch (m_numPlayers) {
				case 9: {
					P3Cycles3();
					P3CyclesCheck();
					P3Cycles3EndAndReturn(ncycles);
				}
				case 15: {
					P3Cycles5();
					P3CyclesCheck();
					P3Cycles5EndAndReturn(ncycles);
				}
				case 21: {
					P3Cycles7();
					P3CyclesCheck();
					P3Cycles7EndAndReturn(ncycles);
				}
				case 27: {
					P3Cycles9();
					P3CyclesCheck();
					P3Cycles9EndAndReturn(ncycles);
				}
			}
		}
	}
	return 0;
}
CC bool alldata::matrixStat(ctchar* table, int nr, bool *pNeedOutput)
{
	if (m_groupSize > 3)
		return true;
	// program returns false if there is not full cycle
	bool ret = true;
	const auto nc = m_numPlayers;
	const auto ncr = pNeedOutput ? MAX_CYCLE_SETS : 1;

	memset(&m_TrCyclesAll, 0, sizeof(m_TrCyclesAll));

	if (m_p1f && !param(t_u1f) && m_groupSize == 3 && !pNeedOutput)
	{
		m_p1f_counter++;
		if (!(m_p1f_counter % 10000000))
			printTable("LR", result(nr - 1), 1, nc, m_groupSize);
		if (!(m_p1f_counter % param(t_p1f_counter)) && nr > 2)
			return true;
	}
	for (int m = nr - 1; m > 0; m--) // start from last row to have option to exit loop if we run it for new row only
	{
		auto* rowm = table + m * nc;
		auto* rowi = table;
		int iend = m;
		for (int i = 0; i < iend; i++, rowi += nc)
		{
#if 1
			if (ncr == 1)
			{
				const auto ncycles = u1fGetCycleLength(ncr, rowi, rowm, result(i), result(m));
				// in case of incorrect cycle length u1fGetCycleLength reports only one cycle (with error)
				if (ncycles == 1 && !cycleLengthOk(m_TrCycles.length[0]))
				{
					tchar* fp = m_TrCycles.fullPath + m_TrCycles.start[0] * 2;
					int fpLength = m_TrCycles.length[0] * 2;
					adjustPlayerPosition(fp, fpLength, nr);
					ret = false;
					break;
				}
				if (cyclesNotOk(ncr, ncycles, m_TrCycles.length))
				{
					ret = false;
					break;
				}
			}
			else
#endif
			{
				int ncycle = u1fGetCycleLength(ncr, rowi, rowm, result(i), result(m));
#if 0
				printTable("ri", result(i), 1, m_numPlayers, m_groupSize);
				printTable("rm", result(m), 1, m_numPlayers, m_groupSize);
				printTable("fullPath", m_TrCycles.fullPath, 1, m_numPlayers * 2, m_groupSize);
				printTable("length", m_TrCycles.length, 1, MAX_CYCLES_PER_SET, 0);
				printTable("start", m_TrCycles.start, 1, MAX_CYCLES_PER_SET, 0);
#endif
			}
		}
		if (ncr == 1)
			return ret; // check new row only
	}
	if (pNeedOutput)
		*pNeedOutput = true;
	return ret;
}
char *alldata::matrixStatOutput(char* str, int maxStr) const 
{
	char* pOut = str;
	SPRINTFS(pOut, str, maxStr, "Cycles:");
	const auto retVal = pOut;
	for (int i = 0; i < MAX_CYCLE_SETS && m_TrCyclesAll[i].counter; i++)
	{
		SPRINTFS(pOut, str, maxStr, "%d(", m_TrCyclesAll[i].counter);
		auto pCycles = m_TrCyclesAll[i].length;
		for (int j = 0; j < MAX_CYCLE_SETS && pCycles[j]; j++)
		{
			if (j)
				SPRINTFS(pOut, str, maxStr, ":");
			SPRINTFS(pOut, str, maxStr, "%d", pCycles[j]);
		}
		SPRINTFS(pOut, str, maxStr, ") ");
	}
	return retVal;
}
CC int CChecklLink::p1fCheck(int nr, const tchar* newRow) const
{
	// checks newRow to be p1f or requested u1f (by checkin newRow with existing rows)
	// returns rejected row number (0:nr-1), or -1(no rejections)
	const auto nc = m_numPlayers;
	ASSERT(nc & 1);

	const int nrMinus1 = nr - 1;
	auto* t1 = p1ftable(nrMinus1);
	CUDA_PRINTF("*** p1fCheck: before p1fSetTableRow  t1 = %p  nr = %d  nc = %d\n", t1, nr, nc);
	// already calculated p1fSetTableRow(t1, newRow);
	CUDA_PRINTF("*** p1fCheck: after p1fSetTableRow  nr = %d\n", nr);
	if (param(t_u1f))
		return u1fCheckFunc(nr, t1);	

	auto* rowi = p1ftable();
	tchar t2[MAX_PLAYER_NUMBER];
	CUDA_PRINTF("*** p1fCheck: befor switch  nc = %d\n", nc);
	switch (nc) {
#if MAX_PLAYER_NUMBER > 10
		case 16:
		{
			for (int i = 0; i < nrMinus1; i++, rowi += nc)
			{
				Stat_p1fCheck("p1f(all)", 0, true);
				tchar k = 0, k1;
				memcpy(t2, rowi, nc);
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				Stat_p1fCheck("not rej", 11, true);
			}
			return -1;
		}
		case 14:
		{
			for (int i = 0; i < nrMinus1; i++, rowi += nc)
			{
				Stat_p1fCheck("p1f(all)", 0, true);
				tchar k = 0, k1;
				memcpy(t2, rowi, nc);
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				Stat_p1fCheck("not rej", 11, true);
			}
			return -1;
		}
		case 12:
		{
			for (int i = 0; i < nrMinus1; i++, rowi += nc)
			{
				Stat_p1fCheck("p1f(all)", 0, true);
				tchar k = 0, k1;
				memcpy(t2, rowi, nc);
				P1F();
				P1F();
				P1F();
				P1F();
				P1F();
				Stat_p1fCheck("not rej", 11, true);
			}
			return -1;
		}
#endif
		default:
		{
			for (int i = 0; i < nrMinus1; i++, rowi += nc)
			{
				Stat_p1fCheck("p1f(all)", 0, true);
				tchar k = 0, k1;
				memcpy(t2, rowi, nc);
				for (int j = 2; j < nc; j += 2)
				{
					P1F();
				}
				Stat_p1fCheck("not rej", 11, true);
			}
			return -1;
		}
	}
	return -1;
}
CC bool CChecklLink::cyclesNotOk(int ncr, int ncycles, tchar* length)
{
	if (ncr != 1)
		return false;
	auto pntr = m_param->u1f[0];
	if (!pntr)
		return ncycles == 1 && length[0] != m_numPlayers;
	const auto ngrp = pntr[0];
	pntr++;
	for (int i = 0; i < ngrp; i++) {
		if (!MEMCMP(pntr, length, ncycles))
			return false;
		pntr += MAX_UNIFOM_CONF_LENGTH;
	}
	return true;
}

CC bool CChecklLink::cycleLengthOk(tchar length)
{
	auto pntr = m_param->u1f[0];
	if (!pntr)
		return length == m_numPlayers;
	const auto ngrp = pntr[0];
	pntr++;
	tchar symb;
	for (int i = 0; i < ngrp; i++) {
		int k = 0;
		while (symb = pntr[k++])
			if (symb == length)
				return true;
		pntr += MAX_UNIFOM_CONF_LENGTH;
	}
	return false;
}
CC int CChecklLink::u1fCheckFunc(const int nr, const tchar* rowm) const
{
	const auto u1f = m_param->u1f[0] + 1;
	auto nu = *(u1f - 1);
	auto u1ft = u1f + nu * MAX_UNIFOM_CONF_LENGTH;
	memset(u1ft, 0, nu);

	const auto nc = m_numPlayers;
	for (auto m = nr; --m > 0; rowm -= nc)
	{
		auto* rowi = p1ftable();
		for (int i = 0; i < m; i++, rowi += nc)
		{
			Stat_p1fCheck("u1f(all)", 0, true);
			auto u1f_j = u1f;
			for (int j = 0; j < nu; j++)
			{
				if (u1fCheck(u1f_j, MAX_UNIFOM_CONF_LENGTH, rowi, i, nc, rowm) < 0)
				{
					u1ft[j] = unset;
					goto nextRow;
				}
				if (u1ft[j] != unset)
					u1ft[j] = m;

				u1f_j += MAX_UNIFOM_CONF_LENGTH;
			}
			return m;
		nextRow:
			continue;
		}
		if (nr != nc - 1)
			break;
	}
	/**/
	if (nr == nc - 1)
	{
		while (nu--)
			if (u1ft[nu] != unset)
				return u1ft[nu];
	}
	/**/
	return -1;
}
	
CC int p1fCheckGroups(int iv, int ic, int nc, ctchar* lnk, ctchar* v)
{
	auto a0 = v[ic - 1];
	auto* la0 = lnk + a0 * nc;
	while (iv < nc)
	{
		auto b0 = iv;
		auto* lb0 = lnk + b0 * nc;
		for (int i = 0; i < ic - 1; i += 2)
		{
			auto a1 = v[i];
			auto b1 = v[i + 1];
			auto a0a1 = *(la0 + a1);
			auto b0b1 = *(lb0 + b1);

			//Stat_p1fCheck("p1f(all)", 0, true);

			if (a0a1 != unset && a0a1 == b0b1)
				goto nextPlayer;
			auto b0a1 = *(lb0 + a1);
			auto a0b1 = *(la0 + b1);
			if (b0a1 != unset && b0a1 == a0b1)
				goto nextPlayer;
#define UseP1FCheck6 0 // if enabled then 50%-100% more time needed
#if UseP1FCheck6
			auto* la1 = lnk + a1 * nc;
			auto* lb1 = lnk + b1 * nc;
			for (int j = i + 2; j < ic - 1; j += 2)
			{
				// a1:b1 a2:b2 a0:b0
				auto a2 = v[j];
				auto b2 = v[j + 1];
				auto b1a2 = *(lb1 + a2);
				auto b0b2 = *(lb0 + b2);
				auto b1b2 = *(lb1 + b2);
				auto b0a2 = *(lb0 + a2);
				P1FCheck6(a0a1, b1a2, b0b2, b1b2, b0a2);
				auto a0a2 = *(la0 + a2);
				auto a0b2 = *(la0 + b2);
				P1FCheck6(b0a1, b1a2, a0b2, b1b2, a0a2);
				auto a1b2 = *(la1 + b2);
				P1FCheck6(a0a2, a1b2, b0b1, b1b2, b0a1);
				P1FCheck6(b0a2, a1b2, a0b1, b1b2, a0a1);
				auto a1a2 = *(la1 + a2);
				P1FCheck6(a0b1, a1a2, b0b2, a1b2, b0a2);
				P1FCheck6(b0b1, a1a2, a0b2, a1b2, a0a2);
				P1FCheck6(a0b2, a1a2, b0b1, b1a2, b0a1);
				P1FCheck6(b0b2, b1a2, a0a1, a1a2, a0b1);
			}
#endif
		}
		break;
#if UseP1FCheck6
	nextPlayer6: iv++; 
		//Stat_p1fCheck("rej6", 2, true);
		break; // only one correction
#endif
	nextPlayer: iv++;
		//Stat_p1fCheck("rej4", 3, true);
		break; // only one correction
	}
	return iv;
}
CC ctchar* alldata::expected2ndRow3p1f(int iSet) const
{
	static tchar _expected2ndRow3p1f_9[9] =
	{ 0,  3,  6,   1,  4,  7,   2,  5,  8 };
	static tchar _expected2ndRow3p1f_15[2*15] = {
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 13,   8, 11, 14, //only this row is pure 3P1F
	  0,  3,  6,   1,  4,  9,   2,  7, 12,   5, 10, 13,   8, 11, 14
	};
	static tchar _expected2ndRow3p1f_21[] = {
	  0,  3,  6,   1,  4,  7,   2,  5,  8,   9, 12, 15,  10, 13, 18,  11, 16, 19,  14, 17, 20, // Use3P1F_21_669_912
#if AllowNotP1FRowsFor3P1F
	  0,  3,  6,   1,  4,  7,   2,  5,  8,   9, 12, 15,  10, 13, 18,  11, 16, 19,  14, 17, 20,
	  0,  3,  6,   1,  4,  7,   2,  5,  9,   8, 10, 12,  11, 15, 18,  13, 16, 19,  14, 17, 20,
	  0,  3,  6,   1,  4,  7,   2,  5,  9,   8, 12, 15,  10, 13, 18,  11, 16, 19,  14, 17, 20,
#endif
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 13,   8, 15, 18,  11, 16, 19,  14, 17, 20, // all cycles 21, has tr to r1, r2
#if Any2RowsConvertToFirst2 == 0
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 15,   8, 11, 18,  13, 16, 19,  14, 17, 20, // all cycles 21
#endif
#if AllowNotP1FRowsFor3P1F
#if Any2RowsConvertToFirst2 == 0
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 15,   8, 13, 18,  11, 16, 19,  14, 17, 20,
#endif
	  0,  3,  6,   1,  4,  9,   2,  7, 10,   5, 12, 15,   8, 13, 18,  11, 16, 19,  14, 17, 20,
	  0,  3,  6,   1,  4,  9,   2,  7, 12,   5, 10, 15,   8, 13, 18,  11, 16, 19,  14, 17, 20,
	  0,  3,  6,   1,  4,  9,   2,  7, 12,   5, 10, 15,   8, 16, 18,  11, 13, 19,  14, 17, 20,
	  0,  3,  6,   1,  4,  9,   2,  7, 12,   5, 13, 15,   8, 16, 18,  10, 14, 19,  11, 17, 20,
	  0,  3,  6,   1,  4,  9,   2,  7, 12,   5, 15, 18,   8, 16, 19,  10, 13, 17,  11, 14, 20,
#endif
	  0,  3,  6,   1,  9, 12,   2, 15, 18,   4, 10, 16,   5, 13, 19,   7, 11, 20,   8, 14, 17, // all cycles 21, has tr to r1, r2
	};
	static tchar _expected2ndRow3p1f_27[] = { 
		// below: one 3U1F {9,9,9} and three pure 3P1F second rows.
	  0,  3,  6,   1,  4,  7,   2,  5,  8,   9, 12, 15,  10, 13, 16,  11, 14, 17,  18, 21, 24,  19, 22, 25,  20, 23, 26, // Use3P1F_27_999
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 13,   8, 15, 18,  11, 16, 19,  14, 21, 24,  17, 22, 25,  20, 23, 26,
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 13,   8, 15, 18,  11, 21, 24,  14, 22, 25,  16, 19, 23,  17, 20, 26,
	  0,  3,  6,   1,  4,  7,   2,  9, 12,   5, 10, 15,   8, 11, 18,  13, 16, 19,  14, 21, 24,  17, 22, 25,  20, 23, 26
	};

	static tchar* _expected2ndRow3p1f[4] = { 
		_expected2ndRow3p1f_9, _expected2ndRow3p1f_15, _expected2ndRow3p1f_21, _expected2ndRow3p1f_27 
	};

	static int _expectedNumSets3p1f[] = {
		sizeof(_expected2ndRow3p1f_9),
		sizeof(_expected2ndRow3p1f_15),
		sizeof(_expected2ndRow3p1f_21),
		sizeof(_expected2ndRow3p1f_27)
	};

	static int iexp = 0;
	static ctchar* pExpected2ndRow3p1f = NULL;
	if (iSet < 0) {
		const auto nc = m_numPlayers;
		const int ind = (nc - 9) / 6;
		if (ind < countof(_expected2ndRow3p1f)) {
			if (iSet == -1) {
				iexp = (nc >= 9 && nc <= 27 && (nc % 3) == 0) ? _expectedNumSets3p1f[ind]/nc - (nc >= 21 ? 1 : 0) : 0;
				pExpected2ndRow3p1f = _expected2ndRow3p1f[ind] + (nc >= 21? nc : 0);
			}
			else {
				iexp = 1;
				pExpected2ndRow3p1f -= nc;
			}
		}
		return NULL;
	}

	return (iexp > iSet) ? pExpected2ndRow3p1f + m_numPlayers * iSet : NULL;
}
CC int alldata::p1fCheck2ndRow() const
{
	if (m_groupSize > 3)
		return 0;

	const auto nc = m_numPlayers;
	const tchar* p2ndRow = result(1);
	if (m_groupSize == 3)
	{
		if (nc < 9 || nc > 27)
			return 0;
#if GenerateSecondRowsFor3P1F
		return 0;
#endif
		ctchar* p;
		int is = 0;
		while ((p = expected2ndRow3p1f(is++)))
		{
			const int icmp = MEMCMP(p2ndRow, p, nc);
			if (icmp <= 0)
				return icmp;
		}
		return 1;
	}

	static char expectedSecondRow[7][MAX_PLAYER_NUMBER] = {
	   { 0,   2,    1,   3},
	   { 0,   2,    1,   4,    3,   5},
	   { 0,   2,    1,   4,    3,   6,    5,   7},
	   { 0,   2,    1,   4,    3,   6,    5,   8,    7,   9},
	   { 0,   2,    1,   4,    3,   6,    5,   8,    7,  10,    9,  11 },
	   { 0,   2,    1,   4,    3,   6,    5,   8,    7,  10,    9,  12,   11,  13 },
	   { 0,   2,    1,   4,    3,   6,    5,   8,    7,  10,    9,  12,   11,  14,   13,  15 } // p1f
//     { 0,   2,    1,   3,    4,   6,    5,   7,    8,  10,    9,  11,   12,  14,   13,  15 } // 4444
	};

	if (nc < 4 || nc > 16)
		return 0;
	return MEMCMP(p2ndRow, expectedSecondRow[(nc - 4) / 2], nc);
}

CC void alldata::p1fCheckStartMatrix(int nr) const
{
	for (int i = 0; i < nr; i++)
	{
		const auto irow = p1fCheck(i+1, result(i));
		CUDA_PRINTF("*** p1fCheck DONE for i = %d  irow = %d\n", i, irow);
		ASSERT(irow >= 0, 
			printfRed("*** Error in input 'Start matrix' - not p1f matrix (rows (%d, %d) are not p1f), Exit\n", irow, i);
			printTable("Incorrect 'Start p1f matrix'", result(), nr, m_numPlayers, m_groupSize, 0, true);
			myExit(1);
		)
	}
}
CC int alldata::getAllV(tchar* allv, int maxv, tchar ir1, tchar ir2, tchar* pt2) const
{
	// Get up to maxv sets of "common" values from rows ir1, ir2.
	// Each value in one set of "common" values present only in one group of row ir1 and in one group of row ir2.
	tchar* t2 = pt2 ? pt2 : p1ftable(ir2);
	int nc = m_numPlayers;
	int gn = m_nGroups;
	tchar* res1 = result(ir1);
	tchar us[MAX_PLAYER_NUMBER];
	tchar t2d3[MAX_PLAYER_NUMBER];
	tchar v[MAX_GROUP_NUMBER];
	for (int it2 = 0; it2 < nc; it2++)
		t2d3[it2] = t2[it2] / 3; // used in macros below
	memset(us, 0, MAX_PLAYER_NUMBER);
	int nv = 0;
	switch (nc) {
		case 9: {
			P3Cycles3();
			memcpy(allv + nv * gn, v, gn);
			nv++;
			if (nv >= maxv)
				return nv;
			P3Cycles3EndAndReturn(nv);
		}
		case 15: {
			P3Cycles5();
			memcpy(allv + nv * gn, v, gn);
			/**
			if (ir1 == 5 && ir2 == 6)
			{
				printTable("r", result(0), 7, 15, 3);
				printTable("v", v, 1, 5, 3);
				if (v[0] > v[1])
					nv = nv;
			}**/
			nv++;
			if (nv >= maxv)
				return nv;
			P3Cycles5EndAndReturn(nv);
		}
		case 21: {
			P3Cycles7();
			memcpy(allv + nv * gn, v, gn);
			nv++;
			if (nv >= maxv)
				return nv;
			P3Cycles7EndAndReturn(nv);
		}
		case 27: {
			P3Cycles9();
			memcpy(allv + nv * gn, v, gn);
			nv++;
			if (nv >= maxv)
				return nv;
			P3Cycles9EndAndReturn(nv);
		}
	}
	return nv;
}