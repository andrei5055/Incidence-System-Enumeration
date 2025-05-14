#include "TripleSys.h"
#include "p1fCheck.h"
CC void CChecklLink::u1fSetTableRow(tchar* ro, const tchar* ri) const
{
	switch (m_groupSize) {
	case 3:
		for (int i = 0; i < m_numPlayers; i += 3)
		{
			U1FT3(i);
		}
		break;
	case 2:
		if (completeGraph()) {
			for (int i = 0; i < m_numPlayers; i += 2)
			{
				U1FT(i);
			}
			break;
		}
	default:
		for (int i = 0; i < m_numPlayers; i++)
		{
			ro[ri[i]] = i;
		}
		break;
	}

}
CC bool p1fCheck2(ctchar* u1fCycles, ctchar* neighborsi, ctchar* neighborsj, int nc)
{
	unsigned int cyclesBitsDef = 0;
	unsigned int checked = 0;
	int ncycles = 0;
	if (u1fCycles) {
		while (ncycles < MAX_CYCLES_PER_SET && u1fCycles[1 + ncycles])
			cyclesBitsDef |= 1 << u1fCycles[1 + ncycles++];
	}
	else {
		cyclesBitsDef = 1 << nc;
		ncycles = 1;
	}
	tchar k = 0;
	for (tchar m = 0; m < nc; m++)
	{
		if (!(checked & (1 << m))) {
			k = m;
			for (int i = 2; i <= nc; i += 2)
			{
				if ((k = neighborsj[neighborsi[k]]) == m) {
					if (!(cyclesBitsDef & (1 << i)))
						return false;
					if (!(--ncycles))
						return true;
					break;
				}
				checked |= 1 << k;
			}
		}
	}
	return false;
}

CC void alldata::sortCycles(tchar* length, tchar* start, int ncycles) const
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
CC int alldata::collectCyclesAndPath(TrCycles* trcAll, TrCycles* trc, eCollectionMode collectionMode) const
{
	int ncr = MAX_CYCLE_SETS;
	int iLength = sizeof(TrCycles);
	const int ncc = MAX_CYCLES_PER_SET;
	int j;
	for (j = 0; j < ncr; j++)
	{
		if (trcAll[j].counter == 0)
		{
			memcpy(&trcAll[j], trc, iLength);
			trcAll[j].counter = 1;
			break;
		}

		// add cycle to array of all cycles
		//ASSERT(!bWithoutPath); // to allow use of bWithoutPath=false you need to change code to create correct m_TrCyclesFirst2Rows
		int ic = (collectionMode == eSameSetsTogether || trc->length[0] == m_numPlayers) ? (MEMCMP(trc->length, trcAll[j].length, ncc)) :
			(MEMCMP(trc, &trcAll[j], iLength - sizeof(trc->counter)));
		 switch (ic) {
		case -1:
			if (j >= ncr) {
				ASSERT(1);
				EXIT_(1);
			}
			for (int i = ncr - 2; i >= j; i--) 
				if (trcAll[i].counter)
					memcpy(&trcAll[i + 1], &trcAll[i], iLength);
			memcpy(&trcAll[j], trc, iLength);
			trcAll[j].counter = 1;
			break;
		case 0: trcAll[j].counter++; break;
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
	maxPlayer = maxPlayer / m_groupSize * m_groupSize + m_groupSize - 1;
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
CC int alldata::getCyclesAndPathFromNeighbors(TrCycles* trc, ctchar* tt1, ctchar* tt2, ctchar* tt3, ctchar* tt4, eCheckForErrors checkErrors) const
{
	// calculate cycle(s) between two rows for group size 2 or 3.
	// return number of cycles calculated, 0 if one of the cycle not selected or -1 if full cycle set not selected.
	//        if 0 trc->ncycles = 1 and trc->length[0] equal to incorrect cycle length
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
		ASSERT(k0 >= nc);
		if (tt2tmp[k0] != unset && k0 != tt2tmp[k0]) // not used before and not common value
		{
			trc->start[ncycles] = ip / 2;
			k1 = tt2tmp[k = k0];
			int i = 0;
			for (; i < nc; i += 2)
			{
				if (k == unset || k1 == unset)
					break;
				ASSERT(ip >= nc*2);
				trc->fullPath[ip] = k;
				if (m_groupSize == 3)
					trc->fullPath[++ip] = tt3[k];
				ASSERT(ip >= nc*2 - 1);
				trc->fullPath[ip + 1] = tt1[k];
				ip += 2;
				ASSERT(k >= nc*2);
				tt2tmp[k] = tt2tmp[k1] = unset;
				k = tt2tmp[k1 = tt1[k]];

				ASSERT(ip >= nc*2);
				trc->fullPath[ip] = k1;
				if (m_groupSize == 3)
					trc->fullPath[++ip] = tt4[k1];
				ASSERT(ip >= nc*2 - 1);
				trc->fullPath[ip + 1] = tt2[k1];
				ip += 2;
			}
			tchar length = (m_groupSize == 3) ? i + i / 2 : i;

			if (checkErrors == eCheckErrors && !cycleLengthOk(length))
			{
				// in check mode (ncr == 1) report only one (first) cycle with incorrect length
				// 
				trc->start[0] = trc->start[ncycles];
				trc->length[0] = length;
				trc->ncycles = 1;
				return 0;
			}
			trc->length[ncycles++] = length;
		}
	}
	if (ncycles == 0)
		return -1;
	trc->ncycles = ncycles;
	if (ncycles > 0) {
		sortCycles(trc->length, trc->start, ncycles);
		if (cyclesNotOk(ncycles, trc->length, checkErrors))
				return -1;
		}

	if (m_TrCyclesCollection)
		collectCyclesAndPath(m_TrCyclesCollection, trc, m_collectionMode);

	return ncycles;
}
CC void getTT14ForG3(tchar* tt1, tchar* tt2, tchar* tt3, tchar* tt4, ctchar* v, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, int gn)
{
	res1 -= 2;
	res2 -= 2;
	for (int i = 0; i < gn; i++)
	{
		const auto iv = v[i];
		tt2[iv] = tt4[iv] = tt1[iv] = tt3[iv] = iv;
		auto res = res1 + t1[iv];
		//tt1[*res] = unset;
		switch (t1[iv] % 3) {
		case 0: res += 3;  
		case 2: tt1[tt1[*res] = *(res + 1)] = *res; tt3[*res] = tt3[*(res + 1)] = iv; break;
		case 1: res++;  tt1[tt1[*res] = *(res + 2)] = *res; tt3[*res] = tt3[*(res + 2)] = iv; break;
		}
		//const auto ig = t1[iv];
		//printf("i=%d iv=%d ig=%d %d %d %d\n", i, iv, ig, tt1[res1[ig / 3 * 3 + 2]], tt1[res1[ig / 3 * 3 + 3]], tt1[res1[ig / 3 * 3 + 4]]);
		res = res2 + t2[iv];
		//tt2[*res] = unset;
		switch (t2[iv] % 3) {
		case 0: res += 3;  
		case 2: tt2[tt2[*res] = *(res + 1)] = *res; tt4[*res] = tt4[*(res + 1)] = iv; break;
		case 1: res++;  tt2[tt2[*res] = *(res + 2)] = *res; tt4[*res] = tt4[*(res + 2)] = iv; break;
		}
		//ig = t1[iv];
		//printf("i=%d iv=%d ig=%d %d %d %d\n", i, iv, ig, tt2[res2[ig / 3 * 3 + 2]], 
		// tt2[res2[ig / 3 * 3 + 3]], tt2[res2[ig / 3 * 3 + 4]]);
	}
}
CC int alldata::p3Cycles(TrCycles* trc, ctchar* t1, ctchar* t2, ctchar* v, ctchar* res1, ctchar* res2,
	eCheckForErrors checkErrors) const
{
	ASSERT(m_groupSize != 3);
	tchar tt1[MAX_PLAYER_NUMBER], tt2[MAX_PLAYER_NUMBER];
	tchar tt3[MAX_PLAYER_NUMBER], tt4[MAX_PLAYER_NUMBER];

	getTT14ForG3(tt1, tt2, tt3, tt4, v, t1, t2, res1, res2, m_nGroups);
	int iret;
	iret = getCyclesAndPathFromNeighbors(trc, tt1, tt2, tt3, tt4, checkErrors);
#if 0
	printf("\niret=%b", iret);
	printTable("ri", res1, 1, m_numPlayers, m_groupSize);
	printTable("rm", res2, 1, m_numPlayers, m_groupSize);
	printTable("v", v, 1, m_nGroups, 0);
	printTable("fullPath", trc->fullPath, 1, m_numPlayers * 2, m_groupSize);
	printTable("length", trc->length, 1, MAX_CYCLES_PER_SET, 0);
	printTable("start", trc->start, 1, MAX_CYCLES_PER_SET, 0);
#endif
	if (iret > 0 && m_TrCyclesCollection)
		collectCyclesAndPath(m_TrCyclesCollection, trc, m_collectionMode);
	return iret;
}
CC int alldata::u1fGetCycleLength(TrCycles* trc, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, 
	eCheckForErrors checkErrors) const
{
	// returns > 0  ok, 0 - one of the cycle not in any list, -1 - cycle set not in list
	// calculate cycle(s) length for rows res1, res2.
	// t1, t2 - precalculated arrays with 
	// for group size = 2: neighbor for each player (for example t1[7] - neighbor of player 7 in row res1)
	// for group size = 3: position of player (for example t2[3] - position of player 3 in row res2)
	if (!completeGraph())
		return u1fGetCycleLengthCBMP(trc, t1, t2, res1, res2, checkErrors);
	else {
		switch (m_groupSize) {
		case 2:
			return getCyclesAndPathFromNeighbors(trc, t1, t2, NULL, NULL, checkErrors);
		case 3:
			return u1fGetCycleLength3(trc, t1, t2, res1, res2, checkErrors);
		}
	}
	return -1;
}
CC int alldata::u1fGetCycleLength3(TrCycles * trc, ctchar * t1, ctchar * t2, ctchar * res1, ctchar * res2, 
	eCheckForErrors checkErrors) const
{
	if (!completeGraph())
	{
		auto v0 = getV0();
		ctchar* r = res1;
		for (int i = 0, k = 0; i < m_numPlayers; i += m_groupSize, k++) {
			for (int j = 0; j < m_groupSize; j++) {
				v0[m_groupSizeRemainder[r[i + j]] * m_nGroups + k] = r[i + j];
			}
		}
		ctchar* vtr = v0; // Common Values
		for (int i = 0; i < m_groupSize; i++, vtr += m_nGroups)
		{
			auto ncycles = p3Cycles(trc, t1, t2, vtr, res1, res2, checkErrors);
			if (cyclesNotOk(ncycles, trc->length, checkErrors))
				return ncycles ? -1 : 0;
		}
		return 1;
	}
	tchar us[MAX_PLAYER_NUMBER];
	tchar v[MAX_GROUP_NUMBER];
	tchar t2d3[MAX_PLAYER_NUMBER];
	for (int it2 = m_numPlayers; it2--;)
		t2d3[it2] = t2[it2] / 3; // used in macros below
	//printTable("t1", t1, 1, m_numPlayers, 3);
	//printTable("t2", t2, 1, m_numPlayers, 3);
	int ncycles = 0;
	memset(us, 0, m_numPlayers);
	switch (m_numPlayers) {
		case 9: {
			P3Cycles3();
			P3CyclesCheck(trc);
			P3Cycles3EndAndReturn(ncycles);
		}
		case 15: {
			P3Cycles5();
			P3CyclesCheck(trc);
			P3Cycles5EndAndReturn(ncycles);
		}
		case 21: {
			P3Cycles7();
			P3CyclesCheck(trc);
			P3Cycles7EndAndReturn(ncycles);
		}
		case 27: {
			P3Cycles9();
			P3CyclesCheck(trc);
			P3Cycles9EndAndReturn(ncycles);
		}
		default: 
			ASSERT(1); EXIT_(1);
	}

	return ncycles;
}

CC bool alldata::matrixStat(ctchar* table, int nr, bool *pNeedOutput)
{
	/**???
	if (m_groupSize > 3 && !pNeedOutput)
		return true;**/
	if (!pNeedOutput && m_allowUndefinedCycles)
		return true;
	bool ret = true;
	const auto nc = m_numPlayers;
	const eCheckForErrors checkErrors = pNeedOutput ? eNoErrorCheck : eCheckErrors;
	if (pNeedOutput) {
		memset(m_TrCyclesAll, 0, sizeof(TrCycles) * MAX_CYCLE_SETS);
		m_TrCyclesCollection = m_TrCyclesAll;
		m_collectionMode = eSameSetsTogether;
	}
	else if (m_use2RowsCanonization /** && !param(t_u1f) **/ && m_groupSize == 3)
	{
		//if (!(m_p1f_counter % 10000000))
		//	printTable("LR", result(nr - 1), 1, nc, m_groupSize);
		if (pNeedOutput && param(t_p1f_counter) && (++m_p1f_counter >= param(t_p1f_counter)) && nr > 2) {
			m_p1f_counter = 0;
			if (!((this->*m_pCheckFunc)(nr, nr - 1))) {
#if !USE_CUDA
				if (m_bPrint)
					reportCurrentMatrix();
#endif
				return false;
			}
		}
	}
	for (int m = nr - 1; m > 0; m--) // start from last row to have option to exit loop if we run it for new row only
	{
		auto* rowm = table + m * nc;
		auto* rowi = table;
		int iend = m;
		for (int i = 0; i < iend; i++, rowi += nc)
		{
			if (checkErrors == eCheckErrors)
			{
				const auto ncycles = u1fGetCycleLength(&m_TrCycles, rowi, rowm, result(i), result(m), eCheckErrors);
				// in case of incorrect one cycle length u1fGetCycleLength reports only one cycle (with error)
				if (ncycles == 0)
				{
					tchar* fp = m_TrCycles.fullPath + m_TrCycles.start[0] * 2;
					int fpLength = m_TrCycles.length[0] * 2;
					adjustPlayerPosition(fp, fpLength, nr);
					ret = false;
					break;
				}
				if (ncycles == -1)
				{
					ret = false;
					break;
				}
			}
			else
			{
				int ncycles = u1fGetCycleLength(&m_TrCycles, rowi, rowm, result(i), result(m), eNoErrorCheck);
#if 0
				if (nr == m_numDaysResult) {
				printf("\nRows %d:", i);
				printTable("", result(i), 1, m_numPlayers, m_groupSize);
				printf("Rows %d:", m);
				printTable("", result(m), 1, m_numPlayers, m_groupSize);
				}
#endif
			}
		}
		if (checkErrors == eCheckErrors)
			return ret; // check new row only
	}
	if (pNeedOutput) {
		*pNeedOutput = true;
		m_TrCyclesCollection = NULL;
		m_collectionMode = eNoCollection;
	}
	return ret;
}
char *alldata::matrixStatOutput(char* str, int maxStr, TrCycles* trs) const
{
	char* pOut = str;
	SPRINTFS(pOut, str, maxStr, "Cycles:");
	const auto retVal = pOut;
	for (int i = 0; i < MAX_CYCLE_SETS && trs[i].counter; i++)
	{
		SPRINTFS(pOut, str, maxStr, "%d(", trs[i].counter);
		auto pCycles = trs[i].length;
		for (int j = 0; j < MAX_CYCLES_PER_SET && pCycles[j]; j++)
		{
			if (j)
				SPRINTFS(pOut, str, maxStr, ":");
			SPRINTFS(pOut, str, maxStr, "%d", pCycles[j]);
		}
		SPRINTFS(pOut, str, maxStr, ") ");
	}
	return retVal;
}

CC bool CChecklLink::cyclesNotOk(int ncycles, tchar* length, eCheckForErrors eCheck) const
{
	if (eCheck == eNoErrorCheck)
		return false;
	if (ncycles <= 0)
		return true;
	if (m_allowUndefinedCycles)
		return false;
	auto pntr = m_param->u1fCycles[0];
	if (!pntr)
		return ncycles == 1 && length[0] != m_numPlayers;
	const auto ngrp = pntr[0];
	pntr++;
	for (int i = 0; i < ngrp; i++) {
		if (!MEMCMP(pntr, length, ncycles))
			return false;
		pntr += MAX_CYCLES_PER_SET;
	}
	return true;
}

CC bool CChecklLink::cycleLengthOk(tchar length) const
{
	if (m_allowUndefinedCycles)
		return true;
	auto pntr = m_param->u1fCycles[0];

	if (!pntr || (pntr[0] == 1 && pntr[1] == m_numPlayers))
		return length == m_numPlayers;
	const auto ngrp = pntr[0];
	pntr++;
	tchar symb;
	for (int i = 0; i < ngrp; i++) {
		int k = 0;
		while (symb = pntr[k++])
			if (symb == length)
				return true;
		pntr += MAX_CYCLES_PER_SET;
	}
	return false;
}
	
CC tchar checkForUnexpectedCycle(ctchar iv, ctchar ic, ctchar nc, ctchar* lnk, ctchar* v)
{
	// v - array with values already in current row
	// iv - candidate for v[ic]
	// 0 v1  v2 v3  v4 v5 ...
	// ... a1 b1 ... a0 b0 - 2 groups from v (b0 = candidate)
	// sw checks that there are two groups in some of prior: [a1 a0...b1 b0] or [a1 b0...b1 a0]
	// if yes, then there is a loop with length 4
	auto a0 = v[ic - 1];
	auto* la0 = lnk + a0 * nc;
	auto b0 = iv;
	auto* lb0 = lnk + b0 * nc;
	for (tchar i = 0; i < ic - 1; i += 2)
	{
		auto a1 = v[i];
		auto b1 = v[i + 1];
		auto a0a1 = *(la0 + a1); // from link table for a0,a1 - day when group used, or unset
		if (a0a1 != unset) {
			auto b0b1 = *(lb0 + b1);
			if (a0a1 == b0b1)
				return iv + 1;
		}
		auto b0a1 = *(lb0 + a1);
		if (b0a1 != unset) {
			auto a0b1 = *(la0 + b1);
			if (b0a1 == a0b1)
				return iv + 1;
		}
#define UseU1FCheck6 0 // if enabled then 50%-100% more time needed
#if UseU1FCheck6
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
			U1FCheck6(a0a1, b1a2, b0b2, b1b2, b0a2);
			auto a0a2 = *(la0 + a2);
			auto a0b2 = *(la0 + b2);
			U1FCheck6(b0a1, b1a2, a0b2, b1b2, a0a2);
			auto a1b2 = *(la1 + b2);
			U1FCheck6(a0a2, a1b2, b0b1, b1b2, b0a1);
			U1FCheck6(b0a2, a1b2, a0b1, b1b2, a0a1);
			auto a1a2 = *(la1 + a2);
			U1FCheck6(a0b1, a1a2, b0b2, a1b2, b0a2);
			U1FCheck6(b0b1, a1a2, a0b2, a1b2, a0a2);
			U1FCheck6(a0b2, a1a2, b0b1, b1a2, b0a1);
			U1FCheck6(b0b2, b1a2, a0a1, a1a2, a0b1);
		}
#endif
	}
	return iv;
#if UseU1FCheck6
	nextPlayer6 : return iv + 1;
#endif
}
CC int alldata::p1fCheck2ndRow() const
{
	if (m_groupSize > 3)
		return 0;

	if (m_createSecondRow)
		return 0;

	const auto nc = m_numPlayers;
	const tchar* p2ndRow = result(1);
	if (1)//m_groupSize == 3)
	{
		ctchar* p;
		for (int is = 0; is < m_pSecondRowsDB->numObjects(); is++)
		{
			p = m_pSecondRowsDB->getObject(is);
			const int icmp = MEMCMP(p2ndRow, p, nc);
			if (icmp <= 0)
				return icmp;
		}
		return 1;
	}

	return MEMCMP(p2ndRow, m_pSecondRowsDB->getObject(0), nc); // only one second row for group size = 2
}

CC void alldata::p1fCheckStartMatrix(int nr) 
{
	bool bCBMP = !completeGraph();
	auto u1fPntr = sysParam()->u1fCycles[0];
	for (int i = 1; i < nr; i++)
	{
		TrCycles trCycles;
		int iret;
		if (bCBMP)
			iret = getCyclesAndPathCBMP(&trCycles, neighbors(0), neighbors(i), result(0), result(i), 0, eCheckErrors) > 0;
		else
			iret = getCyclesAndPathFromNeighbors(&trCycles, neighbors(0), neighbors(i), NULL, NULL, eCheckErrors) > 0;
#if 0
		if (iret) {
			if ((!u1fPntr && trCycles.ncycles != 1) || (u1fPntr && MEMCMP(u1fPntr+1, trCycles.length, trCycles.ncycles)))
				iret = 0;
		}
#endif
		CUDA_PRINTF("*** p1fCheck DONE for rows = 0,%d  iret = %d\n", i, iret);
		ASSERT(iret <= 0, 
			printfRed("*** Error in input 'Start matrix' - rows (0, %d) are not with requested cycles), Exit\n", i);
			printTable("Incorrect 'Start matrix'", result(), nr, m_numPlayers, m_groupSize);
			myExit(1);
		)
	}
}
CC int alldata::getAllV0(tchar* allv, int maxv, ctchar* t2, ctchar* res1) const
{
	// Get up to maxv sets of "common" values from rows ir1, ir2.
	// Each value in one set of "common" values present only in one group of row ir1 and in one group of row ir2.
	int nc = m_numPlayers;
	int gn = m_nGroups;
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
		default: ASSERT(1); EXIT_(1);
	}
	return nv;
}
CC int alldata::getAllV(tchar* allv, int maxv, ctchar* neighbors1, ctchar* result0) const
{
	if (completeGraph()) {
		const auto nv = getAllV0(allv, maxv, neighbors1, result0);
		return nv;
	}
	// creat (groupSize x m_nGroups) array of values of row ir1 (sorted by color)
	ctchar* r = result0;
	for (int i = 0, k = 0; i < m_numPlayers; i+= m_groupSize, k++) {
		for (int j = 0; j < m_groupSize; j++) {
			allv[m_groupSizeRemainder[r[i + j]] * m_nGroups + k] = r[i + j];
		}
	}
	return m_groupSize;
}
CC int alldata::getAllV(tchar* allv, int maxv, int ir0, int ir1) const
{
	return getAllV(allv, maxv, neighbors(ir1), result(ir0));
}
void alldata::cyclesFor2Rows(TrCycles* trcAll, TrCycles* trc, ctchar* neighbors0, ctchar* neighbors1, ctchar* result0, ctchar* result1,
	eCollectionMode collectionMode)
{
	memset(trc, 0, sizeof(TrCycles));
	m_collectionMode = collectionMode;
	m_TrCyclesCollection = trcAll;
	if (trcAll)
		memset(trcAll, 0, sizeof(TrCycles) * MAX_CYCLE_SETS);
	if (!completeGraph())
		u1fGetCycleLengthCBMP(trc, neighbors0, neighbors1, result0, result1, eCheckErrors);
	else {
		switch (m_groupSize) {
			case 2:
				getCyclesAndPathFromNeighbors(trc, neighbors0, neighbors1, NULL, NULL, eCheckErrors);
				break;
			case 3: {
				auto v0 = getV0();
				const int nv0 = getAllV(v0, m_maxCommonVSets, neighbors1, result0); // ! neighbor1 can be different than neighors(1)
				ASSERT(!nv0);
				ctchar* vtr = v0; // Common Values
				for (int i = 0; i < nv0; i++, vtr += m_nGroups)
					p3Cycles(&m_TrCycles, neighbors0, neighbors1, vtr, result0, result1, eCheckErrors);
				break;
			}
		}
	}
	m_TrCyclesCollection = NULL;
}
