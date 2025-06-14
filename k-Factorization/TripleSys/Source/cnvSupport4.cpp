#include "TripleSys.h"
CC bool alldata::cyclesOfTwoRowsOk(TrCycles* trc) const
{
	// check that all requested cycles are present
	if (trc[0].counter == 0)
		return false;
	ctchar* u1fPntr = sysParam()->u1fCycles[0];
	if (m_allowUndefinedCycles) {
		if (!u1fPntr)
			return true;  //all cycles are welcome */ trc[0].length[0] == m_numPlayers && trc[1].counter == 0;
		else {
			// check that selected in params cycles are present
			const int nc = *u1fPntr;
			u1fPntr++;
			for (int itr1 = 0; itr1 < nc; itr1++, u1fPntr += MAX_CYCLES_PER_SET) {
				for (int itr0 = 0; itr0 < MAX_CYCLE_SETS; itr0++)
				{
					if (trc[itr0].counter == 0)
						return false;
					if (!MEMCMP(trc[itr0].length, u1fPntr, MAX_CYCLES_PER_SET))
						//return true;
						break;
					if (itr0 == MAX_CYCLE_SETS)
						return false;
				}
			}
			return true;
		}
	}
	// check that all and only requested cycles are present
	if (!u1fPntr) {
		return trc[0].length[0] == m_numPlayers && (MAX_CYCLE_SETS <= 1 || trc[1].counter == 0);
		//return true;
	}
	else {
		const int nc = *u1fPntr;
		if (nc > MAX_CYCLE_SETS) {
			ASSERT(1);
			EXIT_(1);
		}
		u1fPntr++;
		int itr;
		for (itr = 0; itr < nc; itr++, u1fPntr += MAX_CYCLES_PER_SET)
		{
			if (trc[itr].counter == 0)
				return false;
			// warning! cycles defined in params must be sorted
			if (MEMCMP(trc[itr].length, u1fPntr, MAX_CYCLES_PER_SET))
				return false;
		}
		return itr == MAX_CYCLE_SETS || trc[itr].counter == 0;
	}
	return true;
}

CC int alldata::getCyclesAndPathCBMP(TrCycles* trc, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2, int istart,
	eCheckForErrors checkErrors) const
{
	// calculate cycle(s) between two cbmp-Graph rows.
	// return number of cycles calculated, 0 if one of the cycle not selected or -1 if full cycle set not selected.
	//        if 0 trc->ncycles = 1 and trc->length[0] equal to incorrect cycle length
	// res1, res2 - matrix rows
	// t1, t2 - rows values positions

	const auto nc = m_numPlayers;
	const int ncc = MAX_CYCLES_PER_SET;
	tchar res1tmp[MAX_PLAYER_NUMBER];
	int ncycles = 0;
	memset(trc, 0, sizeof(TrCycles));
	memset(res1tmp, unset, nc);
	tchar* pst = m_groups + istart * m_groupSize;
	int vp = 0;
	tchar v = 0;
	tchar k = 0, k0 = 0, k1 = 0, ip = 0;
	const tchar iGroupSizeM1 = m_groupSize - 1;
	for (; k0 < nc && ncycles < ncc; k0 += m_groupSize)
	{
		if (res1tmp[k0] == unset) // not used before 
		{
			trc->start[ncycles] = ip / 2;
			k = k0;
			int i = 0;
			for (; i < nc; i += m_groupSize)
			{
				k1 = k / m_groupSize * m_groupSize;
				if (i && res1tmp[k1] != unset)
					break;
				res1tmp[k1] = 0;
				for (int j = 0; j < m_groupSize; j++) {
					v = res1[k1 + j];
					vp = pst[m_groupSizeRemainder[v]];
					trc->fullPath[ip + vp] = v;
				}
				ip += m_groupSize;
				k = t2[trc->fullPath[ip - 1]];
				k1 = k / m_groupSize * m_groupSize;
				for (int j = 0; j < m_groupSize; j++) {
					v = res2[k1 + j];
					//vp = pst[m_groupSize - m_groupSizeRemainder[v] - 1] - is;
					//vp = pst[m_groupSizeRemainder[v]];
					vp = iGroupSizeM1 - pst[m_groupSizeRemainder[v]];
					trc->fullPath[ip + vp] = v;
				}
				ip += m_groupSize;
				k = t1[trc->fullPath[ip - 1]];
			}
			tchar length = i;

			if (checkErrors == eCheckErrors && !cycleLengthOk(length))
			{
				trc->start[0] = trc->start[ncycles];
				trc->length[0] = length;
				trc->ncycles = 1;
				return 0;
			}
			trc->length[ncycles++] = length;
			/**
			printTable("r1", res1, 1, nc, m_groupSize);
			printTable("r2", res2, 1, nc, m_groupSize);
			printf("p%d:", ncycles - 1);
			printTableColor("", trc->fullPath + trc->start[ncycles - 1] * 2, 1, trc->length[ncycles - 1] * 2, m_groupSize);
			**/
			if (*(trc->fullPath + trc->start[ncycles - 1] * 2) !=
				*(trc->fullPath + trc->start[ncycles - 1] * 2 + trc->length[ncycles - 1] * 2 - 1)) {
				//*(trc->fullPath + trc->start[ncycles - 1] * 2 + trc->length[ncycles - 1] * 2 - m_groupSize)) {
				memset(trc, 0, sizeof(TrCycles));
				return -1;
			}
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
	return ncycles;
}
CC int alldata::u1fGetCycleLengthCBMP(TrCycles* trc, ctchar* t1, ctchar* t2, ctchar* res1, ctchar* res2,
	eCheckForErrors checkErrors) const
{
	auto const n = m_maxCommonVSets;

	for (int i = 0; i < n; i++) {
		auto ncycles = getCyclesAndPathCBMP(trc, t1, t2, res1, res2, i, checkErrors);
		if (cyclesNotOk(ncycles, trc->length, checkErrors))
			return ncycles ? -1 : 0;
		if (m_TrCyclesCollection) {
			collectCyclesAndPath(m_TrCyclesCollection, trc);
#if 0 // print information about each pair cycles and path
			if (iDay == m_numDaysResult && checkErrors == eNoErrorCheck) {
				printTable("\nCycles length", trc->length, 1, ncycles, 0);
				for (int j = 0; j < ncycles; j++)
					printTable("Full Path", trc->fullPath + trc->start[j] * 2, 1, trc->length[j] * 2, m_groupSize);
			}
#endif
		}
	}
	return 1;
}


