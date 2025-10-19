#include "TripleSys.h"
CC bool alldata::cyclesOfTwoRowsOk(TrCycles* trc) const
{
	// check that all requested cycles are present
	if (trc[0].counter == 0)
		return false;
	ctchar* u1fPntr = sysParam()->u1fCycles[0];
	if (m_allowUndefinedCycles) {
		if (!u1fPntr)
			return true;  //all cycles are welcome or trc[0].length[0] == m_numPlayers && trc[1].counter == 0;
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
						return true;
						//break;
					if (itr0 == MAX_CYCLE_SETS)
						return false;
				}
			}
		}
		return true;
	}
	// check that only requested cycles are present
	if (!u1fPntr) {
		ASSERT(MAX_CYCLE_SETS < 2);
		return trc[0].length[0] == m_numPlayers && trc[1].counter == 0;
		//return true;
	}
	else {
		const int nc = *u1fPntr;
		if (nc > MAX_CYCLE_SETS) {
			ASSERT(1);
			EXIT_(1);
		}
		u1fPntr++;
		for (int j = 0; j < MAX_CYCLE_SETS; j++)
		{
			if (trc[j].counter == 0)
				break;
			int i = 0;
			for (auto u1fPntr1 = u1fPntr; i < nc; i++, u1fPntr1 += MAX_CYCLES_PER_SET) {
				// warning! cycles defined in params must be sorted
				if (MEMCMP(trc[j].length, u1fPntr1, MAX_CYCLES_PER_SET) == 0)
					break;
			}
			if (i == nc)
				return false;
		}
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
	tchar usedGroups[MAX_GROUP_NUMBER];
	int ncycles = 0;
	resetTrCycles(trc);
	memset(usedGroups, unset, sizeof(usedGroups));
	tchar* pst = m_groups + istart * m_groupSize;
	int vp = 0;
	tchar v = 0;
	ctchar cin = 0, cout = 1;
	tchar k = 0, k0 = 0, k1 = 0, ip = 0;
	for (; k0 < m_nGroups && ncycles < ncc; k0++)
	{
		if (usedGroups[k0] == unset) // not used before 
		{
			trc->start[ncycles] = ip / 2;
			k = k0;
			int i = 0;
			for (; i < nc; i += m_groupSize)
			{
				k1 = k * m_groupSize;
				if (usedGroups[k] != unset)
					break;
				usedGroups[k] = 0;
				for (int j = 0; j < m_groupSize; j++) {
					v = res1[k1 + j];
					vp = pst[m_groupSizeRemainder[v]];
					trc->fullPath[ip + vp] = v;
				}
				k1 = t2[trc->fullPath[ip + cout]];
				k1 = k1 - m_groupSizeRemainder[k1];
				ip += m_groupSize;
				//k = t2[trc->fullPath[ip - 1]];
				for (int j = 0; j < m_groupSize; j++) {
					v = res2[k1 + j];
					//vp = m_groupSizeRemainder[v];
					//vp = iGroupSizeM1 - pst[m_groupSizeRemainder[v]];
					vp = pst[m_groupSizeRemainder[v]];
					trc->fullPath[ip + vp] = v;
				}
				k = t1[trc->fullPath[ip + cin]] / m_groupSize;
				ip += m_groupSize;
				//k = t1[trc->fullPath[ip - 1]];
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
			// check that cycle start is equal to cycle end
			if (*(trc->fullPath + trc->start[ncycles - 1] * 2 + cin) !=
				*(trc->fullPath + trc->start[ncycles - 1] * 2 + trc->length[ncycles - 1] * 2 - m_groupSize + cin)) {
				ASSERT(1);
				EXIT_(1);
				resetTrCycles(trc);
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
		auto const ncycles = getCyclesAndPathCBMP(trc, t1, t2, res1, res2, i, checkErrors);
		if (cyclesNotOk(ncycles, trc->length, checkErrors))
			return ncycles ? -1 : 0;
		if (m_TrCyclesCollection) {
			// Collect passes with different cycles (one pass per each cycle)
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


