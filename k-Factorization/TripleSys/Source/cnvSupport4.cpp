#include "TripleSys.h"
CC bool alldata::cyclesOfTwoRowsOk(TrCycles* trc) const
{
	// check that all requested cycles are present
	ctchar* u1fPntr = sysParam()->u1fCycles[0];
	if (allowUndefinedCycles()) {
		if (!u1fPntr)
			return true; // all cycles are welcome //trc[0].length[0] == m_numPlayers && trc[1].counter == 0;
		else {
			
			// check that selected in params cycles are present
			const int nc = *u1fPntr;
			u1fPntr++;
			for (int itr1 = 0; itr1 < nc; itr1++, u1fPntr += MAX_CYCLES_PER_SET) {
				for (int itr0 = 0; itr0 < MAX_3PF_SETS; itr0++)
				{
					if (trc[itr0].counter == 0)
						return false;
					if (!MEMCMP(trc[itr0].length, u1fPntr, MAX_CYCLES_PER_SET))
						//return true;
						break;
				}
			}
			return true;
		}
	}
	// check that all and only requested cycles are present
	if (!u1fPntr) {
		return trc[0].length[0] == m_numPlayers && !trc[1].counter;
		//return true;
	}
	else {
		const int nc = *u1fPntr;
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
		return trc[itr].counter == 0;
	}
	return true;
}