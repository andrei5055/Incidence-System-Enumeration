#include "TripleSys.h"

CC bool alldata::create3U1FTr(tchar* tr, TrCycles* trCycles01, TrCycles* trCycles, ctchar* dir, ctchar* offset, ctchar* start, int iPrint)
{
	m_numCycles = 0;
	if (MEMCMP(trCycles01->length, trCycles->length, MAX_CYCLES_PER_SET) != 0)
		return false;
	int ntr = 0;
	memset(tr, unset, m_numPlayers);
	//printTable("t1", trCycles->fullPath, 1, m_numPlayers * 2, 3);
	//printTable("t2", trCycles01->fullPath, 1, m_numPlayers * 2, 3);
	for (int i = 0; i < trCycles->ncycles; i++)
	{
		m_numCycles++;
		int cycleLength = trCycles->length[i] * 2;
		int ioff = offset[i] * 2;
		int ip = dir[i] ? ioff + 2 : ioff;
		int istep = dir[i] ? -1 : 1;
		for (int j = 0; j < cycleLength; j++, ip += istep)
		{
			ip = ip >= cycleLength ? ip - cycleLength : ip < 0 ? ip + cycleLength : ip;
			int i0 = trCycles->fullPath[ip + start[i] * 2];
			int i1 = trCycles01->fullPath[j + trCycles01->start[i] * 2];
			if (tr[i0] != i1)
			{
				if (tr[i0] != unset)
				{
					m_numCycles = 0;
					return false;
				}
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