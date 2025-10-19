#include "TripleSys.h"

CC bool alldata::createU1FTr(tchar* tr, const TrCycles* trCycles01, const TrCycles* trCycles, ctchar* dir, ctchar* offset, ctchar* start, int iPrint)
{
	m_numCycles = 0;
	if (MEMCMP(trCycles01->length, trCycles->length, MAX_CYCLES_PER_SET) != 0)
		return false;
	int ntr = 0;
	memset(tr, unset, m_numPlayers);
	//printTable("t1", trCycles->fullPath, 1, m_numPlayers * 2, m_groupSize);
	//printTable("t2", trCycles01->fullPath, 1, m_numPlayers * 2, m_groupSize);
	for (int i = 0; i < trCycles->ncycles; i++)
	{
		m_numCycles++;
		int cycleLength = trCycles->length[i] * 2;
		int ioff = offset[i] * 2;
		int ip = dir[i] ? ioff + m_groupSize - 1 : ioff;
		int istep = dir[i] ? -1 : 1;
		auto* fp = trCycles->fullPath + start[i] * 2;
		auto* fp01 = trCycles01->fullPath + trCycles01->start[i] * 2;
		for (int j = 0; j < cycleLength; j++, ip += istep)
		{
			if (ip >= cycleLength)
				ip -= cycleLength;
			else if (ip < 0) 
				ip += cycleLength;
			const auto i0 = fp[ip];
			const auto i1 = fp01[j];
			if (tr[i0] != i1)
			{
				if (tr[i0] != unset)
				{
					m_numCycles = 0;
					return false;
				}
				ASSERT_IF(i0 < 0 || i0 >= m_numPlayers);
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
