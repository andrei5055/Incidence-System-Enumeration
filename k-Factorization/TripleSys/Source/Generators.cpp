
#include "kOrbits.h"

template<>
int Generators<unsigned char>::testNestedGroups(const CGroupInfo* pElemGroup, CGroupInfo* pRowGroup, int rowMin, CKOrbits* pKOrb) const {
	const auto pntr = (const alldata*)pElemGroup;
	if (pntr->param(t_nestedGroups) > 1)
		rowMin = 2;
	else
		if (!pRowGroup && !pKOrb)
			return -1;

	tchar ts[MAX_PLAYER_NUMBER];
	ctchar* mi = pntr->result();
	CGroupInfo* pRowGroupOut = NULL;
	const auto groupOrder = pElemGroup->numObjects();
	const auto rowMax = pntr->numDaysResult();
	for (int j = rowMin; j <= rowMax; j++) {
		if (j == rowMax && !(pRowGroupOut = pRowGroup) && !pKOrb)
			break;

		const auto len = pElemGroup->lenObject() * j;
		for (int i = 0; i < groupOrder; i++) {
			const auto pElemPerm = pElemGroup->getObject(i);
			pntr->kmSortMatrixForReorderedPlayers(mi, j, pElemPerm, ts, false);
			if (MEMCMP(mi, pntr->transformedMatrix(), len))
				return j;

			if (pRowGroupOut)
				pRowGroupOut->updateRepo(ts);

			if (pKOrb)
				pKOrb->updateKSetGroup(mi, pElemPerm, ts);
		}
	}
	return rowMin == rowMax ? -1 : 0;
}