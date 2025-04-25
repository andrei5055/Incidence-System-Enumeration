#include "Table.h"

static void reportNestedGroupCheckResult(int retVal, bool outToScreen) {
	if (retVal > 0) {
		printfRed("Nested groups check failed on row %d\n", retVal);
		myExit(1);
	}
	else {
		if (outToScreen && !retVal)
			printfGreen("Nested groups check OK\n");
	}
}

void Generators::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	// Calculate the orbits and a minimal generating set 
	// of the permutation group under its action on the element set
	setGroupOrder(1);
	releaseAllObjects();
	
	// Adding orbits:
	auto* pOrb = getNextObject();
	for (int i = 0; i < groupDegree(); i++)
		pOrb[i] = i;

	// ...  and trivial permutation:
	addObject(pOrb);
	m_lenStab = groupDegree();
	const auto groupOrder = pElemGroup->numObjects();
	for (int i = 1; i < groupOrder; i++) {
		const auto *c = pElemGroup->getObject(i);
		addAutomorphism(groupDegree(), c, pOrb, true, false, true);
	}

	updateGroupOrder(groupDegree(), pOrb);
	printTable(getObject(0), false, outToScreen, numObjects(), "");

	if (checkNestedGroups) {
		const auto retVal = testNestedGroups(pElemGroup);
		reportNestedGroupCheckResult(retVal, outToScreen);
	}
}

void Generators::savePermutation(ctchar degree, ctchar* permRow, tchar* pOrbits, bool rowPermut, bool savePermut) {
	if (m_lenStab <= stabilizerLength())
		return;

	addObject(permRow);
	m_lenStab = stabilizerLength();
}

int Generators::testNestedGroups(const CGroupInfo* pElemGroup, CGroupInfo* pRowGroup, int rowMin) const {
	const auto pntr = (alldata*)pElemGroup;
	if (pntr->param(t_nestedGroups) > 1)
		rowMin = 2;
	else
		if (!pRowGroup)
			return -1;

	tchar ts[MAX_PLAYER_NUMBER];
	ctchar* mi = pntr->result();
	CGroupInfo* pRowGroupOut = NULL;
	const auto groupOrder = pElemGroup->numObjects();
	const auto rowMax = pntr->numDaysResult();
	for (int j = rowMin; j <= rowMax; j++) {
		if (j == rowMax && !(pRowGroupOut = pRowGroup))
			break;

		const auto len = pElemGroup->lenObject() * j;
		for (int i = 0; i < groupOrder; i++) {
			pntr->kmSortMatrixForReorderedPlayers(mi, j, pElemGroup->getObject(i), ts);
			if (MEMCMP(mi, pntr->transformedMatrix(), len))
				return j;

			if (pRowGroupOut)
				pRowGroupOut->updateGroup(ts);
		}
	}
	return rowMin == rowMax? -1 : 0;
}

void RowGenerators::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	if (!m_pRowGroup)
		m_pRowGroup = new CGroupInfo(lenObject(), 10);
	else
		m_pRowGroup->releaseAllObjects();

	char errBuf[48], *pErr = NULL;
	const auto retVal = testNestedGroups(pElemGroup, m_pRowGroup, lenObject());
	if (retVal > 0)
		snprintf(pErr = errBuf, sizeof(errBuf), "Nested groups check failed on row %d\n", retVal);

	if (m_outGroupMask & 4) {
		m_sName = std::format("\n{}Orbits and generators of Aut(M) acting on matrix rows, "
			"|Aut(R)| = {}", (pErr ? pErr : ""), m_pRowGroup->orderOfGroup());
		Generators::makeGroupOutput(m_pRowGroup, outToScreen, false);

		// To avoid writing this error (if any) to the file twice
		pErr = NULL; 
	}
	if (m_outGroupMask & 8) {
		m_sName = std::format("\n{}Aut(M) acting on matrix rows, "
			"|Aut(R)| = {}", (pErr ? pErr : ""), m_pRowGroup->orderOfGroup());
		COutGroupHandle::makeGroupOutput(m_pRowGroup, outToScreen, false);
	}

	reportNestedGroupCheckResult(retVal, outToScreen);
}
