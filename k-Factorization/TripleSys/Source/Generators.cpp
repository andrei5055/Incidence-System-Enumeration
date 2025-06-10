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

template<typename T>
void Generators<T>::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	// Calculate the orbits and a minimal generating set 
	// of the permutation group under its action on the element set
	this->setGroupOrder(1);
	this->releaseAllObjects();
	
	// Adding orbits:
	auto* pOrb = this->getNextObject();
	for (int i = 0; i < groupDegree(); i++)
		pOrb[i] = i;

	// ...  and trivial permutation:
	this->addObject(pOrb);
	m_lenStab = groupDegree();
	const auto groupOrder = pElemGroup->numObjects();
	for (int i = 1; i < groupOrder; i++) {
		const auto *c = pElemGroup->getObject(i);
		this->addAutomorphism(groupDegree(), c, pOrb, true, false, true);
	}

	this->updateGroupOrder(groupDegree(), pOrb);
	this->printTable(this->getObject(0), false, outToScreen, this->numObjects(), "");

	if (checkNestedGroups) {
		const auto retVal = testNestedGroups(pElemGroup);
		reportNestedGroupCheckResult(retVal, outToScreen);
	}
}

template<typename T>
int Generators<T>::testNestedGroups(const CGroupInfo* pElemGroup, CGroupInfo* pRowGroup, int rowMin, CKOrbits *pKOrb) const {
	const auto pntr = (alldata*)pElemGroup;
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
			pntr->kmSortMatrixForReorderedPlayers(mi, j, pElemGroup->getObject(i), ts, false, pKOrb);
			if (MEMCMP(mi, pntr->transformedMatrix(), len))
				return j;

			if (pRowGroupOut)
				pRowGroupOut->updateGroup(ts);
		}
	}
	return rowMin == rowMax? -1 : 0;
}

RowGenerators::RowGenerators(uint outGroupMask, int rowNumb, int groupSize)
	: Generators(outGroupMask, "", rowNumb, groupSize), m_pRowGroup(NULL) {
	m_outMask = 4;
	m_sActionOn = "matrix rows, |Aut(R)|";
}

void RowGenerators::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
	if (!m_pRowGroup)
		m_pRowGroup = new CGroupInfo(lenObject(), 10);
	else
		m_pRowGroup->releaseAllObjects();

	char errBuf[48], *pErr = NULL;
	const auto retVal = createGroup(pElemGroup);
	if (retVal > 0)
		snprintf(pErr = errBuf, sizeof(errBuf), "Nested groups check failed on row %d\n", retVal);

	if (m_outGroupMask & m_outMask) {
		m_sName = std::format("\n{}Orbits and generators of Aut(M) acting on {} = {}",
			(pErr ? pErr : ""), m_sActionOn, m_pRowGroup->orderOfGroup());
		Generators::makeGroupOutput(m_pRowGroup, outToScreen, false);

		// To avoid writing this error (if any) to the file twice
		pErr = NULL; 
	}
	if (m_outGroupMask & (m_outMask << 1)) {
		m_sName = std::format("\n{}Aut(M) acting on {} = {}",
			(pErr ? pErr : ""), m_sActionOn, m_pRowGroup->orderOfGroup());
		COutGroupHandle::makeGroupOutput(m_pRowGroup, outToScreen, false);
	}

	reportNestedGroupCheckResult(retVal, outToScreen);
}
