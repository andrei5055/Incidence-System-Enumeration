#include "Table.h"

void Generators::makeGroupOutput(const CGroupInfo* pGroup, bool outToScreen) {
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
	const auto groupOrder = pGroup->numObjects();
	for (int i = 1; i < groupOrder; i++) {
		const auto *c = pGroup->getObject(i);
		addAutomorphism(groupDegree(), c, pOrb, true, false, true);
	}

	updateGroupOrder(groupDegree(), pOrb);
	printTable(getObject(0), false, outToScreen, numObjects(), "");
}

void Generators::savePermutation(ctchar degree, ctchar* permRow, tchar* pOrbits, bool rowPermut, bool savePermut) {
	if (m_lenStab <= stabilizerLength())
		return;

	addObject(permRow);
	m_lenStab = stabilizerLength();
}

void RowGenerators::makeGroupOutput(const CGroupInfo* pGroup, bool outToScreen) {
	if (!m_pRowGroup)
		m_pRowGroup = new CGroupInfo(lenObject(), 10);

	const auto len = pGroup->lenObject() * lenObject();
	const auto groupOrder = pGroup->numObjects();
	const auto pntr = (alldata*)pGroup;
	ctchar* mi = pntr->result();
	tchar ts[MAX_PLAYER_NUMBER];
	for (int i = 0; i < groupOrder; i++) {
		pntr->kmSortMatrixForReorderedPlayers(mi, lenObject(), pGroup->getObject(i), ts);
		ASSERT(MEMCMP(mi, pntr->transformedMatrix(), len) != 0);
		m_pRowGroup->updateGroup(ts);
	}

	if (m_outGroupMask & 4) {
		m_sName = std::format("\nOrbits and generators of Aut(M) acting on matrix rows, "
			"|Aut(R)| = {}", m_pRowGroup->orderOfGroup());
		Generators::makeGroupOutput(m_pRowGroup, outToScreen);
	}
	if (m_outGroupMask & 8) {
		m_sName = std::format("\nAut(M) acting on matrix rows, "
			"|Aut(R)| = {}", m_pRowGroup->orderOfGroup());
		COutGroupHandle::makeGroupOutput(m_pRowGroup, outToScreen);
	}
}
