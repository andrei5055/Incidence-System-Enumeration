#include "Table.h"

void Generators::outputGenerators(CGroupInfo* pGroup, bool outToScreen) {
	// Calculate the orbits and a minimal generating set 
	// of the permutation group under its action on the element set
	setGroupOrder(1);
	releaseAllObjects();
	
	// Adding orbits:
	auto* pOrb = getNextObject();
	for (int i = 0; i < groupDegree(); i++)
		pOrb[i] = i;

	//   and trivial permutation:
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

