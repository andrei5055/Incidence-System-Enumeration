#include "TripleSys.h"

int compare_fn(const void* pA, const void* pB) {
	return *(tchar*)pA - *(tchar*)pB;
}

void alldata::checkCommonValues()
{
	if (m_groupSize != 3 || numDaysResult() != 2)
		return;
	tchar v[MAX_GROUP_NUMBER * MAX_3PF_SETS];
	const int nv = getAllV(v, MAX_3PF_SETS, 0, 1);
	for (int i = 0; i < nv; i++)
	{
		std::qsort(v + i * m_nGroups, m_nGroups, 1, compare_fn);
	}
	tchar* vk = v;
	for (int k = 0; k < nv; k++, vk += m_nGroups)
	{
		if (*vk == unset)
			continue;
		tchar* vi = vk;
		for (int i = k + 1; i < nv; i++)
		{
			if (*(vi += m_nGroups) == unset)
				continue;
			tchar vo[MAX_GROUP_NUMBER];
			for (int j = 0; j < orderOfGroup(); j++)
			{
				auto* cmpTr = getObjAddr(j);
				kmTranslate(vo, vi, cmpTr, m_nGroups);
				std::qsort(vo, m_nGroups, 1, compare_fn);
				if (memcmp(vo, vk, m_nGroups) == 0)
				{
					*vi = unset;
					//printfYellow("v%2d of %2d is equal to v%2d (tr%d)\n", i, nv, k, j);
					goto nexti;
				}
			}
		nexti: continue;
		}
		printf("v%2d of %2d canonical for all %d tr:", k, nv, orderOfGroup());
		printTable("", vk, 1, m_nGroups, 0);
	}
	printf("\nEnd of comparision of %d sets of common values.\n\n", nv);

	checkCommonValues(v, nv);
}

void alldata::checkCommonValues(ctchar* pBaseValues, int numSets) {
	auto* pOrbits = new tchar[numSets + m_nGroups];
	auto* pSet = pOrbits + numSets;
	memset(pOrbits, unset, numSets);
	auto* pCurrSet = pBaseValues;
	for (int i = 0; i < numSets; i++, pCurrSet += m_nGroups) {
		int j = orderOfGroup() - 1;
		auto* pPerm = getObjAddr(0);
		while (j--) {
			kmTranslate(pSet, pCurrSet, pPerm += m_numPlayers, m_nGroups);
			qsort(pSet, m_nGroups, sizeof(tchar), compare_fn);
			if (memcmp(pCurrSet, pSet, m_nGroups) > 0)
				break;
		}

		if (j < 0)
			pOrbits[i] = i;
	}

	printfRed("\nNon-isomorphic sets of base elements are:\n");
	char buffer[256];
	int n = 0;
	for (int i = 0; i < numSets; i++) {
		if (pOrbits[i] == i) {
			n++;
			char* pbuf = buffer;
			SPRINTFD(pbuf, buffer, "%3d: ", i);
			auto* pCurrSet = pBaseValues + i * m_nGroups;
			for (int k = 0; k < m_nGroups; k++)
				SPRINTFD(pbuf, buffer, " %3d ", pCurrSet[k]);

			printfRed("%s\n", buffer);
		}
	}

	printfRed("\nOnly %d out of %d sets of orbits are non-isomorphic\n", n, numSets);
	delete[] pOrbits;
}

CC bool CycleMapping::nextPerm(tchar* pPerm, int lenPerm) const {
	int j, i = lenPerm;
	auto val = pPerm[--i];
	while (i-- && pPerm[i] > val)
		val = pPerm[i];

	if (i < 0)
		return false;

	val = pPerm[j = i];
	while (++j < lenPerm && pPerm[j] > val);

	j--;
	SWAP(pPerm[i], pPerm[j]);
	j = lenPerm;
	while (++i < --j)
		SWAP(pPerm[i], pPerm[j]);

	return true;
}

CC ctchar* CycleMapping::InitCycleMapping(ctchar * const pLenCycles, ctchar* pStartCycles, int nCycles, int lenGroup, ctchar ** const pDirection, ctchar** const pStartCycleOut) {
	// Generate identical permutations for all set of cycles with same length.
	m_totalVariants = 1;
	m_currVariant = m_numCycleGroups = 0;
	int j, i = 0;
	auto* permCycles = m_permCycles = m_pNumCycles + nCycles;
	while (i < nCycles) {
		auto len = pLenCycles[j = i];
		while (++i < nCycles && pLenCycles[i] == len);
		
		len /= lenGroup;
		const auto jMax = j = i - j;
		m_pNumCycles[m_numCycleGroups++] = jMax;
		while (j--) {
			m_totalVariants *= len;
			permCycles[j] = j;
		}

		permCycles += jMax;
	}

	// reset element indices and directions.
	memset(m_pIdxElem = m_permCycles + m_nCycles, 0, 2 * (m_nCycles = nCycles));
	*pDirection = m_pDirection = m_pIdxElem + m_nCycles;
	
	// Adjusting the cycle lengths
	m_pLenCycles = m_pDirection + m_nCycles;
	for (int i = nCycles; i--;)
		m_pLenCycles[i] = pLenCycles[i] / lenGroup;

	memcpy(m_pStartCyclesIn = m_pLenCycles + m_nCycles, pStartCycles, nCycles);
	memcpy(m_pStartCyclesOut = m_pStartCyclesIn + m_nCycles, pStartCycles, nCycles);
	*pStartCycleOut = m_pStartCyclesOut;
	m_nLenGroup = lenGroup;
	return m_pIdxElem;
}

CC bool CycleMapping::ProceedToNextMapping() {
	// Attempt to construct the next set of directions 
	// for the current indices of the starting elements 
	// and the current ordering of the cycles.
	if (!m_bCMP_Graph) {
		auto i = m_nCycles;
		while (i-- && m_pDirection[i])
			m_pDirection[i] = 0;

		if (i >= 0) {
			m_pDirection[i] = 1;
			return true;			// managed to construct directions
		}
	}
	else {
		const auto val = m_pDirection[0];
		memset(m_pDirection, 1 - val, m_nCycles);
		if (!val)
			return true;
	}

	// Attempt to construct next set of indices
	if (++m_currVariant < m_totalVariants) {
		auto currVariant = m_currVariant;
		int i = 0;
		while (true) {
			m_pIdxElem[i] = currVariant % m_pLenCycles[i];
			if (m_pIdxElem[i]) {
				m_pIdxElem[i] *= m_nLenGroup;
				return true;	// managed to construct indices
			}

			currVariant /= m_pLenCycles[i++];
		}
	}
	else {
		m_currVariant = 0;
		memset(m_pIdxElem, 0, m_nCycles);
	}

	// Attempt to construct next permutation of cycles
	auto* permCycles = m_permCycles;
	auto* pStartCycleOut = m_pStartCyclesOut;
	auto* pStartCycleIn = m_pStartCyclesIn;
	for (int i = 0; i < m_numCycleGroups; i++) {
		const auto lenGroup = m_pNumCycles[i];
		if (nextPerm(permCycles, lenGroup)) {
			for (int j = 0; j < lenGroup; j++)
				pStartCycleOut[j] = pStartCycleIn[permCycles[j]];

			return true;          // managed to construct next permutation of cycles
		}
		
		auto j = lenGroup;
		while (j--)
			permCycles[j] = j;

		//memcpy(pStartCycleOut, pStartCycleIn, lenGroup);
		//MEMCMP(pStartCycleOut, pStartCycleIn, lenGroup);
		permCycles += lenGroup;
		pStartCycleIn += lenGroup;
		pStartCycleOut += lenGroup;
	}

	return false;   // there are no more directions/indices/permut_cycle combinations
}
