//
//  VariableMapping.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "VariableMapping.h"

template class CVariableMapping<VECTOR_ELEMENT_TYPE>;
template class CDualLambda<VECTOR_ELEMENT_TYPE>;

template<class T>
bool CVariableMapping<T>::findVariable(T varID, uint *pIdx)
{
	// Try to find variable in the ordered array of variable IDs 
	while (*pIdx < getMapPosition()) {
		const auto curVarId = *(getMappingPntr() + *pIdx);
		if (curVarId > varID)
			return false;

		*pIdx += 2;
		if (curVarId == varID)
			return true;
	}

	return false;
}


#if USE_EXRA_EQUATIONS

template<class T>
int CVariableMapping<T>::addVarValue(T varID, T varValue, const T * const *pBoundaryMapping, bool from)
{
	// Let's search for the variable varID. Perhaps, we already added it before
	const T *pMapping = NULL;
	if (pBoundaryMapping) {
		if (from)
			pMapping = findMapping(varID, pBoundaryMapping[0], getMappingPntr());
		else
			pMapping = findMapping(varID, getLastMapping(), pBoundaryMapping[0]);
	}

	if (pMapping) {
		if (*(pMapping + 1) != varValue)
			return -2;  // current values is not which was previously defined

		return pMapping  < pBoundaryMapping[1]? 1 : 0;
	}

	printDefinedVariable(varID, varValue);
	addMapping(varID, varValue);
	return 2;
}

template<class T>
CVariableMapping<T> *CVariableMapping<T>::addVariablMinValue(const T *pMinValues, CVariable **pVar, const T *pMaxValues, CVariableMapping<T> **ppResMaxList, CVariableMapping<T> *pResList) const
{
	// Actually this method will be called only for instances of CDualLambda class ...
	CVariableMapping<T> *pResMaxList = *ppResMaxList = NULL;
	auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		const auto varIdx = *(pTo -= 2) - 1;    // ... and it's why we calculate varIdx this way

		// Skip all lambda-variables which 
		//    a) have no minimal values or 
		//    b) are not in the system of extra equations
		if (!*(pMinValues + varIdx) || !*(pVar + varIdx))
			continue;

		if (*(pMinValues + varIdx) == *(pMaxValues + varIdx)) {
			if (!pResMaxList)
				*ppResMaxList = pResMaxList = new CVariableMapping(nElement());
			
			pResMaxList->addMapping(varIdx, *(pMinValues + varIdx));
		}
		else {

			if (!pResList)
				pResList = new CVariableMapping(nElement());

			pResList->addMapping(varIdx, *(pMinValues + varIdx));
		}
	}

	return pResList;
}

template<class T>
CVariableMapping<T> *CDualLambda<T>::addDefinedVariables(CVariableMapping<T> *pResMapping, const T *pResult, CVariable **pVar, const CVariableMapping<T> *pVariation) const
{
	const auto pLastVar = pVariation->getLastMapping();
	auto pCurrVar = pVariation->getMappingPntr();

	// Define pLastMapping which will be used as a boundary for the search of previously added mappings in pResMapping->findMapping(...) call
	const auto pLastMapping = pResMapping && !pResMapping->isEmpty() ? pResMapping->getLastMapping() : NULL;
	const VECTOR_ELEMENT_TYPE *pMapBoundaries[] = { pLastMapping, NULL };
	const auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		const auto varIdx = *(pTo -= 2) - 1;
		if (!*(pVar + varIdx))
			continue;

		// Variable is used in the system of extra equations
		while (pCurrVar < pLastVar && *pCurrVar > varIdx)
			pCurrVar += 2;

		if (pCurrVar < pLastVar) {
			if (*pCurrVar == varIdx) {
				// We found variable in the list of undefined variables
				pCurrVar += 2;
				continue;
			}
		}

		// The variable is NOT in the list of undefined variables
		// It means that its value already defined and, perhaps, should be added to pResMapping
		pResMapping->addVarValue(varIdx, *(pResult + varIdx), pMapBoundaries);
	}

	return pResMapping;
}

template<class T>
CVariableMapping<T> *CVariableMapping<T>::addDefinedVariables(CVariableMapping<T> *pResMapping, const T *pResult, CVariable **pVar, const CVariableMapping<T> *pVariation) const
{
	const auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		if (!*(pVar + *(pTo -= 2)))
			continue;

		// Variable is used in the system of extra equations
		pResMapping->addMapping(*pTo, *(pResult + *pTo));
	}

	return pResMapping;
}

#endif