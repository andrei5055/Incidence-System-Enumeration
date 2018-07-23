#include "InSysRowEquation.h"
#include "EquSystem.h"

template class CInSysRowEquation<VECTOR_ELEMENT_TYPE>;

#if USE_EXRA_EQUATIONS
template<class T>
int CInSysRowEquation<T>::resolveExtraMapping(const CVariableMapping<T> *pVarValue, const T *pMaxVal, T *pResult, CVariableMapping<T> *pVarValueRemove, CVariableMapping<T> *pVarValueOut, size_t shift) const
{
	// Function for adjustment of partialy constructed solution by the values
	// of variables obtained from the system of extra equations
	// Parameters are:
	// pVarValue - mapping of variable indices to their values

	// The array of indices which stores minimal index for {t_singleNoLambda, t_singleLambda, t_dual} 
	// we need to start comparison for current variable. These indicies could be changed in varMapping(i)->findVariable.
	uint idx[t_dual + 1];
	memset(idx, 0, sizeof(idx));

	int lambda = 0;
	const auto *pSingleVar = varMapping(0);
	const auto *pSingleVarLambda = varMapping(1);
	const auto *pPair = pVarValue->getMappingPntr() + shift - 2;
	const auto *pPairLast = pVarValue->getLastMapping();
	while ((pPair += 2) < pPairLast) {
		int i = -1; 
		// When we call CVariableMapping::findVariable, perhaps, we need to adjust variable ID when i == t_dual
		// NOTE: We should do it only when the variable with the ID *pPair is the "lambda"-variable.
		//       For now (Apr. 9, 2016) we could be here only for t-design and because of that all variables are 
		//       the "lambda"-variable
		//       If it' not true, change flag to 0
		const int flag = 1;
		while (++i <= t_dual && !varMapping(i)->findVariable(*pPair + flag && i == t_dual ? 1 : 0, idx + i));
		if (i > t_dual)
			continue;

		const int diff = *(pPair + 1) - *(pResult + *pPair);
		switch (i) {
			case t_singleNoLambda:
			case t_singleLambda:
				// The values for current variable obtained from right part AND
				// from extra equations should be the same
				if (diff)
					return -1;   // they are not

				continue;
		}
		
		int index;
		// Because in CDualLambda we keep index of second variable from the pair,
		// we need revert the condition to be checked here
		if (flag == 1) {
			pVarValueRemove->removeMapping(*pPair);
			lambda += diff;
			index = 1;
		}
		else {
			// this part was not debugged and we should not be here when flag = 1
			index = -1;
		}

		// New variable value cannot be bigger than it's maximum
		if (*(pMaxVal + *pPair) < *(pPair + 1))
			return -1;

		*(pResult + *pPair) = *(pPair + 1);

		// The value of another variable from the pair is defined
		// Let's keep it and try to use this value to simplify the system of extra equations
		const VECTOR_ELEMENT_TYPE newVar = *pPair + index;

		// Set variable value for another variable from the pair
		// and add this variable to the mapping for output (which will be used on the next iteration)
		pVarValueOut->addMapping(newVar, *(pResult + newVar) -= diff);
	}

	return lambda;
}

template<class T>
size_t CInSysRowEquation<T>::excludeVariables(CVariableMapping<T> *pVarValue, int adj) const
{ 
	if (pVarValue->isEmpty())
		return 0;

	return equSystem()->excludeVariables(pVarValue, adj, adj > -2);
}

template<class T>
CVariableMapping<T> *CInSysRowEquation<T>::addDefinedVariables(CVariableMapping<T> *pResMapping, T *pResult, CVariableMapping<T> *pVariation) const
{
	auto **pVar = equSystem()->varPntr();
	// NOTE: For t-design the variables from varMapping(t_singleNoLambda) could not appear in extra equations
	for (int i = t_DesigneEnum() ? t_singleLambda : t_singleNoLambda; i < t_dual; i++) {
		const auto pVarMapping = varMapping(i);
		if (pVarMapping->isEmpty())
			continue;

		pResMapping = pVarMapping->addDefinedVariables(pResMapping, pResult, pVar);
	}

	// Add variables from varMapping(t_dual) which were defined 	
	// either as 0's OR as their max values by currently used right part 
	// These variable should be NOT in pVariation
	return varMapping(t_dual)->addDefinedVariables(pResMapping, pResult, pVar, pVariation);
}

template<class T>
CVariableMapping<T> *CInSysRowEquation<T>::AdjustExtraEquationRightParts(CVariableMapping<T> *pVarList, bool addVarr)
{
	int adj;
	CVariableMapping<T> *pResMaxList = NULL;
	if (!pVarList) {
		// Create list of variable which 
		//    a) have minimal values > 0 AND 
		//    b) are included in the extra equations
		static int uuu = 0; uuu++;
		pVarList = varMapping(t_dual)->addVariablMinValue(variableMinValPntr(), equSystem()->varPntr(), variableMaxValPntr(), &pResMaxList, pVarList);

		if (pResMaxList)
			excludeVariables(pResMaxList, 1);

		if (!pVarList)
			return NULL;

		adj = 0;
	}
	else {
		// To revert the variable's values in the list of variable minimal values
		// We need to add variables, and not just their values variables 
		// when we are adding the variables, the forcibly is defined by the right parts of the equations
		adj = addVarr? -2 : -1;
	}

	// Adjust the right parts of the extra equations
	excludeVariables(pVarList, adj);
	return pVarList;
}

#endif