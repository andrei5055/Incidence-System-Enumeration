#include "EquSystem.h"

#if USE_EXRA_EQUATIONS
#include "VariableMapping.h"

CEquSystem::CEquSystem(size_t nVar, size_t nRow, int t) : CForcibleCol(nVar), m_t(t)
{
	m_ppVariable = new CVariable *[nVar<<1];
	m_pVarDefined = new char[nVar];
	m_pEquDef = new CEquation *[nVar];
	resetEquNumb();

	// For 3 design it will be enough and we won't realloc that memory
	allocateMemoryForEquations(nRow);
	m_pEquArray = new CEquArray();
}

CEquSystem::~CEquSystem()	
{ 
	delete[] varPntr();
	delete[] varDefinedPtr();
	delete[] getDefEquationPntr();

	releaseVariables();
	for (auto i = equNumbMax(); i--;)
		delete equation(i);

	delete[] equPntr();
	delete[] equIndx();
	delete equArray();
}

void CEquSystem::releaseVariables()
{
	for (auto i = equNumb(); i--;)
		equation(i)->releaseEquation();

	resetEquNumb();
}

void CEquSystem::allocateMemoryForEquations(size_t nEqu)
{
	setEquNumbMax(nEqu);
	m_pEquPntr = new CEquation *[equNumbMax()];
	m_pEquIndx = new VECTOR_ELEMENT_TYPE[equNumbMax()];
	for (auto i = equNumb(); i < nEqu; i++)
		m_pEquPntr[i] = new CEquation();
}

void CEquSystem::resetVarPtr(size_t nVar)					
{
	releaseVariables();
	memset(varPntr(), 0, (m_nVar = nVar) * sizeof(*varPntr())); 
}

void CEquSystem::addVariable(CVariable *pVar, size_t numVar)
{
	auto **ppVar = varPntr() + numVar;
	if (*ppVar)
		pVar->setSameVarNextEqu(*(ppVar + nVar()));
	else
		*ppVar = pVar;

	*(ppVar + nVar()) = pVar;
}

void CEquSystem::closeVarloop() const
{
	// Make the loop of all appearance of given variable in different equations
	for (auto i = nVar(); i--;) {
		auto pVar = variable(i);
		if (pVar) {
			auto pVarNextEqu = variable(i + nVar());
			pVar->setSameVarNextEqu(pVarNextEqu);
			while (pVarNextEqu->sameVarNextEqu() != pVar)
				pVarNextEqu = pVarNextEqu->sameVarNextEqu();

			pVar->setSameVarPrevEqu(pVarNextEqu);
		}
	}
}

void CEquSystem::addEquation(CVariable *pFirstVar, size_t nVar, size_t lambda)
{
	if (equNumb() == equNumbMax()) {
		auto *pTmp = equPntr();
		auto *pIdx = equIndx();
		allocateMemoryForEquations(equNumbMax() << 1);
		memcpy(equPntr(), pTmp, equNumb() * sizeof(*pTmp));
		memcpy(equIndx(), pIdx, equNumb() * sizeof(*pIdx));
		delete[] pTmp;
		delete[] pIdx;
	}

	const auto num = equNumb();
	equPntr()[num]->initEquation(pFirstVar, nVar, lambda);
	equIndx()[num] = num;
	setEquNumb(num + 1);
}

CVariableMapping *CEquSystem::solveExtraEquations()
{
	static int hhh; ++hhh;
	CVariableMapping *pVarValue = NULL;
	// Loop over all constructed equations

	auto iMax = equNumb();
	size_t i = 0;
	while (i < iMax) {
		CEquation *pEquation = equation(i);
		if (pEquation->numbVar() != 1) {
			i++;
			continue;
		}

		// We found trivial equation
		const auto pVariable = pEquation->firstVariable();
		const auto varValue = pEquation->rightPart();

		if (!pVarValue)
			pVarValue = new CVariableMapping(iMax - i);

		pVarValue->addMapping(pVariable->varIndex(), varValue);

		// Get next equation with the current variable
		// (we need to do it BEFORE releasing this variable)
		auto pVariableNext = pVariable->sameVarNextEqu();
		// Release current equation and copy last equation of the system to the current slot
		pEquation->releaseEquation();
		if (i < --iMax)
			setEquation(equation(iMax), i);

		CVariable *pVariableTmp, *pVarTmp, *pTmp;
		while (pVariableNext != pVariable) {
			pVariableNext = (pVariableTmp = pVariableNext)->sameVarNextEqu();

			CEquation *pCurrEquation = equation(pVariableTmp);
			if (pCurrEquation->rightPart() < varValue) {
				delete pVarValue;
				return NULL;    // no solution
			}

			auto *pVarNext = pVariableTmp->nextVarSameEqu();
			if (pCurrEquation->rightPart() == varValue) {
				// The values of all remaining variables of current equation are 0's
				while (pVarNext != pVariableTmp) {
					pVarNext = (pVarTmp = pVarNext)->nextVarSameEqu();
					// Map variable with it's value 0
					pVarValue->addMapping(pVarTmp->varIndex(), 0);
					// ... and remove this variable from remaining equations
					auto *pNxt = pVarTmp->sameVarNextEqu();
					while (pNxt != pVarTmp) {
						pNxt = (pTmp = pNxt)->sameVarNextEqu();
						pTmp->prevVarSameEqu()->linkVariable(pTmp->nextVarSameEqu());
						delete pTmp;
					}
					delete pVarTmp;
				}

				// Current equation already solved and we need to delete it and move the last equation to its place
				setEquation(equation(--iMax), pCurrEquation);
			}
			else {
				if (pVarNext == pVariableTmp) {
					delete pVarValue;
					return NULL;	// Only one variable in this equation
				}

				// Adjust right part
				pCurrEquation->adjustRightPart(varValue);
				// Remove variable
				pVariableTmp->prevVarSameEqu()->linkVariable(pVarNext);
				delete pVariableTmp;
			}
		}
	}

	setEquNumb(iMax);
	return pVarValue;
}

void CEquSystem::setEquation(CEquation *pEqu, size_t idx)
{
	*(equPntr() + idx) = pEqu;

	// Change equIndx to be able find this equation by it's initial index  
	const auto pVar = pEqu->firstVariable();
	*(equIndx() + pVar->equIdx()) = idx;
}

void CEquSystem::setEquation(CEquation *pEquFrom, const CEquation *pEquTo)
{
	const auto pVar = pEquTo->firstVariable();
	setEquation(pEquFrom, *(equIndx() + pVar->equIdx()));
}

int CEquSystem::excludeVariables(CVariableMapping *pVarValue, int adj, bool exclude, size_t numVar, const CEquArray *pEquations, bool addToAllEquation)
{
	// the return value is the number of added mappings 
	auto iMax = equNumb();
	if (!iMax)
		return 0;   // Empty system of equations

	size_t idx;

	// The number of variables in the main loop we need to process, BEFORE we will need to go to the next pEquUsed 
	// (Will be used only when pEquations != NULL)
	size_t nVarToProcess = 0;

	// Iterating in different orders over the elements from pVarValue
	int step = -2;
	int jTo = 0;
	int jFrom = pVarValue->nElement() << 1;

	CEquation *pEquUsed = NULL;
	if (!exclude) {
		if (pEquations) {
			pEquUsed = pEquations->GetAt(idx = pEquations->GetSize() - 1);
			nVarToProcess = pEquUsed->numDefinedVar();
		}
	} 
	else {
		equArray()->RemoveAll(); // reset the array of used equations
		pVarValue->setMapBoundaries(numVar ? numVar : pVarValue->nElement());
	}

	const bool setDefFlag = adj > 0;
	bool addVariable = false;
	int varDefined = exclude ? 1 : 0;
	int nAddedVarTotal = 0;
	auto pMapping = pVarValue->getMappingPntr();
	bool checkFlg = false;
	int retVal;
	size_t i = 0;
/*
	if (exclude)
	{
		step = -2;
		jTo = 0;
		jFrom = pVarValue->nElement() << 1;
	} else {
		step = 2;
		jTo = (pVarValue->nElement() << 1) - 2;
		jFrom = numVar? jTo - (numVar << 1) : -2;
	}
*/
	bool reOrderMapping = false;

	for (int j = jFrom; j != jTo && (!numVar || i < numVar); i++) {
		j += step;
		const auto varIdx = *(pMapping + j);
		CVariable *pVar = variable(varIdx);
		if (!pVar)
			continue;		// Variable is not in the system
		
		if (setDefFlag)		// we are here to exclude variables and not just to adjust right parts
			setVarDefined(varIdx, varDefined);	// Mark the variable as defined

		// Value of variable - varMin:
		const auto varVal = *(pMapping + j + 1);

		// For all equations, containing pVar variable
		auto pVarCurr = pVar;
		do {
			CEquation *pEquCurr = equation(pVarCurr);
			CVariable *pRemovedVar = NULL;
			if (exclude) {
				if (!pEquCurr->solved()) {
					const size_t nTotal = pVarValue->nElement();
					printAddEquationVariable(pEquCurr, varIdx, varVal, false);
					/*const auto */retVal = pEquCurr->removeVariable(varVal, varDefinedPtr(), pVarValue, getDefEquationPntr(), adj);
					if (retVal == -2) {
						// Some equation could not be solved
						if (nAddedVarTotal) {
							// Some variables got their values here
							// We need to add these variables to the equations from which their values were defined
							excludeVariables(pVarValue, 1, false, nAddedVarTotal, equArray());
							pVarValue->removeLastMapping(nAddedVarTotal);
						}

						// We need to run this loop over the equations for current varIdx 
						// and add this variable to all equations where it was just removed
						// Define the boundaries for the loop:
						auto *pTmp = pVar->sameVarPrevEqu();
						pVar = pVarCurr;
						pVarCurr = pTmp;
						exclude = false;
						addVariable = true;
					} else
					if (retVal >= 0) {
						const size_t nAddedVar = pVarValue->nElement() - nTotal;
						if (nAddedVar) {
							nAddedVarTotal += nAddedVar;
							equArray()->Add(pEquCurr);
						}
					}
					else 
						if (setDefFlag && retVal == -1 && !numVar) {
						// We excluding the "lambda" variables from the system AND their exclusion DOES NOT eliminate current equation
						// It means that current equation still contains two or more variables (it could happend for t-design, when t >= 3)
						// Let's exclude this variable from the list of variables included in its equation.
						(pRemovedVar = pVarCurr)->excludeVariableFromEquation();

						// To be able to restore the equation, we need to keep some information for pVarCurr here ????
						// Actually, it's NOT TRUE. We need to exlude this variable only once for ALL right parts. 
						// Therefor, no need to restore equation
					}
				} else {
					// The equation is marked as "solved"
					// let's check if varIdx is in the list of indices
					if (pEquCurr->isVarIndexInIndexArray(varIdx)) {
						// If we are here, then the value of current variable is exactly 
						// the same as it would be defined from current equation. 
						// For correct rollback, we should re-order the group of variables in pMapping
						// Since we need to do that outside of this loop, let's raise flag for that
						reOrderMapping = true;
					}
				}
			}
			else {
				// when addToAllEquation is true, we add carrent variable to all equations, where from it was extracted before
				//		addToAllEquation is false, we add it only to the equation, where from it was defined 
				if (addToAllEquation || pEquCurr->solved() || pEquCurr->numDefinedVar()) {
					if (!pEquUsed || pEquCurr == pEquUsed) {
						// Is it the equation where from the current variable was defined?
						if (!checkFlg || *(getDefEquationPntr() + varIdx) == pEquCurr) {
							printAddEquationVariable(pEquCurr, varIdx, varVal, true);
							pEquCurr->addVariable(varVal);
							pVarCurr->restoreVariableInEquation();
							resetDefEquation(varIdx);

							if (!addToAllEquation) {
								// When we are here because we adding variables only to the equations 
								// where from the variables were defined, we need to decrease numDefinedVar
								// (perhaps, several variables were defined from the same equation (as 0's)
								pEquCurr->decreaseNumDefinedVar();

							}

							if (pEquUsed && !--nVarToProcess) {
								// Try to find current equation in the list of used equations
								// If it's there, it should be the last one
								if (idx--)
									pEquUsed = pEquations->GetAt(idx); // new pointer to the equation to search for 
								else
									return 0; // no more used equations, we can exit from this method
							}
						}
					}
				}
			}

			pVarCurr = pVarCurr->sameVarNextEqu();
//			delete pRemovedVar;
		} while (pVarCurr != pVar);

		if (reOrderMapping) {

			// We need to reoder elements in the mapping and make current variable stored undex index jFrom + step
			int kNext;
			size_t i1 = 0;
			VECTOR_ELEMENT_TYPE *pMap = (VECTOR_ELEMENT_TYPE *)pMapping;
			for (int k = j; i1 < i; i1++, k = kNext) {
				*(pMap + k) = *(pMap + (kNext = k - step));
				*(pMap + k + 1) = *(pMap + kNext + 1);
			}

			// kNext should be equal to (jFrom + step) at that point
			*(pMap + kNext) = varIdx;
			*(pMap + kNext + 1) = varVal;
		}

		if (addVariable) {
			addVariable = false;    // we need to do it only once when we are in that loop
			
			// We don't need to do it for current variable because we could be here only when retVal == -2
			// and in that case we did not change the equation for that variable
#if VAR_1
			if (i)
#endif								
			{
				// Exclude all previously defined variables if we have them
				excludeVariables(pVarValue, 1, false, i);
					// Mark as undefined current variable (its value was just excluded from right parts)
/* version 0
				do {
					// Remove it from the equation where its value was defined, it was not done yet
					// ANDREI   ????????
					CEquation *pEqu = getDefEquation(*(pMapping + j));
					if (pEqu)
						pEqu->addVariable(*(pMapping + j + 1));
				} while ((j += 2) < jMax);

/*	
//					setVarDefined(varIdx, 0);
				// Exclude all previously defined variables
				excludeVariables(pVarValue, 1, false, i);

				// Remove all variables scheduled for exclusion from the list:
				pVarValue->removeLastMapping(numVar);
				*/
			}

			// Everyone of remaining scheduled variables should be excluded from 
			// the extra equation where it was defined
			exclude = false;
			varDefined = 0;
			nAddedVarTotal = -1;
			checkFlg = true;
			pEquUsed = NULL;
			if (numVar > 1) {
/*				pEquations = equArray(); 
				idx = pEquations->GetSize()
				pEquUsed = pEquations->GetAt(idx = pEquations->GetSize() - 1);
				nVarToProcess = pEquUsed->numDefinedVar();
*/
				// To do the correct roll-back, we need to revert the order of variables in the mapping
				step = 2;
				jTo = jFrom;
				j = jTo - (numVar << 1) - 2;
/*
				int jLast = jMax - (numVar << 1);
				VECTOR_ELEMENT_TYPE *pMap = (VECTOR_ELEMENT_TYPE*)pMapping;
				for (j = jMax; (j -= 2) > jLast; jLast += 2) {
					for (int k = 0; k < 2; k++) {
						const auto tmp = *(pMap + j + k);
						*(pMap + j + k) = *(pMap + jLast + k);
						*(pMap + jLast + k) = tmp;
					}
				} */
			}
			else
				return -1;

			//j = jMax;  // we need to re-initiate this loop and run it with adding variables
			i = -1;    // because we do have i++ in the loop
//			return -1;
		}
	}

	return nAddedVarTotal;
}

int CEquSystem::excludeVariable(VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE varValue)
{
	CVariableMapping varMapping(varIdx, varValue);
	return excludeVariables(&varMapping);
}

int CEquSystem::excludeVariables(VECTOR_ELEMENT_TYPE *pRes, VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE varValueDelta, const VECTOR_ELEMENT_TYPE *pVarMinVal, CVariableMapping *varMapping)
{
	if (varValueDelta) {
		pRes += varIdx;
		// Add the mapping for first variable
		const VECTOR_ELEMENT_TYPE varVal = *pRes += varValueDelta;
		// Since we adjusted the right parts of extra equations, by pVarMinVal[varIdx], 
		// we need to extract this value from varVal
		varMapping->addMapping(varIdx, varVal - pVarMinVal[varIdx]);

		// Add the mapping for second (with the next index) variable 
		if (getT() > 2) {
			*(pRes + 1) -= varValueDelta;
			if (!variable(varIdx))  // this variable is not in the system
				return 0;
		} else {
			// NOTE: What should be done with pVarMinVal[varIdx] ???
			varMapping->addMapping(varIdx + 1, *(pRes + 1) -= varValueDelta);
		}
	}

	// Exclude both from the system
	return excludeVariables(varMapping, 1, true, getT() > 2? 1 : 2);
}

bool CEquSystem::isSolved() const
{
	for (auto i = equNumb(); i--;)
		if (!equation(i)->solved())
			return false;

	return true;
}

void CEquSystem::resetVariables(size_t nVar)
{ 
	memset(varDefinedPtr(), 0, nVar * sizeof(*varDefinedPtr()));
	memset(getDefEquationPntr(), 0, nVar * sizeof(*getDefEquationPntr()));
}

#if 0
CForcibleCol::CForcibleCol(size_t nVar, size_t nRow) : CColOrbInfoArray(nVar)
{
	m_pForcibleColNumb = new size_t[nRow];
	memset(forcibleColNumb(), nRow, sizeof(forcibleColNumb()[0]));
}

CForcibleCol::~CForcibleCol()
{
	delete[] forcibleColNumb();
}

void CForcibleCol::addForciblyConstructedColOrbit(CColOrbit *pColOrbit, S nPart, size_t rowNumb)
{
	CColOrbInfo *pCurr, *pPrev = NULL;
	// Prepare the loop over the previously defined instances of CColOrbInfo
	const auto lenArray = GetSize();

	size_t i = 0;
	for (; i < lenArray; i++) {
		pCurr = GetElement(i);
		if (pColOrbit < pCurr->colOrb())
			break;

		pPrev = pCurr;
	}

	// Activate next instance of  CColOrbInfo:
	incNumb();
	GetElement(lenArray)->Init(pColOrbit, pPrev, i < lenArray ? pCurr : NULL);

	// Increase the number of forsibly constructed columns for current row
	forcibleColNumb()[rowNumb]++;
}

void CForcibleCol::removeForciblyConstructedColOrbit(size_t rowNumb)
{

}
#endif
#endif