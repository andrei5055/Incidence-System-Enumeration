//
//  InSysSolver.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__InSysSolver__
#define __BIBD_Mac__InSysSolver__

#include "VariableMapping.h"

Class2Def(CRowSolution);
class CEquSystem;

#if USE_EXRA_EQUATIONS
#define CONST
#else
#define CONST	const
#endif

Class2Def(CInSysSolver) : public CVariableMapping<S>
{
public:
	CK CInSysSolver(size_t len, int t);
	CK ~CInSysSolver();
	CK S *findAllSolutionsForLambda(S *pResult, S lambdaToSplit) CONST;
    CK void initSolver(RowSolutionPntr pntr, const S *pVarMinValPntr);
    CK virtual bool isValidSolution(const S *pSol) const				{ return true; }
protected:
    CK inline RowSolutionPntr rowSolution() const						{ return m_pRowSolution; }
#if USE_EXRA_EQUATIONS
	CK inline void setEquSystem(CEquSystem *pEquSystem)					{ m_pEquSystem = pEquSystem; }
	inline void resetExtra()											{ getVarValues()->resetMapping(); resetAddedVarStack(); }
	inline CVariableMapping *getVarValues() const						{ return m_pVarValues; }
	void addVarIdxToStack(size_t nAddedVarIdx, size_t varShift);
#endif
private:
	CK bool findDiffIndex(S &lambdaToSplit, S *pResult, int *pMapIdx) CONST;
    CK int splitLambda(S &lambdaToSplit, S *pResult, int mapIdx = 0) CONST;
    CK inline void setRowSolution(RowSolutionPntr pntr)					{ m_pRowSolution = pntr; }
    CK inline void setVarMinValPntr(const S *pntr)						{ m_pVarMinVal = pntr; }
    CK inline const S *varMinValPntr() const							{ return m_pVarMinVal; }
	inline bool tDesignEnum() const										{ return m_t > 2; }
#if USE_EXRA_EQUATIONS
	inline CEquSystem *equSystem() const								{ return m_pEquSystem; }
	inline void setVarValues(CVariableMapping *pntr)					{ m_pVarValues = pntr; }
	inline void setAddedVarStack(CVariableMapping *pntr)				{ m_pAddedVarStack = pntr; }
	inline CVariableMapping *addedVarStack() const						{ return m_pAddedVarStack; }
	inline void resetAddedVarStack()									{ addedVarStack()->resetMapping(); }
	inline size_t getVarShift() const									{ return *(addedVarStack()->getLastMapping() + 1); }
	size_t getVarIdxFromStack();

	int assignVariableValues(T *pResult, size_t nVar, const T *pVariation, int lambdaToSplit) const;
	size_t undoAssignVariableValues(T *pResul, size_t nVar, size_t adj = 1, size_t n = -1, int j = 0, bool backStep = true) const;
	int removeVariableFromExtraEquations(int varIdx, T *pResult, T variation, int &lambdaToSplit) CONST;
	size_t undoAssignVariableValue(T idxVar, char *pVarDefined, T *pResult) const;
	int unassignVariable(int nVar, T *pResult, size_t adj = 1, bool adjustSecond = false, int varIdx = -1, bool excludeVar = true, bool undoAssignNeeded = true);

	CEquSystem *m_pEquSystem;
	CVariableMapping *m_pVarValues;			// Mapping of variable indices and their values
	CVariableMapping *m_pAddedVarStack;
#endif

    const S *m_pVarMinVal;
	RowSolutionPntr m_pRowSolution;
	const int m_t;
};

FClass2(CInSysSolver)::CInSysSolver(size_t len, int t) : CVariableMapping<S>(len), m_t(t)
{
#if USE_EXRA_EQUATIONS
	setVarValues(new CVariableMapping(len));
	setAddedVarStack(new CVariableMapping(len));
#endif
}

FClass2(CInSysSolver)::~CInSysSolver()
{
#if USE_EXRA_EQUATIONS
	delete getVarValues();
#endif
}

FClass2(CInSysSolver, void)::initSolver(RowSolutionPntr pntr, const S *pVarMinValPntr)
{
	setRowSolution(pntr);
	setVarMinValPntr(pVarMinValPntr);
}

FClass2(CInSysSolver, S *)::findAllSolutionsForLambda(S *pResult, S lambdaToSplit) CONST
{
	// lambdaToSplit - part of current lambda to be splited between lambda-variables
	if (!lambdaToSplit) {
		pResult = rowSolution()->copySolution(this);
		printResults(pResult, lambdaToSplit, 0, -1);
		return pResult;
	}

	bool firstPath = true;
	int varMapIdx = 0;
	do {
		varMapIdx = splitLambda(lambdaToSplit, pResult, varMapIdx);

		if (varMapIdx < 0) {
			if (firstPath) {
				// For regular case it means that there are no place for lambdaToSplit units
				// It could not be the case when we use additional equation. All lambdaToSplit units were placed somewhere,
				// but at least one extra equation is not satisfied. 

#if USE_EXRA_EQUATIONS
				if (!equSystem()->equNumb())
#endif
					return NULL;
			}
		}
		else {
			firstPath = false;
			pResult = rowSolution()->copySolution(this);
		}

	} while (findDiffIndex(lambdaToSplit, pResult, &varMapIdx));

	return pResult;
}

#if USE_EXRA_EQUATIONS
#define FLAG_SHIFT		10   // value, which allows to use negative parameters mapIdx
// and pass it to CInSysSolver::splitLambda
#endif

FClass2(CInSysSolver,  int)::splitLambda(S &lambdaToSplit, S *pResult, int mapIdx) CONST
{
	const auto *pToLast = this->getLastMapping();
#if USE_EXRA_EQUATIONS
	auto *pVarDefined = equSystem()->varDefinedPtr();
	bool changeVar = mapIdx >= 0; // Do we need change first variable here?
	const auto *pTo = getMappingPntr() + (changeVar ? mapIdx : -mapIdx - FLAG_SHIFT) - 2;
	printResults(pResult, lambdaToSplit, -1, *(pTo + 2));

#else
	const auto *pTo = this->getMappingPntr() + mapIdx - 2;
#endif

	VECTOR_ELEMENT_TYPE variation;

	while ((pTo += 2) < pToLast) {
#if USE_EXRA_EQUATIONS
		// This part should be changed
		if (*(pVarDefined + *pTo))
			continue;	// Skip variables defined from the system of extra equations

		if (changeVar) {
			// Variation cannot be bigger than any right part of any extra-equation 
			// where current variable is included
			// !!!!!!!!!!!!!!!!!!!

			// When lambdaToSplit = 0, we need only to check external equations !!!!!!!
#endif
			variation = *(pTo + 1);
			if (variation > lambdaToSplit) {
				variation = lambdaToSplit;
				lambdaToSplit = 0;
			}
			else
				lambdaToSplit -= variation;
#if USE_EXRA_EQUATIONS
		}
		else {
			changeVar = true;
			variation = 0;
		}

		if (equSystem()->equNumb()) {
			const size_t nTotal = getVarValues()->nElement();
			if (removeVariableFromExtraEquations(*pTo, pResult, variation, lambdaToSplit) < 0) {
				lambdaToSplit += unassignVariable(1, pResult, 0, 0, -1, true);
				return -1;
			}

			// Adding number of just defined variables to the stack
			const size_t nAddedVar = getVarValues()->nElement() - nTotal + (variation ? 0 : 1);
			addVarIdxToStack(nAddedVar, pTo - getMappingPntr() + 2);
		}
		else {
			*(pResult + *pTo) += variation;
			*(pResult + *pTo + 1) -= variation;
		}
#else
		*(pResult + *pTo) += variation;
		*(pResult + *pTo + 1) -= variation;
#endif

		if (lambdaToSplit)
			continue;

#if USE_EXRA_EQUATIONS
			// Nothing to split, we need to make sure that all extra equations are solved 
		if (!equSystem()->isSolved())
			return -1;
#endif
		return (int)(pTo - this->getMappingPntr());
	}

	return -1;
}

FClass2(CInSysSolver, bool)::findDiffIndex(S &lambdaToSplit, S *pResult, int *pMapIdx) CONST
{
#if USE_EXRA_EQUATIONS
	printResults(pResult, lambdaToSplit, *pMapIdx, -1);
	if (equSystem()->equNumb()) {
		int j;
		size_t varIdx, nVar = 0;
		auto *pVarValue = getVarValues()->getLastMapping();
		while (!addedVarStack()->isEmpty()) {	// while stack is not empty
			// Increase the number of variables we need to remove from stack
			nVar += j = getVarIdxFromStack();
			varIdx = *(pVarValue - (j <<= 1));

			assert(varIdx < equSystem()->nVar());

			if (*(varMinValPntr() + varIdx) < *(pResult + varIdx)) {
				// Found first possibility to decrease the value of leading variable of the group
				lambdaToSplit += unassignVariable(nVar, pResult, 0, false, varIdx);

				if (j == 2) {
					*pMapIdx = getVarShift();					// leading variable is the only one element in the group
					addedVarStack()->restoreLastMapping();		// we should restore last element of the stack
					// ... and exclude from extra equations
					// It easy  to prove that in that case our algorithm could not define the values of some new variables
					equSystem()->excludeVariables(pResult, varIdx, 0, varMinValPntr(), getVarValues());
				} else
					*pMapIdx = -((int)(getVarShift() - 2) + FLAG_SHIFT);
				// addedVarStack()->restoreLastMapping();		// we should restore last element of the stack

				return true;
			}
			else {
				nVar += 0;
			}

			pVarValue -= j;
		}   

		lambdaToSplit += unassignVariable(nVar, pResult, 0, false, -1);
		return false;
	}
#else
	lambdaToSplit = 0;
#endif

	const auto *pFirst = this->getMappingPntr();
    const auto *pCurr  = this->getLastMapping();
    S *pCurrVar;
    
    // Counts total amount of splitted units AFTER pTo and define if it's possible to split one more unit in the same area

    bool spaceNotFound = true;
    while ((pCurr -= 2) >= pFirst) {
		int currMinVal, varIdx, val;
        if (!(val = *(pCurrVar = pResult + (varIdx = *pCurr))) ||
            val == (currMinVal = *(varMinValPntr() + varIdx))) {
            // If val = 0 || val == curMinVal, we don't need to adjust lambdaToSlit or *pCurrVar
            // But, definitely, we found space to one unit
            spaceNotFound = false;
            continue;
        }
        
        if (!spaceNotFound) {
            *pCurrVar -= 1;
            *(pCurrVar + 1) += 1;
            *pMapIdx = (int)(pCurr - pFirst + 2);
			lambdaToSplit++;
            return true;
        }
        
        spaceNotFound = val == currMinVal + *(pCurr+1);
        lambdaToSplit += (val -= currMinVal);
        *pCurrVar = currMinVal;
        *(pCurrVar + 1) += val;
    }

    return false;
}


#endif /* defined(__BIBD_Mac__InSysSolver__) */
