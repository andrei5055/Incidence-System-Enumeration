//
//  InSysSolver.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "InSysSolver.h"
#include "RowSolution.h"

void CInSysSolver::initSolver(CRowSolution *pntr, const VECTOR_ELEMENT_TYPE *pVarMinValPntr)
{
    setRowSolution(pntr);
    setVarMinValPntr(pVarMinValPntr);
}

int CInSysSolver::splitLambda(int lambdaToSplit, VECTOR_ELEMENT_TYPE *pResult, int mapIdx) const
{
    const auto *pTo = getMappingPntr();
    const auto *pToLast = pTo + getMapPosition();
    pTo += mapIdx;
    while (pTo < pToLast) {
        const auto variation = *(pTo + 1);
        if (variation >= lambdaToSplit) {
            *(pResult + *pTo) += lambdaToSplit;
            *(pResult + *pTo + 1) -= lambdaToSplit;
            return  (int)(pTo - getMappingPntr());
        }
        
        *(pResult + *pTo) += variation;
        *(pResult + *pTo + 1) -= variation;
        lambdaToSplit -= variation;
        pTo += 2;
    }
    
    return -1;
}

int CInSysSolver::findDiffIndex(VECTOR_ELEMENT_TYPE *pResult, int *pMapIdx) const
{
    const auto *pFirst = getMappingPntr();
    const auto *pCurr  = getLastMappping();
    VECTOR_ELEMENT_TYPE *pCurrVar;
    
    // Counts total amount of splitted units AFTER pTo and define if it's possible to split one more unit in the same area
    int lambdaToSplit = 0;
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
            return lambdaToSplit + 1;
        }
        
        spaceNotFound = val == currMinVal + *(pCurr+1);
        lambdaToSplit += (val -= currMinVal);
        *pCurrVar = currMinVal;
        *(pCurrVar + 1) += val;
    }
    
    return -1;
}

VECTOR_ELEMENT_TYPE *CInSysSolver::findAllSolutionsForLambda(VECTOR_ELEMENT_TYPE *pResult, int lambdaToSplit) const
{
    // lambdaToSplit - part of current lambda to be splited between lambda-variables
    if (!lambdaToSplit)
        return rowSolution()->copySolution(this);
    
    // Find first solution for current lambdaToSplit
    int varMapIdx = splitLambda(lambdaToSplit, pResult);
    if (varMapIdx < 0)
        return NULL;
    
    // First (lexicographicaly smallest) solution is constructed
    while (true) {
        pResult = rowSolution()->copySolution(this);
        lambdaToSplit = findDiffIndex(pResult, &varMapIdx);
        if (lambdaToSplit < 0)
            break;
        
        varMapIdx = splitLambda(lambdaToSplit, pResult, varMapIdx);
    }

    return pResult;
}
