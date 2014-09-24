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

class CRowSolution;

class CInSysSolver : public CVariableMapping
{
public:
    CInSysSolver(size_t len) : CVariableMapping(len)                    {}
    VECTOR_ELEMENT_TYPE *findAllSolutionsForLambda(VECTOR_ELEMENT_TYPE *pResult, int lambdaToSplit) const;
    void initSolver(CRowSolution *pntr, const VECTOR_ELEMENT_TYPE *pVarMinValPntr);
    virtual bool isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const { return true; }
private:
    int findDiffIndex(VECTOR_ELEMENT_TYPE *pResult, int *pMapIdx) const;
    int splitLambda(int lambdaToSplit, VECTOR_ELEMENT_TYPE *pResult, int mapIdx = 0) const;
    inline void setRowSolution(CRowSolution *pntr)                      { m_pRowSolution = pntr; }
    inline CRowSolution *rowSolution() const                            { return m_pRowSolution; }
    inline void setVarMinValPntr(const VECTOR_ELEMENT_TYPE *pntr)       { m_pVarMinVal = pntr; }
    inline const VECTOR_ELEMENT_TYPE *varMinValPntr() const             { return m_pVarMinVal; }
    
    const VECTOR_ELEMENT_TYPE *m_pVarMinVal;
    CRowSolution *m_pRowSolution;
};

#endif /* defined(__BIBD_Mac__InSysSolver__) */
