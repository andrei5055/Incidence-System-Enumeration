//
//  VariableMapping.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef BIBD_Mac_VariableMapping_h
#define BIBD_Mac_VariableMapping_h

#include "DataTypes.h"

// Class supporting the mapping of variable and right parts of equations
// for trivial equations used for the enumeration of the incidence systems
class CVariableMapping : public CMapping<VECTOR_ELEMENT_TYPE>
{
public:
    CVariableMapping(size_t len) : CMapping<VECTOR_ELEMENT_TYPE>(len) {}
    virtual int resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const;
};

class CSingleLambda : public CVariableMapping
{
public:
    CSingleLambda(size_t len) : CVariableMapping(len)                 {}
    virtual int resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const;
};

class CDualLambda : public CVariableMapping
{
public:
    CDualLambda(size_t len) : CVariableMapping(len)                   {}
    virtual int resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const;
};

#endif
