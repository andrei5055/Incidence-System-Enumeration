//
//  VariableMapping.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "VariableMapping.h"

int CVariableMapping::resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const
{
    const VECTOR_ELEMENT_TYPE *pTo = getMappingPntr() + getMapPosition();
    while (pTo > getMappingPntr()) {
        pTo -= 2;
        *(pResult + *pTo) = *(pRightPart + *(pTo + 1));
    }
    
    return 0;
}

int CSingleLambda::resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const
{
    int lambda = 0;
    const VECTOR_ELEMENT_TYPE *pTo = getMappingPntr() + getMapPosition();
    while (pTo > getMappingPntr()) {
        pTo -= 2;
        lambda += *(pResult + *pTo) = *(pRightPart + *(pTo + 1));
    }
    
    return lambda;
}

int CDualLambda::resolveMapping(const VECTOR_ELEMENT_TYPE *pRightPart, const VECTOR_ELEMENT_TYPE *pResultMax, VECTOR_ELEMENT_TYPE *pResult, CVariableMapping *pVariation) const
{
    // In this function we also construct pVariation - the list of variables
    // with their limits to be added to their minimal values
    int lambda = 0;
    const VECTOR_ELEMENT_TYPE *pTo = getMappingPntr() + getMapPosition();
    while (pTo > getMappingPntr()) {
        pTo -= 2;
        // Index of the current "lambda" variable
        const VECTOR_ELEMENT_TYPE idx = *pTo - 1;
        // Define variation as the maximal value for current "lambda" variable
        VECTOR_ELEMENT_TYPE variation = *(pResultMax + idx);
        // Take value which should be splited between two variables
        VECTOR_ELEMENT_TYPE val = *(pRightPart + *(pTo+1));
        if (variation > val)
            variation = val;
        
        if (val > *(pResultMax + idx + 1)) {
            // This value is bigger then maximum for "non-lambda" variable
            // Set the value of "non-lambda" variable to its maximum
            // and set value of "lambda" variable to remaining part of val
            // (this vill be the minimal value of "lambda" variable)
			val -= (*(pResult + idx + 1) = *(pResultMax + idx + 1));
            lambda += *(pResult + idx) = val;
            if (variation == val)
				continue;

            if (variation < val)
				return -1;

            variation -= val;
        } else {
            *(pResult + idx + 1) = val;
            *(pResult + idx) = 0;
            if (!val)       // Right part is 0
                continue;   // both variables will be 0's
        }

        pVariation->addMapping(idx, variation);
    }
    
    return lambda;
}

