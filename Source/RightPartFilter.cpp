//
//  RightPartFilter.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/23/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "RightPartFilter.h"

VECTOR_ELEMENT_TYPE *CRightPartFilter::getRightPart(const VECTOR_ELEMENT_TYPE *pRightSide, const VECTOR_ELEMENT_TYPE *pVarMaxVal, size_t lenRightPart, VECTOR_ELEMENT_TYPE *pRighPartMem) const
{
    if (!numFilter())
        return (VECTOR_ELEMENT_TYPE *)pRightSide;
    
    int j, n;
    const FilterData *pFilterData = getFilterData();
    for (uint i = j = n = 0; i < numFilter(); i++, pFilterData++) {
        // Copy all nonfiltered right parts
        const VECTOR_ELEMENT_TYPE len = pFilterData->equIdx() - j;
        if (len) {
            memcpy(pRighPartMem + n, pRightSide + j, len * sizeof(pRighPartMem[0]));
            n += len;
            j += len;
        }
        
        const VECTOR_ELEMENT_TYPE rightSideVal = *(pRightSide+j);
		const int val = pFilterData->val();
        const int diff = val - rightSideVal;
        if (pFilterData->varIdx() == -1) {
            if (diff)
                return NULL;
            
            j++;        // we need to skip current right part
            continue;
        }
        
        if (!val) {
            if (*(pVarMaxVal + pFilterData->varIdx()) < rightSideVal)
                return NULL;        // we cannot use all units of current right part for corresponding variable

            continue;
        }

        if (diff > 0)
            return NULL;            // not enough units in the current right part
                
        j++;
        *(pRighPartMem + n++) = -diff;
    }
    
    if ((lenRightPart -= j) > 0)
        memcpy(pRighPartMem + n, pRightSide + j, lenRightPart * sizeof(pRighPartMem[0]));
        
    return pRighPartMem;
}