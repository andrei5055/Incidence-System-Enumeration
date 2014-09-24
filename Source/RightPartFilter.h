//
//  RightPartFilter.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/23/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__RightPartFilter__
#define __BIBD_Mac__RightPartFilter__

#include "VariableMapping.h"

typedef struct {
    int m_equIdx;
    size_t m_varIdx;
    VECTOR_ELEMENT_TYPE m_val;
    inline void setFilterData(int eIdx, size_t vIdx, VECTOR_ELEMENT_TYPE v)
                                                { m_equIdx = eIdx; m_varIdx = vIdx; m_val = v; }
    inline int equIdx() const                   { return m_equIdx; }
    inline size_t varIdx() const                { return m_varIdx; }
    inline VECTOR_ELEMENT_TYPE val() const      { return m_val; }
} FilterData;

class CRightPartFilter
{
public:
    CRightPartFilter(size_t len)                { m_pFilterData = new FilterData [len]; }
    ~CRightPartFilter()                         { delete [] getFilterData(); }
    void addFilter(int equIdx, VECTOR_ELEMENT_TYPE val, size_t varIdx = -1)
                                                { getFilterData(m_nFilter++)->setFilterData(equIdx, varIdx, val); }
    void reset()                                { m_nFilter = 0; }
    VECTOR_ELEMENT_TYPE *getRightPart(const VECTOR_ELEMENT_TYPE *pRightSide, const VECTOR_ELEMENT_TYPE *pVarMaxVal, size_t lenRightPart, VECTOR_ELEMENT_TYPE *pRighPartMem) const;
protected:
private:
    inline uint numFilter() const                       { return m_nFilter; }
    inline FilterData *getFilterData(int i = 0) const   { return m_pFilterData + i; }

    FilterData *m_pFilterData;
    uint m_nFilter;
};

#endif /* defined(__BIBD_Mac__RightPartFilter__) */
