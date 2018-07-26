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
    CK inline void setFilterData(int eIdx, size_t vIdx, VECTOR_ELEMENT_TYPE v)
                                                { m_equIdx = eIdx; m_varIdx = vIdx; m_val = v; }
    CK inline int equIdx() const                { return m_equIdx; }
    CK inline size_t varIdx() const             { return m_varIdx; }
    CK inline VECTOR_ELEMENT_TYPE val() const					{ return m_val; }
} FilterData;

template<class T>
class CRightPartFilter
{
public:
    CK CRightPartFilter(size_t len)             { m_pFilterData = new FilterData [len]; }
    CK ~CRightPartFilter()                      { delete [] getFilterData(); }
    CK inline void addFilter(int equIdx, T val, size_t varIdx = (size_t)-1)
                                                { getFilterData(m_nFilter++)->setFilterData(equIdx, varIdx, val); }
    CK inline void reset()                      { m_nFilter = 0; }
    CK T *getRightPart(const T *pRightSide, const T *pVarMaxVal, 
		size_t lenRightPart, T *pRighPartMem) const {
		if (!numFilter())
			return (T *)pRightSide;

		int j, n;
		const FilterData *pFilterData = getFilterData();
		for (uint i = j = n = 0; i < numFilter(); i++, pFilterData++) {
			// Copy all nonfiltered right parts
			const auto len = pFilterData->equIdx() - j;
			if (len) {
				memcpy(pRighPartMem + n, pRightSide + j, len * sizeof(pRighPartMem[0]));
				n += len;
				j += len;
			}

			const auto rightSideVal = *(pRightSide + j);
			const int val = pFilterData->val();
			const int diff = val - rightSideVal;
			if (pFilterData->varIdx() == (size_t)-1) {
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

protected:
private:
    CK inline uint numFilter() const            { return m_nFilter; }
    CK inline FilterData *getFilterData(int i = 0) const	{ return m_pFilterData + i; }

    FilterData *m_pFilterData;
    uint m_nFilter;
};

#endif /* defined(__BIBD_Mac__RightPartFilter__) */
