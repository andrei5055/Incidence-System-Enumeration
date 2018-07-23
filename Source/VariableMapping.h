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

typedef enum {
	t_singleNoLambda,
	t_singleLambda,
	t_dual
} t_varType;

class CVariable;

// Class supporting the mapping of variable and right parts of equations
// for trivial equations used for the enumeration of the incidence systems
template <class T>
class CVariableMapping : public CMapping<T>
{
public:
	CK CVariableMapping(size_t len = 1) : CMapping<T>(len)	{ resetMapping(); }
	CK CVariableMapping(T varIdx, T varValue, size_t len = 1) : CMapping<T>(varIdx, varValue, len) {}
    CK virtual int resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping *pVariation) const;
	bool findVariable(T varID, uint *pIdx);
#if USE_EXRA_EQUATIONS
	virtual CVariableMapping *addDefinedVariables(CVariableMapping *pResMapping, const T *pResult, CVariable **pVar, const CVariableMapping *pVariation = NULL) const;
	virtual CVariableMapping *addVariablMinValue(const T *pMinValues, CVariable **pVar, const T *pMaxValues, CVariableMapping **ppResMaxList, CVariableMapping *pResList = NULL) const;
	int addVarValue(T varID, T varValue, const T * const *pLastMapping = NULL, bool from = true);
	void setMapBoundaries(size_t nVar)						{ m_nMapBoundary[0] = (m_nMapBoundary[1] = getLastMapping()) - (nVar << 1); }
	const T *mapBoundary(int i) const						{ return m_nMapBoundary[i]; }
	const T *const *mapBoundaries() const					{ return m_nMapBoundary; }
private:
	const T *m_nMapBoundary[2];  // the mapping which we will use to start the search of newly added variables
#endif
};

template <class T>
class CSingleLambda : public CVariableMapping<T>
{
public:
    CK CSingleLambda(size_t len) : CVariableMapping<T>(len)      {}
    CK virtual int resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping *pVariation) const;
};

template <class T>
class CDualLambda : public CVariableMapping<T>
{
public:
    CK CDualLambda(size_t len) : CVariableMapping<T>(len)        {}
    CK virtual int resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping *pVariation) const;
#if USE_EXRA_EQUATIONS
	virtual CVariableMapping *addDefinedVariables(CVariableMapping *pResMapping, const T *pResult, CVariable **pVar, const CVariableMapping *pVariation = NULL) const;
#endif
};

template<class T>
int CVariableMapping<T>::resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping<T> *pVariation) const
{
	const auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		pTo -= 2;
		*(pResult + *pTo) = *(pRightPart + *(pTo + 1));
	}

	return 0;
}

template<class T>
int CSingleLambda<T>::resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping<T> *pVariation) const
{
	int lambda = 0;
	const auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		pTo -= 2;
		lambda += *(pResult + *pTo) = *(pRightPart + *(pTo + 1));
	}

	return lambda;
}

template<class T>
int CDualLambda<T>::resolveMapping(const T *pRightPart, const T *pResultMax, T *pResult, CVariableMapping<T> *pVariation) const
{
	// In this function we also construct pVariation - the list of variables
	// with their limits to be added to their minimal values
	int lambda = 0;
	const auto *pTo = getLastMapping();
	while (pTo > getMappingPntr()) {
		pTo -= 2;
		// Index of the current "lambda" variable
		const auto idx = *pTo - 1;
		// Define variation as the maximal value for current "lambda" variable
		auto variation = *(pResultMax + idx);
		// Take value which should be splited between two variables
		auto val = *(pRightPart + *(pTo + 1));
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
		}
		else {
			*(pResult + idx + 1) = val;
			*(pResult + idx) = 0;
			if (!val)       // Right part is 0
				continue;   // both variables will be 0's
		}

		pVariation->addMapping(idx, variation);
	}

	return lambda;
}

#endif
