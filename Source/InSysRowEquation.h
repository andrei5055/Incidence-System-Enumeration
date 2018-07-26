#pragma once
#include "Vector.h"

//#define CRowEquation CSimpleArray

class CEquSystem;

template <class T>
class CInSysRowEquation : public CSimpleArray<T>
{
public:
	CK CInSysRowEquation(size_t len, bool tDesignEnum);
	CK virtual ~CInSysRowEquation();
	CK inline void setVariableMaxVal(size_t varIdx, T maxVal)		{ *(variableMaxValPntr() + varIdx) = maxVal; }
	CK inline void setVariableMaxLimit(size_t varIdx, T maxVal)		{ *(variableMaxLimitPntr() + varIdx) = maxVal; }
	inline T variableMaxVal(uint varIdx) const						{ return *(variableMaxValPntr() + varIdx); }
	CK inline void resetMappings()                                  { for (int i = 0; i <= t_dual; i++) varMapping(i)->resetMapping(); }

	CK void addVarMapping(int idx, T to, T from)					{ varMapping(idx)->addMapping(to, from); }
	CK int resolveTrivialEquations(const T *pRightPart, T *pResult, size_t nVar, CVariableMapping<T> *pVariation) const;
	CK inline T *variableMinValPntr() const						{ return CSimpleArray<T>::elementPntr(); }
	CK inline T *variableMaxValPntr() const							{ return CSimpleArray<T>::elementPntr() + memShift(); }
	CK inline T *variableMaxLimitPntr() const						{ return CSimpleArray<T>::elementPntr() + (memShift() << 1); }
	CK inline void setEquSystem(CEquSystem *pEquSystem)				{ m_pEquSystem = pEquSystem; }
#if USE_EXRA_EQUATIONS
	size_t excludeVariables(CVariableMapping<T> *pVarValue, int adj = 1) const;
	CVariableMapping *addDefinedVariables(CVariableMapping<T> *pResMapping, T *pResult, CVariableMapping<T> *pVariation) const;
	int resolveExtraMapping(const CVariableMapping<T> *pVarValue, const T *pMaxVal, T *pResult, CVariableMapping<T> *pRemove, CVariableMapping<T> *pOut, size_t shift) const;
	CVariableMapping<T> *AdjustExtraEquationRightParts(CVariableMapping<T> *pVarList = NULL, bool addVar = false);
#endif
private:
	CK inline size_t memShift() const								{ return m_memShift; }
	CK CVariableMapping<T> *varMapping(int idx = 0) const			{ return m_pVarMapping[idx]; }
	inline CEquSystem *equSystem() const							{ return m_pEquSystem; }
	inline bool t_DesigneEnum() const								{ return m_bTDesignEnum; }

	CVariableMapping<T> **m_pVarMapping;
	CEquSystem *m_pEquSystem;
	const size_t m_memShift;
	const bool m_bTDesignEnum;
};

template<class T>
CInSysRowEquation<T>::CInSysRowEquation(size_t len, bool tDesignEnum) : m_memShift(len), m_bTDesignEnum(tDesignEnum), CSimpleArray<T>(3 * len)
{
	m_pVarMapping = new CVariableMapping<T> *[3];
	m_pVarMapping[t_singleNoLambda] = new CVariableMapping<T>(len);
	m_pVarMapping[t_singleLambda] = new CSingleLambda<T>(len);
	m_pVarMapping[t_dual] = new CDualLambda<T>(len);
}

template<class T>
CInSysRowEquation<T>::~CInSysRowEquation()
{
	for (int i = 0; i <= t_dual; i++)
		delete varMapping(i);

	delete[] m_pVarMapping;
}

template<class T>
int CInSysRowEquation<T>::resolveTrivialEquations(const T *pRightPart, T *pResult, size_t nVar, CVariableMapping<T> *pVariation) const
{
	const auto pResultMax = variableMaxValPntr();
	int lambda = varMapping(t_dual)->resolveMapping(pRightPart, pResultMax, pResult, pVariation);
	if (lambda < 0)
		return -1;

	for (int i = 0; i < t_dual; i++)
		lambda += varMapping(i)->resolveMapping(pRightPart, pResultMax, pResult, NULL);

	// Save minimal values
	memcpy(variableMinValPntr(), pResult, sizeof(*variableMinValPntr()) * nVar);
	return lambda;
}

