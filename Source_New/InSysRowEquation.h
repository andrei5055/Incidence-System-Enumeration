#pragma once
#include "Vector.h"

class CEquSystem;

Class1Def(CInSysRowEquation) : public Class1(CSimpleArray)
{
public:
	CK CInSysRowEquation(size_t len, bool tDesignEnum) : CSimpleArray<S>(3 * len) {
		InitRowEquation(len, tDesignEnum, false);
	}
	CK CInSysRowEquation()                                      	{}
	CK virtual ~CInSysRowEquation();
	void InitRowEquation(size_t len, bool tDesignEnum, bool initArray = true);
	CK inline void setVariableMaxVal(S varIdx, S maxVal)			{ *(variableMaxValPntr() + varIdx) = maxVal; }
	CK inline void setVariableMaxLimit(S varIdx, S maxVal)			{ *(variableMaxLimitPntr() + varIdx) = maxVal; }
	inline S variableMaxVal(uint varIdx) const						{ return *(variableMaxValPntr() + varIdx); }
	CK inline void resetMappings()                                  { for (int i = 0; i <= t_dual; i++) varMapping(i)->resetMapping(); }

	CK void addVarMapping(int idx, S to, S from)					{ varMapping(idx)->addMapping(to, from); }
	CK int resolveTrivialEquations(const S *pRightPart, S *pResult, size_t nVar, CVariableMapping<S> *pVariation) const;
	CK inline S *variableMinValPntr() const							{ return CSimpleArray<S>::elementPntr(); }
	CK inline S *variableMaxValPntr() const							{ return CSimpleArray<S>::elementPntr() + memShift(); }
	CK inline S *variableMaxLimitPntr() const						{ return CSimpleArray<S>::elementPntr() + (memShift() << 1); }
	CK inline void setEquSystem(CEquSystem *pEquSystem)				{ m_pEquSystem = pEquSystem; }
#if USE_EXRA_EQUATIONS
	size_t excludeVariables(CVariableMapping<S> *pVarValue, int adj = 1) const;
	CVariableMapping *addDefinedVariables(CVariableMapping<S> *pResMapping, S *pResult, CVariableMapping<T> *pVariation) const;
	int resolveExtraMapping(const CVariableMapping<S> *pVarValue, const S *pMaxVal, S *pResult, CVariableMapping<S> *pRemove, CVariableMapping<S> *pOut, size_t shift) const;
	CVariableMapping<S> *AdjustExtraEquationRightParts(CVariableMapping<S> *pVarList = NULL, bool addVar = false);
#endif
private:
	CK inline auto memShift() const									{ return m_memShift; }
	CK auto varMapping(int idx = 0) const							{ return m_pVarMapping[idx]; }
	inline auto *equSystem() const									{ return m_pEquSystem; }
	inline bool t_DesigneEnum() const								{ return m_bTDesignEnum; }

	CVariableMapping<S> **m_pVarMapping;
	CEquSystem *m_pEquSystem;
	size_t m_memShift;
	bool m_bTDesignEnum;
};

FClass1(CInSysRowEquation, void)::InitRowEquation(size_t len, bool tDesignEnum, bool initArray) {
	m_memShift = len;
	m_bTDesignEnum = tDesignEnum;
	m_pVarMapping = new CVariableMapping<S> *[3];
	m_pVarMapping[t_singleNoLambda] = new CVariableMapping<S>(len);
	m_pVarMapping[t_singleLambda] = new CSingleLambda<S>(len);
	m_pVarMapping[t_dual] = new CDualLambda<S>(len);
	if (initArray)
		Init(3 * len, new S[3 * len]);
}

FClass1(CInSysRowEquation)::~CInSysRowEquation()
{
	for (int i = 0; i <= t_dual; i++)
		delete varMapping(i);

	delete[] m_pVarMapping;
}


FClass1(CInSysRowEquation, int)::resolveTrivialEquations(const S *pRightPart, S *pResult, size_t nVar, CVariableMapping<S> *pVariation) const
{
	const auto pResultMax = variableMaxValPntr();
	int lambda = varMapping(t_dual)->resolveMapping(pRightPart, pResultMax, pResult, pVariation);
	if (lambda < 0)
		return -1;

	for (int i = 0; i < t_dual; i++)
		lambda += varMapping(i)->resolveMapping(pRightPart, pResultMax, pResult, NULL);

	// Save minimal values
	memcpy(variableMinValPntr(), pResult, sizeof(S) * nVar);
	return lambda;
}

