#pragma once
#include <vector>
#include "VariableMapping.h"

typedef std::vector<const VECTOR_ELEMENT_TYPE *> CSolutionStorage;

class CVector
{
 public: 
 	CVector(size_t size = 0, CArrayOfVectorElements *pElement = NULL);
	virtual ~CVector();
	inline CArrayOfVectorElements *getCoord() const		{ return m_pCoord; }
	inline void AddElement(VECTOR_ELEMENT_TYPE val)		{ getCoord()->Add(val); }
	inline VECTOR_ELEMENT_TYPE GetAt(size_t idx) const	{ return getCoord()->GetAt(idx); }
	inline size_t GetSize() const						{ return getCoord()->GetSize(); }
	inline VECTOR_ELEMENT_TYPE *GetData() const			{ return getCoord()->GetData(); }
    inline void IncreaseVectorSize(size_t len)          { getCoord()->SetSize(GetSize() + len, -1); }
protected:
	void InitCoordArray(CArrayOfVectorElements *pCoord, bool freeFlag = true);
    inline void SetAt(size_t idx, VECTOR_ELEMENT_TYPE val)	{ getCoord()->SetAt(idx, val); }
 private:
	bool m_Owner;
	CArrayOfVectorElements *m_pCoord;
};

typedef CSimpleArray<VECTOR_ELEMENT_TYPE> CRowEquation;

class CInSysRowEquation : public CRowEquation
{
public:
    CInSysRowEquation(size_t len);
    virtual ~CInSysRowEquation();
    inline void setVariableMaxVal(size_t varIdx, VECTOR_ELEMENT_TYPE maxVal)
                                                                    { *(variableMaxValPntr() + varIdx) = maxVal; }
    inline void setVariableMaxLimit(size_t varIdx, VECTOR_ELEMENT_TYPE maxVal)
                                                                    { *(variableMaxLimitPntr() + varIdx) = maxVal; }
    inline VECTOR_ELEMENT_TYPE variableMaxVal(uint varIdx) const    { return *(variableMaxValPntr() + varIdx); }
    void resetMappings()                                            { for (int i = 0; i < 3; i++) varMapping(i)->resetMapping(); }

    void addVarMapping(int idx, VECTOR_ELEMENT_TYPE to, VECTOR_ELEMENT_TYPE from)
                                                                    { varMapping(idx)->addMapping(to, from); }
    int resolveTrivialEquations(const VECTOR_ELEMENT_TYPE *pRightPart, VECTOR_ELEMENT_TYPE *pResult, size_t nVar, CVariableMapping *pVariation) const;
    inline VECTOR_ELEMENT_TYPE *variableMinValPntr() const          { return CSimpleArray<VECTOR_ELEMENT_TYPE>::elementPntr(); }
    inline VECTOR_ELEMENT_TYPE *variableMaxValPntr() const          { return CSimpleArray<VECTOR_ELEMENT_TYPE>::elementPntr() + memShift(); }
    inline VECTOR_ELEMENT_TYPE *variableMaxLimitPntr() const        { return CSimpleArray<VECTOR_ELEMENT_TYPE>::elementPntr() + (memShift()<<1); }
private:
    inline const size_t memShift() const                            { return m_memShift; }
    CVariableMapping *varMapping(int idx=0) const                   { return m_pVarMapping[idx]; }
    
    CVariableMapping **m_pVarMapping;
    const size_t m_memShift;
};