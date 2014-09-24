#include "Variable.h"

class CEquation
{
public:
	CEquation()												{}
	void initEquation(CVariable *pntr, VECTOR_ELEMENT_TYPE nVar, VECTOR_ELEMENT_TYPE rPart);
	void releaseEquation();
	inline void setSolved(bool val = true)					{ m_bSolved = val; }
	inline bool solved() const								{ return m_bSolved; }
	inline CVariable *firstVariable() const					{ return m_pVar; }
	inline VECTOR_ELEMENT_TYPE numbVar() const				{ return m_nVar; }
	inline VECTOR_ELEMENT_TYPE rightPart() const			{ return m_nRightPart; }
	inline void adjustRightPart(VECTOR_ELEMENT_TYPE val)    { setRightPart(rightPart() - val); }
private: 
	inline void setNumbVar(VECTOR_ELEMENT_TYPE val)			{ m_nVar = val; }
	inline void setFirstVariable(CVariable *pVar)			{ m_pVar = pVar; }
	inline void setRightPart(VECTOR_ELEMENT_TYPE val)		{ m_nRightPart = val; }

	CVariable *m_pVar;
	VECTOR_ELEMENT_TYPE m_nVar;
	VECTOR_ELEMENT_TYPE m_nRightPart;
	bool m_bSolved;
};