#include "DataTypes.h"

class CVariable
{
public:
	CVariable(VECTOR_ELEMENT_TYPE equIdx, VECTOR_ELEMENT_TYPE varIdx, CVariable *pntr = NULL);
	inline CVariable *nextVarSameEqu() const			{ return m_pNextVarSameEqu; }
	inline CVariable *prevVarSameEqu() const			{ return m_pPrevVarSameEqu; }
	inline CVariable *sameVarNextEqu() const			{ return m_pSameVarNextEqu; }
	inline void setNextVarSameEqu(CVariable *pntr)		{ m_pNextVarSameEqu = pntr; }
	inline void setPrevVarSameEqu(CVariable *pntr)		{ m_pPrevVarSameEqu = pntr; }
	inline void setSameVarNextEqu(CVariable *pntr)		{ m_pSameVarNextEqu = pntr; }
	inline const VECTOR_ELEMENT_TYPE varIndex() const	{ return m_nVarIdx; }
	inline VECTOR_ELEMENT_TYPE equIdx()	const			{ return m_nEquIdx; }
	void linkVariable(CVariable *pVar);
private:

	CVariable *m_pNextVarSameEqu; // Pointer to the entrance of the next variable into same equestion
	CVariable *m_pPrevVarSameEqu; // Pointer to the entrance of the prev variable into same equestion
	CVariable *m_pSameVarNextEqu; // Pointer to the entrance of the same variable into next equestion
	const VECTOR_ELEMENT_TYPE m_nEquIdx;
	const VECTOR_ELEMENT_TYPE m_nVarIdx;
};
