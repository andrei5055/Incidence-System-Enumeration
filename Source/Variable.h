#pragma once

#include "DataTypes.h"

template<typename TYPE, size_t IDX>
class CLinks
{
public:
	inline TYPE next(int idx) const						{ return m_pNext[idx]; }
	inline TYPE prev(int idx) const						{ return m_pPrev[idx]; }
	inline void setNext(TYPE pntr, int idx)				{ m_pNext[idx] = pntr; }
	inline void setPrev(TYPE pntr, int idx)				{ m_pPrev[idx] = pntr; }
private:
	TYPE m_pNext[IDX];
	TYPE m_pPrev[IDX];
};

typedef enum {
	t_sameEqu,
	t_sameVar
} t_sameType;

class CVariable : public CLinks<CVariable *, 2>
{
public:
	CVariable(VECTOR_ELEMENT_TYPE equIdx, VECTOR_ELEMENT_TYPE varIdx, CVariable *pntr = NULL);
#define nextVarSameEqu()			next(t_sameEqu) // Pointer to the entrance of the next variable into same equation
#define prevVarSameEqu()			prev(t_sameEqu)	// Pointer to the entrance of the prev variable into same equation
#define sameVarNextEqu()			next(t_sameVar)	// Pointer to the entrance of the same variable into next equation
#define sameVarPrevEqu()			prev(t_sameVar)
#define setNextVarSameEqu(pntr)		setNext(pntr, t_sameEqu)
#define setPrevVarSameEqu(pntr)		setPrev(pntr, t_sameEqu)
#define setSameVarNextEqu(pntr)		setNext(pntr, t_sameVar)
#define setSameVarPrevEqu(pntr)		setPrev(pntr, t_sameVar)
	inline VECTOR_ELEMENT_TYPE varIndex() const			{ return m_nVarIdx; }
	inline VECTOR_ELEMENT_TYPE equIdx()	const			{ return m_nEquIdx; }
	void linkVariable(CVariable *pVar);
	void excludeVariableFromEquation();
	void restoreVariableInEquation();
private:

	const VECTOR_ELEMENT_TYPE m_nEquIdx;
	const VECTOR_ELEMENT_TYPE m_nVarIdx;
};
