#include "Variable.h"

CVariable::CVariable(VECTOR_ELEMENT_TYPE equIdx, VECTOR_ELEMENT_TYPE varIdx, CVariable *pntr) : m_nEquIdx(equIdx), m_nVarIdx(varIdx)
{
	if (pntr)
		linkVariable(pntr);
	else
		setNextVarSameEqu(this);

	setSameVarNextEqu(NULL);
}

void CVariable::linkVariable(CVariable *pVar)
{
	// NOTE: We use oposite direction for linked list of variable in the equation.
	// For instance, for the equation
	//
	//  x2 + x5 + x7 + x9 = XXX
	//
	// for x5: x2 will be "next" and x7 will be "prev"
	// for x2: x9 will be "next" and x5 will be "prev"

	setNextVarSameEqu(pVar);
	pVar->setPrevVarSameEqu(this);
}

void CVariable::excludeVariableFromEquation()
{ 
	// Here we are using m_pPrev[1] to keep the pointer which will be used 
	// when we will need to add current variable to the equation
	// setPrev(prevVarSameEqu(), t_sameVar);

	prevVarSameEqu()->linkVariable(nextVarSameEqu()); 
}

void CVariable::restoreVariableInEquation()
{
//	auto pntr = prev(t_sameVar);
	auto pPrev = prevVarSameEqu();
	auto pNext = nextVarSameEqu();
	if (pPrev != pNext) {
		pPrev->setNextVarSameEqu(this);
		pNext->setPrevVarSameEqu(this);
	}
}


