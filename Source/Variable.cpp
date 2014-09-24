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
	setNextVarSameEqu(pVar);
	pVar->setPrevVarSameEqu(this);
}
