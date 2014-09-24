#include "Equation.h"

void CEquation::initEquation(CVariable *pntr, VECTOR_ELEMENT_TYPE nVar, VECTOR_ELEMENT_TYPE rPart)
{
	setFirstVariable(pntr);
	setNumbVar(nVar);
	setRightPart(rPart);
	setSolved(false);
}

void CEquation::releaseEquation()
{
	auto pFirstVar = firstVariable();
	auto pVarTmp = pFirstVar;
	do {
		auto pNext = pVarTmp->nextVarSameEqu();
		delete pVarTmp;
		pVarTmp = pNext;
	} while (pVarTmp != pFirstVar);
}