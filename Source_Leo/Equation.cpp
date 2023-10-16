#include "Equation.h"
#include "VariableMapping.h"

#if USE_EXRA_EQUATIONS

CIndexArray::CIndexArray()
{
	// Allocate the array to store the indices of the variables defined last from that equation
	m_numDefinedVarMax = 1;
	m_idxVarLast = NULL;
	AlocateArrayForIndices();
	resetIndexArray();
}

CIndexArray::~CIndexArray()
{
	delete[] m_idxVarLast;
}

void CIndexArray::addVarIndex(VECTOR_ELEMENT_TYPE varIdx)
{
	if (m_numDefinedVar == m_numDefinedVarMax)
		AlocateArrayForIndices();

	m_idxVarLast[m_numDefinedVar++] = varIdx;
}

bool CIndexArray::isVarIndexInIndexArray(VECTOR_ELEMENT_TYPE idx) const
{
	for (size_t i = m_numDefinedVar; i--;)
		if (m_idxVarLast[i] == idx)
			return true;

	return false;
}

void CIndexArray::AlocateArrayForIndices()
{
	VECTOR_ELEMENT_TYPE * pIdxVarTmp = m_idxVarLast;
	m_idxVarLast = new VECTOR_ELEMENT_TYPE[m_numDefinedVarMax <<= 1];
	if (!pIdxVarTmp)
		return;

	memcpy(m_idxVarLast, pIdxVarTmp, (m_numDefinedVarMax >> 1) * sizeof(m_idxVarLast[0]));
	delete[] pIdxVarTmp;
}

void CEquation::initEquation(CVariable *pntr, VECTOR_ELEMENT_TYPE nVar, VECTOR_ELEMENT_TYPE rPart)
{
	setFirstVariable(pntr);
	setNumbVar(nVar);
	setRightPart(rPart);
	setSolved(false);
	setNumDefinedVar(0);
	assignPrintResultInfo(myID, this);
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

int CEquation::removeVariable(const VECTOR_ELEMENT_TYPE valVar, char *pVarDefined, CVariableMapping *pVarValue, CEquation **pDefVar, int adj)
{
	if (adj > 0) {
		if (rightPart() < valVar) {
			// Right part of current equation is less than the value of current variable
			// ???				if (pVarValueOut)
			// ???					pVarValueOut->resetMapping();

			return -2;
		}
		adjustNumbVar(adj);
	}

	const size_t rightPartValue = adjustRightPart(valVar, adj >= 0? 1 : -1);
	if (numbVar() != 1 && rightPartValue)
		return -1;

	const auto *pSavedLast = pVarValue->getLastMapping();
	// Current equation either
	//    a) contains exactly two variables OR
	//    b) has right part equal to 0
	// In both cases the remaining variables are defined 
	// and they values are equal to current right part
	setNumDefinedVar(numbVar());
	resetIndexArray();

	CVariable *pCurrVar = firstVariable();
	for (int i = 0; i < numbVar(); i++) {
		while (*(pVarDefined + (pCurrVar = pCurrVar->nextVarSameEqu())->varIndex())) {
			if (pCurrVar == firstVariable()) {			
				pVarValue->removeLastMapping(i);// Roll back just added variables
				return -2;
			}
		}

		const VECTOR_ELEMENT_TYPE varIdx = pCurrVar->varIndex();
		switch (pVarValue->addVarValue(varIdx, rightPartValue, pVarValue->mapBoundaries(), false)) {
			case 2: *(pDefVar + varIdx) = this;	// Memorize the equation where from the value of variable was defined
					continue;
												// Variable was already added with the same varValue
			case 1:	addVarIndex(varIdx);		// and it was added into CURRENT set of defined variables
			case 0:	continue;					//     it was added in a NEXT set of defined variable	
		}

		pVarValue->removeLastMapping(i);		// Roll back just added variables
		adjustNumbVar(-1);						// Restore the number of variables
		adjustRightPart(valVar, -1);			// and the right part value
		return -2;
	}

	// Regardless of what we have (a) or (b), for correct reconstruction of the equation 
	// we have to set the number of current equation variable to 1
	setNumbVar(1);
	setSolved();
	return 1;
}


void CEquation::addVariable(VECTOR_ELEMENT_TYPE valVar)
{
	if (!solved()) {
		if (valVar)
			adjustRightPart(-valVar);

		adjustNumbVar(-1);
	}
	else {
		setSolved(false);
		setNumDefinedVar(0);
	}
}

void CEquation::addVarIndex(VECTOR_ELEMENT_TYPE varIdx)
{
	if (!m_pIndexArray)
		m_pIndexArray = new CIndexArray();

	m_pIndexArray->addVarIndex(varIdx);
}

void CEquation::resetIndexArray() const
{
	if (m_pIndexArray)
		m_pIndexArray->resetIndexArray();
}

#endif
