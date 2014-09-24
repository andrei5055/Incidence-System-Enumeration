#include "Vector.h"

#ifdef _MSC_VER
#ifndef _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#endif
#include <crtdbg.h>
#endif
#if defined(_MSC_VER) && defined(_DEBUG)
#define new new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CVector::CVector(size_t size, CArrayOfVectorElements *pCoord)
{
	if (pCoord)	{
		InitCoordArray(pCoord);
		return;
	}

	m_pCoord = new CArrayOfVectorElements();
	getCoord()->SetSize(size);
	m_Owner = true;
}

CVector::~CVector()
{
	if (m_Owner)
		delete getCoord();
}

void CVector::InitCoordArray(CArrayOfVectorElements *pCoord, bool freeFlag)
{
	if (freeFlag && m_Owner)
		delete getCoord();

	m_pCoord = pCoord;
	m_Owner = false;
}

CInSysRowEquation::CInSysRowEquation(size_t len) : m_memShift(len), CRowEquation(3*len)
{
    m_pVarMapping = new CVariableMapping *[3];
    m_pVarMapping[0] = new CVariableMapping(len);
    m_pVarMapping[1] = new CSingleLambda(len);
    m_pVarMapping[2] = new CDualLambda(len);
}

CInSysRowEquation::~CInSysRowEquation()
{
    for (int i = 0; i < 3; i++)
        delete varMapping(i);
    
    delete [] m_pVarMapping;
}

int CInSysRowEquation::resolveTrivialEquations(const VECTOR_ELEMENT_TYPE *pRightPart, VECTOR_ELEMENT_TYPE *pResult, size_t nVar, CVariableMapping *pVariation) const
{
    const VECTOR_ELEMENT_TYPE *pResultMax = variableMaxValPntr();
    int lambda = varMapping(2)->resolveMapping(pRightPart, pResultMax, pResult, pVariation);
	if (lambda < 0)
		return -1;

    for (int i = 0; i < 2; i++)
		lambda += varMapping(i)->resolveMapping(pRightPart, pResultMax, pResult, pVariation);

    // Save minimal values
    memcpy(variableMinValPntr(), pResult, sizeof(*variableMinValPntr()) * nVar);
    return lambda;
}
