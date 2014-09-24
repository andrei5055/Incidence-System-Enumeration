#include "matrix.h"

int CColOrbitCS::m_maxElement = 0;

CColOrbit::CColOrbit()
{ 
#if !WAIT_THREADS
	setNext(NULL);
#endif
}

CColOrbit *CColOrbit::InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit *pColOrbit, int idx)
{
    CColOrbit *pColOrbitNext = (CColOrbit *)((char *)this + lenFragm * colOrbitLen);
    Init(lenFragm, pColOrbitNext);
    InitEntryCntrs(pColOrbit, idx);
    return pColOrbitNext;
}

void CColOrbit::clone(const CColOrbit *pColOrb)
{
	setLenght(pColOrb->lenght());
	InitEntryCntrs(pColOrb);
	CColOrbit *pNext = pColOrb->next();
	if (pNext) {
		CColOrbit *pNxt = (CColOrbit *)((char *)this + ((const char *)pNext - (const char *)pColOrb));
		pNxt->clone(pNext);
		setNext(pNxt);
	} else
		setNext(NULL);
}

void CColOrbitCS::InitOrbit(int maxElement)
{
	 m_pEntryCntrs = new unsigned int [maxElement];
	 memset(getEntryCntrs(), 0, maxElement * sizeof(getEntryCntrs()[0]));
}

void CColOrbitCS::InitEntryCntrs(const CColOrbit *pParent, int idx)
{ 
	memcpy(getEntryCntrs(), pParent->getEntryCntrs(), maxElement() * sizeof(getEntryCntrs()[0]));
	++*(getEntryCntrs() + idx);
}


