#pragma once
#include "Vector.h"

//template<class T> class CRowSolution;

IClass1Def(S, ColOrbit) {
 public:
	CC CColOrbit() {
#if !WAIT_THREADS
		setNext(NULL);
#endif
		assignOrbID(myOrbID);
	}

	CC virtual ~CColOrbit()							{}   
 	CC void Init(S length, CColOrbit *pNext = NULL)	{ setLength(length); setNext(pNext); }
	CC inline S length() const						{ return m_nLength; }
	CC inline CColOrbit *next() const				{ return m_pNext; }
	CC inline void setNext(CColOrbit *pntr)			{ m_pNext = pntr; }
	CC virtual void InitEntryCntrs(const CColOrbit *pParent, int idx = 0) = 0;
	CC virtual unsigned int *getEntryCntrs() const	{ return NULL; }
	CC virtual int columnWeight() const				{ return 0; }
    CC CColOrbit *InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit *pColOrbit, int idx);
	CC void clone(const CColOrbit *pColOrb);
	CC inline void setLength(S len)					{ m_nLength = len; }
private:
	S m_nLength;
	IClass1(S, ColOrbit) *m_pNext;
 protected:
	 MY_ORB_ID
};

TClass1(S, ColOrbit, ColOrbPntr)::InitOrbit(int lenFragm, size_t colOrbitLen, const ColOrbPntr pColOrbit, int idx)
{
	auto pColOrbitNext = (ColOrbPntr)((char *)this + lenFragm * colOrbitLen);
	Init(lenFragm, pColOrbitNext);
	InitEntryCntrs(pColOrbit, idx);
	return pColOrbitNext;
}

TClass1(S, ColOrbit, void)::clone(const ColOrbPntr pColOrb)
{
	setLength(pColOrb->length());
	InitEntryCntrs(pColOrb);
	auto pNext = pColOrb->next();
	if (pNext) {
		auto pNxt = (ColOrbPntr)((char *)this + ((const char *)pNext - (const char *)pColOrb));
		pNxt->clone(pNext);
		setNext(pNxt);
	}
	else
		setNext(NULL);
}

IClass1Def(S, ColOrbitIS) : public IClass1(S, ColOrbit)
{// Orbits for Incidence Systems
 public:
	 CC CColOrbitIS()													{ setColumnWeight(0); } 
	 CC ~CColOrbitIS()													{} 
	 virtual CC void InitEntryCntrs(const ColOrbPntr pParent, int idx = 0)	{ setColumnWeight(pParent->columnWeight() + idx); }
	 virtual CC int columnWeight() const								{ return m_nColumnWeight; }
 protected:

 private:
	 virtual CC void setColumnWeight(int val)							{ m_nColumnWeight = val; }
	 int m_nColumnWeight;
};

IClass1Def(S, ColOrbitCS) : public IClass1(S, ColOrbit)
{// Orbits for Colored Incidence Systems
 public:
	CC ~CColOrbitCS()													 { delete [] getEntryCntrs(); }
	CC virtual void InitEntryCntrs(const ColOrbPntr pParent, int idx) {
		memcpy(getEntryCntrs(), pParent->getEntryCntrs(), maxElement() * sizeof(getEntryCntrs()[0]));
		++*(getEntryCntrs() + idx);
	}

	void InitOrbit(int maxElement) {
		m_pEntryCntrs = new unsigned int[maxElement];
		memset(getEntryCntrs(), 0, maxElement * sizeof(getEntryCntrs()[0]));
	}

	inline static void setMaxElement(int maxElement)					{ m_maxElement = maxElement; }
	CC inline unsigned int *getEntryCntrs() const						{ return m_pEntryCntrs; }
 
 protected:
 private:
	CC inline static int maxElement()									{ return m_maxElement; }
	unsigned int *m_pEntryCntrs;
	static int m_maxElement;
};
