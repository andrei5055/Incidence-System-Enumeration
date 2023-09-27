#pragma once
#include "Vector.h"

Class1Def(CColOrbit) {
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
	CC inline auto next() const						{ return m_pNext; }
	CC inline void setNext(CColOrbit *pntr)			{ m_pNext = pntr; }
	CC virtual void InitEntryCntrs(const CColOrbit *pParent, int idx = 0) = 0;
	CC virtual unsigned int *getEntryCntrs() const	{ return NULL; }
	CC virtual int columnWeight() const				{ return 0; }
	CC CColOrbit *InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit *pColOrbit, int idx);
	CC void clone(const CColOrbit *pColOrb);
	CC inline void setLength(S len)					{ m_nLength = len; }
private:
	S m_nLength;
	ColOrbPntr m_pNext;
 protected:
	 MY_ORB_ID
};

FClass1(CColOrbit, ColOrbPntr)::InitOrbit(int lenFragm, size_t colOrbitLen, const ColOrbPntr pColOrbit, int idx)
{
	auto pColOrbitNext = (ColOrbPntr)((char *)this + lenFragm * colOrbitLen);
	Init(lenFragm, pColOrbitNext);
	InitEntryCntrs(pColOrbit, idx);
	return pColOrbitNext;
}

FClass1(CColOrbit, void)::clone(const ColOrbPntr pColOrb)
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

Class1Def(CColOrbitIS) : public Class1(CColOrbit)
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

Class1Def(CColOrbitCS) : public Class1(CColOrbit)
{// Orbits for Colored Incidence Systems
 public:
	CC ~CColOrbitCS()													 { delete [] getEntryCntrs(); }
	CC virtual void InitEntryCntrs(const ColOrbPntr pParent, int idx) {
		memcpy(getEntryCntrs(), pParent->getEntryCntrs(), maxElement() * sizeof(getEntryCntrs()[0]));
		++*(getEntryCntrs() + idx);
	}

	void InitOrbit(int maxElement) {
		setMaxElement(maxElement);
		m_pEntryCntrs = new unsigned int[maxElement];
		memset(getEntryCntrs(), 0, maxElement * sizeof(getEntryCntrs()[0]));
	}

	CC inline unsigned int *getEntryCntrs() const						{ return m_pEntryCntrs; }
 
 protected:
 private:
	inline void setMaxElement(int maxElement)							{ m_maxElement = maxElement; }
	CC inline int maxElement()											{ return m_maxElement; }
	unsigned int *m_pEntryCntrs = NULL;
	int m_maxElement = -1;
};

class CRank {
public:
	CRank(int nRows, int rank = 2) : m_rank(rank)	{
		if (nRows <= 8)  // We could use sizeof(void *) as the index in m_pShift array 
			nRows = 9;

		m_pShift = new size_t[nRows];
		m_pShift[0] = 0;
		for (int i = 1; i < nRows; i++)
			m_pShift[i] = rank + m_pShift[i - 1];
	}
	~CRank()                                                            { delete [] m_pShift; }
	CC inline auto rank() const											{ return m_rank; }
//	CK inline auto shiftToUnforcedOrbit(S nRow) const					{ return m_pShift[nRow]; }
protected:
	size_t* m_pShift;
private:
	int m_rank;
};