#pragma once
#include "Vector.h"

template<class T> class CRowSolution;

template<class T> class CColOrbit //: public CArray<int, int>
{
 public:
	CC CColOrbit() {
#if !WAIT_THREADS
		setNext(NULL);
#endif
		assignOrbID(myOrbID);
	}

	virtual CC ~CColOrbit()								{}   
 	void CC Init(T length, CColOrbit *pNext = NULL)		{ setLength(length); setNext(pNext); }
	inline CC T length() const							{ return m_nLength; }
	inline CC CColOrbit *next() const					{ return m_pNext; }
	inline CC void setNext(CColOrbit *pntr)             { m_pNext = pntr; }
	virtual CC void InitEntryCntrs(const CColOrbit *pParent, int idx = 0) = 0;
	virtual unsigned int *getEntryCntrs() const			{ return NULL; }
	virtual CC int colomnWeight() const					{ return 0; }
    CC CColOrbit *InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit *pColOrbit, int idx);
	CC void clone(const CColOrbit *pColOrb);
	inline CC void setLength(T len)						{ m_nLength = len; }
private:
	T m_nLength;
	CColOrbit<T> *m_pNext;
 protected:
	 MY_ORB_ID
};

template<class T>
CColOrbit<T> *CColOrbit<T>::InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit<T> *pColOrbit, int idx)
{
	CColOrbit<T> *pColOrbitNext = (CColOrbit<T> *)((char *)this + lenFragm * colOrbitLen);
	Init(lenFragm, pColOrbitNext);
	InitEntryCntrs(pColOrbit, idx);
	return pColOrbitNext;
}

template<class T>
void CColOrbit<T>::clone(const CColOrbit<T> *pColOrb)
{
	setLength(pColOrb->length());
	InitEntryCntrs(pColOrb);
	CColOrbit<T> *pNext = pColOrb->next();
	if (pNext) {
		CColOrbit<T> *pNxt = (CColOrbit<T> *)((char *)this + ((const char *)pNext - (const char *)pColOrb));
		pNxt->clone(pNext);
		setNext(pNxt);
	}
	else
		setNext(NULL);
}

template<class T>class CColOrbitIS : public CColOrbit<T>
{// Orbits for Incidence Systems
 public:
	 CC CColOrbitIS()													{ setColomnWeight(0); } 
	 CC ~CColOrbitIS()													{} 
	 virtual CC void InitEntryCntrs(const CColOrbit<T> *pParent, int idx = 0)	{ setColomnWeight(pParent->colomnWeight() + idx); }
	 virtual CC int colomnWeight() const								{ return m_nColomnWeight; }
 protected:

 private:
	 virtual CC void setColomnWeight(int val)							{ m_nColomnWeight = val; }
	 int m_nColomnWeight;
};

#ifndef USE_CUDA		// NOT yet implemented for GPU
template<class T>class CColOrbitCS : public CColOrbit<T>
{// Orbits for Colored Incidence Systems
 public:
	~CColOrbitCS()														{ delete [] getEntryCntrs(); }
	virtual CC void InitEntryCntrs(const CColOrbit<T> *pParent, int idx) {
		memcpy(getEntryCntrs(), pParent->getEntryCntrs(), maxElement() * sizeof(getEntryCntrs()[0]));
		++*(getEntryCntrs() + idx);
	}

	void InitOrbit(int maxElement) {
		m_pEntryCntrs = new unsigned int[maxElement];
		memset(getEntryCntrs(), 0, maxElement * sizeof(getEntryCntrs()[0]));
	}

	inline static void setMaxElement(int maxElement)					{ m_maxElement = maxElement; }
	inline unsigned int *getEntryCntrs() const							{ return m_pEntryCntrs; }
 
 protected:
 private:
	inline static int maxElement()										{ return m_maxElement; }
	unsigned int *m_pEntryCntrs;
	static int m_maxElement;
};
#endif
