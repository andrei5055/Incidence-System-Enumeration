#include "Vector.h"

class CRowSolution;

class CColOrbit //: public CArray<int, int>
{
 public:
	CColOrbit();
	virtual ~CColOrbit()								{}   
 	void Init(size_t lenght, CColOrbit *pNext = NULL)	{ setLenght(lenght); setNext(pNext); }
	inline size_t lenght() const						{ return m_nLenght; }
	inline CColOrbit *next() const						{ return m_pNext; }
	inline void setNext(CColOrbit *pntr)                { m_pNext = pntr; }
	virtual void InitEntryCntrs(const CColOrbit *pParent, int idx = 0) = 0;
	virtual unsigned int *getEntryCntrs() const			{ return NULL; }
	virtual int colomnWeight() const					{ return 0; }
    CColOrbit *InitOrbit(int lenFragm, size_t colOrbitLen, const CColOrbit *pColOrbit, int idx);
	void clone(const CColOrbit *pColOrb);
 private:
	inline void setLenght(size_t len)					{ m_nLenght = len; }
	size_t m_nLenght;
	CColOrbit *m_pNext;
 protected:
};

class CColOrbitIS : public CColOrbit
{// Orbits for Incidence Sistems
 public:
	 CColOrbitIS()														{ setColomnWeight(0); } 
	 ~CColOrbitIS()														{} 
	 virtual void InitEntryCntrs(const CColOrbit *pParent, int idx = 0)	{ setColomnWeight(pParent->colomnWeight() + idx); }
	 virtual int colomnWeight() const									{ return m_nColomnWeight; }
 protected:

 private:
	 virtual void setColomnWeight(int val)								{ m_nColomnWeight = val; }
	 int m_nColomnWeight;
};

class CColOrbitCS : public CColOrbit
{// Orbits for Colored Incidence Sistems
 public:
	~CColOrbitCS()														{ delete [] getEntryCntrs(); }
	virtual void InitEntryCntrs(const CColOrbit *pParent, int idx);
	void InitOrbit(int maxElement);
	inline static void setMaxElement(int maxElement)					{ m_maxElement = maxElement; }
	inline unsigned int *getEntryCntrs() const							{ return m_pEntryCntrs; }
 
 protected:
 private:
	inline static int maxElement()										{ return m_maxElement; }
	unsigned int *m_pEntryCntrs;
	static int m_maxElement;
};
