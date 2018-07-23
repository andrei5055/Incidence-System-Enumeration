#pragma once
#include <vector>
#include "VariableMapping.h"

typedef std::vector<const VECTOR_ELEMENT_TYPE *> CSolutionStorage;

template<class T>
class CVector
{
 public: 
 	CK CVector(size_t size = 0, CArrayOfVectorElements *pElement = NULL) {
		if (pElement) {
			InitCoordArray(pElement);
			return;
		}

		m_pCoord = new CArrayOfVectorElements();
		getCoord()->SetSize(size);
		m_Owner = true;
	}

	CK virtual ~CVector()								{ releaseCoord(); }
	CK inline void AddElement(T val)					{ getCoord()->Add(val); }
	CK inline T GetAt(size_t idx) const					{ return getCoord()->GetAt(idx); }
	CK inline size_t GetSize() const					{ return getCoord()->GetSize(); }
	CK inline T *GetData() const						{ return getCoord()->GetData(); }
    CK inline void IncreaseVectorSize(size_t len)       { getCoord()->SetSize(GetSize() + len, -1); }
protected:
	CK inline CArrayOfVectorElements *getCoord() const	{ return m_pCoord; }
	CK void InitCoordArray(CArrayOfVectorElements *pCoord, bool freeFlag = true) {
		if (freeFlag)
			releaseCoord();

		m_pCoord = pCoord;
		m_Owner = false;
	}

    CK inline void SetAt(size_t idx, T val)				{ getCoord()->SetAt(idx, val); }
 private:
	CK void releaseCoord()								{ if (m_Owner) delete getCoord(); }
	bool m_Owner;
	CArrayOfVectorElements *m_pCoord;
};

