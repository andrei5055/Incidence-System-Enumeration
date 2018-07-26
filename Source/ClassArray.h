/*
 *  ClassArray.h
 *  BlockDesign
 *
 *  Created by Andrei Ivanov on 12/16/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#if !defined(__ClassArr_H__)
#define __ClassArr_H__

#include <assert.h>
#include <memory>
#include <new>   //only supports Win32 and Mac
#include <cstring>
#include "DataTypes.h"

template<class TYPE, class ARG_TYPE>
class CClassArray
{
public:
	// Construction & destruction
	CC CClassArray()								{ m_pData = NULL; m_nSize = m_nMaxSize = m_nGrowBy = 0; }
	CC ~CClassArray();
	
	// Attributes
	CC inline size_t GetSize() const				{ return m_nSize; }
	CK inline size_t GetUpperBound() const			{ return m_nSize-1; }
	CC void SetSize(size_t nNewSize, int nGrowBy = -1);
	
	// Operations
	// Clean up
	CK void FreeExtra();
	CK void RemoveAll()								{ SetSize(0, -1); }
	
	// Accessing elements
	CC TYPE GetAt(size_t nIndex) const				{ return m_pData[nIndex]; }
	CK void SetAt(size_t nIndex, ARG_TYPE element)	{ m_pData[nIndex] = element; }
	CK TYPE& ElementAt(size_t nIndex)				{ return m_pData[nIndex]; }
	
	// Direct Access to the element data (may return NULL)
	CK const TYPE* GetData() const					{ return (const TYPE*)m_pData; }
	CK TYPE* GetData()								{ return (TYPE*)m_pData; }
	
	// Potentially growing the array
	CC void SetAtGrow(size_t nIndex, ARG_TYPE newElement);
	CC size_t Add(ARG_TYPE newElement)				{ const size_t nIndex = m_nSize; SetAtGrow(nIndex, newElement); return nIndex; }
	CK int Append(const CClassArray& src);
	CK void Copy(const CClassArray& src);
	
	// overloaded operator helpers
	CK TYPE operator[](int nIndex) const			{ return GetAt(nIndex); }
	CK TYPE& operator[](int nIndex)					{ return ElementAt(nIndex); }
	
	// Operations that move elements around
	CC void InsertAt(size_t nIndex, ARG_TYPE newElement, int nCount = 1);
	CK void RemoveAt(size_t nIndex, int nCount = 1);
	CK void InsertAt(size_t nStartIndex, CClassArray* pNewArray);
	
	// Implementation
protected:
	CC void ConstructElements(TYPE* pElements, size_t nCount);
	CC void DestructElements(TYPE* pElements, size_t nCount);
	CK void CopyElements(TYPE* pDest, const TYPE* pSrc, int nCount);
	
	TYPE* m_pData;   // the actual array of data
	size_t m_nSize;  // # of elements (upperBound - 1)
	size_t m_nMaxSize;  // max allocated
	int m_nGrowBy;   // grow amount
};

/////////////////////////////////////////////////////////////////////////////
// CArray<TYPE, ARG_TYPE> out-of-line functions

template<class TYPE, class ARG_TYPE>
CClassArray<TYPE, ARG_TYPE>::~CClassArray()
{
	if (m_pData != NULL)
	{
		DestructElements(m_pData, m_nSize);
		delete[] (char*)m_pData;
	}
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::SetSize(size_t nNewSize, int nGrowBy)
{	
	if (nGrowBy != -1)
		m_nGrowBy = nGrowBy;  // set new size
	
	if (nNewSize == 0)
	{
		// shrink to nothing
		if (m_pData != NULL)
		{
			DestructElements(m_pData, m_nSize);
			delete[] (char*)m_pData;
			m_pData = NULL;
		}
		m_nSize = m_nMaxSize = 0;
	}
	else if (m_pData == NULL)
	{
		// create one with exact size
#ifdef SIZE_T_MAX
		assert(nNewSize <= SIZE_T_MAX/sizeof(TYPE));    // no overflow
#endif
		m_pData = (TYPE*) new char[nNewSize * sizeof(TYPE)];
		ConstructElements(m_pData, nNewSize);
		m_nSize = m_nMaxSize = nNewSize;
	}
	else if (nNewSize <= m_nMaxSize)
	{
		// it fits
		if (nNewSize > m_nSize)
		{
			// initialize the new elements
			ConstructElements(&m_pData[m_nSize], nNewSize-m_nSize);
		}
		else if (m_nSize > nNewSize)
		{
			// destroy the old elements
			DestructElements(&m_pData[nNewSize], m_nSize-nNewSize);
		}
		m_nSize = nNewSize;
	}
	else
	{
		// otherwise, grow array
		size_t nGrowBy = m_nGrowBy;
		if (nGrowBy == 0)
		{
			// heuristically determine growth when nGrowBy == 0
			//  (this avoids heap fragmentation in many situations)
			nGrowBy = m_nSize / 8;
			nGrowBy = (nGrowBy < 4) ? 4 : ((nGrowBy > 1024) ? 1024 : nGrowBy);
		}
		size_t nNewMax;
		if (nNewSize < m_nMaxSize + nGrowBy)
			nNewMax = m_nMaxSize + nGrowBy;  // granularity
		else
			nNewMax = nNewSize;  // no slush
		
		assert(nNewMax >= m_nMaxSize);  // no wrap around
#ifdef SIZE_T_MAX
		assert(nNewMax <= SIZE_T_MAX/sizeof(TYPE)); // no overflow
#endif
		TYPE* pNewData = (TYPE*) new char[nNewMax * sizeof(TYPE)];
		
		// copy new data from old
		memcpy(pNewData, m_pData, m_nSize * sizeof(TYPE));
		
		// construct remaining elements
		assert(nNewSize > m_nSize);
		ConstructElements(&pNewData[m_nSize], nNewSize-m_nSize);
		
		// get rid of old stuff (note: no destructors called)
		delete[] (char*)m_pData;
		m_pData = pNewData;
		m_nSize = nNewSize;
		m_nMaxSize = nNewMax;
	}
}

template<class TYPE, class ARG_TYPE>
int CClassArray<TYPE, ARG_TYPE>::Append(const CClassArray& src)
{
	assert(this != &src);   // cannot append to itself
	
	int nOldSize = m_nSize;
	SetSize(m_nSize + src.m_nSize);
	CopyElements(m_pData + nOldSize, src.m_pData, src.m_nSize);
	return nOldSize;
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::Copy(const CClassArray& src)
{
	assert(this != &src);   // cannot append to itself
	
	SetSize(src.m_nSize);
	CopyElements(m_pData, src.m_pData, src.m_nSize);
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::FreeExtra()
{
	if (m_nSize != m_nMaxSize)
	{
		// shrink to desired size
#ifdef SIZE_T_MAX
		assert(m_nSize <= SIZE_T_MAX/sizeof(TYPE)); // no overflow
#endif
		TYPE* pNewData = NULL;
		if (m_nSize != 0)
		{
			pNewData = (TYPE*) new char[m_nSize * sizeof(TYPE)];
			// copy new data from old
			memcpy(pNewData, m_pData, m_nSize * sizeof(TYPE));
		}
		
		// get rid of old stuff (note: no destructors called)
		delete[] (char*)m_pData;
		m_pData = pNewData;
		m_nMaxSize = m_nSize;
	}
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::SetAtGrow(size_t nIndex, ARG_TYPE newElement)
{
	if (nIndex >= m_nSize)
		SetSize(nIndex+1, -1);

	m_pData[nIndex] = newElement;
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::InsertAt(size_t nIndex, ARG_TYPE newElement, int nCount /*=1*/)
{
	assert(nCount > 0);     // zero or negative size not allowed
	
	if (nIndex >= m_nSize)
	{
		// adding after the end of the array
		SetSize(nIndex + nCount, -1);   // grow so nIndex is valid
	}
	else
	{
		// inserting in the middle of the array
		size_t nOldSize = m_nSize;
		SetSize(m_nSize + nCount, -1);  // grow it to new size
		// destroy intial data before copying over it
		DestructElements(&m_pData[nOldSize], nCount);
		// shift old data up to fill gap
#define USE_MEM_MOVE	0
#if  USE_MEM_MOVE
		memmove(&m_pData[nIndex+nCount], &m_pData[nIndex],
				(nOldSize-nIndex) * sizeof(TYPE));
#else
		for (size_t i = GetSize(); i-- > nIndex + nCount;)
			memcpy(m_pData + i, m_pData + i - nCount, sizeof(TYPE));
#endif
		// re-init slots we copied from
		ConstructElements(&m_pData[nIndex], nCount);
	}
	
	// insert new value in the gap
	assert(nIndex + nCount <= m_nSize);
	while (nCount--)
		m_pData[nIndex++] = newElement;
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::RemoveAt(size_t nIndex, int nCount)
{
	assert(nCount >= 0);
	assert(nIndex + nCount <= m_nSize);
	
	// just remove a range
	int nMoveCount = m_nSize - (nIndex + nCount);
	DestructElements(&m_pData[nIndex], nCount);
	if (nMoveCount)
		memmove(&m_pData[nIndex], &m_pData[nIndex + nCount],
				nMoveCount * sizeof(TYPE));
	m_nSize -= nCount;
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::InsertAt(size_t nStartIndex, CClassArray* pNewArray)
{
	assert(pNewArray != NULL);
	
	if (pNewArray->GetSize() > 0)
	{
		InsertAt(nStartIndex, pNewArray->GetAt(0), pNewArray->GetSize());
		for (int i = 0; i < pNewArray->GetSize(); i++)
			SetAt(nStartIndex + i, pNewArray->GetAt(i));
	}
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::ConstructElements(TYPE* pElements, size_t nCount)
{
	// first do bit-wise zero initialization
	std::memset((void*)pElements, 0, nCount * sizeof(TYPE));
	// then call the constructor(s)
	for (; nCount--; pElements++)
		::new((void*)pElements) TYPE;
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::DestructElements(TYPE* pElements, size_t nCount)
{
	// call the destructor(s)
	for (; nCount--; pElements++)
		pElements->~TYPE();
}

template<class TYPE, class ARG_TYPE>
void CClassArray<TYPE, ARG_TYPE>::CopyElements(TYPE* pDest, const TYPE* pSrc, int nCount)
{
	// default is element-copy using assignment
	while (nCount--)
		*pDest++ = *pSrc++;
}

#endif
