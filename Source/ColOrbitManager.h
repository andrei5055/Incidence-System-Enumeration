//
//  OrbitManager.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/22/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__OrbitManager__
#define __BIBD_Mac__OrbitManager__

#include <iostream>
#include "DataTypes.h"
#include "ColOrbits.h"
#ifdef USE_CUDA
#include "cuda_runtime.h"
#endif

template<class T> class CColOrbitManager
{
public:
	CC CColOrbitManager(int rank, T nRows, T nCol)		{ InitiateColOrbitManager(rank, nRows, nCol); }
	CK CColOrbitManager()								{}
	CC ~CColOrbitManager()								{ ReleaseColOrbitManager(); }
	CC void InitiateColOrbitManager(int matrRank, T nRows, T nCol, void *pMem = NULL);
	CC void ReleaseColOrbitManager();
	CC inline size_t colOrbitLen() const				{ return m_nColOrbLen; }
	CC inline CColOrbit<T> **colOrbitPntr() const		{ return m_ppOrb; }
	CC inline T colNumb() const							{ return m_nCol; }
	CC inline CColOrbit<T> **colOrbits() const			{ return m_pColOrb; }
    CC inline T currentRowNumb() const					{ return m_nCurrRow; }
	CC inline CColOrbit<T> **colOrbitsIni() const		{ return m_pColOrbIni; }
	CC void initiateColOrbits(size_t nRow, bool using_IS_enumerator, const CColOrbitManager<T> *pMaster = NULL, void *pMem = NULL);
	CK void copyColOrbitInfo(const CColOrbitManager<T> *pColOrb, T nRow);
	CC void restoreColOrbitInfo(T nRow, const size_t *pColOrbInfo) const;
    CC void closeColOrbits() const;
protected:
	CK inline CColOrbit<T> *colOrbitIni(size_t n) const	{ return *(colOrbitsIni() + n); }
	CK virtual CColOrbit<T> **unforcedOrbits(size_t n) const { return NULL; }
	CK inline void setColOrbitCurr(CColOrbit<T> *pntr)  { m_pColOrb[currentRowNumb()] = pntr; }
	CK inline void resetColOrbitCurr()                  { setColOrbitCurr(*(colOrbitPntr() + colNumb() * currentRowNumb())); }
	CK void resetUnforcedColOrb();
	CK virtual void resetFirstUnforcedRow()				{}
	CK void addForciblyConstructedColOrbit(CColOrbit<T> *pColOrbit, int n);
	CK inline CColOrbit<T> **currUnforcedOrbPtr() const { return m_ppUnforcedColOrbCurr; }
	CK inline void setCurrUnforcedOrbPtr(size_t nRow)	{ m_ppUnforcedColOrbCurr = unforcedColOrbPntr() + unfColIdx(nRow); }
    CC inline void setCurrentRowNumb(T n)				{ m_nCurrRow = n; }
	CC inline CColOrbit<T> **unforcedColOrbPntr() const	{ return m_ppUnforcedColOrb; }
	CK inline CColOrbit<T> *colOrbit(size_t idx) const	{ return m_pColOrb[idx]; }
    CC inline int rankMatr() const                      { return m_nRank; }
	CC inline size_t rowMaster() const					{ return m_nRowMaster; }
private:
	CK inline size_t unfColIdx(size_t r, int idx = 0) const{ return r * rankMatr() + idx; }
	CC inline void setColOrbitLen(size_t len)			{ m_nColOrbLen = len; }
	CC inline void setRowMaster(size_t val)				{ m_nRowMaster = val; }

	T m_nCurrRow;
    int m_nRank;
    T m_nCol;
	CColOrbit<T> **m_ppUnforcedColOrb;
	CColOrbit<T> **m_ppUnforcedColOrbCurr;
    CColOrbit<T> **m_pColOrb;
    CColOrbit<T> **m_ppOrb;
	CColOrbit<T> **m_pColOrbIni;
	size_t m_nColOrbLen;
	size_t m_nRowMaster;
};

template<class T>
void CColOrbitManager<T>::InitiateColOrbitManager(int matrRank, T nRows, T nCol, void *pMem)
{
	m_ppOrb = NULL;
	m_nRank = matrRank;
	m_nCol = nCol;
	setCurrentRowNumb(0);
	if (!pMem) {
		m_pColOrb = new CColOrbit<T> *[2 * nRows];
		const auto len = nRows * rankMatr();
		m_ppUnforcedColOrb = new CColOrbit<T> *[len];
		memset(unforcedColOrbPntr(), 0, len * sizeof(unforcedColOrbPntr()[0]));
	} else {
		// Using previously allocated memory
		m_pColOrb = (CColOrbit<T> **)pMem;
		m_ppUnforcedColOrb = NULL;
	}

	m_pColOrbIni = m_pColOrb + nRows;
}

template<class T>
void CColOrbitManager<T>::ReleaseColOrbitManager()
{
	delete[] unforcedColOrbPntr();
	if (unforcedColOrbPntr())
		delete[] m_pColOrb;
}

template<class T>
void CColOrbitManager<T>::initiateColOrbits(size_t nRow, bool using_IS_enumerator, const CColOrbitManager<T> *pMaster, void *pMem)
{
	// Number of CColOrbits taken from pMaster
	setRowMaster(pMaster ? pMaster->currentRowNumb() + 1 : 0);
	const size_t fromMaster = WAIT_THREADS ? rowMaster() * colNumb() : 0;
	const size_t nCol_2 = nRow * colNumb();
#ifndef USE_CUDA		// NOT yet implemented for GPU
	const size_t lenColOrbitElement = using_IS_enumerator? sizeof(CColOrbitIS<T>) : sizeof(CColOrbitCS<T>);
#else
	const size_t lenColOrbitElement = sizeof(CColOrbitIS<T>);
#endif
	setColOrbitLen(lenColOrbitElement);

	if (pMem) {
		m_ppOrb = (CColOrbit<T> **)pMem;
		pMem = (char *)pMem + nCol_2 * sizeof(CColOrbit<T> *);
	} else
		m_ppOrb = new CColOrbit<T> *[nCol_2];


	if (!using_IS_enumerator) {
#ifndef USE_CUDA		// NOT yet implemented for GPU
		const int maxElement = rankMatr() - 1;
		auto pColOrbitsCS = pMem ? (CColOrbitCS<T> *)pMem : new CColOrbitCS<T>[nCol_2 - fromMaster];
		CColOrbitCS<T>::setMaxElement(maxElement);
		for (auto i = nCol_2; i-- > fromMaster;) {
			pColOrbitsCS[i].InitOrbit(maxElement);
			m_ppOrb[i] = pColOrbitsCS + i - fromMaster;
		}
#endif
	}
	else {
		auto pColOrbitsIS = pMem? (CColOrbitIS<T> *)pMem : new CColOrbitIS<T>[nCol_2 - fromMaster];
		for (auto i = nCol_2; i-- > fromMaster;)
			m_ppOrb[i] = pColOrbitsIS + i - fromMaster;
	}

	auto i = nRow;
	auto iMin = WAIT_THREADS ? rowMaster() : 0;
	while (i-- > iMin)
		m_pColOrbIni[i] = m_pColOrb[i] = m_ppOrb[i * colNumb()];

	if (rowMaster() > 0) {
#if WAIT_THREADS
		const size_t len = rowMaster() * sizeof(m_pColOrbIni[0]);
		memcpy(m_pColOrbIni, pMaster->colOrbitsIni(), len);
		memcpy(m_pColOrb, pMaster->colOrbits(), len);
		memcpy(unforcedColOrbPntr(), pMaster->unforcedColOrbPntr(), rowMaster() * 2 * sizeof(*unforcedColOrbPntr()));
#else
		// When we are not waiting for all threads to finish on rowMaster() level, 
		// we need to reconstruct the information regarding master's column orbits here 
		// because it could be changed when master will continue its calculation

		for (size_t i = 0; i < rowMaster(); i++) {
			m_pColOrbIni[i]->clone(pMaster->colOrbitsIni()[i]);
			const size_t idx = ((char *)pMaster->colOrbits()[i] - (char *)pMaster->colOrbitsIni()[i]) / lenColOrbitElement;
			if (idx)
				m_pColOrb[i] = (CColOrbit<T> *)((char *)m_pColOrb[i] + idx * lenColOrbitElement);

			m_pColOrb[i]->clone(pMaster->colOrbits()[i]);

			for (int j = 0; j < 2; j++) {
				const size_t idx = (i << 1) + j;
				if (pMaster->unforcedColOrbPntr()[idx]) {
					const size_t shift = (char *)pMaster->unforcedColOrbPntr()[idx] - (char *)pMaster->colOrbitsIni()[i];
					unforcedColOrbPntr()[idx] = (CColOrbit<T> *)((char *)colOrbitsIni()[i] + shift);
					unforcedColOrbPntr()[idx]->clone(pMaster->unforcedColOrbPntr()[idx]);
				}
			}
		}
#endif
	}
	else
		m_pColOrb[0]->Init(colNumb());
}

template<class T>
void CColOrbitManager<T>::copyColOrbitInfo(const CColOrbitManager<T> *pColOrb, T nRow)
{
	// Function copy existing information about the orbits of columns  
	// into new structures, which could be used on GPU
	const auto colOrbit = pColOrb->colOrbits();
	const auto colOrbitIni = pColOrb->colOrbitsIni();
	while (nRow--) {
		const auto pColOrbitIni = colOrbitIni[nRow];
		auto pColOrbit = colOrbit[nRow];
		if (!pColOrbit) {
			colOrbits()[nRow] = NULL;
			continue;
		}

		const auto pColOrbitIniTo = (char *)colOrbitsIni()[nRow];
		auto pColOrbitToNext = (CColOrbit<T> *)(pColOrbitIniTo + ((char *)pColOrbit - (char *)pColOrbitIni));
		colOrbits()[nRow] = pColOrbitToNext;
		while (pColOrbit) {
			auto pColOrbitTo = pColOrbitToNext;
			const auto len = pColOrbit->length();
			pColOrbit = pColOrbit->next();
			pColOrbitToNext = pColOrbit ? (CColOrbit<T> *)(pColOrbitIniTo + ((char *)pColOrbit - (char *)pColOrbitIni)) : NULL;
			pColOrbitTo->Init(len, pColOrbitToNext);
		}
	}
}

template<class T>
void CColOrbitManager<T>::restoreColOrbitInfo(T nRow, const size_t *pColOrbInfo) const
{
	// Function copy existing information about the orbits of columns  
	// into new structures, which could be used on GPU
	const auto colOrbit = colOrbits();
	const auto colOrbitIni = colOrbitsIni();
	while (nRow--) {
		auto iMax = *pColOrbInfo++;
		if (iMax == UINT64_MAX) {
			colOrbits()[nRow] = NULL;
			continue;
		}

		const auto pColOrbitIniTo = (char *)colOrbitsIni()[nRow];
		auto pColOrbitToNext = (CColOrbit<T> *)(pColOrbitIniTo + *pColOrbInfo++);
		colOrbits()[nRow] = pColOrbitToNext;
		for (int i = 1; i <= iMax; i++) {
			auto pColOrbitTo = pColOrbitToNext;
			auto len = *pColOrbInfo++;
			pColOrbitToNext = i < iMax ? (CColOrbit<T> *)(pColOrbitIniTo + *pColOrbInfo++) : NULL;
			pColOrbitTo->Init((T)len, pColOrbitToNext);
		}
	}
}

template<class T>
void CColOrbitManager<T>::closeColOrbits() const
{
	delete[] colOrbitsIni()[WAIT_THREADS ? rowMaster() : 0];
	delete[] colOrbitPntr();
}

template<class T>
void CColOrbitManager<T>::addForciblyConstructedColOrbit(CColOrbit<T> *pColOrbit, int n)
{
	auto ppTmp = currUnforcedOrbPtr() + n;
	pColOrbit->setNext(*ppTmp);
	*ppTmp = pColOrbit;
}

template<class T>
void CColOrbitManager<T>::resetUnforcedColOrb()
{
	memset(unforcedColOrbPntr() + unfColIdx(currentRowNumb()), 0, rankMatr() * sizeof(*unforcedColOrbPntr()));
	resetFirstUnforcedRow();
}

#endif /* defined(__BIBD_Mac__OrbitManager__) */
