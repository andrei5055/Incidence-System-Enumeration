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

Class1Def(CColOrbitManager)
{
public:
	CC CColOrbitManager(int rank, S nRows, S nCol)		{ InitiateColOrbitManager(rank, nRows, nCol); }
	CK CColOrbitManager()								{}
	CC ~CColOrbitManager()								{ ReleaseColOrbitManager(); }
	CC void InitiateColOrbitManager(int matrRank, S nRows, S nCol, void *pMem = NULL);
	CC void ReleaseColOrbitManager();
	CC inline auto colOrbitLen() const					{ return m_nColOrbLen; }
	CC inline auto colOrbitPntr() const					{ return m_ppOrb; }
	CC inline S colNumb() const							{ return m_nCol; }
	CC inline auto colOrbits() const					{ return m_pColOrb; }
	CC inline S currentRowNumb() const					{ return m_nCurrRow; }
	CC inline auto colOrbitsIni() const					{ return m_pColOrbIni; }
	CC inline int rankMatr() const						{ return m_nRank; }
	CC void initiateColOrbits(S nRows, S firstRow, const Class1(BlockGroupDescr) *pGroupDesct, bool using_IS_enumerator, const Class1(CColOrbitManager) *pMaster = NULL, void *pMem = NULL);
	CK void copyColOrbitInfo(const Class1(CColOrbitManager) *pColOrb, S nRow);
	CC void restoreColOrbitInfo(S nRow, const size_t *pColOrbInfo) const;
	CC void closeColOrbits() const;
protected:
	CK inline auto colOrbitIni(size_t n) const			{ return *(colOrbitsIni() + n); }
	CK virtual ColOrbPntr *unforcedOrbits(size_t n) const { return NULL; }
	CK inline void setColOrbitCurr(ColOrbPntr pntr)		{ m_pColOrb[currentRowNumb()] = pntr; }
	CK inline void resetColOrbitCurr()                  { setColOrbitCurr(*(colOrbitPntr() + colNumb() * currentRowNumb())); }
	CK void resetUnforcedColOrb();
	CK virtual void resetFirstUnforcedRow()				{}
	CK void addForciblyConstructedColOrbit(ColOrbPntr pColOrbit, int n);
	CK inline ColOrbPntr* currUnforcedOrbPtr() const	{ return m_ppUnforcedColOrbCurr; }
	CK inline void setCurrUnforcedOrbPtr(size_t nRow)	{ m_ppUnforcedColOrbCurr = unforcedColOrbPntr() + unfColIdx(nRow); }
	CC inline void setCurrentRowNumb(S n)				{ m_nCurrRow = n; }
	CC inline ColOrbPntr *unforcedColOrbPntr() const	{ return m_ppUnforcedColOrb; }
	CK inline ColOrbPntr colOrbit(S idx) const			{ return m_pColOrb[idx]; }
	CC inline size_t rowMaster() const					{ return m_nRowMaster; }
private:
	CK inline size_t unfColIdx(size_t r, int idx = 0) const{ return r * rankMatr() + idx; }
	CC inline void setColOrbitLen(size_t len)			{ m_nColOrbLen = len; }
	CC inline void setRowMaster(size_t val)				{ m_nRowMaster = val; }

	S m_nCurrRow;
    int m_nRank;
    S m_nCol;
	ColOrbPntr *m_ppUnforcedColOrb;
	ColOrbPntr *m_ppUnforcedColOrbCurr;
	ColOrbPntr *m_pColOrb;
	ColOrbPntr *m_ppOrb;
	ColOrbPntr *m_pColOrbIni;
	size_t m_nColOrbLen;
	size_t m_nRowMaster;
	bool m_IS_enumerator;
};

FClass1(CColOrbitManager, void)::InitiateColOrbitManager(int matrRank, S nRows, S nCol, void *pMem)
{
	m_ppOrb = NULL;
	m_nRank = matrRank;
	m_nCol = nCol;
	if (!pMem) {
//		???
		m_pColOrb = new ColOrbPntr[2 * nRows];
		const auto len = nRows * rankMatr();
		m_ppUnforcedColOrb = new ColOrbPntr[len];
		memset(unforcedColOrbPntr(), 0, len * sizeof(unforcedColOrbPntr()[0]));
	} else {
		// Using previously allocated memory
		m_pColOrb = (ColOrbPntr *)pMem;
		m_ppUnforcedColOrb = NULL;
	}

	m_pColOrbIni = m_pColOrb + nRows;
}

FClass1(CColOrbitManager, void)::ReleaseColOrbitManager()
{
	delete[] unforcedColOrbPntr();
	if (unforcedColOrbPntr())
		delete[] m_pColOrb;
}

FClass1(CColOrbitManager, void)::initiateColOrbits(S nRows, S firstRow, const Class1(BlockGroupDescr) *pGroupDesct, bool using_IS_enumerator, const  Class1(CColOrbitManager) *pMaster, void *pMem)
{
	m_IS_enumerator = using_IS_enumerator;
	// Number of CColOrbits taken from pMaster
	setRowMaster(pMaster ? pMaster->currentRowNumb() + 1 : 0);

	// In order not to re-create orbits on the way back, we will store them in different places.
	// Therefore, we must have number of rows as a multiplier for fromMaster and nCol_2.
	const size_t fromMaster = WAIT_THREADS ? rowMaster() * colNumb() : 0;
	const size_t nCol_2 = nRows * colNumb();
	const size_t numParts = pGroupDesct ? pGroupDesct->numParts() : 1;
#ifndef USE_CUDA		// NOT yet implemented for GPU
	const size_t lenColOrbitElement = using_IS_enumerator? sizeof(Class1(CColOrbitIS)) : sizeof(Class1(CColOrbitCS));
#else
	const size_t lenColOrbitElement = sizeof(CColOrbitIS<S>);
#endif
	setColOrbitLen(lenColOrbitElement);

	if (pMem) {
		m_ppOrb = (ColOrbPntr *)pMem;
		pMem = (char *)pMem + nCol_2 * sizeof(ColOrbPntr);
	} else
		m_ppOrb = new ColOrbPntr[nCol_2];


	const int maxElement = rankMatr();
	if (!using_IS_enumerator) {
#ifndef USE_CUDA		// NOT yet implemented for GPU
		auto pColOrbitsCS = pMem ? (Class1(CColOrbitCS) *)pMem : new  Class1(CColOrbitCS)[nCol_2 - fromMaster];
		for (auto i = nCol_2; i-- > fromMaster;) {
			pColOrbitsCS[i].InitOrbit(maxElement - 1);
			m_ppOrb[i] = pColOrbitsCS + i - fromMaster;
		}
#endif
	}
	else {
		auto pColOrbitsIS = pMem? (Class1(CColOrbitIS) *)pMem : new  Class1(CColOrbitIS)[nCol_2 - fromMaster];
		for (auto i = nCol_2; i-- > fromMaster;)
			m_ppOrb[i] = pColOrbitsIS + i - fromMaster;
	}

	auto i = nRows;
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
				m_pColOrb[i] = (ColOrbPntr)((char *)m_pColOrb[i] + idx * lenColOrbitElement);

			m_pColOrb[i]->clone(pMaster->colOrbits()[i]);
			for (int j = 0; j < maxElement; j++) {
				const size_t idx = i * maxElement + j;
				if (pMaster->unforcedColOrbPntr()[idx]) {
					const size_t shift = (char *)pMaster->unforcedColOrbPntr()[idx] - (char *)pMaster->colOrbitsIni()[i];
					unforcedColOrbPntr()[idx] = (ColOrbPntr)((char *)colOrbitsIni()[i] + shift);
					unforcedColOrbPntr()[idx]->clone(pMaster->unforcedColOrbPntr()[idx]);
				}
			}
		}
#endif
	}
	else {
		if (numParts > 1) {
			//			for (S i = 0; i < pMatrix->numParts(); i++)
			//	m_pColOrb[firstRow + i]->Init(m_nPartInfo->element((i << 1) + 1));
		}
		else
			m_pColOrb[firstRow]->Init(colNumb());
	}
}

FClass1(CColOrbitManager, void)::copyColOrbitInfo(const  Class1(CColOrbitManager) *pColOrb, S nRow)
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
		auto pColOrbitToNext = (ColOrbPntr)(pColOrbitIniTo + ((char *)pColOrbit - (char *)pColOrbitIni));
		colOrbits()[nRow] = pColOrbitToNext;
		while (pColOrbit) {
			auto pColOrbitTo = pColOrbitToNext;
			const auto len = pColOrbit->length();
			pColOrbit = pColOrbit->next();
			pColOrbitToNext = pColOrbit ? (ColOrbPntr)(pColOrbitIniTo + ((char *)pColOrbit - (char *)pColOrbitIni)) : NULL;
			pColOrbitTo->Init(len, pColOrbitToNext);
		}
	}
}

FClass1(CColOrbitManager, void)::restoreColOrbitInfo(S nRow, const size_t *pColOrbInfo) const
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
		auto pColOrbitToNext = (ColOrbPntr)(pColOrbitIniTo + *pColOrbInfo++);
		colOrbits()[nRow] = pColOrbitToNext;
		for (int i = 1; i <= iMax; i++) {
			auto pColOrbitTo = pColOrbitToNext;
			auto len = *pColOrbInfo++;
			pColOrbitToNext = i < iMax ? (ColOrbPntr)(pColOrbitIniTo + *pColOrbInfo++) : NULL;
			pColOrbitTo->Init(len, pColOrbitToNext);
		}
	}
}

FClass1(CColOrbitManager, void)::closeColOrbits() const
{
	auto pntr = colOrbitsIni()[WAIT_THREADS ? rowMaster() : 0];
	if (m_IS_enumerator)
		delete[] (Class1(CColOrbitIS) *)pntr;
	else
		delete[] (Class1(CColOrbitCS) *)pntr;

	delete[] colOrbitPntr();
}

FClass1(CColOrbitManager, void)::addForciblyConstructedColOrbit(ColOrbPntr pColOrbit, int n)
{
	auto ppTmp = currUnforcedOrbPtr() + n;
	pColOrbit->setNext(*ppTmp);
	*ppTmp = pColOrbit;
}

FClass1(CColOrbitManager, void)::resetUnforcedColOrb()
{
	memset(unforcedColOrbPntr() + unfColIdx(currentRowNumb()), 0, rankMatr() * sizeof(*unforcedColOrbPntr()));
	resetFirstUnforcedRow();
}

#endif /* defined(__BIBD_Mac__OrbitManager__) */
