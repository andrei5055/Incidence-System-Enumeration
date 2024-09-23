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

Class1Def(CColOrbitManager) : public CRank
{
public:
	CC inline CColOrbitManager(uint rank, S nRows, S nCol, S nParts = 1) : CRank(nRows, rank) {
		InitiateColOrbitManager(rank, nRows, nCol, nParts);
	}
	CK CColOrbitManager()								{ m_ppColOrb = NULL; }
	CC ~CColOrbitManager()								{ closeColOrbits();
														  ReleaseColOrbitManager();
														}
	CC void InitiateColOrbitManager(uint matrRank, S nRows, S nCol, S numParts = 1, void *pMem = NULL);
	CC void ReleaseColOrbitManager();
	CC inline auto colOrbitLen() const					{ return m_nColOrbLen; }
	CC inline auto *colOrbitPntr() const				{ return m_ppOrb; }
	CC inline S colNumb() const							{ return m_nCol; }
	CC inline auto *colOrbits(S iPart = 0) const		{ return m_ppColOrb[iPart]; }
	CC inline S currentRowNumb() const					{ return m_nCurrRow; }
	CC inline auto colOrbitsIni(S iPart = 0) const		{ return m_ppColOrbIni[iPart]; }
	CC void initiateColOrbits(S nRows, S firstRow = 0, const Class1(BlockGroupDescr) *pGroupDesct = NULL, bool using_IS_enumerator = false,
		                      int use_master_solutions = 0, const Class1(CColOrbitManager) *pMaster = NULL, void *pMem = NULL);
	CK void copyColOrbitInfo(const Class1(CColOrbitManager) *pColOrb, S nRow);
	CC void restoreColOrbitInfo(S nRow, const size_t *pColOrbInfo) const;
	CC void closeColOrbits();
	CC inline auto colOrbit(S idx, S idxPart = 0) const				{ return m_ppColOrb[idxPart][idx]; }
	CC inline auto colOrbitIni(S nRow, S idxPart = 0) const			{ return *(colOrbitsIni(idxPart) + nRow); }
protected:
	CK inline void setColOrbitCurr(ColOrbPntr pntr, S idxPart)		{ setColOrbit(pntr, currentRowNumb(), idxPart); }
	CK inline void resetUnforcedColOrb(S idxPart, S nRow = 0)		{
		memset(pntr2UnforcedColOrb(nRow? nRow : currentRowNumb(), idxPart), 0, rank()*sizeof(*unforcedColOrbPntr()));
	}
	CK inline ColOrbPntr* pntr2UnforcedColOrb(S nRow, S nPart) const { return unforcedColOrbPntr(nPart) + m_pShift[nRow]; }
	CK void addForciblyConstructedColOrbit(ColOrbPntr pColOrbit, S nParts, S n);
	CK inline auto *currUnforcedOrbPtr(S nPart) const				{ return m_ppUnforcedColOrbCurr[nPart]; }
	CK inline void setCurrUnforcedOrbPtr(S nRow, S nPart)			{ m_ppUnforcedColOrbCurr[nPart] = pntr2UnforcedColOrb(nRow, nPart); }
	CC inline void setCurrentRowNumb(S n)							{ m_nCurrRow = n; }
	CC inline auto *unforcedColOrbPntr(S idxPart = 0) const			{ return m_ppUnforcedColOrb[idxPart]; }
	CC inline void setRowMaster(S val)								{ m_nRowMaster = val; }
	CC inline auto rowMaster() const								{ return m_nRowMaster; }
private:
	CK inline void setColOrbit(ColOrbPntr pntr, S idx, S idxPart)	{ m_ppColOrb[idxPart][idx] = pntr; }
	CC inline void setColOrbitLen(size_t len)						{ m_nColOrbLen = len; }

	S m_nCurrRow;
	uint m_nShiftMult;
    S m_nCol;
	ColOrbPntr **m_ppColOrb;
	ColOrbPntr **m_ppUnforcedColOrb;
	ColOrbPntr **m_ppUnforcedColOrbCurr;
	ColOrbPntr *m_ppOrb;
	ColOrbPntr **m_ppColOrbIni;
	size_t m_nColOrbLen;
	S m_nRowMaster = 0;
	bool m_IS_enumerator;
	int m_Memoryidx = -1;      // when >= 0, it's index where allocated memory is stored
};

FClass1(CColOrbitManager, void)::InitiateColOrbitManager(uint matrRank, S nRows, S nCol, S nParts, void *pMem)
{
	m_ppOrb = NULL;
	m_nCol = nCol;
	m_nShiftMult = matrRank * nRows;
	m_ppColOrb = new ColOrbPntr *[4 * nParts];
	m_ppUnforcedColOrbCurr = (m_ppUnforcedColOrb = (m_ppColOrbIni = m_ppColOrb + nParts) + nParts) + nParts;
	const auto nElem = nParts * nRows;
	if (!pMem) {
		// We need two sets of pointers for each part of Block Design 
		// one sets is accessed by m_pColOrb, the other one by m_pColOrbIni
		m_ppColOrb[0] = new ColOrbPntr[2 * nElem];
		const auto len = nElem * rank();
		m_ppUnforcedColOrb[0] = new ColOrbPntr[len];
		memset(m_ppUnforcedColOrb[0], 0, len * sizeof(*m_ppUnforcedColOrb[0]));
	} else {
		// Using previously allocated memory
		m_ppColOrb[0] = (ColOrbPntr *)pMem;
		memset(m_ppUnforcedColOrb, 0, nParts * sizeof(m_ppUnforcedColOrb[0]));
	}
	
	m_ppColOrbIni[0] = m_ppColOrb[0] + nRows;
	for (S i = 1; i < nParts; i++) {
		m_ppColOrb[i] = m_ppColOrbIni[i - 1] + nRows;
		m_ppColOrbIni[i] = m_ppColOrb[i] + nRows;
		m_ppUnforcedColOrb[i] = m_ppUnforcedColOrb[i - 1] + nRows * rank();
	}
}

FClass1(CColOrbitManager, void)::ReleaseColOrbitManager()
{
	if (unforcedColOrbPntr()) {
		delete[] unforcedColOrbPntr();
		delete[] m_ppColOrb[0];
	}

	delete[] m_ppColOrb;
}

FClass1(CColOrbitManager, void)::initiateColOrbits(S nRows, S firstRow, const Class1(BlockGroupDescr) *pGroupDescr, bool using_IS_enumerator,
												   int use_master_solutions, const Class1(CColOrbitManager) *pMaster, void *pMem)
{
	m_IS_enumerator = using_IS_enumerator;
	// Number of CColOrbits taken from pMaster
	if (!pMaster || !rowMaster())
	    setRowMaster(pMaster ? pMaster->currentRowNumb() + use_master_solutions : 0);

	// In order not to re-create orbits on the way back, we will store them in different places.
	// Therefore, we must have number of rows as a multiplier for fromMaster and nCol_2.
	const auto nCol_2 = nRows * colNumb();
	auto numParts = pGroupDescr ? pGroupDescr->numParts() : 1;
#ifndef USE_CUDA		// NOT yet implemented for GPU
	const auto lenColOrbitElement = using_IS_enumerator? sizeof(Class1(CColOrbitIS)) : sizeof(Class1(CColOrbitCS));
#else
	const auto lenColOrbitElement = sizeof(CColOrbitIS<S>);
#endif
	setColOrbitLen(lenColOrbitElement);

	if (pMem) {
		m_ppOrb = (ColOrbPntr *)pMem;
		pMem = (char *)pMem + nCol_2 * sizeof(ColOrbPntr);
	}
	else
		m_ppOrb = new ColOrbPntr[nCol_2];

	const int maxElement = rank();
	const auto fromMaster = WAIT_THREADS ? (rowMaster() + 1 - use_master_solutions) * colNumb() : 0;
	if (!using_IS_enumerator) {
#ifndef USE_CUDA		// NOT yet implemented for GPU
		auto pColOrbitsCS = pMem ? (Class1(CColOrbitCS) *)pMem : new Class1(CColOrbitCS)[nCol_2 - fromMaster];
		for (auto i = nCol_2; i-- > fromMaster;) {
			pColOrbitsCS[i].InitOrbit(maxElement - 1);
			m_ppOrb[i] = pColOrbitsCS + i - fromMaster;
		}
#endif
	}
	else {
		auto pColOrbitsIS = pMem? (Class1(CColOrbitIS) *)pMem : new Class1(CColOrbitIS)[nCol_2 - fromMaster];
		for (auto i = nCol_2; i-- > fromMaster;)
			m_ppOrb[i] = pColOrbitsIS + i - fromMaster;
	}

	m_Memoryidx = pMem ? -1 : fromMaster;

	const auto iMin = WAIT_THREADS ? rowMaster() + 1 - use_master_solutions : 0;
	const auto iMax = rowMaster() + 1 - use_master_solutions;
	for (S j = 0; j < numParts; j++) {
		const auto shift = pGroupDescr? pGroupDescr->getShift(j) : 0;
		auto i = nRows;
		while (i-- > iMin)
			m_ppColOrbIni[j][i] = m_ppColOrb[j][i] = m_ppOrb[i * colNumb() + shift];
	}

	if (rowMaster() > 0) {
#if WAIT_THREADS
		const size_t len = rowMaster() * sizeof(m_ppColOrbIni[0][0]);
		memcpy(m_ppColOrbIni[0], pMaster->colOrbitsIni(), len);
		memcpy(m_pColOrb, pMaster->colOrbits(), len);
		memcpy(unforcedColOrbPntr(0), pMaster->unforcedColOrbPntr(0), rowMaster() * 2 * sizeof(*unforcedColOrbPntr(0)));
#else
		// When we are not waiting for all threads to finish on rowMaster() level, 
		// we need to reconstruct the information regarding master's column orbits here 
		// because it could be changed when master will continue its calculation
		ColOrbPntr pColOrb;
		for (S iPart = 0; iPart < numParts; iPart++) {
			for (S i = 0; i < iMax; i++) {
				const auto* pOrbIni = pMaster->colOrbitsIni(iPart)[i];
				const auto* pOrb = pMaster->colOrbits(iPart)[i];
				m_ppColOrbIni[iPart][i]->clone(pOrbIni);
				pColOrb = m_ppColOrb[iPart][i] = (ColOrbPntr)((char*)m_ppColOrb[iPart][i] + ((char*)pOrb - (char*)pOrbIni));
				pColOrb->clone(pOrb);

				const size_t idx = i * maxElement;
				for (int j = maxElement; j--;) {
					const auto* pUnforcedOrb = pMaster->unforcedColOrbPntr(iPart)[idx+j];
					if (!pUnforcedOrb)
						continue;

					pColOrb = unforcedColOrbPntr(iPart)[idx+j] = (ColOrbPntr)((char*)colOrbitsIni(iPart)[i] + ((char*)pUnforcedOrb - (char*)pOrbIni));
					pColOrb->clone(pUnforcedOrb);
				}
			}
		}
#endif
	}
	else {
		m_ppColOrb[0][0]->Init(colNumb());
		if (numParts > 1) {
			// Initiating the leading column orbits of all block
			S len(0);
			while (numParts--) {
				const auto shift = pGroupDescr->GetPartInfo(numParts, &len);
				m_ppColOrb[numParts][firstRow]->Init(len);
			}
		}
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
		const auto iMax = *pColOrbInfo++;
		if (iMax == UINT64_MAX) {
			colOrbits()[nRow] = NULL;
			continue;
		}

		const auto pColOrbitIniTo = (char *)colOrbitsIni()[nRow];
		auto pColOrbitToNext = (ColOrbPntr)(pColOrbitIniTo + *pColOrbInfo++);
		colOrbits()[nRow] = pColOrbitToNext;
		for (int i = 1; i <= iMax; i++) {
			auto pColOrbitTo = pColOrbitToNext;
			const auto len = *pColOrbInfo++;
			pColOrbitToNext = i < iMax ? (ColOrbPntr)(pColOrbitIniTo + *pColOrbInfo++) : NULL;
			pColOrbitTo->Init(len, pColOrbitToNext);
		}
	}
}

FClass1(CColOrbitManager, void)::closeColOrbits()
{
	if (colOrbitPntr() && m_Memoryidx >= 0) {
		delete[] colOrbitPntr()[m_Memoryidx];
		delete[] colOrbitPntr();
		m_ppOrb = NULL;
	}
}

FClass1(CColOrbitManager, void)::addForciblyConstructedColOrbit(ColOrbPntr pColOrbit, S nPart, S n)
{
	auto *ppTmp = currUnforcedOrbPtr(nPart) + n;
	pColOrbit->setNext(*ppTmp);
#if TEST
	// Check if the loop was created. It helped to investigate some bug.
	if (*ppTmp && (*ppTmp)->next() == pColOrbit)
		*ppTmp = pColOrbit;
#endif
	*ppTmp = pColOrbit;
}

#endif /* defined(__BIBD_Mac__OrbitManager__) */
