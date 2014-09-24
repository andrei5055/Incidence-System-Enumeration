//
//  OrbitManager.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/22/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "ColOrbitManager.h"
#include "ColOrbits.h"

CColOrbitManager::CColOrbitManager(int rankMatr, size_t nRows, size_t nCol) : m_nRank(rankMatr), m_nCol(nCol)
{
    setCurrentRowNumb(0);
    m_pColOrb = new CColOrbit *[2 * nRows];
	m_pColOrbIni = m_pColOrb + nRows;
    const auto len = nRows * rank();
	m_ppUnforcedColOrb = new CColOrbit *[len];
	memset(unforcedColOrbPntr(), 0, len * sizeof(unforcedColOrbPntr()[0]));
}

CColOrbitManager::~CColOrbitManager()
{
    delete [] unforcedColOrbPntr();
    delete [] m_pColOrb;
}

size_t CColOrbitManager::initiateColOrbits(size_t nRow, bool using_IS_enumerator, const CColOrbitManager *pMaster)
{
	// Number of CColOrbits taken from pMaster
	setRowMaster(pMaster ? pMaster->currentRowNumb() + 1 : 0);
	const size_t fromMaster = WAIT_THREADS? rowMaster() * colNumb() : 0;
	const size_t nCol_2 = nRow * colNumb();
	m_ppOrb = new CColOrbit *[nCol_2];
	size_t lenColOrbitElement;
	if (!using_IS_enumerator) {
		CColOrbitCS *pColOrbitsCS = new CColOrbitCS[nCol_2 - fromMaster];
		lenColOrbitElement = sizeof(pColOrbitsCS[0]);
		const int maxElement = rank() - 1;
		CColOrbitCS::setMaxElement(maxElement);
		for (auto i = nCol_2; i-- > fromMaster;) {
			pColOrbitsCS[i].InitOrbit(maxElement);
			m_ppOrb[i] = pColOrbitsCS + i - fromMaster;
		}
	}
	else {
		CColOrbitIS *pColOrbitsIS = new CColOrbitIS[nCol_2 - fromMaster];
		lenColOrbitElement = sizeof(pColOrbitsIS[0]);
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
				m_pColOrb[i] = (CColOrbit *)((char *)m_pColOrb[i] + idx * lenColOrbitElement);

			m_pColOrb[i]->clone(pMaster->colOrbits()[i]);

			for (int j = 0; j < 2; j++) {
				const size_t idx = (i << 1) + j;
				if (pMaster->unforcedColOrbPntr()[idx]) {
					const size_t shift = (char *)pMaster->unforcedColOrbPntr()[idx] - (char *)pMaster->colOrbitsIni()[i];
					unforcedColOrbPntr()[idx] = (CColOrbit *)((char *)colOrbitsIni()[i] + shift);
					unforcedColOrbPntr()[idx]->clone(pMaster->unforcedColOrbPntr()[idx]);
				}
			}
		}
#endif
	} else
		m_pColOrb[0]->Init(colNumb());

	setColOrbitLen(lenColOrbitElement);	
    return lenColOrbitElement;
}

void CColOrbitManager::closeColOrbits()
{
	delete[] colOrbitsIni()[WAIT_THREADS? rowMaster() : 0];
	delete[] colOrbitPntr();
}

void CColOrbitManager::addForciblyConstructedColOrbit(CColOrbit *pColOrbit, int n)
{
    CColOrbit **ppTmp = currUnforcedOrbPtr() + n;
    pColOrbit->setNext(*ppTmp);
    *ppTmp = pColOrbit;
}

void CColOrbitManager::resetUnforcedColOrb()
{ 
	memset(unforcedColOrbPntr() + unfColIdx(currentRowNumb()), 0, rank() * sizeof(*unforcedColOrbPntr())); 
	resetFirstUnforcedRow();
}
