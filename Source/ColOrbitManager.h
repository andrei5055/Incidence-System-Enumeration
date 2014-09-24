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

class CColOrbit;

class CColOrbitManager
{
public:
    CColOrbitManager(int rank, size_t nRows,size_t nCol);
    ~CColOrbitManager();
	inline size_t colOrbitLen() const					{ return m_nColOrbLen; }
	inline CColOrbit **colOrbitPntr() const				{ return m_ppOrb; }
	inline const size_t colNumb() const					{ return m_nCol; }
	inline CColOrbit **colOrbits() const				{ return m_pColOrb; }
	inline FILE *outFile() const                        { return m_pFile; }
    inline size_t currentRowNumb() const                { return m_nCurrRow; }
	inline CColOrbit **colOrbitsIni() const				{ return m_pColOrbIni; }
protected:
	inline CColOrbit *colOrbitIni(size_t nRow) const	{ return *(colOrbitsIni() + nRow); }
	virtual CColOrbit **unforcedOrbits(size_t nRow)	const { return NULL; }
    size_t initiateColOrbits(size_t nRow, bool using_IS_enumerator, const CColOrbitManager *pMaster = NULL);
    void closeColOrbits();
    inline void setColOrbitCurr(CColOrbit *pntr)        { m_pColOrb[currentRowNumb()] = pntr; }
    inline void resetColOrbitCurr()                     { setColOrbitCurr(*(colOrbitPntr() + colNumb() * currentRowNumb())); }
	void resetUnforcedColOrb();
	virtual void resetFirstUnforcedRow()				{}
    void addForciblyConstructedColOrbit(CColOrbit *pColOrbit, int n);
	inline CColOrbit **currUnforcedOrbPtr() const       { return m_ppUnforcedColOrbCurr; }
	inline void setCurrUnforcedOrbPtr(size_t nRow)	    { m_ppUnforcedColOrbCurr = unforcedColOrbPntr() + unfColIdx(nRow); }
    inline void setCurrentRowNumb(size_t nRow)          { m_nCurrRow = nRow; }
	inline CColOrbit **unforcedColOrbPntr() const		{ return m_ppUnforcedColOrb; }
	inline CColOrbit *colOrbit(size_t idx) const		{ return m_pColOrb[idx]; }
    inline int rank() const                             { return m_nRank; }
	inline void setOutFile(FILE *file)                  { m_pFile = file; }
	inline size_t rowMaster() const						{ return m_nRowMaster; }
private:
	inline size_t unfColIdx(size_t r, int idx = 0) const{ return r * rank() + idx; }
	inline void setColOrbitLen(size_t len)				{ m_nColOrbLen = len; }
	inline void setRowMaster(size_t val)				{ m_nRowMaster = val; }

    size_t m_nCurrRow;
    const int m_nRank;
    const size_t m_nCol;
	CColOrbit **m_ppUnforcedColOrb;
	CColOrbit **m_ppUnforcedColOrbCurr;
    CColOrbit **m_pColOrb;
    CColOrbit **m_ppOrb;
	CColOrbit **m_pColOrbIni;
	size_t m_nColOrbLen;
	size_t m_nRowMaster;
	FILE *m_pFile;
};

#endif /* defined(__BIBD_Mac__OrbitManager__) */
