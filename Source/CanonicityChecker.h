//
//  CanonicityChecker.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__CanonicityChecker__
#define __BIBD_Mac__CanonicityChecker__
#pragma once

#include "DataTypes.h"
#include "PermutStorage.h"

typedef CSimpleArray<size_t>   CPermut;
typedef CContainer<size_t>     CColNumbStorage;

class CEnumerator;

typedef enum {
	t_saveRowPermutations	= 1,
	t_saveRowToChange		= (1 << 1)
} t_canonOutInfo;

class CCanonicityChecker : CPermut
{
public:
    CCanonicityChecker(size_t nRow, size_t nCol, int rank);
    ~CCanonicityChecker();
	bool TestCanonicity(size_t nRowMax, const CEnumerator *pEnum, int outInfo = 0, size_t *pRowOut = NULL, CRowSolution *pRowSolution = NULL);
    void outputAutomorphismInfo(FILE *file) const;
	inline uint groupOrder() const					{ return m_nGroupOrder; }
	inline size_t numColOrb() const					{ return m_nNumColOrb; }
	inline void setGroupOrder(uint val)				{ m_nGroupOrder = val; }
	inline size_t lenPerm() const					{ return permStorage()->lenPerm(); }
	inline void adjustGenerators(int *idx, size_t len) { permStorage()->adjustGenerators(idx, len); }
	inline size_t constructGroup()					{ return permStorage()->constructGroup(); }
	inline size_t findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon)
													{ return permStorage()->findSolutionIndex(pFirst, idx, pMem, pCanonIdx, nCanon); }
	inline CPermutStorage *permStorage() const		{ return m_pPermutStorage; }
	inline VECTOR_ELEMENT_TYPE *improvedSolution() const { return m_pImprovedSol; }
	bool groupIsTransitive() const;
private:
    void init(size_t nRow, bool savePerm);
    size_t next_permutation(size_t idx = -1);
    void addAutomorphism(bool rowPermut = true);
    int checkColOrbit(size_t orbLen, size_t nColCurr, const MATRIX_ELEMENT_TYPE *pRow, const MATRIX_ELEMENT_TYPE *pRowPerm) const;
    inline size_t *permRow() const                  { return m_pPermutRow->elementPntr(); }
    inline size_t *permCol() const                  { return m_pPermutCol->elementPntr(); }
    inline void setStabilizerLength(size_t len)     { m_nStabLength = len; }
    inline size_t stabilizerLength() const          { return m_nStabLength; }
    inline void setStabilizerLengthAut(size_t len)  { m_nStabLengthAut = len; }
    inline size_t stabilizerLengthAut() const       { return m_nStabLengthAut; }
	inline void setNumRow(size_t nRow)				{ m_nNumRow = nRow; }
	inline size_t numRow() const					{ return m_nNumRow; }
    inline size_t numCol() const                    { return m_pPermutCol->numElement(); }
#define orbits()	elementPntr()
	void updateGroupOrder();
    inline CColNumbStorage **colNumbStorage() const { return m_nColNumbStorage; }
    inline const int rank() const                   { return m_rank; }
    inline CCounter<int> *counter() const           { return m_pCounter; }
	inline void setColIndex(size_t *pntr)			{ m_pColIndex = pntr; }
	inline size_t *colIndex() const					{ return m_pColIndex; }
	inline void setNumColOrb(size_t val)			{ m_nNumColOrb = val; }
    void revert(size_t i);
	size_t getLenPermutCol(size_t **permCol);
	size_t *constructColIndex(const CColOrbit *pColOrbit, const CColOrbit *pColOrbitIni, size_t colOrbLen);
	inline void setPermStorage(CPermutStorage *p)	{ m_pPermutStorage = p; }
	void rowToChange(size_t *pRowOut, size_t nRow) const;
	void reconstructSolution(const CColOrbit *pColOrbitStart, const CColOrbit *pColOrbit, size_t colOrbLen, const CColOrbit *pColOrbitIni, const MATRIX_ELEMENT_TYPE *pRowPerm, const VECTOR_ELEMENT_TYPE *pRowSolution, size_t solutionSize);
#if USE_STRONG_CANONICITY
	inline void setSolutionStorage(CSolutionStorage *p) { m_pSolutionStorage = p; }
	inline CSolutionStorage *solutionStorage() const { return m_pSolutionStorage; }

	CSolutionStorage *m_pSolutionStorage;
#else
#define setSolutionStorage(x)
#define solutionStorage()		(CSolutionStorage *)NULL
#endif

    size_t m_nStabLength;
    size_t m_nStabLengthAut;
    const int m_rank;
    CPermut *m_pPermutRow;
    CPermut *m_pPermutCol;
    CColNumbStorage **m_nColNumbStorage;
	CPermutStorage *m_pPermutStorage;
    CCounter<int> *m_pCounter;
	size_t *m_pColIndex;
	VECTOR_ELEMENT_TYPE *m_pImprovedSol;

	uint m_nGroupOrder;
	size_t m_nNumRow;
	size_t m_nNumColOrb;
};

#define initPermStorage()			permStorage()->initPermutStorage()
#define savePerm(x, y)				permStorage()->savePermut(x, y)
#define outputPerm(x, y, z)			permStorage()->outputPerm(x, y, z)
#define outputPermutations(x, y)	permStorage()->outputPermutations(x, y)

#endif /* defined(__BIBD_Mac__CanonicityChecker__) */
