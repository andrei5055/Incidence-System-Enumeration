//
//  RowSolution.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 1/30/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__RowSolution__
#define __BIBD_Mac__RowSolution__

#include "Vector.h"

class CMatrix;
class CColOrbit;
class CCanonicityChecker;
class CInSysSolver;

typedef CArray<uchar, uchar> CArrayOfCanonFlags;

class CSolutionPerm : public CArraySolutionPerm
{
public:
	CSolutionPerm()										{ m_CanonFlgs = new CArrayOfCanonFlags(); }
	~CSolutionPerm()									{ delete canonFlgs(); }
	PERMUT_ELEMENT_TYPE *initSorting(size_t num, uchar **pntr = NULL);
	inline uchar *canonFlags() const					{ return canonFlgs()->GetData(); }
private:
	inline CArrayOfCanonFlags *canonFlgs() const		{ return m_CanonFlgs; }

	CArrayOfCanonFlags *m_CanonFlgs;
};

class CRowSolution : public CVector
{
public:
 	CRowSolution(size_t size = 0, uint nVect = 1, CArrayOfVectorElements *pCoord = NULL);
    ~CRowSolution();
    inline const VECTOR_ELEMENT_TYPE *firstSolution() const		{ return getCoord()->GetData(); }
	inline const VECTOR_ELEMENT_TYPE *solution(size_t i) const  { return firstSolution() + variantIndex(i) * solutionSize(); }
	inline const VECTOR_ELEMENT_TYPE *currSolution() const		{ return firstSolution() + variantIndex() * solutionSize(); }
	inline PERMUT_ELEMENT_TYPE nextSolutionIndex()				{ return ++m_nSolutionIndex; }
    inline void prevSolutionIndex()								{ --m_nSolutionIndex; }
    VECTOR_ELEMENT_TYPE *newSolution();
    VECTOR_ELEMENT_TYPE *copySolution(const CInSysSolver *pSysSolver);
	CRowSolution *NextSolution(bool useCanonGroup = false);
	void removeNoncanonicalSolutions(size_t startIndex) const;
	size_t moveNoncanonicalSolutions(const VECTOR_ELEMENT_TYPE *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage = NULL, size_t *pSolIdx = NULL);
	bool isValidSolution(size_t idx) const;
	void InitSolutions(size_t size = 0, uint nVect = 1, CArrayOfVectorElements *pCoord = NULL);
	inline size_t solutionSize() const							{ return m_Size; }
	inline void setSolutionSize(size_t size)					{ m_Size = size; }
    CRowSolution *getSolution();
    bool findFirstValidSolution(const VECTOR_ELEMENT_TYPE *pMax, const VECTOR_ELEMENT_TYPE *pMin = NULL);
	bool checkChoosenSolution(const CColOrbit *pColOrbit, size_t nRowToBuild, size_t kMin);
	void sortSolutions(CCanonicityChecker *pCanonChecker = NULL);
    void printSolutions(FILE *file, bool markNextUsed = false) const;
	inline PERMUT_ELEMENT_TYPE solutionIndex() const			{ return m_nSolutionIndex; }
	inline void setSolutionIndex(PERMUT_ELEMENT_TYPE val)		{ m_nSolutionIndex = val; }
	inline CSolutionPerm *solutionPerm() const					{ return m_pSolutionPerm; }
    void resetSolution();
	inline size_t numSolutions() const							{ return m_nNumSolutions; }
	inline bool isLastSolution() const							{ return solutionIndex() + 1 == numSolutions(); }
	inline void setLenOrbitOfSolution(size_t len)               { m_nLenSolOrb = len; }
private:
	void sortSolutionByGroup(CCanonicityChecker *pCanonChecker);
	inline void setSolutionPerm(CSolutionPerm *perm)			{ m_pSolutionPerm = perm; }
	inline void setNumSolutions(size_t val)						{ m_nNumSolutions = val; }
	inline PERMUT_ELEMENT_TYPE variantIndex() const				{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(solutionIndex()) : solutionIndex(); }
	inline PERMUT_ELEMENT_TYPE variantIndex(size_t i) const		{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(i) : i; }
    void printRow(FILE *file = NULL, PERMUT_ELEMENT_TYPE *pPerm = NULL, const VECTOR_ELEMENT_TYPE *pSol = NULL) const;
	int setSolutionFlags(char *buffer, size_t lenBuf, int solIdx) const;
	inline PERMUT_ELEMENT_TYPE *initSorting(uchar **pntr = NULL){ return solutionPerm()->initSorting(numSolutions(), pntr); }
	size_t findSolution(const VECTOR_ELEMENT_TYPE *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const;
	inline size_t lenOrbitOfSolution() const                    { return m_nLenSolOrb; }
    inline const VECTOR_ELEMENT_TYPE *lastSolution() const      { return currSolution() - solutionSize(); }
#if USE_THREADS || MY_QUICK_SORT
	void quickSort(PERMUT_ELEMENT_TYPE *arr, long left, long right) const;
	int compareVectors(const PERMUT_ELEMENT_TYPE idx, const VECTOR_ELEMENT_TYPE *pSecnd) const;
#endif

	size_t m_Size;
	size_t m_nNumSolutions;
	PERMUT_ELEMENT_TYPE m_nSolutionIndex;
	CSolutionPerm *m_pSolutionPerm;
	size_t m_nLenSolOrb;
};

#endif /* defined(__BIBD_Mac__RowSolution__) */
