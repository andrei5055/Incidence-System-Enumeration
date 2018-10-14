//
//  RowSolution.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 1/30/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__RowSolution__
#define __BIBD_Mac__RowSolution__

#include "InSysRowEquation.h"
#include "ColOrbits.h"
//class CMatrix;
//class CColOrbit;
template<class T> class CCanonicityChecker;
template<class T> class CInSysSolver;

typedef CArray<uchar, uchar> CArrayOfCanonFlags;

class CSolutionPerm : public CArraySolutionPerm
{
public:
	CK CSolutionPerm()									{ m_CanonFlgs = new CArrayOfCanonFlags(); }
	CK ~CSolutionPerm()									{ delete canonFlgs(); }
	CK PERMUT_ELEMENT_TYPE *initSorting(size_t num, uchar **pntr = NULL) {
															SetSize(num);
															canonFlgs()->SetSize(num);
															if (pntr)
																*pntr = canonFlags();

															return GetData();
														}
	CK inline uchar *canonFlags() const					{ return canonFlgs()->GetData(); }
private:
	CK inline CArrayOfCanonFlags *canonFlgs() const		{ return m_CanonFlgs; }
public:
	CArrayOfCanonFlags *m_CanonFlgs;
};

template<class T>
class CRowSolution : public CVector<T>
{
public:
 	CK CRowSolution(size_t size = 0, uint nVect = 1, CArrayOfVectorElements *pCoord = NULL) : CVector<T>(size * nVect, pCoord) {
																	setSolutionPerm(NULL);
																	InitSolutions(size, nVect, pCoord);
																}
    CK ~CRowSolution()											{ delete solutionPerm(); }
    CK inline const T *firstSolution() const					{ return this->GetData(); }
	CK inline const T *solution(PERMUT_ELEMENT_TYPE i) const	{ return firstSolution() + variantIndex(i) * solutionSize(); }
	CK inline const T *currSolution() const						{ return firstSolution() + variantIndex() * solutionSize(); }
	CK inline PERMUT_ELEMENT_TYPE nextSolutionIndex()			{ return ++m_nSolutionIndex; }
    CK inline void prevSolutionIndex()							{ --m_nSolutionIndex; }
    CK T *newSolution();
    CK T *copySolution(const CInSysSolver<T> *pSysSolver);
	CK CRowSolution *NextSolution(bool useCanonGroup = false);
	void removeNoncanonicalSolutions(size_t startIndex) const;
	size_t moveNoncanonicalSolutions(const T *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage = NULL, size_t *pSolIdx = NULL);
	CK inline bool isValidSolution(size_t idx) const			{ const uchar *pCanonFlags = solutionPerm()->canonFlags();
																  return pCanonFlags ? *(pCanonFlags + idx) != 0xff : true; }
	CK void InitSolutions(size_t size = 0, uint nVect = 1, CArrayOfVectorElements *pCoord = NULL);
	CK inline size_t solutionSize() const						{ return m_Size; }
	CK inline void setSolutionSize(size_t size)					{ m_Size = size; }
    CK CRowSolution *getSolution();
    CK bool findFirstValidSolution(const T *pMax, const T *pMin = NULL);
	CK bool checkChoosenSolution(const CColOrbit<T> *pColOrbit, size_t nRowToBuild, size_t kMin);
	CK void sortSolutions(CCanonicityChecker<T> *pCanonChecker = NULL);
    void printSolutions(FILE *file, bool markNextUsed = false) const;
	CK inline PERMUT_ELEMENT_TYPE solutionIndex() const			{ return m_nSolutionIndex; }
	CK inline void setSolutionIndex(PERMUT_ELEMENT_TYPE val)	{ m_nSolutionIndex = val; }
	CK inline CSolutionPerm *solutionPerm() const				{ return m_pSolutionPerm; }
    CK void resetSolution()										{ setSolutionIndex(0); solutionPerm()->RemoveAll(); }
	CK inline size_t numSolutions() const						{ return m_nNumSolutions; }
	inline bool isLastSolution() const							{ return solutionIndex() + 1 == numSolutions(); }
	inline void setLenOrbitOfSolution(size_t len)               { m_nLenSolOrb = len; }
private:
	CK void sortSolutionByGroup(CCanonicityChecker<T> *pCanonChecker);
	CK inline void setSolutionPerm(CSolutionPerm *perm)			{ m_pSolutionPerm = perm; }
	CK inline void setNumSolutions(size_t val)					{ m_nNumSolutions = val; }
	CK inline PERMUT_ELEMENT_TYPE variantIndex() const			{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(solutionIndex()) : solutionIndex(); }
	CK inline PERMUT_ELEMENT_TYPE variantIndex(PERMUT_ELEMENT_TYPE i) const	{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(i) : i; }
    void printRow(FILE *file = NULL, PERMUT_ELEMENT_TYPE *pPerm = NULL, const T *pSol = NULL) const;
	size_t setSolutionFlags(char *buffer, size_t lenBuf, size_t solIdx) const;
	CK inline PERMUT_ELEMENT_TYPE *initSorting(uchar **pntr = NULL){ return solutionPerm()->initSorting(numSolutions(), pntr); }
	size_t findSolution(const T *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const;
	inline size_t lenOrbitOfSolution() const                    { return m_nLenSolOrb; }
    CK inline T *lastSolution() const							{ return (T *)currSolution() - solutionSize(); }
#if USE_THREADS || MY_QUICK_SORT
	CK void quickSort(PERMUT_ELEMENT_TYPE *arr, long left, long right) const;
	CK int compareVectors(const PERMUT_ELEMENT_TYPE idx, const VECTOR_ELEMENT_TYPE *pSecnd) const;
#endif

	size_t m_Size;
	size_t m_nNumSolutions;
	PERMUT_ELEMENT_TYPE m_nSolutionIndex;
	CSolutionPerm *m_pSolutionPerm;
	size_t m_nLenSolOrb;
};

#define USE_PERM    1   // Should be 1. Version for 0 has a bug

template<class T>
void CRowSolution<T>::InitSolutions(size_t size, unsigned int nVect, CArrayOfVectorElements *pCoord)
{
	setSolutionIndex(0);
	setSolutionSize(size);					// vector length
	delete solutionPerm();
	setSolutionPerm(new CSolutionPerm());
	setNumSolutions(nVect);
}

template<class T>
T *CRowSolution<T>::newSolution()
{
	const auto numSolution = solutionIndex() + 1;
	if (this->GetSize() < numSolution * solutionSize())
		this->IncreaseVectorSize(solutionSize());

	auto pntr = currSolution();
	nextSolutionIndex();
	return (T *)pntr;
}

template<class T>
CRowSolution<T> *CRowSolution<T>::getSolution()
{
	if (!solutionIndex())
		return NULL;

	setNumSolutions(solutionIndex() - 1);
	// Reset the index of current solution
	setSolutionIndex((PERMUT_ELEMENT_TYPE)-1);
	return this;
}

template<class T>
T *CRowSolution<T>::copySolution(const CInSysSolver<T> *pSysSolver)
{
	auto pSolution = lastSolution();
	if (!pSysSolver->isValidSolution(pSolution))
		return pSolution;

	pSolution = newSolution();
	memcpy(pSolution, pSolution - solutionSize(), solutionSize() * sizeof(pSolution[0]));
	return pSolution;
}

template<class T>
CRowSolution<T> *CRowSolution<T>::NextSolution(bool useCanonGroup)
{
	uchar *pCanonFlags = useCanonGroup ? solutionPerm()->canonFlags() : NULL;
	if (!pCanonFlags)
		return nextSolutionIndex() < numSolutions() ? this : NULL;

	size_t solIndex;
	while ((solIndex = nextSolutionIndex()) < numSolutions()) {
		if (*(pCanonFlags + solIndex) == 1)
			return this;
	}

	return NULL;
}

template<class T>
bool CRowSolution<T>::findFirstValidSolution(const T *pMax, const T *pMin)
{
#if USE_PERM
	const auto pFirst = firstSolution();
	const size_t nSolutions = numSolutions();
	CSolutionPerm *pPerm = solutionPerm();
	size_t n, idx = 0;
	for (auto i = solutionSize(); i--;) {
		// Current min and max values we have to test
		const T minVal = pMin ? *(pMin + i) : 1;
		const T maxVal = *(pMax + i);
		const T currentVal = *(pFirst + pPerm->GetAt(idx) * solutionSize() + i);
		if (currentVal >= minVal) {
			// No need to check minVal anymore
			// We need to make sure that before pCurrentPoz there is at least one value, which is less than maxVal
			if (currentVal < maxVal)
				continue;

			// If we are here, then *pCurrentPoz == maxVal
			n = idx;
			while (n-- > 0 && *(pFirst + pPerm->GetAt(n) * solutionSize() + i) == maxVal);

			if (n != SIZE_MAX)
				continue;       // no need to change first valid solution

								// Need to go to the "bigger" solutions and find first which would be < maxVal
			n = idx;
			while (++n < nSolutions && *(pFirst + pPerm->GetAt(n) * solutionSize() + i) == maxVal);
		}
		else {

			n = idx;
			while (n-- > 0 && *(pFirst + pPerm->GetAt(n) * solutionSize() + i) < minVal);

			if (n != SIZE_MAX)
				continue;       // no need to change first valid solution

								// If we are here, then *pCurrentPoz == maxVal
								// Need to go to "bigger" solutions and find first which would be >= minVal
			n = idx;
			while (++n < nSolutions && *(pFirst + pPerm->GetAt(n) * solutionSize() + i) < minVal);
		}

		if (n >= nSolutions) {
			// Found coordinate, with all values < minVal OR == maxVal
			return false;
		}

		idx = n;
	}

	setSolutionIndex(idx - 1);
#else
	const T *pFirst = firstSolution();
	const T *pLast = pFirst + solutionSize() * numSolutions();
	const T *pCurrentPoz = pFirst + solutionSize();
	for (auto i = solutionSize(); i--;) {
		pCurrentPoz--;
		const T *pTmp;

		// Current min and max values we have to test
		const T minVal = pMin ? *(pMin + i) : 1;
		const T maxVal = *(pMax + i);
		if (*pCurrentPoz >= minVal) {
			// No need to check minVal anymore
			// We need to make sure that before pCurrentPoz there is at least one value, which is less than maxVal
			if (*pCurrentPoz < maxVal)
				continue;

			// If we are here, then *pCurrentPoz == maxVal
			pTmp = pCurrentPoz;
			while ((pTmp -= solutionSize()) >= pFirst && *pTmp == maxVal);

			if (pTmp >= pFirst)
				continue;       // no need to change first valid solution

								// Need to go to the "bigger" solutions and find first which would be < maxVal
			pTmp = pCurrentPoz;
			while ((pTmp += solutionSize()) < pLast && *pTmp == maxVal);
		}
		else {

			pTmp = pCurrentPoz;
			while ((pTmp -= solutionSize()) >= pFirst && *pTmp < minVal);

			if (pTmp >= pFirst)
				continue;       // no need to change first valid solution

								// If we are here, then *pCurrentPoz == maxVal
								// Need to go to "bigger" solutions and find first which would be >= minVal
			pTmp = pCurrentPoz;
			while ((pTmp += solutionSize()) < pLast && *pTmp < minVal);
		}

		if (pTmp >= pLast) {
			// Found coordinate, with all values < minVal OR == maxVal
			return false;
		}

		pCurrentPoz = pTmp;
	}

	setSolutionIndex(uint(pCurrentPoz - pFirst + 1) / solutionSize() - 1);
#endif

	return true;
}

template<class T>
void CRowSolution<T>::sortSolutions(CCanonicityChecker<T> *pCanonChecker)
{
	if (!this || !numSolutions() || !solutionSize())
		return;

	uchar *pCanonFlags;
	PERMUT_ELEMENT_TYPE *pPerm = initSorting(&pCanonFlags);
	for (PERMUT_ELEMENT_TYPE i = 0; i < numSolutions(); i++)
		*(pPerm + i) = i;

	if (numSolutions() > 1) {
#if USE_THREADS || MY_QUICK_SORT
		// When we use threads we cannot use qsort since in our implementation 
		// qsort will use global variables - pntrSolution and sizeSolution
		quickSort(pPerm, 0, static_cast<long>(numSolutions() - 1));
#else
		extern size_t sizeSolution;
		extern const T* pntrSolution;
		sizeSolution = solutionSize();
		pntrSolution = firstSolution();
		int compareSolutions(const void *p1, const void *p2);
		qsort(pPerm, numSolutions(), sizeof(pPerm[0]), compareSolutions);
#endif
	}

	if (pCanonChecker)
		sortSolutionByGroup(pCanonChecker);
	else
		memset(pCanonFlags, 1, numSolutions());
}

#if USE_THREADS || MY_QUICK_SORT
template<class T>
int CRowSolution<T>::compareVectors(const PERMUT_ELEMENT_TYPE idx, const VECTOR_ELEMENT_TYPE *pSecnd) const
{
	const VECTOR_ELEMENT_TYPE *pFirst = firstSolution() + idx * solutionSize();
	for (size_t i = 0; i < solutionSize(); i++) {
		if (*(pFirst + i) > *(pSecnd + i))
			return 1;

		if (*(pFirst + i) < *(pSecnd + i))
			return -1;
	}

	return 0;
}

template<class T>
void CRowSolution<T>::quickSort(PERMUT_ELEMENT_TYPE *arr, long left, long right) const
{
	long i = left, j = right;
	const auto pivotIdx = (left + right) >> 1;
	const VECTOR_ELEMENT_TYPE *pivot = firstSolution() + arr[pivotIdx] * solutionSize();

	/* partition */
	while (i <= j) {
		while (i != pivotIdx && compareVectors(arr[i], pivot) == -1)
			i++;

		while (j != pivotIdx && compareVectors(arr[j], pivot) == 1)
			j--;

		if (i <= j) {
			const auto tmp = arr[i];
			arr[i++] = arr[j];
			arr[j--] = tmp;
		}
	}

	/* recursion */
	if (left < j)
		quickSort(arr, left, j);

	if (i < right)
		quickSort(arr, i, right);
}
#endif


template<class T>
void CRowSolution<T>::sortSolutionByGroup(CCanonicityChecker<T> *pCanonChecker)
{
	pCanonChecker->constructGroup();
	// Since we removed the indices corresponding to the forcibly constructed colOrbits 
	// from the generators of the group, the order of just constructe group could be less than |Aut(D)|

	uchar *pCanonFlags;
	PERMUT_ELEMENT_TYPE *pPerm = initSorting(&pCanonFlags);
	// Suppose that all solutions are not canonical
	memset(pCanonFlags, 0, numSolutions() * sizeof(pCanonFlags[0]));

	VECTOR_ELEMENT_TYPE buffer[256];
	const size_t lenMem = solutionSize() << 1;
	VECTOR_ELEMENT_TYPE *pMem = lenMem <= countof(buffer) ? buffer : new VECTOR_ELEMENT_TYPE[lenMem];
	size_t canonIdx[256];
	size_t *pCanonIdx = numSolutions() <= countof(canonIdx) ? canonIdx : new size_t[numSolutions()];
	int nCanon = 0;

	// Loop for all solutions, starting from the biggest ones
	const VECTOR_ELEMENT_TYPE *pFirst = firstSolution();
	for (PERMUT_ELEMENT_TYPE i = numSolutions(); i--;) {
		int nCanonPrev = nCanon;
		const PERMUT_ELEMENT_TYPE idxPerm = *(pPerm + i);
		const size_t idx = pCanonChecker->findSolutionIndex(pFirst, idxPerm, pMem, pCanonIdx, nCanon);
		if (idxPerm == idx) {
			if (nCanonPrev != nCanon)
				*(pCanonFlags + i) = 1; // The solution we just processed is the canonical one

			continue;
		}

		// Current solution is NOT the canonical one			
		// and it's canonical is NOT the smallest one
		PERMUT_ELEMENT_TYPE last;
		if (idx != SIZE_MAX) {
			// We need to find idx in pPerm array
			last = i;
			while (*(pPerm + last + 1) != idx)
				last++;
		}
		else
			setNumSolutions(last = numSolutions() - 1);

		for (PERMUT_ELEMENT_TYPE j = i; j < last; j++) {
			*(pPerm + j) = *(pPerm + j + 1);
			*(pCanonFlags + j) = *(pCanonFlags + j + 1);
		}

		if (idx != SIZE_MAX) {
			*(pPerm + last) = idxPerm;
			*(pCanonFlags + last) = 0;
		}
	}

	if (pMem != buffer)
		delete[] pMem;

	if (pCanonIdx != canonIdx)
		delete[] pCanonIdx;
}

template<class T>
bool CRowSolution<T>::checkChoosenSolution(const CColOrbit<T> *pColOrbit, size_t nRowToBuild, size_t kMin)
{
	size_t idx = solutionIndex() + 1;
	for (uint i = 0; i < solutionSize(); i++, pColOrbit = pColOrbit->next()) {
		auto minVal = pColOrbit->length();
		const size_t limitWeight = minVal * (kMin - pColOrbit->columnWeight());
		// To start, we need to set the solution index to -1
		setSolutionIndex((PERMUT_ELEMENT_TYPE)-1);
		do {
			while (minVal && nextSolutionIndex() <= idx) {
				const auto curVal = *(currSolution() + i);
				if (minVal > curVal)
					minVal = curVal;
			}

			if (!minVal || minVal * nRowToBuild <= limitWeight)
				break;

			prevSolutionIndex();
		} while (++idx < numSolutions());

		if (idx == numSolutions())
			return false;
	}

	setSolutionIndex(idx - 1);
	assert(idx < m_nNumSolutions + 1);
	return true;
}

#endif /* defined(__BIBD_Mac__RowSolution__) */
