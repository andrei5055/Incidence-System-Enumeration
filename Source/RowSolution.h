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
#include "Sorter.h"

Class2Def(CInSysSolver);
Class2Def(CPermutStorage);

typedef CArray<uchar, uchar> CArrayOfCanonFlags;


#define USING_FORMER_CANON_FLAG 0	// For now we are using t_formerCanon_solution as t_not_canon_solution
#define HARD_REMOVE	0				// Remove solution instead of marking them as t_invalid_as_righ_part

enum solutionType {
	t_not_canon_solution   = 0,
	t_canon_solution       = 1 << 0,
#if USING_FORMER_CANON_FLAG
	t_formerCanon_solution = 1 << 1,
#endif
#if !HARD_REMOVE
	t_invalid_as_righ_part = 0xff,
#endif
};

#if !USING_FORMER_CANON_FLAG
#define t_formerCanon_solution   t_not_canon_solution
#endif

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
	CSolutionPerm& operator = (const CSolutionPerm& src) {
		RemoveAll();
		Append(src);
		*m_CanonFlgs = *src.canonFlgs();
		return *this;
	}
private:
	CK inline CArrayOfCanonFlags *canonFlgs() const		{ return m_CanonFlgs; }
public:
	CArrayOfCanonFlags *m_CanonFlgs;
};

Class2Def(CRowSolution) : public CVector<S>, CSorter<unsigned char>
{
public:
	CK CRowSolution(T length = 0, size_t nVect = 1, CArrayOfVectorElements *pCoord = NULL) :
		                                                        CVector<S>(length * nVect, pCoord),
																CSorter(length * sizeof(T)) {
																	setSolutionPerm(NULL);
																	InitSolutions(length, nVect, pCoord);
																}
	CK ~CRowSolution()											{ delete solutionPerm(); }
	CRowSolution& operator = (const CRowSolution& src);
	CK inline const auto *firstSolution() const					{ return this->GetData(); }
	CK inline const auto *solution(PERMUT_ELEMENT_TYPE i) const	{ return firstSolution() + variantIndex(i) * solutionLength(); }
	CK inline const auto *currSolution() const					{ return firstSolution() + variantIndex() * solutionLength(); }
	CK inline PERMUT_ELEMENT_TYPE nextSolutionIndex()			{ return ++m_nSolutionIndex; }
	CK inline void prevSolutionIndex()							{ --m_nSolutionIndex; }
	CK inline auto allSolutionChecked()                         { return nextSolutionIndex() >= numSolutions(); }
	CK T *newSolution();
	CK T *copySolution(const InSysSolverPntr pSysSolver, T λ);
	CK CRowSolution *NextSolution(bool useCanonGroup = false);
	void removeNoncanonicalSolutions(size_t startIndex);
	size_t moveNoncanonicalSolutions(const T *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage = NULL, size_t *pSolIdx = NULL);
#if !HARD_REMOVE
	CK inline bool validSolution(size_t idx) const				{ const uchar *pCanonFlags = solutionPerm()->canonFlags();
																  return pCanonFlags ? *(pCanonFlags + idx) != t_invalid_as_righ_part : true; }
#endif
	CK void InitSolutions(T size = 0, size_t nVect = 1, CArrayOfVectorElements *pCoord = NULL, PERMUT_ELEMENT_TYPE lastIdx = 0);
	CK inline auto solutionLength() const						{ return m_Length; }
	CK inline void setSolutionLength(T length)					{ setRecordLength((m_Length = length)*sizeof(T)); }
	CK CRowSolution *getSolution();
	CK bool findFirstValidSolution(const T *pMax, const T *pMin = NULL);
	CK bool checkChoosenSolution(const CColOrbit<S> *pColOrbit, T nRowToBuild, T kMin);
	CK void sortSolutions(bool doSorting, PermutStoragePntr pPermStorage);
	CK inline auto numRemainingSolutions() const				{ return numSolutions() - solutionIndex(); }
#if PRINT_SOLUTIONS
	void printSolutions(FILE *file, bool markNextUsed, T nRow, T nPortion, bool addPortionNumb = false) const;
#endif
	CK inline auto solutionIndex() const						{ return m_nSolutionIndex; }
	CK inline void setSolutionIndex(PERMUT_ELEMENT_TYPE val)	{ m_nSolutionIndex = val; }
	CK inline void resetSolutionIndex()							{ setSolutionIndex(static_cast<PERMUT_ELEMENT_TYPE>(-1)); }
	CK inline auto *solutionPerm() const						{ return m_pSolutionPerm; }
	CK void resetSolution()										{ setSolutionIndex(0); solutionPerm()->RemoveAll(); }
	CK inline auto numSolutions() const							{ return m_nNumSolutions; }
	inline bool isLastSolution() const							{ return numRemainingSolutions() == 1; }
	inline void setLenOrbitOfSolution(size_t len)               { m_nLenSolOrb = len; }
	CK inline void saveSolutionIndex()							{ m_nSavedSolutionIndex = solutionIndex(); }
	CK inline void restoreSolutionIndex()						{ setSolutionIndex(m_nSavedSolutionIndex); }
	CK inline void makeDummySolution()							{ setNumSolutions(1);}
	CK inline PERMUT_ELEMENT_TYPE variantIndex() const			{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(solutionIndex()) : solutionIndex(); }
private:
	CK void sortSolutionsByGroup(PermutStoragePntr pPermutStorage);
	CK inline void setSolutionPerm(CSolutionPerm *perm)			{ m_pSolutionPerm = perm; }
	CK inline void setNumSolutions(size_t val)					{ m_nNumSolutions = val; }
	CK inline PERMUT_ELEMENT_TYPE variantIndex(PERMUT_ELEMENT_TYPE i) const	{ return solutionPerm()->GetData() ? solutionPerm()->GetAt(i) : i; }
#if PRINT_SOLUTIONS
	void printRow(FILE *file = NULL, PERMUT_ELEMENT_TYPE *pPerm = NULL, const S *pSol = NULL) const;
	void outPortionInfo(char *buffer, char* pBuf, size_t lenBuf, FILE* file, T nPortion, bool addPortionNumb) const;
	size_t setSolutionFlags(char *buffer, size_t lenBuf, size_t solIdx) const;
#endif
	CK inline PERMUT_ELEMENT_TYPE *initSorting(uchar **pntr = NULL){ return solutionPerm()->initSorting(numSolutions(), pntr); }
	size_t findSolution(const T *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const;
	inline auto lenOrbitOfSolution() const						{ return m_nLenSolOrb; }
	CK inline T *lastSolution() const							{ return (T *)currSolution() - solutionLength(); }

	T m_Length;
	size_t m_nNumSolutions;
	PERMUT_ELEMENT_TYPE m_nSolutionIndex;
	PERMUT_ELEMENT_TYPE m_nSavedSolutionIndex;
	CSolutionPerm *m_pSolutionPerm;
	size_t m_nLenSolOrb;
};

#define USE_PERM    1   // Should be 1. Version for 0 has a bug

FClass2(CRowSolution, CRowSolution<T,S>&)::operator = (const CRowSolution<T,S>& src) {
	// This function is used only for initiation of row solutions used by threads
	// We need to find maximum index of solution which is potentially could be used as a right part
	const auto iMax = src.solutionIndex();
	const auto* pPerm = src.solutionPerm()->GetData();
	// iMax could be zero for i-th (i > 1) part of CombBIBD
	// nSols is the number of solutions used for current thread
	const PERMUT_ELEMENT_TYPE nSols = iMax ? iMax + 1 : src.numSolutions();
	PERMUT_ELEMENT_TYPE i = nSols;
	PERMUT_ELEMENT_TYPE maxIdx = pPerm[--i];
	while (i--) {
		if (maxIdx < pPerm[i])
			maxIdx = pPerm[i];
	}

	InitSolutions(src.solutionLength(), nSols, src.getCoord(), maxIdx);
	*solutionPerm() = *src.solutionPerm();
	setSolutionIndex(iMax);
	return *this;
}

FClass2(CRowSolution, T *)::newSolution() {
	const auto numSolution = solutionIndex() + 1;
	if (this->GetSize() < numSolution * solutionLength())
		this->IncreaseVectorSize(solutionLength());

	auto pntr = currSolution();
	nextSolutionIndex();
	return (S *)pntr;
}

FClass2(CRowSolution, RowSolutionPntr)::getSolution() {
	if (!solutionIndex())
		return NULL;

	setNumSolutions(solutionIndex() - 1);
	// Reset the index of current solution
	resetSolutionIndex();
	return this;
}

FClass2(CRowSolution, T *)::copySolution(const InSysSolverPntr pSysSolver, T λ)
{
	auto *pSolution = lastSolution();
	if (pSysSolver->isValidSolution(pSolution, λ)) {
		pSolution = newSolution();
		memcpy(pSolution, pSolution - solutionLength(), solutionLength() * sizeof(pSolution[0]));
	}
	return pSolution;
}

FClass2(CRowSolution, RowSolutionPntr)::NextSolution(bool useCanonGroup) {
	uchar *pCanonFlags = useCanonGroup ? solutionPerm()->canonFlags() : NULL;
	if (!pCanonFlags)
		return allSolutionChecked() ? NULL : this;

	while (!allSolutionChecked()) {
		if (*(pCanonFlags + solutionIndex()) == 1)
			return this;
	}

	return NULL;
}

FClass2(CRowSolution, bool)::findFirstValidSolution(const T *pMax, const T *pMin) {
#if USE_PERM
	const auto pFirst = firstSolution();
	const auto nSolutions = numSolutions();
	auto *pPerm = solutionPerm();
	PERMUT_ELEMENT_TYPE n, idx = 0;
	for (auto i = solutionLength(); i--;) {
		// Current min and max values we have to test
		const auto minVal = pMin ? *(pMin + i) : 1;
		const S currentVal = *(pFirst + pPerm->GetAt(idx) * solutionLength() + i);
		if (currentVal >= minVal) {
			const auto maxVal = *(pMax + i);
			// No need to check minVal anymore
			// We need to make sure that before pCurrentPoz there is at least one value, which is less than maxVal
			if (currentVal < maxVal)
				continue;

			// If we are here, then *pCurrentPoz == maxVal
			n = idx;
			while (n-- > 0 && *(pFirst + pPerm->GetAt(n) * solutionLength() + i) == maxVal);

			if (n != PERMUT_ELEMENT_MAX)
				continue;       // no need to change first valid solution

								// Need to go to the "bigger" solutions and find first which would be < maxVal
			n = idx;
			while (++n < nSolutions && *(pFirst + pPerm->GetAt(n) * solutionLength() + i) == maxVal);
		}
		else {

			n = idx;
			while (n-- > 0 && *(pFirst + pPerm->GetAt(n) * solutionLength() + i) < minVal);

			if (n != PERMUT_ELEMENT_MAX)
				continue;       // no need to change first valid solution

								// If we are here, then *pCurrentPoz == maxVal
								// Need to go to "bigger" solutions and find first which would be >= minVal
			n = idx;
			while (++n < nSolutions && *(pFirst + pPerm->GetAt(n) * solutionLength() + i) < minVal);
		}

		if (n >= nSolutions) {
			// Found coordinate, with all values < minVal OR == maxVal
			return false;
		}

		idx = n;
	}

	setSolutionIndex(idx - 1);
#else
	const S *pFirst = firstSolution();
	const S *pLast = pFirst + solutionLength() * numSolutions();
	const S *pCurrentPoz = pFirst + solutionLength();
	for (auto i = solutionLength(); i--;) {
		pCurrentPoz--;
		const S *pTmp;

		// Current min and max values we have to test
		const auto minVal = pMin ? *(pMin + i) : 1;
		if (*pCurrentPoz >= minVal) {
			const auto maxVal = *(pMax + i);
			// No need to check minVal anymore
			// We need to make sure that before pCurrentPoz there is at least one value, which is less than maxVal
			if (*pCurrentPoz < maxVal)
				continue;

			// If we are here, then *pCurrentPoz == maxVal
			pTmp = pCurrentPoz;
			while ((pTmp -= solutionLength()) >= pFirst && *pTmp == maxVal);

			if (pTmp >= pFirst)
				continue;       // no need to change first valid solution

								// Need to go to the "bigger" solutions and find first which would be < maxVal
			pTmp = pCurrentPoz;
			while ((pTmp += solutionLength()) < pLast && *pTmp == maxVal);
		}
		else {

			pTmp = pCurrentPoz;
			while ((pTmp -= solutionLength()) >= pFirst && *pTmp < minVal);

			if (pTmp >= pFirst)
				continue;       // no need to change first valid solution

								// If we are here, then *pCurrentPoz == maxVal
								// Need to go to "bigger" solutions and find first which would be >= minVal
			pTmp = pCurrentPoz;
			while ((pTmp += solutionLength()) < pLast && *pTmp < minVal);
		}

		if (pTmp >= pLast) {
			// Found coordinate, with all values < minVal OR == maxVal
			return false;
		}

		pCurrentPoz = pTmp;
	}

	setSolutionIndex(uint(pCurrentPoz - pFirst + 1) / solutionLength() - 1);
#endif

	return true;
}

FClass2(CRowSolution, void)::sortSolutionsByGroup(PermutStoragePntr pPermStorage)
{
	pPermStorage->constructGroup();
	// Since we removed the indices corresponding to the forcibly constructed colOrbits 
	// from the generators of the group, the order of just constructed group could be less than |Aut(D)|
	uchar *pCanonFlags;
	auto *pPerm = initSorting(&pCanonFlags);
	// Suppose that all solutions are not canonical
	memset(pCanonFlags, 0, numSolutions() * sizeof(pCanonFlags[0]));

	S buffer[256];
	const auto lenMem = solutionLength() << 1;
	auto *pMem = lenMem <= countof(buffer) ? buffer : new S[lenMem];
	size_t canonIdx[256];
	auto pCanonIdx = numSolutions() <= countof(canonIdx) ? canonIdx : new size_t[numSolutions()];
	int nCanon = 0;

	// Loop for all solutions, starting from the biggest ones
	const auto *pFirst = firstSolution();
	for (auto i = numSolutions(); i--;) {
		const auto nCanonPrev = nCanon;
		const auto idxPerm = *(pPerm + i);
		const auto idx = pPermStorage->findSolutionIndex(pFirst, idxPerm, pMem, pCanonIdx, nCanon);
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

		for (auto j = i; j < last; j++) {
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

FClass2(CRowSolution, bool)::checkChoosenSolution(const CColOrbit<S> *pColOrbit, T nRowToBuild, T kMin) {
	size_t idx = solutionIndex() + 1;
	for (uint i = 0; i < solutionLength(); i++, pColOrbit = pColOrbit->next()) {
		size_t minVal = pColOrbit->length();
		const size_t limitWeight = minVal * (kMin - pColOrbit->columnWeight());
		// To start, we need to set the solution index to -1
		resetSolutionIndex();
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
