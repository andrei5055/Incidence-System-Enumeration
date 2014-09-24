//
//  RowSolution.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 1/30/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "matrix.h"
#include "RowSolution.h"
#include "CanonicityChecker.h"
#include "InSysSolver.h"
#ifdef _MSC_VER
#ifndef _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#endif
#include <crtdbg.h>
#endif
#if defined(_MSC_VER) && defined(_DEBUG)
#define new new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#define USE_PERM 1

CRowSolution::CRowSolution(size_t size, unsigned int nVect, CArrayOfVectorElements *pCoord) : CVector(size * nVect, pCoord)
{
	setSolutionPerm(NULL);
	InitSolutions(size, nVect, pCoord);
}

CRowSolution::~CRowSolution()
{
    delete solutionPerm();
}

void CRowSolution::InitSolutions(size_t size, unsigned int nVect, CArrayOfVectorElements *pCoord)
{
	setSolutionIndex(0);
	setSolutionSize(size);					// vector length
	delete solutionPerm();
	setSolutionPerm(new CSolutionPerm());
	setNumSolutions(nVect);
}

VECTOR_ELEMENT_TYPE *CRowSolution::newSolution()
{
	const auto numSolution = solutionIndex() + 1;
	if (GetSize() < numSolution * solutionSize())
		IncreaseVectorSize(solutionSize());

    VECTOR_ELEMENT_TYPE *pntr = (VECTOR_ELEMENT_TYPE *)currSolution();
    nextSolutionIndex();
    return pntr;
}

VECTOR_ELEMENT_TYPE *CRowSolution::copySolution(const CInSysSolver *pSysSolver)
{
    VECTOR_ELEMENT_TYPE *pSolution = (VECTOR_ELEMENT_TYPE *)lastSolution();
    if (!pSysSolver->isValidSolution(pSolution))
        return pSolution;
        
    pSolution = newSolution();
    memcpy(pSolution, pSolution - solutionSize(), solutionSize() * sizeof(pSolution[0]));
    return pSolution;
}

CRowSolution *CRowSolution::getSolution()
{
    if (!solutionIndex())
        return NULL;
    
    setNumSolutions(solutionIndex() - 1);
    // Reset the index of current solution
    setSolutionIndex(-1);
    return this;
}

#if USE_THREADS || MY_QUICK_SORT
int CRowSolution::compareVectors(const PERMUT_ELEMENT_TYPE idx, const VECTOR_ELEMENT_TYPE *pSecnd) const
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

void CRowSolution::quickSort(PERMUT_ELEMENT_TYPE *arr, long left, long right) const
{
	long i = left, j = right;
	const long pivotIdx = (left + right) >> 1;
	const VECTOR_ELEMENT_TYPE *pivot = firstSolution() + arr[pivotIdx] * solutionSize();

	/* partition */
	while (i <= j) {
		while (i != pivotIdx && compareVectors(arr[i], pivot) == -1)
			i++;

		while (j != pivotIdx && compareVectors(arr[j], pivot) == 1)
			j--;

		if (i <= j) {
			const PERMUT_ELEMENT_TYPE tmp = arr[i];
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

#else
#if USE_PERM
const VECTOR_ELEMENT_TYPE *pntrSolution;
#endif
size_t sizeSolution;

int compareSolutions (const void *p1, const void *p2)
{
#if USE_PERM
    const VECTOR_ELEMENT_TYPE *pFirst = pntrSolution + *(PERMUT_ELEMENT_TYPE *)p1 * sizeSolution;
    const VECTOR_ELEMENT_TYPE *pSecnd = pntrSolution + *(PERMUT_ELEMENT_TYPE *)p2 * sizeSolution;
#else
    const VECTOR_ELEMENT_TYPE *pFirst = (VECTOR_ELEMENT_TYPE *)p1;
    const VECTOR_ELEMENT_TYPE *pSecnd = (VECTOR_ELEMENT_TYPE *)p2;
#endif
    for (size_t i = 0; i < sizeSolution; i++) {
        if (*(pFirst+i) > *(pSecnd+i))
            return 1;
        
        if (*(pFirst+i) < *(pSecnd+i))
            return -1;
    }

    return 0;
}
#endif

void CRowSolution::sortSolutions(CCanonicityChecker *pCanonChecker)
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
		quickSort(pPerm, 0, numSolutions() - 1);
#else
		sizeSolution = solutionSize();
		pntrSolution = firstSolution();
		qsort(pPerm, numSolutions(), sizeof(pPerm[0]), compareSolutions);
#endif
	}

	if (pCanonChecker)
		sortSolutionByGroup(pCanonChecker);
	else
		memset(pCanonFlags, 1, numSolutions());
}

void CRowSolution::resetSolution()
{
    setSolutionIndex(0);
	solutionPerm()->RemoveAll();
}

CRowSolution *CRowSolution::NextSolution(bool useCanonGroup)
{ 
	uchar *pCanonFlags = useCanonGroup? solutionPerm()->canonFlags() : NULL;
	if (!pCanonFlags)
		return nextSolutionIndex() < numSolutions()? this : NULL;

	size_t solIndex;
	while ((solIndex = nextSolutionIndex()) < numSolutions()) {
		if (*(pCanonFlags + solIndex) == 1)
			return this;
	}

	return NULL;
}

void CRowSolution::removeNoncanonicalSolutions(size_t startIndex) const
{
	const CSolutionPerm *pSolPerm = solutionPerm();
	uchar *pCanonFlags = pSolPerm->canonFlags();
	// The forward or backward loop over all solutions isomorphic to the just tested solutions 
	if (startIndex > 0) {
		size_t solIndex = startIndex;
		do {
			*(pCanonFlags + solIndex) = 0xff;  // invalidate solution as a right part
		} while (!*(pCanonFlags + ++solIndex));
		*(pCanonFlags + solIndex) = 0xff;
	}
	else {
		size_t solIndex = solutionIndex();
		do {
			*(pCanonFlags + solIndex) = 0xff;  // invalidate solution as a right part
		} while (solIndex-- > 0 && !*(pCanonFlags + solIndex));
	}
}

size_t CRowSolution::findSolution(const VECTOR_ELEMENT_TYPE *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const
{
	const uchar *pCanonFlags = pSolPerm->canonFlags();
	const size_t len = solutionSize() * sizeof(*pSolution);
	while (++i < iMax && memcmp(pSolution, firstSolution() + pSolPerm->GetAt(i) * solutionSize(), len)) {
		if (*(pCanonFlags + i) == 1) {
			if (lastCanonIdx == SIZE_MAX)
				*pNextSolutionIdx = i;

			lastCanonIdx = i;
		}
	}

	return i;
}

size_t findStartingIndex(size_t from, const uchar *pCanonFlags)
{
	while (from-- > 0 && !*(pCanonFlags + from));
	return from;
}

size_t findCanonIndex(size_t idx, const uchar *pCanonFlags)
{
	while (!*(pCanonFlags + idx)) idx++;
	return idx;
}

size_t CRowSolution::moveNoncanonicalSolutions(const VECTOR_ELEMENT_TYPE *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage, size_t *pSolIdx)
{
	// For improved solution find the canonical solution:
	CSolutionPerm *pSolPerm = solutionPerm();
	uchar *pCanonFlags = pSolPerm->canonFlags();
	PERMUT_ELEMENT_TYPE *pPerm = pSolPerm->GetData();
    size_t len, nextSolutionIdx, lastCanonIdx = SIZE_MAX;
	size_t from;
	bool found = false;
	size_t i, curCanon = solutionIndex();
	const bool firstCall = !pSolutionStorage || pSolutionStorage->empty();
	if (!firstCall) {
		// Current solution was moved
		// Since we keep solutionIndex() as the previous index, we need to use "++curCanon" here
		from = findStartingIndex(++curCanon, pCanonFlags);
		i = findSolution(pSolution, from, startIndex, pSolPerm, lastCanonIdx, &nextSolutionIdx);
		found = i < startIndex;
		if (!found) { 
			// Solution was not found BEFORE the new interval of the current solution group
			// Let's try to find it AFTER this interval
			if (!lenOrbitOfSolution()) {
				// The length of the current solution group was not defined yet
				// When we are here, it could happend only when this group was merged with the next one
				// during previous call of current function
				curCanon--;
				setLenOrbitOfSolution(curCanon - findStartingIndex(curCanon, pCanonFlags));
				
				// In that case startIndex is not define. At this point startIndex should be 0 and
				//  solutionIndex() was not changed when we was here for thr first time
				*(pCanonFlags + curCanon) = 0;
				startIndex = curCanon - lenOrbitOfSolution() + 1;
			} else
				curCanon = startIndex + lenOrbitOfSolution() - 1;

			lastCanonIdx = SIZE_MAX;
		}
	}

	if (!found) {
		// We did not found the improved solution
		i = findSolution(pSolution, curCanon, numSolutions(), pSolPerm, lastCanonIdx, &nextSolutionIdx);
	}

	if (i < numSolutions() && pCanonFlags[i] != 0xff) {
		if (pSolutionStorage)
			pSolutionStorage->push_back(firstSolution() + pSolPerm->GetAt(i) * solutionSize());

		// We found solution and it was not marked as deleted
		bool changeSolutionIndex;
		if (!(changeSolutionIndex = firstCall)) {
			if (found) {
				// The whole orbit needs to be moved after startIndex + lenOrbitOfSolution position
				if (lastCanonIdx == SIZE_MAX) {
					// The orbit to be moved is the current one
					changeSolutionIndex = true;
					nextSolutionIdx = findCanonIndex(curCanon + 1, pCanonFlags);
				}
				else {
					// Going forward till next canonical solution
					curCanon = findCanonIndex(i, pCanonFlags);
				}

				lastCanonIdx = startIndex + lenOrbitOfSolution() - 1;
			}
			else {
				if (lastCanonIdx == SIZE_MAX)   // startIndex already in that orbit
					return startIndex;            // nothing to do

				curCanon = findCanonIndex(startIndex + lenOrbitOfSolution(), pCanonFlags);
				if (lastCanonIdx == curCanon) {
					// Two consequitive orbits need to be merged
#if PRINT_SOLUTIONS
					*(pCanonFlags + curCanon) = 0;
#endif
					return startIndex;
				}
			}
		}

		// Move all isomorphic solutions of the current solution into the new group:
		if (lastCanonIdx != SIZE_MAX) {
			// Find index of the first solution from current solution group
            from = findStartingIndex(curCanon, pCanonFlags);
			len = curCanon - from++;

			PERMUT_ELEMENT_TYPE buffer[256], *pTmp = buffer;
			if (len > countof(buffer))
				pTmp = new PERMUT_ELEMENT_TYPE[len];

			startIndex = lastCanonIdx - len + 1;
			lastCanonIdx -= curCanon++;

			// Moving lastCanonIdx elements starting from curCanon
			// Move the permutation's part
			memcpy(pTmp, pPerm + from, len * sizeof(pPerm[0]));
			memcpy(pPerm + from, pPerm + curCanon, lastCanonIdx * sizeof(pPerm[0]));
			memcpy(pPerm + startIndex, pTmp, len * sizeof(pPerm[0]));
			
			// Move the canonical flags's part
			memcpy(pCanonFlags + from, pCanonFlags + curCanon, lastCanonIdx * sizeof(pCanonFlags[0]));
			memset(pCanonFlags + startIndex, 0, len * sizeof(pCanonFlags[0]));
			
			if (firstCall)
				setLenOrbitOfSolution(len);
			else
			if (i < startIndex)		// When we add orbit which is LESS than the canonical solution of tested 
				startIndex -= len;	// solution associated on previous step, we need to adjust startIndex		

			// If our solution index was moved: 
			if (pSolIdx && *pSolIdx >= curCanon  && *pSolIdx < curCanon + lastCanonIdx)
				*pSolIdx -= len;  // adjust solution index

			// Move current solution index to the solution which precedes next canonical one
			if (changeSolutionIndex)
				setSolutionIndex(nextSolutionIdx - len - 1);

			if (pTmp != buffer)
				delete[] pTmp;
		}
	} else {
		// New group for current solution and its group was not found OR
		// this group is equivalent to already deleted solution
		// It means that it was removed (not costructed) on previous levels
		// and we don't need to use this group
		removeNoncanonicalSolutions(startIndex);
		return -1;
	}

	return startIndex;
}

bool CRowSolution::isValidSolution(size_t idx) const
{
	const uchar *pCanonFlags = solutionPerm()->canonFlags();
	return pCanonFlags? *(pCanonFlags + idx) != 0xff : true;
}

bool CRowSolution::findFirstValidSolution(const VECTOR_ELEMENT_TYPE *pMax, const VECTOR_ELEMENT_TYPE *pMin)
{
#if USE_PERM
    const VECTOR_ELEMENT_TYPE *pFirst = firstSolution();
	const size_t nSolutions = numSolutions();
    CSolutionPerm *pPerm = solutionPerm();
    size_t n, idx = 0;
    for (auto i = solutionSize(); i--;) {
        // Current min and max values we have to test
        const VECTOR_ELEMENT_TYPE minVal = pMin? *(pMin + i) : 1;
        const VECTOR_ELEMENT_TYPE maxVal = *(pMax + i);    
		const VECTOR_ELEMENT_TYPE currentVal = *(pFirst + pPerm->GetAt(idx) * solutionSize() + i);
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
        } else {
            
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
    const VECTOR_ELEMENT_TYPE *pFirst = firstSolution();
    const VECTOR_ELEMENT_TYPE *pLast = pFirst + solutionSize() * numSolutions();
    const VECTOR_ELEMENT_TYPE *pCurrentPoz = pFirst + solutionSize();
    for (int i = solutionSize(); i--;) {
        pCurrentPoz--;
        const VECTOR_ELEMENT_TYPE *pTmp;
        
        // Current min and max values we have to test
        const VECTOR_ELEMENT_TYPE minVal = pMin? *(pMin + i) : 1;
        const VECTOR_ELEMENT_TYPE maxVal = *(pMax + i);
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
        } else {
            
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

bool CRowSolution::checkChoosenSolution(const CColOrbit *pColOrbit, size_t nRowToBuild, size_t kMin)
{
	size_t idx = solutionIndex() + 1;
	for (uint i = 0; i < solutionSize(); i++, pColOrbit = pColOrbit->next()) {
		auto minVal = pColOrbit->lenght();
		const size_t limitWeight = minVal * (kMin - pColOrbit->colomnWeight());
		// To start, we need to set the solution index to -1
		setSolutionIndex(-1);
		do {
			while (minVal && nextSolutionIndex() <= idx) {
				const VECTOR_ELEMENT_TYPE curVal = *(currSolution() + i);
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

void CRowSolution::printRow(FILE *file, PERMUT_ELEMENT_TYPE *pPerm, const VECTOR_ELEMENT_TYPE *pSol) const
{
    char buffer[512], *pBuf = buffer;
    for (size_t j = 0; j < numSolutions(); j++) {
        size_t idx = pPerm? *(pPerm + j) : j;
        if (pSol)
            idx = *(pSol + idx * solutionSize());
        
        pBuf += sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "%2lu", idx % 100);
        if (pBuf >= buffer + sizeof(buffer) - 2)
			outString(pBuf = buffer, file);
    }
    
	sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "\n");
	outString(buffer, file);
}

int CRowSolution::setSolutionFlags(char *buffer, size_t lenBuf, int solIdx) const
{
	memset(buffer, ' ', lenBuf);
	const uchar *pCanonFlags = solutionPerm()->canonFlags();
	if (!pCanonFlags)
		return solIdx += lenBuf / 2;

	for (size_t i = 0; i < lenBuf; i += 2, solIdx++) {
		if (*(pCanonFlags + solIdx) == 1)
			buffer[i + 1] = '!';
		else
		if (!isValidSolution(solIdx))
			buffer[i + 1] = '-';
	}

	return solIdx;
}

void CRowSolution::printSolutions(FILE *file, bool markNextUsed) const
{
	if (!solutionSize() || !numSolutions())
        return;
    
	MUTEX_LOCK(out_mutex);
	char buffer[2048];
	if (numSolutions() >= sizeof(buffer) / 2) {
        if (markNextUsed) {
			SPRINTF(buffer, "Using solution # %lu out of %lu\n", solutionIndex(), numSolutions());
			outString(buffer, file);
		}
	}
	else {
		if (markNextUsed) {
			const size_t len2 = sizeof(buffer) << 1;
			const size_t len = solutionIndex() * 2;
			const size_t nLoops = len / len2;
			int idx = 0;
			for (size_t j = 0; j < nLoops; j++) {
				idx = setSolutionFlags(buffer, sizeof(buffer), idx);
				buffer[sizeof(buffer)-1] = '\0';
				outString(buffer, file);
			}

			const size_t lastLen = 2 * (numSolutions() - idx);
			setSolutionFlags(buffer, lastLen, idx);
			buffer[len % len2 + 1] = '*';
			strcpy_s(buffer + lastLen, sizeof(buffer)-lastLen, "\n");
			outString(buffer, file);
		}

		printRow(file);

		PERMUT_ELEMENT_TYPE *pPerm = solutionPerm() ? solutionPerm()->GetData() : NULL;
		if (pPerm)
			printRow(file, pPerm);

		const VECTOR_ELEMENT_TYPE *pSolution = firstSolution();
		for (unsigned int i = 0; i < solutionSize(); i++)
			printRow(file, pPerm, pSolution + i);
	}

	MUTEX_UNLOCK(out_mutex);
}

void CRowSolution::sortSolutionByGroup(CCanonicityChecker *pCanonChecker)
{
	pCanonChecker->constructGroup();
	// Since we removed the indeces corresponding to the forcibly constructed colOrbits 
	// from the generators of the group, the order of just constructe group could be less than |Aut(D)|

	uchar *pCanonFlags;
	PERMUT_ELEMENT_TYPE *pPerm = initSorting(&pCanonFlags);
	// Suppose that all solutions are not canonical
	memset(pCanonFlags, 0, numSolutions() * sizeof(pCanonFlags[0]));

	VECTOR_ELEMENT_TYPE buffer[256];
	const size_t lenMem = solutionSize() << 1;
	VECTOR_ELEMENT_TYPE *pMem = lenMem <= countof(buffer) ? buffer : new VECTOR_ELEMENT_TYPE [lenMem];
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
		} else
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
		delete [] pMem;

	if (pCanonIdx != canonIdx)
		delete[] pCanonIdx;
}

PERMUT_ELEMENT_TYPE *CSolutionPerm::initSorting(size_t num, uchar **pntr)
{
	SetSize(num); 
	canonFlgs()->SetSize(num);
	if (pntr)
		*pntr = canonFlags();

	return GetData();
}
