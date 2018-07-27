//
//  RowSolution.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 1/30/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "RowSolution.h"
#include "CanonicityChecker.h"
#include "InSysSolver.h"
#include "ColOrbits.h"

template class CRowSolution<MATRIX_ELEMENT_TYPE>;

#if !USE_THREADS && !MY_QUICK_SORT
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

template<class T>
void CRowSolution<T>::removeNoncanonicalSolutions(size_t startIndex) const
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

template<class T>
size_t CRowSolution<T>::findSolution(const T *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const
{
	const uchar *pCanonFlags = pSolPerm->canonFlags();
	const size_t len = solutionSize() * sizeof(*pSolution);
	while (++i < iMax && MEMCMP(pSolution, firstSolution() + pSolPerm->GetAt(i) * solutionSize(), len)) {
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

template<class T>
size_t CRowSolution<T>::moveNoncanonicalSolutions(const T *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage, size_t *pSolIdx)
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
		return (size_t)-1;
	}

	return startIndex;
}

template<class T>
void CRowSolution<T>::printRow(FILE *file, PERMUT_ELEMENT_TYPE *pPerm, const T *pSol) const
{
    char buffer[512], *pBuf = buffer;
    for (size_t j = 0; j < numSolutions(); j++) {
        size_t idx = pPerm? *(pPerm + j) : j;
        if (pSol)
            idx = *(pSol + idx * solutionSize());
        
        pBuf += sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "%2zd", idx % 100);
        if (pBuf >= buffer + sizeof(buffer) - 2)
			outString(pBuf = buffer, file);
    }
    
	sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "\n");
	outString(buffer, file);
}

template<class T>
size_t CRowSolution<T>::setSolutionFlags(char *buffer, size_t lenBuf, size_t solIdx) const
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

template<class T>
void CRowSolution<T>::printSolutions(FILE *file, bool markNextUsed) const
{
	if (!solutionSize() || !numSolutions())
        return;
    
	MUTEX_LOCK(out_mutex);
	char buffer[2048];
	if (numSolutions() >= sizeof(buffer) / 2) {
        if (markNextUsed) {
			SPRINTF(buffer, "Using solution # %d out of %d\n", solutionIndex(), numSolutions());
			outString(buffer, file);
		}
	}
	else {
		if (markNextUsed) {
			const size_t len2 = sizeof(buffer) << 1;
			const size_t len = solutionIndex() * 2;
			const size_t nLoops = len / len2;
			size_t idx = 0;
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

