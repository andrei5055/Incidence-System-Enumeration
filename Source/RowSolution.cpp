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

template class CRowSolution<TDATA_TYPES>;

#if !USE_THREADS && !USE_MY_QUICK_SORT
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

FClass2(CRowSolution, void)::removeNoncanonicalSolutions(size_t startIndex) {
	auto *pCanonFlags = solutionPerm()->canonFlags();
	// The forward or backward loop over all solutions isomorphic to the just tested solutions 
	if (startIndex > 0) {
		auto solIndex = startIndex;
#if HARD_REMOVE
		const auto iDst = startIndex;
		while (*(pCanonFlags + ++solIndex) != t_canon_solution);
		const size_t nElem = numSolutions() - solIndex;
		moveOrbitSolution(solutionPerm()->GetData(), pCanonFlags, iDst, solIndex, nElem);
		// Total number of sollutions and index of current should be adjusted
		setNumSolutions(iDst + nElem);
		setSolutionIndex(iDst - 1);
#else
		do {
			*(pCanonFlags + solIndex) = t_invalid_as_righ_part;  // invalidate solution as a right part
		} while (!*(pCanonFlags + ++solIndex));
		*(pCanonFlags + solIndex) = t_invalid_as_righ_part;
#endif
	}
	else {
		auto solIndex = solutionIndex();
#if HARD_REMOVE
		const auto iSrc = solIndex + 1;
		while (solIndex-- > 0 && *(pCanonFlags + solIndex) != t_canon_solution);
		const size_t nElem = numSolutions() - iSrc;
		moveOrbitSolution(solutionPerm()->GetData(), pCanonFlags, solIndex + 1, iSrc, nElem);
		// Total number of sollutions and index of current should be adjusted
		setNumSolutions(solIndex + nElem + 1);
		setSolutionIndex(solIndex);
#else
		do {
			*(pCanonFlags + solIndex) = t_invalid_as_righ_part;  // invalidate solution as a right part
		} while (solIndex-- > 0 && !*(pCanonFlags + solIndex));
#endif
	}
}

FClass2(CRowSolution, size_t)::findSolution(const T *pSolution, size_t i, size_t iMax, const CSolutionPerm *pSolPerm, size_t &lastCanonIdx, size_t *pNextSolutionIdx) const
{
	const auto *pCanonFlags = pSolPerm->canonFlags();
	const auto len = solutionLength() * sizeof(*pSolution);
	while (++i < iMax && MEMCMP(pSolution, firstSolution() + pSolPerm->GetAt(i) * solutionLength(), len)) {
		if (*(pCanonFlags + i) == t_canon_solution) {
			if (lastCanonIdx == SIZE_MAX)
				*pNextSolutionIdx = i;

			lastCanonIdx = i;
		}
	}

	return i;
}

static size_t findStartingIndex(size_t from, const uchar *pCanonFlags)
{
	while (from-- > 0 && !*(pCanonFlags + from));
	return from;
}

static size_t findCanonIndex(size_t idx, const uchar *pCanonFlags)
{
	while (!*(pCanonFlags + idx)) idx++;
	return idx;
}

void moveOrbitSolution(PERMUT_ELEMENT_TYPE *pPerm, uchar *pCanonFlags, size_t iDst, size_t iSrc, size_t nElem, size_t len=0, size_t iDstNew=0)
{
	if (!len) {
		// It is NOT necessary to store information addressed by the destination index.
		memcpy(pPerm + iDst, pPerm + iSrc, nElem * sizeof(*pPerm));
		memcpy(pCanonFlags + iDst, pCanonFlags + iSrc, nElem * sizeof(*pCanonFlags));
		return;
	}

	// Allocate memory for temporay storage of information addressed by the destination index
	PERMUT_ELEMENT_TYPE buffer[256], *pTmp = len <= countof(buffer) ? buffer : new PERMUT_ELEMENT_TYPE[len];

	size_t length = len * sizeof(*pPerm);
	memcpy(pTmp, pPerm + iDst, length);			// Moving the permutation part out of the way.
	memcpy(pPerm + iDst, pPerm + iSrc, nElem * sizeof(*pPerm));
	memcpy(pPerm + iDstNew, pTmp, length);       // Moving the permutation's part on its new plase

	// Moving the canonical flags's parts
	length = len * sizeof(pCanonFlags[0]);
	memcpy(pTmp, pCanonFlags + iDst, length);
	memcpy(pCanonFlags + iDst, pCanonFlags + iSrc, nElem * sizeof(*pCanonFlags));
	memcpy(pCanonFlags + iDstNew, pTmp, length);
	if (pTmp != buffer)
		delete[] pTmp;
}

FClass2(CRowSolution, size_t)::moveNoncanonicalSolutions(const T *pSolution, size_t startIndex, CSolutionStorage *pSolutionStorage, size_t *pSolIdx)
{
	// For improved solution find the canonical solution:
	auto *pSolPerm = solutionPerm();
	auto *pCanonFlags = pSolPerm->canonFlags();
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
				//  solutionIndex() was not changed when we was here for the first time
				*(pCanonFlags + curCanon) = t_not_canon_solution;
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

	if (i < numSolutions()
#if !HARD_REMOVE
		&& pCanonFlags[i] != t_invalid_as_righ_part
#endif
		) {
		if (pSolutionStorage)
			pSolutionStorage->push_back(firstSolution() + pSolPerm->GetAt(i) * solutionLength());

		// We found solution and it was not marked as deleted
		bool changeSolutionIndex = firstCall;
		if (!changeSolutionIndex) {
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
					*(pCanonFlags + curCanon) = t_formerCanon_solution;
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

			*(pCanonFlags + curCanon) = t_formerCanon_solution;
			startIndex = lastCanonIdx - len + 1;
			lastCanonIdx -= curCanon++;  // it's a number of elements to move

			moveOrbitSolution(pPerm, pCanonFlags, from, curCanon, lastCanonIdx, len, startIndex);
			if (firstCall)
				setLenOrbitOfSolution(len);
			else
			if (i < startIndex)		// When we add orbit which is LESS than the canonical solution of tested 
				startIndex -= len;	// solution associated on previous step, we need to adjust startIndex		

			// If our solution index was moved: 
			if (pSolIdx && *pSolIdx >= curCanon  && *pSolIdx < curCanon + lastCanonIdx)
				*pSolIdx -= len;  // adjust solution index

			// Move the index of the current solution to the solution that precedes the next canonical
			if (changeSolutionIndex)
				setSolutionIndex(nextSolutionIdx - len - 1);
		}
		else {
			*(pCanonFlags + curCanon) = t_formerCanon_solution;
		}
	} else {
		// New group for current solution and its group was not found OR
		// this group is equivalent to already deleted solution
		// It means that it was removed (not costructed) on previous levels
		// and we don't need to use this group
		removeNoncanonicalSolutions(startIndex);
//		void moveOrbitSolution(PERMUT_ELEMENT_TYPE* pPerm, uchar* pCanonFlags, size_t iDst, size_t iSrc, size_t nElem, size_t len, size_t iDstNew)

		return (size_t)-1;
	}

	return startIndex;
}

#if PRINT_SOLUTIONS
FClass2(CRowSolution, void)::printRow(FILE *file, PERMUT_ELEMENT_TYPE *pPerm, const S *pSol) const
{
    char buffer[512], *pBuf = buffer;
    for (size_t j = 0; j < numSolutions(); j++) {
        size_t idx = pPerm? *(pPerm + j) : j;
        if (pSol)
            idx = *(pSol + idx * solutionLength());
        
        pBuf += sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "%2zd", idx % 100);
        if (pBuf >= buffer + sizeof(buffer) - 2)
			outString(pBuf = buffer, file);
    }
    
	sprintf_s(pBuf, sizeof(buffer) - (pBuf - buffer), "\n");
	outString(buffer, file);
}

FClass2(CRowSolution, size_t)::setSolutionFlags(char *buffer, size_t lenBuf, size_t solIdx) const
{
	memset(buffer, ' ', lenBuf);
	const uchar *pCanonFlags = solutionPerm()->canonFlags();
	if (!pCanonFlags)
		return solIdx += lenBuf / 2;

	for (size_t i = 0; i < lenBuf; i += 2, solIdx++) {
		switch (*(pCanonFlags + solIdx)) {
		case t_canon_solution:		 buffer[i + 1] = '!'; break;
#if USING_FORMER_CANON_FLAG
		case t_formerCanon_solution: buffer[i + 1] = '#'; break;
#endif
#if !HARD_REMOVE
		default: if (!isValidSolution(solIdx))
					buffer[i + 1] = '-';
#endif
		}
	}

	return solIdx;
}

FClass2(CRowSolution, void)::printSolutions(FILE *file, bool markNextUsed, T nRow, T nPortion, bool addPortionNumb) const
{
	if (!solutionLength() || !numSolutions())
        return;
    
	char buffer[2048], *pBuf = buffer;
	const auto lenBuf = countof(buffer);
	if (markNextUsed)
		pBuf += SNPRINTF(pBuf, lenBuf, "\nRow #%2d: the solution # %zd out of %zd will be used", nRow, solutionIndex() + 1, numSolutions());
	else
		pBuf += SNPRINTF(pBuf, lenBuf, "\nRow #%2d: %zd solution%s constructed", nRow, numSolutions(), numSolutions() > 1? "s were" : " was");


	if (addPortionNumb || nPortion)
		SNPRINTF(pBuf, lenBuf - (pBuf - buffer), " for portion %d\n", nPortion);
	else
		SNPRINTF(pBuf, lenBuf - (pBuf - buffer), "\n");

	outString(buffer, file);

	if (numSolutions() >= sizeof(buffer) / 2)
		return;   // The buffer is not large enough for output

	if (markNextUsed) {
		const auto len2 = sizeof(buffer) << 1;
		const auto len = solutionIndex() << 1;
		const auto nLoops = len / len2;
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
	for (unsigned int i = 0; i < solutionLength(); i++)
		printRow(file, pPerm, pSolution + i);
}
#endif