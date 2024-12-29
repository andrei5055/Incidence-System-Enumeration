#include "TripleSys.h"

#define UseIPX				0 // works faster with 0
#define USE_INTRINSIC		!USE_CUDA

#if USE_INTRINSIC || USE_64_BIT_MASK
#include <immintrin.h> // Header for AVX2 intrinsics

void bitwise_multiply(const long long* a, const long long* b, long long* result, size_t size) {
	size_t i = 0;

	// Process 4 `long long` elements at a time using AVX2
	for (; i + 4 <= size; i += 4) {
		// Load 4 `long long` elements from each array
		__m256i vec_a = _mm256_loadu_si256((__m256i*) & a[i]);
		__m256i vec_b = _mm256_loadu_si256((__m256i*) & b[i]);

		// Perform bitwise AND operation
		__m256i vec_result = _mm256_and_si256(vec_a, vec_b);

		// Store the result back to the result array
		_mm256_storeu_si256((__m256i*) & result[i], vec_result);
	}

	// Handle the remainder elements (if `size` is not a multiple of 4)
	for (; i < size; ++i) {
		result[i] = a[i] & b[i];
	}
}
#endif

CC bool CRowStorage::p1fCheck2(ctchar* neighborsi, ctchar* neighborsj) const {
	uint checked = 0;
	for (tchar m = 0; m < m_numPlayers; m++)
	{
		if (!(checked & (1 << m))) {
			tchar k = m;
			for (int i = 2; i <= m_numPlayers; i += 2)
			{
				if ((k = neighborsj[neighborsi[k]]) == m)
					return i == m_numPlayers;

				checked |= 1 << k;
			}
		}
	}
	return false;
}

CC bool CRowStorage::checkCompatibility(ctchar* neighborsi, const long long* rm, uint idx) const {
	// Let's check if the masks are mutually compatible
	auto* pMask = (const long long*)(m_pMaskStorage->getObject(idx));
	int j = m_lenMask;
	while (j-- && !(rm[j] & pMask[j]));

	return j < 0 && p1fCheck2(neighborsi, getObject(idx) + m_numPlayers);
}

CC bool CRowStorage::maskForCombinedSolutions(tmask* pMaskOut, uint & solIdx) const {
	// Constructing a mask for the "combined" solution, which represents 
	// the combination of solutions for the first two "non-predefined" rows.
	const auto n = m_numRec[1];
	do {
		const auto idx1 = solIdx / n;
		auto* rm1 = (const long long*)m_pMaskStorage->getObject(idx1);
		const auto pNeighbors1 = getObject(idx1) + m_numPlayers;
		if (checkCompatibility(pNeighbors1, rm1, solIdx % n + m_numRecAdj)) {
			memcpy(pMaskOut, getSolutionMask(solIdx % n), m_numSolutionTotalB);
			const auto pCompSol = pMaskOut - (m_numSolutionTotalB >>(SHIFT - 3));
			memcpy(pCompSol, pMaskOut, m_numSolutionTotalB);

			const auto last = m_numSolutionTotal - m_numRecAdj;
			auto first = m_numRecAdj2 - m_numRecAdj;
			const auto lastB = last >> 3;
			while (true) {
				// Skip all bytes/longs equal to 0
				auto firstB = first >> SHIFT;
				while (firstB < lastB && !pCompSol[firstB])
					firstB++;

				if (firstB >= lastB)
					break;
#if USE_64_BIT_MASK
				unsigned long iBit;
				_BitScanForward64(&iBit, pCompSol[firstB]);
#else
				const auto iBit = this->firstOnePosition(pCompSol[firstB]);
#endif
				if ((first = (firstB << SHIFT) + iBit) >= last)
					break;

				pCompSol[firstB] ^= (tmask)1 << iBit;
				if (!checkCompatibility(pNeighbors1, rm1, first + m_numRecAdj))
					RESET_MASK_BIT(pMaskOut, first);     //  reset bit
			}

			return true;
		}
	} while ((solIdx += m_step) < m_lastInFirstSet);

	return false;
}

CC void CRowStorage::generateCompatibilityMasks(tmask* pMaskOut, uint solIdx, uint idx) const {
	auto* rm = (const long long*)m_pMaskStorage->getObject(solIdx);
	const auto pNeighbors = getObject(solIdx) + m_numPlayers;
	do {
		if (checkCompatibility(pNeighbors, rm, idx)) {
			const auto newIdx = idx - m_numRecAdj;
			SET_MASK_BIT(pMaskOut, newIdx);     // 1 - means OK
		}
	} while (++idx < m_numSolutionTotal);
}

CC void CRowStorage::initCompatibilityMasks(ctchar* u1fCycles) {
	m_u1fCycles = u1fCycles;
	for (int i = 1; i < m_numPlayers; i++)
		m_pRowSolutionCntr[i] += m_pRowSolutionCntr[i - 1];

	// Define the number of first long long's we don't need to copy to the next row.
	memset(m_pNumLongs2Skip, 0, m_numPlayers * sizeof(m_pNumLongs2Skip[0]));
	int i = m_numPreconstructedRows;
	m_pNumLongs2Skip[i] = m_pRowSolutionCntr[i] >> 6;
	m_lastInFirstSet = m_numRecAdj = m_pRowSolutionCntr[i];
	if (NEW && sysParam()->val[t_useCombinedSolutions])
		m_lastInFirstSet *= (m_numRec[1] = ((m_numRecAdj2 = m_pRowSolutionCntr[i+1]) - m_numRecAdj));

	while (++i < m_numPlayers)
		m_pNumLongs2Skip[i] = ((m_pRowSolutionCntr[i] - m_numRecAdj) >> 6);

	m_numSolutionTotal = m_pRowSolutionCntr[m_numPlayers - 1];
	delete[] m_fullExcludeTable;
	m_numSolutionTotalB = ((m_numSolutionTotal - m_numRecAdj + 7) / 8 + 7) / 8 * 8;

	auto len = (m_numSolutionTotal - m_numRecAdj) * m_numSolutionTotalB;
	m_fullExcludeTable = new tchar[len];
	memset(m_fullExcludeTable, 0, len);

	delete[] m_pRowSolutionMasksIdx;
	delete[] m_pRowSolutionMasks;
	len = numPlayers();
	m_pRowSolutionMasksIdx = new uint[len];

	m_pRowSolutionMasks = new tmask[len];
	memset(m_pRowSolutionMasks, 0, len * sizeof(m_pRowSolutionMasks[0]));

#if !USE_64_BIT_MASK
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

	auto* pFullIncludeTable = (tmask*)m_fullExcludeTable;
	unsigned int last = 0;
	m_pRowSolutionMasksIdx[0] = 0;
	i = m_numPreconstructedRows;
	const auto shift = m_numSolutionTotalB / sizeof(tmask);
	m_lenMask >>= 3;  // Length of the mask in long long's 

#if 1
	while (i < m_numPlayers) {
		auto first = last;
		last = m_pRowSolutionCntr[i];
		if (m_pRowSolutionMasksIdx && i < m_numPlayers - 1) {
			const auto lastAdj = last - m_numRecAdj;
			m_pRowSolutionMasksIdx[i] = lastAdj >> SHIFT;
			if (i == m_numPreconstructedRows || m_pRowSolutionMasksIdx[i] > m_pRowSolutionMasksIdx[i-1]) {
				const auto rem = REM(lastAdj);
				if (rem)
					m_pRowSolutionMasks[i] = (tmask)(-1) << rem;
			}
			else {
				// We cannot use code for UseSolutionMasks, because now our code is not ready for such cases
				delete[] m_pRowSolutionMasksIdx;
				delete[] m_pRowSolutionMasks;
				m_pRowSolutionMasksIdx = NULL;
				m_pRowSolutionMasks = NULL;
			}
		}

		i++;
#if NEW
		if (!first) {
			// Skip construction of masks for the first set of solutions.
			// The threads will do this latter.
			continue;
		}
#endif
		while (first < last) {
			generateCompatibilityMasks(pFullIncludeTable, first++, last);
			pFullIncludeTable += shift;
		}
	}

	if (m_numRecAdj) {
		for (int i = m_numPreconstructedRows; i < m_numPlayers; i++)
			m_pRowSolutionCntr[i] -= m_numRecAdj;
	}

#else
	// Calculate the number of mutually compatible pairs of solutions
	int cntrs[8]; 
	memset(cntrs, 0, sizeof(cntrs));
	unsigned long long fff = 0;
	int a = 0;
	while (i < m_numPlayers - 1) {
		auto first = last;
		last = m_pRowSolutionCntr[i];
		i++;
		while (first < last) {
			auto* rm = (const long long*)m_pMaskStorage->getObject(first);
			const auto pRow = getObject(first++);
			ASSERT(pRow[1] != i);
			const auto pNeighbors = pRow + m_numPlayers;
			auto idx = last - 1;
			while (++idx < m_pRowSolutionCntr[i]) {
				// Let's check if the masks are mutually compatible
				auto* pMask = (const long long*)(m_pMaskStorage->getObject(idx));
				int j = m_lenMask;
				while (j-- && !(rm[j] & pMask[j]));

				if (j < 0 && p1fCheck2(pNeighbors, getObject(idx) + m_numPlayers)) {
					cntrs[i / 2]++;
					if (first <= m_numRecAdj && idx < m_numRecAdj2)
						a++;
				}
			}
		}
		fff += cntrs[i / 2];
		last = m_pRowSolutionCntr[i++];
	}
#endif
#if !NEW
	delete m_pMaskStorage;
	m_pMaskStorage = NULL;
#endif
}

long long cntr = 0;
CC int CRowUsage::getRow(int iRow, int ipx)
{
#if !USE_CUDA
	//cntr++;
#endif
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows);
	auto last = m_pRowSolutionIdx[iRow + 1] = m_pRowStorage->numRowSolutions(iRow);
	auto& first = m_pRowSolutionIdx[iRow];
	if (iRow == numPreconstructedRows) {
		if (first >= (last = m_pRowStorage->lastInFirstSet()))
			return 0;
	
		if (m_bUseCombinedSolutions) {
			m_bSolutionReady = m_pRowStorage->maskForCombinedSolutions((tmask*)(m_pCompatibleSolutions + m_numSolutionTotalB), first);
			if (!m_bSolutionReady)
				return 0;
		}
		else {
			memset(m_pCompatibleSolutions, 0, m_numSolutionTotalB);
			m_pRowStorage->generateCompatibilityMasks((tmask*)m_pCompatibleSolutions, first, last);
		}

		first += m_step;
		return 1;
	}

	const auto nRow = iRow - numPreconstructedRows - 1;
	if (m_bUseCombinedSolutions && !nRow) {
		if (m_bSolutionReady) {
			m_bSolutionReady = false;
			return 1;
		}

		return 0;
	}

	auto* pCompSol = (tmask*)(m_pCompatibleSolutions + nRow * m_numSolutionTotalB);
#if UseIPX
	ctchar* pPrevSolution = ipx > 0 ? m_pRowStorage->getObject(first - 1) : NULL;
#endif
	const auto lastB = IDX(last);

#if UseSolutionMasks || UseIPX
	while (true) {
#endif
		// Skip all bytes/longs equal to 0
		auto firstB = first >> SHIFT;
		while (firstB < lastB && !pCompSol[firstB])
			firstB++;

		if (firstB >= lastB)
			return 0;

#if USE_64_BIT_MASK
		unsigned long iBit;
		_BitScanForward64(&iBit, *(pCompSol + firstB));
#else
		const auto iBit = m_pRowStorage->firstOnePosition(pCompSol[firstB]);
#endif
		if ((first = (firstB << SHIFT) + iBit) >= last)
			return 0;

		pCompSol[firstB] ^= (tmask) 1 << iBit;
#if UseIPX
		// Previous solution should be different in first ipx bytes
		if (ipx && pPrevSolution && !memcmp(pPrevSolution, m_pRowStorage->getObject(first), ipx + 1)) {
			first++; // We need to try next solution 
			continue;
		}
#endif
		if (iRow < m_nRowMax) {
			// Construct the intersection of compatible solutions only if we will use it.
			const auto numLongs2Skip = m_pRowStorage->numLongs2Skip(iRow);
			auto pPrevA = (const long long*)(pCompSol) + numLongs2Skip;
			const auto shift = m_numSolutionTotalB >> 3;
			auto pToA = (long long*)(pPrevA + shift);
			auto pFromA = m_pRowStorage->getSolutionMask(first) + numLongs2Skip;
			const auto len = shift - numLongs2Skip;
#if USE_INTRINSIC
			bitwise_multiply(pPrevA, pFromA, pToA, len);
#else
			for (auto j = len; j--;)
				pToA[j] = pPrevA[j] & pFromA[j];
#endif

#if UseSolutionMasks
			auto pRowSolutionMasksIdx = m_pRowStorage->rowSolutionMasksIdx();
			if (pRowSolutionMasksIdx) {
				auto* pToASol = (tmask*)(pToA - numLongs2Skip);
				auto pRowSolutionMasks = m_pRowStorage->rowSolutionMasks();

				int i = iRow + 1;
				auto jMax = pRowSolutionMasksIdx[iRow];
				for (; i <= m_nRowMax; i++) {
					auto j = jMax;
					jMax = pRowSolutionMasksIdx[i];

					// Check left, middle and right parts of the solution interval for i-th row
					auto mask = pRowSolutionMasks[i - 1];
					if (mask && (mask & pToASol[j++]))
						continue;  // at least one solution masked by left part of the interval is still valid

					// middle part
					while (j < jMax && !pToASol[j])
						j++;

					if (j < jMax)
						continue;   // at least one solution masked by middle part of the interval is still valid

					// There are no valid solutions with the indices inside 
					// the interval defined by set of long longs
					mask = pRowSolutionMasks[i];
					// If mask != 0, we need to check the right side of the intervals.
					if (!mask || !((~mask) & pToASol[jMax]))
						break;
				}

				if (i <= m_nRowMax) {
					first++;
					continue;
				}

#if UseSolutionCliques
				if (m_pRowStorage->useCliques(iRow)) {
					first++;
					if (ConstructCompatibleSolutionGraph(pToASol, iRow))
						return 2;   // Ready to proceed with the getMatrix2() call.

					continue;
				}
#endif  // UseSolutionCliques
			}
#endif  // UseSolutionMasks
		}
#if UseSolutionMasks || UseIPX
		break;
	}
#endif

	first++;
	return 1;
}
