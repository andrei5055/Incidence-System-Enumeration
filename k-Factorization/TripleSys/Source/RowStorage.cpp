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

CC void CRowStorage::initCompatibilityMasks(ctchar* u1fCycles)
{
	m_u1fCycles = u1fCycles;
	for (int i = 1; i < m_numPlayers; i++)
		m_pRowSolutionCntr[i] += m_pRowSolutionCntr[i - 1];

	// Define the number of first long long's we don't need to copy to the next row.
	memset(m_pNumLongs2Skip, 0, m_numPlayers * sizeof(m_pNumLongs2Skip[0]));
	int i = m_numPreconstructedRows;
	m_pNumLongs2Skip[i] = m_pRowSolutionCntr[i] >> 6;
	m_numRecAdj = NEW? m_pRowSolutionCntr[i] : 0;

	while (++i < m_numPlayers)
		m_pNumLongs2Skip[i] = ((m_pRowSolutionCntr[i] - m_numRecAdj) >> 6);

	const auto& numSolutionTotal = m_pRowSolutionCntr[m_numPlayers - 1];
	delete[] m_fullExcludeTable;
	m_numSolutionTotalB = ((numSolutionTotal - m_numRecAdj + 7) / 8 + 7) / 8 * 8;

	auto len = numSolutionTotal * m_numSolutionTotalB;
	m_fullExcludeTable = new tchar[len];
	memset(m_fullExcludeTable, 0, len);

	delete[] m_pRowSolutionMasksIdx;
	delete[] m_pRowSolutionMasks;
	len = numPlayers();
	m_pRowSolutionMasksIdx = new uint[len];

	m_pRowSolutionMasks = new tmask[len];
	memset(m_pRowSolutionMasks, 0, len * sizeof(m_pRowSolutionMasks[0]));

#if !USE_64_BIT_MASK || !NEW_GET_ROW
	// Filling the lookup table m_FirstOnePosition
	memset(m_FirstOnePosition, 0, sizeof(m_FirstOnePosition));
	for (int i = 2; i < 256; i += 2)
		m_FirstOnePosition[i] = m_FirstOnePosition[i >> 1] + 1;
#endif

	const auto jMax = m_lenMask >> 3;
	auto* pFullIncludeTable = (tmask*)m_fullExcludeTable;
	unsigned int last = 0;
	m_pRowSolutionMasksIdx[0] = 0;
	i = m_numPreconstructedRows;
	const auto shift = m_numSolutionTotalB / sizeof(tmask);
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
		while (first < last) {
			auto* rm = (const long long*)m_pMaskStorage->getObject(first);
			const auto pRow = getObject(first++);
			ASSERT(pRow[1] != i);
			const auto pNeighbors = pRow + m_numPlayers;
			auto idx = last - 1;
			while (++idx < numSolutionTotal) {
				// Let's check if the masks are mutually compatible
				auto* pMask = (const long long*)(m_pMaskStorage->getObject(idx));
				int j = jMax;
				while (j-- && !(rm[j] & pMask[j]));

				if (j < 0 && p1fCheck2(m_u1fCycles, pNeighbors, getObject(idx) + m_numPlayers, m_numPlayers)) {
					// The masks are compatible and the length of the cycle is equal to m_numPlayers
					const auto newIdx = idx - m_numRecAdj;
					SET_MASK_BIT(pFullIncludeTable, newIdx);     // 1 - means OK
				}
			}

			pFullIncludeTable += shift;
		}
	}

	if (m_numRecAdj) {
		for (int i = m_numPreconstructedRows; i < m_numPlayers; i++)
			m_pRowSolutionCntr[i] -= m_numRecAdj;
	}

#else
	int cntrs[8]; 
	memset(cntrs, 0, sizeof(cntrs));
	unsigned long long fff = 0;
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
				int j = jMax;
				while (j-- && !(rm[j] & pMask[j]));

				if (j < 0 && p1fCheck2(m_u1fCycles, pNeighbors, getObject(idx) + m_numPlayers, m_numPlayers))
					cntrs[i / 2]++;
			}
		}
		fff += cntrs[i / 2];
		last = m_pRowSolutionCntr[i++];
	}

#endif
	delete m_pMaskStorage;
	m_pMaskStorage = NULL;
}

long long cntr = 0;
CC int CRowUsage::getRow(int iRow, int ipx)
{
#if !USE_CUDA
	//cntr++;
#endif
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows);
	const auto last = m_pRowSolutionIdx[iRow + 1] = m_pRowStorage->numRowSolutions(iRow);
	auto& first = m_pRowSolutionIdx[iRow];
	const auto numLongs2Skip = m_pRowStorage->numLongs2Skip(iRow);
#if NEW_GET_ROW

	if (iRow == numPreconstructedRows) {
		if (first >= last + m_pRowStorage->numRecAdj())
			return 0;

		const auto shift = NEW? 0 : numLongs2Skip << 3;
		memcpy(m_pCompatibleSolutions + shift, ((const tchar*)m_pRowStorage->getSolutionMask(first)) + shift, m_numSolutionTotalB - shift);
		first += m_step;
		return 1;
	}

	auto* pCompSol = (tmask*)(m_pCompatibleSolutions + (iRow - numPreconstructedRows - 1) * m_numSolutionTotalB);
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
			auto pPrevA = (const long long*)(pCompSol)+numLongs2Skip;
			const auto shift = m_numSolutionTotalB >> 3;
			auto pToA = (long long*)(pPrevA + shift);
			auto pFromA = m_pRowStorage->getSolutionMask(first + m_pRowStorage->numRecAdj()) + numLongs2Skip;
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
#else
	if (iRow == numPreconstructedRows) {
		if (first >= last)
			return 0;
		m_excludeForRow[iRow] = (tchar*)m_pRowStorage->getSolutionMask(first);
		first += m_step;
		return 1;
	}
#if UseIPX
	ctchar* pPrevSolution = ipx > 0 ? m_pRowStorage->getObject(first - 1) : NULL;
	while (true) {
#endif
		while (first < last) {
			unsigned int firstB = first >> 3;
			tchar shift = first & 0x7;
			tchar msk = !shift ? 0xff : ~((1 << shift) - 1);
			for (int j = numPreconstructedRows; j < iRow; j++)
			{
				msk &= (m_excludeForRow[j])[firstB];
				//if (!msk)
				//	break;
			}
			if (!msk)
				first = (first + 8) & 0xfffffff8;
			else {
				first = first + m_pRowStorage->firstOnePosition(msk) - shift;
				break;
			}
		}

		if (first >= last)
			return 0;
#if UseIPX
		// Previous solution should be different in first ipx bytes
		if (!pPrevSolution || memcmp(pPrevSolution, m_pRowStorage->getObject(first), ipx + 1))
			break;

		first++; // We need to try next solution 
	}
#endif

	m_excludeForRow[iRow] = (tchar *)m_pRowStorage->getSolutionMask(first);
#endif

	first++;
	return 1;
}
