#include "TripleSys.h"

#define UseIPX				0 // works faster with 0
#define USE_INTRINSIC		!USE_CUDA

#if USE_INTRINSIC || USE_64_BIT_MASK
#include <immintrin.h> // Header for AVX2 intrinsics

void bitwise_multiply(const ll* a, const ll* b, ll* result, size_t size) {
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

ll cntr = 0;
CC int CRowUsage::getRow(int iRow, int ipx)
{
#if !USE_CUDA
	//cntr++;
#endif
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows);
	uint last;
	auto& first = m_pRowStorage->getSolutionInterval(m_pRowSolutionIdx, iRow, &last);
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

	while (true) {
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

		pCompSol[firstB] ^= (tmask)1 << iBit;
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
			auto pPrevA = (const ll*)(pCompSol)+numLongs2Skip;
			const auto shift = m_numSolutionTotalB >> 3;
			auto pToA = (ll*)(pPrevA + shift);
			auto pFromA = (const ll*)(m_pRowStorage->getSolutionMask(first)) + numLongs2Skip;
#if USE_INTRINSIC
			bitwise_multiply(pPrevA, pFromA, pToA, shift - numLongs2Skip);
#else
			for (auto j = shift - numLongs2Skip; j--;)
				pToA[j] = pPrevA[j] & pFromA[j];
#endif

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
		}

		break;
	}

	first++;
	return 1;
}
