#include "TripleSys.h"

#define UseIPX				0 // works faster with 0
#define USE_INTRINSIC		0 //!USE_CUDA // with 0 we can calculate ptoa on the fly

#if USE_INTRINSIC
#define CalculatePtoAOnTheFly 0 // do not change, must be 0
#else
#define CalculatePtoAOnTheFly 1
#endif

#if USE_INTRINSIC || USE_64_BIT_MASK
#include <immintrin.h> // Header for AVX2 intrinsics

void bitwise_multiply(const ll* a, const ll* b, ll* result, size_t size) {
	size_t i = 0;
	/**
	// Process 4 `long long` elements at a time using AVX2
	for (; i + 4 <= size; i += 4) {
		// Load 4 `long long` elements from each array
		__m256i vec_a = _mm256_loadu_si256((__m256i*) & a[i]);
		__m256i vec_b = _mm256_loadu_si256((__m256i*) & b[i]);

		// Perform bitwise AND operation
		__m256i vec_result = _mm256_and_si256(vec_a, vec_b);

		// Store the result back to the result array
		_mm256_storeu_si256((__m256i*) & result[i], vec_result);

			// Load 4 `long long` elements from each array
			vec_a = _mm256_loadu_si256((__m256i*) & a[i + 4]);
			vec_b = _mm256_loadu_si256((__m256i*) & b[i + 4]);

			// Perform bitwise AND operation
			vec_result = _mm256_and_si256(vec_a, vec_b);

			// Store the result back to the result array
			_mm256_storeu_si256((__m256i*) & result[i + 4], vec_result);
		}
	}
	/**/
	// Handle the remainder elements (if `size` is not a multiple of 4)
	for (; i < size; ++i) {
		result[i] = a[i] & b[i];
	}
}
#endif

CC void CRowUsage::init(int iThread, int numThreads) {
	m_lenMask = m_pRowStorage->initRowUsage(&m_pCompatibleSolutions, &m_bSelectPlayerByMask);
	const auto iRow = m_pRowStorage->numPreconstructedRows();
	if (iThread) {
		// Determine the maximum number of solutions for the first non-predetermined row.  
		m_pRowSolutionIdx[0] = 0;
		uint last = iRow;
		m_pRowStorage->getSolutionInterval(m_pRowSolutionIdx, &last, m_pRowStorage->getPlayersMask());
	}

	m_threadID = m_pRowSolutionIdx[iRow] = iThread;
	m_step = numThreads;
}

#if COUNT_GET_ROW_CALLS && !USE_CUDA
ll cntr = 0;
#define increaseGetRowCounter()  cntr++
#else
#define increaseGetRowCounter()
#endif


#if 0 && !USE_CUDA
#define OUT_SELECTED_ROW_SOLUTION(pRowStorage, iRow, first, last) \
					pRowStorage->outSelectedSolution(iRow, first, last, m_threadID)
#else
#define OUT_SELECTED_ROW_SOLUTION(pRowStorage, iRow, first, last)  
#endif

CC int CRowUsage::getRow(int iRow, int ipx) {
	increaseGetRowCounter();
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows || iRow >= m_pRowStorage->numDaysResult());

	const auto nRow = iRow - numPreconstructedRows - 1;
	const ll availablePlayers = nRow >= 0
		? *((const ll*)(m_pCompatibleSolutions + (nRow + 1) * m_lenMask) - 1)
		: m_pRowStorage->getPlayersMask();

	uint last = iRow;
	auto& first = m_pRowStorage->getSolutionInterval(m_pRowSolutionIdx+last, &last, availablePlayers);
	if (last == UINT_MAX)
		return availablePlayers? 0 : -1;

	if (last > m_pRowStorage->getNumSolution())
		return 0;

	ASSERT(last > m_pRowStorage->getNumSolution());

	if (iRow == numPreconstructedRows) {
		if (first >= last)
			return 0;

		if (m_bUseCombinedSolutions) {
			m_bSolutionReady = m_pRowStorage->maskForCombinedSolutions(m_pCompatibleSolutions + m_lenMask, first);
			if (!m_bSolutionReady)
				return 0;
		}
		else {
			m_pRowStorage->passCompatibilityMask(m_pCompatibleSolutions, first, last);
		}

		OUT_SELECTED_ROW_SOLUTION(m_pRowStorage, iRow, first, last);
		first += m_step;
		return 1;
	}

	if (m_bUseCombinedSolutions && !nRow) {
		if (m_bSolutionReady) {
			m_bSolutionReady = false;
			return 1;
		}

		return 0;
	}

	auto* pCompSol = m_pCompatibleSolutions + nRow * m_lenMask;
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
#define multiplyAll() { const uint j8max = m_lenMask - 7; \
						uint  j = m_pRowStorage->numLongs2Skip(iRow); \
						if (m_lenMask > 7) { \
							for (; j < j8max; j += 8) { \
								pToA[j] = pPrevA[j] & pFromA[j]; \
								pToA[j + 1] = pPrevA[j + 1] & pFromA[j + 1]; \
								pToA[j + 2] = pPrevA[j + 2] & pFromA[j + 2]; \
								pToA[j + 3] = pPrevA[j + 3] & pFromA[j + 3]; \
								pToA[j + 4] = pPrevA[j + 4] & pFromA[j + 4]; \
								pToA[j + 5] = pPrevA[j + 5] & pFromA[j + 5]; \
								pToA[j + 6] = pPrevA[j + 6] & pFromA[j + 6]; \
								pToA[j + 7] = pPrevA[j + 7] & pFromA[j + 7]; \
						}} \
						for (; j < m_lenMask; j++) \
							pToA[j] = pPrevA[j] & pFromA[j]; \
					  }
#define multiplyAllOld()  for (uint j = m_pRowStorage->numLongs2Skip(iRow); j < m_lenMask; j++) \
							pToA[j] = pPrevA[j] & pFromA[j]; 

#if CalculatePtoAOnTheFly
#define ptoa(j) (pToA[j] = pPrevA[j] & pFromA[j])
#define multiply(j, jMax)  while (++j <= jMax) ptoa(j)
#define multiply_1()
#define multiply_2()    multiplyAll()
#else
#define ptoa(j) pToA[j]
#define multiply(j, jMax)  // empty macro
#define multiply_1()	multiplyAll()
#define multiply_2()    
#endif
			auto pPrevA = (const ll*)(pCompSol);
			auto pToA = (ll*)(pCompSol + m_lenMask);
			auto pFromA = (const ll*)(m_pRowStorage->getSolutionMask(first));
#if USE_INTRINSIC
			unsigned int numLongs2Skip = m_pRowStorage->numLongs2Skip(iRow);
			bitwise_multiply(pPrevA + numLongs2Skip, pFromA + numLongs2Skip, pToA + numLongs2Skip, m_lenMask - numLongs2Skip);
#else
			multiply_1();
#endif
			const auto pRowSolutionMasksIdx = m_pRowStorage->rowSolutionMasksIdx();
			if (pRowSolutionMasksIdx) {
				if (!selectPlayerByMask()) {
					// Usially, we should be here only when groupSize == 2 and it's NOT 
					// a complete balanced maltipartite graph case.
#if NEW
					if (!m_pRowStorage->checkSolutionByMask(iRow, pToASol)) {
						first++;
						continue;
					}
#else
					auto pRowSolutionMasks = m_pRowStorage->rowSolutionMasks();
					int i = iRow;
					auto jMax = pRowSolutionMasksIdx[iRow];
					for (; ++i <= m_nRowMax;) {
						auto j = jMax;
						jMax = pRowSolutionMasksIdx[i];

						// Check left, middle and right parts of the solution interval for i-th row
						auto mask = pRowSolutionMasks[i - 1];
						if (mask) {
							if (mask & ptoa(j)) {
								multiply(j, jMax);
								continue;  // at least one solution masked by left part of the interval is still valid
							}
							j++;
						}

						// middle part
						while (j < jMax && !ptoa(j))
							j++;

						if (j < jMax) {
							multiply(j, jMax);
							continue;   // at least one solution masked by middle part of the interval is still valid
						}

						// There are no valid solutions with the indices inside 
						// the interval defined by set of long longs
						mask = pRowSolutionMasks[i];
						// If mask != 0, we need to check the right side of the intervals.
						if (!mask || !((~mask) & ptoa(jMax))) {
							break;
						}
					}

					if (i <= m_nRowMax) {
						first++;
						continue;
					}
#endif

#if UseSolutionCliques
					if (m_pRowStorage->useCliques(iRow)) {
						first++;
						if (ConstructCompatibleSolutionGraph((tmask*)(pToA), iRow))
							return 2;   // Ready to proceed with the getMatrix2() call.

						continue;
					}
#endif  // UseSolutionCliques
				}  else {
					multiply_2();
					break;
				}
			} else
				multiply_2();
		}

		break;
	}

	OUT_SELECTED_ROW_SOLUTION(m_pRowStorage, iRow, first, last);
	first++;
	return 1;
}
