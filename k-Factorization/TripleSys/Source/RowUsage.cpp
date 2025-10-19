#include "TripleSys.h"
#include "OneApp.h"
//#include <ppl.h>
//#include <array>
//#include <sstream>

#define USE_INTRINSIC		!USE_CUDA

#if USE_INTRINSIC
#include <immintrin.h> // Header for AVX2 intrinsics
void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	int i = 0;
	// Process 4 `long long` elements at a time using AVX2
	for (; i + 4 <= n; i += 4) {
		// Load 4 `long long` elements from each array
		__m256i vec_b = _mm256_loadu_si256((__m256i*) & b[i]);
		__m256i vec_c = _mm256_loadu_si256((__m256i*) & c[i]);

		// Perform bitwise AND operation
		__m256i vec_a = _mm256_and_si256(vec_b, vec_c);

		// Store the result back to the result array
		_mm256_storeu_si256((__m256i*) & a[i], vec_a);
	}

	// Handle the remainder elements (if `size` is not a multiple of 4)
	for (; i < n; ++i) {
		a[i] = b[i] & c[i];
	}
}
#else
void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	const int n8max = n & (~7);
	int j = 0;
	for (; j < n8max; j += 8) {
		a[j] = b[j] & c[j];
		a[j + 1] = b[j + 1] & c[j + 1];
		a[j + 2] = b[j + 2] & c[j + 2];
		a[j + 3] = b[j + 3] & c[j + 3];
		a[j + 4] = b[j + 4] & c[j + 4];
		a[j + 5] = b[j + 5] & c[j + 5];
		a[j + 6] = b[j + 6] & c[j + 6];
		a[j + 7] = b[j + 7] & c[j + 7];
	}
	for (; j < n; j++)
		a[j] = b[j] & c[j];
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

#define nc 0//500000 // 100000000;
#if nc
void testLogicalMultiplication(const long long* h_A, const long long* h_B, long long* h_C, long long *pCPU_res, int N, int nRep)
{
#define MULT_FLG 3
#if MULT_FLG & 1
	__declspec(dllimport) void logical_mul(const long long* h_A, const long long* h_B, long long* h_C, int N, int nRep);
	logical_mul(h_A, h_B, h_C, N, nRep);
#endif
#if MULT_FLG & 2

	MultArrays(h_A, h_B, h_C, N, nRep); // from j to <= jMax
#endif

	if (pCPU_res) {
		auto tm = clock();

		for (int k = 0; k < nRep; k++)
			multiplyAll(pCPU_res, h_A, h_B, N); // from j to <= jMax 

		printfYellow("cpu x %d iterations  lenArray = %d: %ld ms\n", nRep, N, clock() - tm);
		const char* cmpRes = memcmp(h_C, pCPU_res, N * sizeof(*h_C)) == 0 ?
			"the same as" : "different than";
		printfYellow("Graphic card data are %s cpu", cmpRes);
	}
	exit(0);

}
#else
#define testLogicalMultiplication(...)
#endif

CC int CRowUsage::getRow(int iRow, int ipx) {
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows || iRow >= m_pRowStorage->numDaysResult());

	const auto nRow = iRow - numPreconstructedRows - 1;
	const ll availablePlayers = nRow >= 0
		? *((const ll*)(m_pCompatibleSolutions + (nRow + 1) * m_lenMask) - 1)
		: m_pRowStorage->getPlayersMask();

	uint last = iRow;
	auto& first = m_pRowStorage->getSolutionInterval(m_pRowSolutionIdx + last, &last, availablePlayers);
	if (last == UINT_MAX)
		return availablePlayers ? 0 : -1;

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
	const auto lastB = IDX(last);

	while (true) {
		// Skip all bytes/longs equal to 0
		auto firstB = first >> SHIFT;
		while (firstB < lastB && !pCompSol[firstB])
			firstB++;

		if (firstB >= lastB)
			return 0;

		unsigned long iBit;
		_BitScanForward64(&iBit, *(pCompSol + firstB));
		if ((first = (firstB << SHIFT) + iBit) >= last)
			return 0;

		pCompSol[firstB] ^= (tmask)1 << iBit;
		if (iRow < m_nRowMax) {
			// Construct the intersection of compatible solutions only if we will use it.
			const auto pPrevA = (const ll*)(pCompSol);
			auto pToA = (ll*)(pCompSol + m_lenMask);
			const auto pFromA = (const ll*)(m_pRowStorage->getSolutionMask(first));
			unsigned int numLongs2Skip = m_pRowStorage->numLongs2Skip(iRow);
			const auto pPrevAStart = pPrevA + numLongs2Skip;
			int jNum = m_lenMask - numLongs2Skip;
			auto pToAStart = pToA + numLongs2Skip;
			const auto pFromAStart = pFromA + numLongs2Skip;
			const auto pRowSolutionMasksIdx = m_pRowStorage->rowSolutionMasksIdx();
			if (pRowSolutionMasksIdx) {
				testLogicalMultiplication(pPrevAStart, pFromAStart, pToAStart, pToAStart, jNum, nc);
				if (!selectPlayerByMask()) {
					// Usually, we should be here only when groupSize == 2 and it's NOT 
					// a complete balanced multipartite graph case.
					auto pRowSolutionMasks = m_pRowStorage->rowSolutionMasks();
					int i = iRow;
					auto jMax = pRowSolutionMasksIdx[iRow];
					//static int ccc[16] = { 0 };
					for (; ++i <= m_nRowMax;) {
						/**
						ccc[i]++;
						if (ccc[i] > 10000000) {
							for (int j = 4; j < 15; j++)
								printf(" %4d:%2d", ccc[j] / 1000, j);
							printf("\n");
							memset(ccc, 0, sizeof(ccc));
						}**/
						auto j = jMax;
						jMax = pRowSolutionMasksIdx[i];
						jNum = jMax - j + 1;
						testLogicalMultiplication(pFromA + j, pPrevA + j, pToA + j, pToA + jMax + 1, jNum, nc);
						multiplyAll(pToA + j, pFromA + j, pPrevA + j, jNum); // from j to <= jMax 

						// Check left, middle and right parts of the solution interval for i-th row
						auto mask = pRowSolutionMasks[i - 1];

						if (mask) {
							if (mask & pToA[j]) {
								continue;  // at least one solution masked by left part of the interval is still valid
							}
							j++;
						}
						// middle part
						while (j < jMax && !pToA[j]) {
							j++;
						}

						if (j < jMax) {
							continue;   // at least one solution masked by middle part of the interval is still valid
						}
						// There are no valid solutions with the indices inside 
						// the interval defined by set of long longs
						mask = pRowSolutionMasks[i];
						// If mask != 0, we need to check the right side of the intervals.

						if (!mask || !((~mask) & pToA[jMax])) {
							break;
						}
					}

					if (i <= m_nRowMax) {
						first++;
						continue;
					}
				}
				else {
					multiplyAll(pToAStart, pPrevAStart, pFromAStart, jNum);
					//???break;
				}
			}
			else
				multiplyAll(pToAStart, pPrevAStart, pFromAStart, jNum);
		}

		break;
	}
	first++;
	return 1;
}
