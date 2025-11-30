#include "TripleSys.h"
#include "OneApp.h"
//#include <ppl.h>
//#include <array>
//#include <sstream>

#define USE_INTRINSIC		!USE_CUDA
#if 0
#include <DirectXMath.h>
using namespace DirectX;
void multiplyAll(ll* aa, const ll* bb, const ll* cc, int n) {
	XMUINT4* a = (XMUINT4*)aa, * b = (XMUINT4*)bb, * c = (XMUINT4*)cc;
	int i = 0;
	n += n;
	// Process 4 `uint` elements at a time 
	for (; i + 4 <= n; i += 4) {  // Load four unsigned integer components from each input array.
		XMVECTOR v1 = XMLoadUInt4(&b[i]);
		XMVECTOR v2 = XMLoadUInt4(&c[i]);

		// Perform the bitwise AND operation on the vectors.
		XMVECTOR v_result = XMVectorAndInt(v1, v2);

		// Store the result back into the output array.
		XMStoreUInt4(&a[i], v_result);
	}
	// Handle the remainder elements (if `size` is not a multiple of 4)
	for (; i < n; ++i) {
		*((uint*)(b + i)) = *((uint *)(b+i)) & *((uint*)(c + i));
	}
}
#elif USE_INTRINSIC
#include <immintrin.h> // Header for AVX2 intrinsics
#include <algorithm>
void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	int i = 0;
#if 0
	static int iAvx512Support = 0;
	switch (iAvx512Support) {
	case 0:
		if (n >= 8) {
			__try {
				__m512i vec1 = _mm512_loadu_si512(&b[0]);
				__m512i vec2 = _mm512_loadu_si512(&c[0]);
				__m512i and_result = _mm512_and_si512(vec1, vec2);
				_mm512_storeu_si512(&a[0], and_result);
				iAvx512Support = 1;
				i += 8;
			}
			__except (GetExceptionCode() == EXCEPTION_ILLEGAL_INSTRUCTION ? EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH) {
				iAvx512Support = -1;
			}
		}
		break;
	case 1:
		for (; i + 8 <= n; i += 8) {
			// Load 8 64-bit integers from arr1 into a 512-bit register (zmm)
			__m512i vec1 = _mm512_loadu_si512(&b[i]);
			// Load 8 64-bit integers from arr2 into a 512-bit register
			__m512i vec2 = _mm512_loadu_si512(&c[i]);

			// Perform the bitwise AND operation on the two registers
			__m512i and_result = _mm512_and_si512(vec1, vec2);

			// Store the result back into the output array
			_mm512_storeu_si512(&a[i], and_result);
		}
		break;
	}
#endif
	// Process 4 `long long` elements at a time using AVX2
	for (; i + 4 <= n; i += 4) {
		// Load 4 `long long` elements from each array
		__m256i vec_b = _mm256_loadu_si256((__m256i*) & b[i]);
#if 1
		if (_mm256_testz_si256(vec_b, vec_b)) {
			_mm256_storeu_si256((__m256i*) & a[i], vec_b);
			continue;
		}
#endif
		__m256i vec_c = _mm256_loadu_si256((__m256i*) & c[i]);

		// Perform bitwise AND operation
		__m256i vec_a = _mm256_and_si256(vec_b, vec_c);

		// Store the result back to the result array
		_mm256_storeu_si256((__m256i*) & a[i], vec_a);
	}

	// Handle the remainder elements (if `size` is not a multiple of 4)
	switch (n - i) {
	case 3: a[i + 2] = b[i + 2] & c[i + 2];
	case 2: a[i + 1] = b[i + 1] & c[i + 1];
	case 1: a[i] = b[i] & c[i];
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
	m_threadID = iThread;
	m_step = numThreads;
	m_lenMask = m_pRowStorage->initRowUsage(&m_pCompatibleSolutions, &m_bSelectPlayerByMask);
	m_pRowSolutionIdx[m_pRowStorage->numPreconstructedRows()] = 0;
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

#if COUNT_GET_ROW_CALLS
ll getRowCallsCalls = 0;
size_t totalWeighChange = 0;
#define incGetRowCalls()	getRowCallsCalls++
#else
#define incGetRowCalls()
#endif

CC int CRowUsage::getRow(int iRow, int ipx) {
	incGetRowCalls();
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT_IF(iRow < numPreconstructedRows || iRow >= m_pRowStorage->numDaysResult());

	const auto nRow = iRow - numPreconstructedRows - 1;
	const ll availablePlayers = nRow >= 0
		? *((const ll*)(m_pCompatibleSolutions + (nRow + 1) * m_lenMask) - 1)
		: m_pRowStorage->getPlayersMask();

	uint last = iRow;
	auto& first = m_pRowStorage->getSolutionInterval(m_pRowSolutionIdx + last, &last, availablePlayers);
	if (!first)
		first += m_threadID;

	if (last == UINT_MAX)
		return availablePlayers ? 0 : -1;

	if (last > m_pRowStorage->getNumSolution())
		return 0;

	ASSERT_IF(last > m_pRowStorage->getNumSolution());

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
#if 0
	static int c[16], cc;
	c[m_threadID]++;
	cc++;
	if (cc >= 10000000) {
		printTable("nr", c, 1, 16, 0);
		cc = 0;
	}
#endif

	while (true) {
		// Skip all bytes/longs equal to 0
		auto firstB = first >> SHIFT;
		while (firstB < lastB && !pCompSol[firstB])
			firstB++;

		if (firstB >= lastB)
			return 0;

		unsigned long iBit;
		//iBit = (unsigned long)_tzcnt_u64(*(pCompSol + firstB));
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
					for (; ++i <= m_nRowMax;) {
						auto j = jMax;
						jMax = pRowSolutionMasksIdx[i];
						jNum = jMax - j + 1;
						testLogicalMultiplication(pFromA + j, pPrevA + j, pToA + j, pToA + jMax + 1, jNum, nc);
						multiplyAll(pToA + j, pPrevA + j, pFromA + j, jNum); // from j to <= jMax 
						/**
						if (i == m_nRowMax - 2)
						{
							static int ic = 0, c = 0;
							for (int k = 3; k <= jNum; k+=4)
							{

								if (pToA[k] || pToA[k - 1] || pToA[k - 2] || pToA[k - 3])
									ic++;
							}
							c += jNum / 4;
							if (c > 1000000) {
								printf("%d ", ic * 100 / c); c = ic = 0;
							}
						}
						*/
						// Check left, middle and right parts of the solution interval for i-th row
						auto mask = pRowSolutionMasks[i - 1];
						/*
						static int c[16], cc;
						c[i]++;
						cc++;
						if ((cc % 100000000) == 0) {
							//printTable("nr", c, 1, 16, 0);
							cc = 0;
						}**/
						if (mask) {
							if (mask & pToA[j]) {
								continue;  // at least one solution masked by left part of the interval is still valid
							}
							j++;
						}
						// middle part
#if 0
						while (j < jMax && !pToA[j]) {
							j++;
						}
						// we do not need "if" below if we use "while" below instead of while above
						if (j < jMax) {
							continue;   // at least one solution masked by middle part of the interval is still valid
						}
#else
						while (j + 4 <= jMax) {
#if 1 // USE_INTRINSIC
							__m256i fourValues = _mm256_loadu_si256((__m256i*) & pToA[j]);
							if (!_mm256_testz_si256(fourValues, fourValues)) goto Cont1;
#else
							if (pToA[j] || pToA[j + 1] || pToA[j + 2] || pToA[j + 3]) 
								goto Cont1;
#endif
							j += 4;
						}
						switch (jMax - j) {
						case 3: if (pToA[j + 2]) goto Cont1;
						case 2: if (pToA[j + 1] ) goto Cont1;
						case 1: if (pToA[j]) goto Cont1;
						}
#endif
						// There are no valid solutions with the indices inside 
						// the interval defined by set of long longs
						mask = pRowSolutionMasks[i];
						// If mask != 0, we need to check the right side of the intervals.

						if (!mask || !((~mask) & pToA[jMax])) {
							break;
						}
					Cont1: continue;
					}

					if (i <= m_nRowMax) {
						first++;
						continue;
					}
				}
				else {
					multiplyAll(pToAStart, pPrevAStart, pFromAStart, jNum);
#if 0
					{
						static int ic = 0, c = 0;
						for (int k = 3; k <= jNum; k += 4)
						{

							if (pPrevAStart[k] || pPrevAStart[k - 1] || pPrevAStart[k - 2] || pPrevAStart[k - 3])
								ic++;
						}
						c += jNum / 4;
						if (c > 100000000) {
							printf("%d:%d ", jNum, ic * 100 / c); c = ic = 0;
						}
					}
#endif
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
