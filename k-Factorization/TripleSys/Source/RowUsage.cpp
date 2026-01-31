#include "TripleSys.h"
#include "Table.h"

#include "OneApp.h"
#define UseIPX				0 // works faster with 0
#if USE_INTRINSIC
#include <immintrin.h> // Header for AVX2 intrinsics
#include <algorithm>
void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
	int i = 0;
	// Unroll to process 8 elements (2 vectors) at a time
	for (; i + 8 <= n; i += 8) {
		__m256i v_b1 = _mm256_loadu_si256((__m256i*) & b[i]);
		__m256i v_b2 = _mm256_loadu_si256((__m256i*) & b[i + 4]);
		__m256i v_c1 = _mm256_loadu_si256((__m256i*) & c[i]);
		_mm256_storeu_si256((__m256i*) & a[i], _mm256_and_si256(v_b1, v_c1));
		__m256i v_c2 = _mm256_loadu_si256((__m256i*) & c[i + 4]);
		_mm256_storeu_si256((__m256i*) & a[i + 4], _mm256_and_si256(v_b2, v_c2));
	}

	// Process remaining 4-element chunks if they exist
	for (; i + 4 <= n; i += 4) {
		__m256i v_b = _mm256_loadu_si256((__m256i*) & b[i]);
		__m256i v_c = _mm256_loadu_si256((__m256i*) & c[i]);
		_mm256_storeu_si256((__m256i*) & a[i], _mm256_and_si256(v_b, v_c));
	}

	// Scalar remainder loop
	for (; i < n; ++i) {
		a[i] = b[i] & c[i];
	}
}
#else
CC void multiplyAll(ll* a, const ll* b, const ll* c, int n) {
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

#ifdef USE_CUDA
#define BIT_SCAN_FORWARD64 BitScanForward64_Cuda
// Dummy function written by AI and never tested
CC void BitScanForward64_Cuda(unsigned long* _Index, unsigned long long _Mask) {
	unsigned long index = 0;
	if (_Mask == 0)
		return;
	while (!(_Mask & 1)) {
		_Mask >>= 1;
		index++;
	}
	*_Index = index;
}
#else
#define BIT_SCAN_FORWARD64  _BitScanForward64
#endif

CC void CRowUsage::init(int iThread, int numThreads) {
	m_threadID = iThread;
	m_step = numThreads;
	m_pRowStorage->initRowUsage(&m_pCompatibleSolutions, &m_bSelectPlayerByMask);
	m_pRowSolutionIdx[m_pRowStorage->numPreconstructedRows()] = 0;
}

#define nRep 0//500000 // 100000000;
#if nRep
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
#define incGetRowCalls()	getRowCallsCalls++
#else
#define incGetRowCalls()
#endif

CC int CRowUsage::getRow(int iRow, int ipx, const alldata* pAllData) {
	incGetRowCalls();
	const auto numPreconstructedRows = m_pCompatMasks->numPreconstructedRows();
	ASSERT_IF(iRow < numPreconstructedRows || iRow >= m_pCompatMasks->numDaysResult());

#if CHECK_GET_ROW
	static int cntr = 0;
#define FILE_NAME "C:\\Users\\andre\\source\\repos\\andrei5055\\Incidence-System-Enumeration\\LogsTestB\\Complete_graphs\\14\\14x13x2\\P0000000001.txt"
	FOPEN_F(f, FILE_NAME, cntr++ ? "a" : "w");
	fprintf(f, "cntr = %3d: iRow = %2d\n", cntr, iRow);
	FCLOSE_F(f);
#endif
	const auto nRow = iRow - numPreconstructedRows - 1;
	const auto lenMask = m_pCompatMasks->lenSolutionMask();
	const ll availablePlayers = nRow >= 0
		? *((const ll*)(m_pCompatibleSolutions + (nRow + 1) * lenMask) - 1)
		: m_pCompatMasks->getPlayersMask();

	uint last = iRow;
	auto& first = m_pCompatMasks->getSolutionInterval(m_pRowSolutionIdx + last, &last, availablePlayers);

	if (last == UINT_MAX)
		return availablePlayers ? 0 : -1;

	if (last > m_pCompatMasks->getNumSolution())
		return 0;

	if (iRow == numPreconstructedRows) {
		if (!first) {
			int mode = -1;
			pAllData->cnvPrecalcRowsCompCheck(mode);
			first = m_threadID;
		}
		if (first >= last)
			return 0;

		if (m_bUseCombinedSolutions) {
			m_bSolutionReady = m_pRowStorage->maskForCombinedSolutions(m_pCompatibleSolutions + lenMask, first);
			if (!m_bSolutionReady)
				return 0;
		}
		else {
			m_pRowStorage->passCompatibilityMask(m_pCompatibleSolutions, first, last, pAllData, (CCompatMasks **)&m_pCompatMasks);
			if (m_pCompatMasks != m_pRowStorage) {
				delete[] m_pCompatibleSolutions;
				m_pCompatibleSolutions = NULL;
				m_pCompatMasks->initRowUsage(&m_pCompatibleSolutions, &m_bSelectPlayerByMask);
				// NOTE; Let's make a trivial mask for now and improve it later
//				memcpy(m_pCompatibleSolutions, m_pCompatMasks->rowsCompatMasks(), m_pCompatMasks->numSolutionTotalB());
				memset(m_pCompatibleSolutions + 1, 0xff, m_pCompatMasks->numSolutionTotalB() - 1);
				m_pCompatibleSolutions[0] = 0xfe;
				m_pRowSolutionIdx[iRow] = 0;     // first on current row to 0 - it will be inreased by m_step
				m_pRowSolutionIdx[iRow + 1] = 1; // index of the solution from the compressed set which will be used first for next row
			}
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

	auto* pCompSol = m_pCompatibleSolutions + nRow * lenMask;
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

		unsigned long iBit;
		//iBit = (unsigned long)_tzcnt_u64(*(pCompSol + firstB));
		BIT_SCAN_FORWARD64(&iBit, *(pCompSol + firstB));
		if ((first = (firstB << SHIFT) + iBit) >= last)
			return 0;

		pCompSol[firstB] ^= (tmask)1 << iBit;
		if (iRow >= m_nRowMax)
			break;

#if		UseIPX
		// Previous solution should be different in first ipx bytes
		if (ipx && nRow > 0 && pPrevSolution && !memcmp(pPrevSolution, m_pRowStorage->getObject(first), ipx + 1)) {
			//printTable("Skipped", m_pRowStorage->getObject(first), 1, 24, 2);
			first++; // We need to try next solution 
			continue;
		}
#endif
		// Construct the intersection of compatible solutions only if we will use it.
		const auto pPrevA = (const ll*)(pCompSol);
		auto pToA = (ll*)(pCompSol + lenMask);
		const auto pFromA = (const ll*)(m_pCompatMasks->getSolutionMask(first));
		unsigned int numLongs2Skip = m_pCompatMasks->numLongs2Skip(iRow);
		const auto pPrevAStart = pPrevA + numLongs2Skip;
		int jNum = lenMask - numLongs2Skip;
		auto pToAStart = pToA + numLongs2Skip;
		const auto pFromAStart = pFromA + numLongs2Skip;
		const auto pRowSolutionMasksIdx = m_pCompatMasks->rowSolutionMasksIdx();
		if (pRowSolutionMasksIdx) {
			testLogicalMultiplication(pPrevAStart, pFromAStart, pToAStart, pToAStart, jNum, nRep);
			if (!selectPlayerByMask()) {
				// Usually, we should be here only when groupSize == 2 and it's NOT 
				// a complete balanced multipartite graph case.
				auto pRowSolutionMasks = m_pCompatMasks->rowSolutionMasks();
				int i = iRow;
				auto jMax = pRowSolutionMasksIdx[iRow];
				for (; ++i <= m_nRowMax;) {
					auto j = jMax;
					jMax = pRowSolutionMasksIdx[i];
					jNum = jMax - j + 1;
					testLogicalMultiplication(pFromA + j, pPrevA + j, pToA + j, pToA + jMax + 1, jNum, nRep);
					multiplyAll(pToA + j, pPrevA + j, pFromA + j, jNum); // from j to <= jMax 

					// Check left, middle and right parts of the solution interval for i-th row
					auto mask = pRowSolutionMasks[i - 1];
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
#if USE_INTRINSIC
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
			}
		}
		else
			multiplyAll(pToAStart, pPrevAStart, pFromAStart, jNum);

		break;
	}

	first++;
#if CHECK_GET_ROW
	FOPEN_F(f1, FILE_NAME, "a");
	fprintf(f1, "first = %2d\n", first);
	FCLOSE_F(f1);

	extern TableAut * pReslt;
	auto pRes = m_pRowStorage->allData()->result();
	getMatrix(pRes, m_pRowStorage->allData()->neighbors(), iRow+1);
	pReslt->printTable(pRes, true, false, iRow+1);
#endif
	return 1;
}
