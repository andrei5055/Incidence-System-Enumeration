#include "TripleSys.h"
#include "Table.h"

#include "OneApp.h"
#define UseIPX				0 // works faster with 0

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
	const auto fullMatrix = !m_pRowStorage->sysParam()->useFeature(t_useCompressedMasks);
	m_pRowStorage->initRowUsage(&m_pCompatibleSolutions, fullMatrix, &m_bSelectPlayerByMask);
	m_pCompatSolutions = m_pCompatibleSolutions;
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

	const auto nRow = iRow - numPreconstructedRows - 1;
	CCompressedMask* pCompatMasksTmp = NULL;
	if (nRow < 0 && m_pCompatMasks != m_pRowStorage) {
		pCompatMasksTmp = (CCompressedMask *)m_pCompatMasks;
		m_pCompatMasks = m_pRowStorage;
	}

	const auto lenMask = m_pCompatMasks->lenSolutionMask();
	ll availablePlayers = nRow >= 0
		? *((const ll*)(m_pCompatSolutions + (nRow + 1) * lenMask) - 1)
		: m_pRowStorage->getPlayersMask();

#if CHECK_GET_ROW
	static ll ggg; ggg++;
	const ll gBnd = 15405721;// 1444331332;// 1685633 //6486548; //
	static int cntr = 0;
	extern TableAut* pReslt;
	const char* fileName = pReslt->outFileName();
	if (ggg > gBnd - 50) {
		FOPEN_F(f, fileName, "a");
		fprintf(f, "cntr = %3d: iRow = %2d  availablePlayers = %llx\n", ++cntr, iRow, selectPlayerByMask() || nRow < 0 ? availablePlayers : -1);
		FCLOSE_F(f);
	}
#endif
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
		if (first >= last) {
			delete pCompatMasksTmp;
			return 0;
		}

		if (m_bUseCombinedSolutions) {
			m_bSolutionReady = m_pRowStorage->maskForCombinedSolutions(m_pCompatibleSolutions + lenMask, first);
			if (!m_bSolutionReady)
				return 0;
		}
		else {
			m_pRowStorage->passCompatibilityMask(m_pCompatibleSolutions, first, last, pAllData, (CCompatMasks**)&m_pCompatMasks, pCompatMasksTmp);
			if (m_pCompatMasks != m_pRowStorage) {
				availablePlayers = m_pCompatibleSolutions[m_pRowStorage->lenSolutionMask() - 1];
				auto pPntr = (tmask **)m_pCompatMasks->compatibleSolutionsPntr();
				m_pCompatMasks->initRowUsage(pPntr, true, &m_bSelectPlayerByMask, &availablePlayers);
				m_pCompatSolutions = m_pCompatMasks->compatibleSolutions();
				m_pRowSolutionIdx[iRow] = first;
				m_pRowSolutionIdx[iRow + 1] = selectPlayerByMask()? 0 : 1;
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

	auto* pCompSol = m_pCompatSolutions + nRow * lenMask;
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
		const auto* const pPrevA = (const ll*)(pCompSol);
		auto pToA = (ll*)(pCompSol + lenMask);
		const auto* const pFromA = (const ll*)(m_pCompatMasks->getSolutionMask(first));
		const auto* const pRowSolutionMasksIdx = m_pCompatMasks->rowSolutionMasksIdx();
		if (pRowSolutionMasksIdx && !selectPlayerByMask()) {
			// Usually, we should be here only when groupSize == 2 and 
			// it's NOT a complete balanced multipartite graph case.
			const auto* const pRowSolutionMasks = m_pCompatMasks->rowSolutionMasks();
			int i = iRow;
			auto jMax = pRowSolutionMasksIdx[i];
			int k = 0;
			for (; ++i <= m_nRowMax;) {
				auto j = jMax;
				jMax = pRowSolutionMasksIdx[i];
				const auto jRead = j + k;
				const auto jNum = jMax - jRead + 1;
				k = 1;
				auto mask = pRowSolutionMasks[i - 1];
				testLogicalMultiplication(pFromA + j, pPrevA + j, pToA + j, pToA + jMax + 1, jNum, nRep);
				ASSERT_IF(jRead + jNum > lenMask);
				multiplyAll(pToA + jRead, pPrevA + jRead, pFromA + jRead, jNum); // from j to <= jMax 

				// Check left, middle and right parts of the solution interval for i-th row
				if (mask) {
					if (mask & pToA[j])
						continue;  // at least one solution masked by left part of the interval is still valid

					j++;
				}

				// middle part
#if 0
				while (j < jMax && !pToA[j])
					j++;
						
				// we do not need "if" below if we use "while" below instead of while above
				if (j < jMax)
					continue;   // at least one solution masked by middle part of the interval is still valid					}
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
				case 3: if (pToA[j + 2]) continue;
				case 2: if (pToA[j + 1]) continue;
				case 1: if (pToA[j]) continue;
				}
#endif
				// There are no valid solutions with the indices inside 
				// the interval defined by set of long longs
				mask = pRowSolutionMasks[i];
				// If mask != 0, we need to check the right side of the intervals.

				if (!mask || !((~mask) & pToA[jMax]))
					break;
					
			Cont1:;
			}

			if (i <= m_nRowMax) {
				first++;
				continue;
			}
		}
		else {
			const auto numLongs2Skip = m_pCompatMasks->numLongs2Skip(iRow);
			auto pToAStart = pToA + numLongs2Skip;
			const auto* const pPrevAStart = pPrevA + numLongs2Skip;
			const auto* const pFromAStart = pFromA + numLongs2Skip;
			testLogicalMultiplication(pPrevAStart, pFromAStart, pToAStart, pToAStart, lenMask - numLongs2Skip, nRep);
			multiplyAll(pToAStart, pPrevAStart, pFromAStart, lenMask - numLongs2Skip);
		}


		break;
	}

	first++;

#if CHECK_GET_ROW	
	if (ggg > gBnd - 50) {
		if (selectPlayerByMask()) {
			FOPEN_F(f1, fileName, "a");
			//fprintf(f1, "first = %2d  availablePlayers = %lld\n", first, availablePlayers);
			availablePlayers = *((const ll*)(m_pCompatSolutions + (nRow + 2) * lenMask) - 1);
			fprintf(f1, "availablePlayers = %llx\n", availablePlayers);
			FCLOSE_F(f1);
		}

		auto pRes = m_pRowStorage->allData()->result();
		getMatrix(pRes, m_pRowStorage->allData()->neighbors(), iRow + 1);
		pReslt->printTable(pRes, true, false, iRow + 1);
	}
#endif
	return 1;
}
