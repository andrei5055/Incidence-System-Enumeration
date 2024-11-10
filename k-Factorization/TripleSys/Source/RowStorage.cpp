#include "TripleSys.h"

#if USE_64_BIT_MASK
#include <bitset>
#endif

int cntr = 0;

CC bool CRowUsage::getRow(int iRow, int ipx) const
{
	const auto numPreconstructedRows = m_pRowStorage->numPreconstructedRows();
	ASSERT(iRow < numPreconstructedRows);
	const auto last = m_pRowSolutionIdx[iRow + 1] = m_pRowStorage->numRowSolutions(iRow);
	auto& first = m_pRowSolutionIdx[iRow];
#if NEW_GET_ROW

	if (iRow == numPreconstructedRows) {
		if (first >= last)
			return false;

		memcpy(m_pCompatibleSolutions, m_pRowStorage->getSolutionMask(first), m_numSolutionTotalB);
		first += m_step;
		return true;
	}

	auto* pCompSol = (tmask*)(m_pCompatibleSolutions + (iRow - numPreconstructedRows - 1) * m_numSolutionTotalB);
	ctchar* pPrevSolution = ipx > 0 ? m_pRowStorage->getObject(first - 1) : NULL;
	const auto lastB = IDX(last);

#define UseIPX 0 // works faster with 0
#if UseSolutionMasks
	while (true)
#endif
	{
#if UseIPX
		while (true) {
#endif
#if USE_64_BIT_MASK
			// Skip all longs equal to 0
			auto firstB = first >> 6;
			unsigned long iBit = 0;
			while (firstB < lastB && !pCompSol[firstB])
				firstB++;

			if (firstB >= lastB || !_BitScanForward64(&iBit, *(pCompSol + firstB))) {
				first = last;
				return false;
			}

			first = (firstB << 6) + iBit;
			pCompSol[firstB] ^= (long long)1 << iBit;
#else
			// Skip all bytes equal to 0
			auto firstB = first >> 3;
			while (firstB < lastB && !pCompSol[firstB])
				firstB++;

			if (firstB >= lastB)
				return false;

			const auto shift = m_pRowStorage->firstOnePosition(pCompSol[firstB]);
			if ((first = (firstB << 3) + shift) >= last)
				return false;

			pCompSol[firstB] ^= 1 << shift;
#endif
#if UseIPX
			// Previous solution should be different in first ipx bytes
			if (!ipx || !pPrevSolution || memcmp(pPrevSolution, m_pRowStorage->getObject(first), ipx + 1))
				break;
			first++; // We need to try next solution 
		}
#endif
		if (iRow < m_nRowMax) {
			// Construct the intersection of compatible solutions only if we will use it.
			const auto numLongs2Skip = m_pRowStorage->numLongs2Skip(iRow);
			auto pPrevA = (const long long*)(pCompSol)+numLongs2Skip;
			auto pToA = (long long*)(pPrevA + (m_numSolutionTotalB >> 3));
			auto pFromA = m_pRowStorage->getSolutionMask(first) + numLongs2Skip;
			for (auto j = (m_numSolutionTotalB >> 3) - numLongs2Skip; j--;)
				pToA[j] = pPrevA[j] & pFromA[j];

#if UseSolutionMasks
			auto pRowSolutionMasksIdx = m_pRowStorage->rowSolutionMasksIdx();
			auto pRowSolutionMasks = m_pRowStorage->rowSolutionMasks();

			auto jMax = pRowSolutionMasksIdx[iRow - 1];
			int i = iRow;
			for (; i < m_nRowMax; i++) {
				auto j = jMax;
				jMax = pRowSolutionMasksIdx[i];

				// Check left, middle and right parts of the solution interval for i - th row
				auto mask = pRowSolutionMasks[i];
				if (mask && (mask & pToA[j - 1]))
					continue;  // at least one solution masked by first left part of the interval is still valid

				// middle part
				while (j < jMax && !pToA[j])
					j++;

				if (j < jMax)
					continue;   // at least one solution masked by middle part of the interval is still valid

				// There are no valid solutions with the indices inside 
				// the interval defined by set of long longs
				// We need to check the left and right sides of the intervals.
				mask = pRowSolutionMasks[i + 1];
				if (!mask || !((~mask) & pToA[jMax]))
					break;
			}

#if 0
			first++;
			if (i >= m_nRowMax)
				return true;
#else
			if (i >= m_nRowMax)
				break;

			first++;
#endif
#endif
		}
		else
			break;
	}
#else
	if (iRow == numPreconstructedRows) {
		if (first >= last)
			return false;
		m_excludeForRow[iRow] = (tchar*)m_pRowStorage->getSolutionMask(first);
		first += m_step;
		return true;
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
			return false;
#if UseIPX
		// Previous solution should be different in first ipx bytes
		if (!pPrevSolution || memcmp(pPrevSolution, m_pRowStorage->getObject(first), ipx + 1))
			break;

		first++; // We need to try next solution 
	}
#endif

	m_excludeForRow[iRow] = (tchar *)m_pRowStorage->getSolutionMask(first);
#endif

#if 0	
	FOPEN_F(f, "aaa.txt", cntr++ ? "a" : "w");
	fprintf(f, "cntr = %2d:  iRow = %2d  first = %4d\n", cntr, iRow, first);

	FCLOSE_F(f);
	if (cntr >= 10)
		cntr += 0;
#else
	cntr++;
#endif

	first++;
	return true;
}
