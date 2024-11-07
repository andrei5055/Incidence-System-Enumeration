#include "TripleSys.h"

#if USE_64_BIT_MASK
#include <bitset>
#endif

CC bool CRowUsage::getRow(int iRow, int ipx)
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
	const auto lastB =
#if USE_64_BIT_MASK
	(last + 63) >> 6;
#else
		(last + 7) >> 3;
#endif
#define UseIPX 0 // works faster with 0
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

	first++;
	return true;
}
