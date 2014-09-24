#include "BIBD_Enumerator.h"
#include "EnumInfo.h"
#ifdef _MSC_VER
#ifndef _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#endif
#include <crtdbg.h>
#endif
#if defined(_MSC_VER) && defined(_DEBUG)
#define new new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
int CBIBD_Enumerator::unforcedElement(const CColOrbit *pOrb, int nRow) const
{
	const size_t diffWeight = getInSys()->GetK() - pOrb->colomnWeight();
	if (diffWeight)
		return diffWeight == rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

bool CBIBD_Enumerator::sortSolutions(CRowSolution *pSolution, size_t idx)
{
	if (currentRowNumb() + 1 == rowNumb())
		return true;        // We will be here when lambda > 1 AND one of the colOrbits was splitted into 2 parts  
	// in previous row. (6,10,5,3,2) is one such example.

	if (!pSolution->numSolutions())
		return false;

	pSolution->sortSolutions(useCanonGroup() ? canonChecker() : NULL);
	if (!pSolution->findFirstValidSolution(inSysRowEquation()->variableMaxLimitPntr(), GetData()))
		return false;

	return checkChoosenSolution(pSolution, currentRowNumb(), idx);
}

bool CBIBD_Enumerator::solutionsForRightSideNeeded(const VECTOR_ELEMENT_TYPE *pRighPart, const VECTOR_ELEMENT_TYPE *pCurrSolution, size_t nRow) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allow to do not consider the solutions for some right parts
	if (nRow == 3) {
		// Condition should eliminate right part with x1 = 0 for some BIBD, and (8, 4, 3) is one of them
		if (!*pRighPart && *pCurrSolution == 1) {
			const auto k = getInSys()->GetK();
			const auto lambda = getInSys()->lambda();
			const auto v = rowNumb();
			if ((k - 2) * lambda == v - 2)
				return false;
		}
	}

	return true;
}

bool CBIBD_Enumerator::checkChoosenSolution(CRowSolution *pCurrSolution, size_t nRow, size_t usedSolIndex) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allow to skip testing of some solutions for some rows
	if (nRow == 3) {
		const CRowSolution *pPrevSolution = rowStuff(nRow - 1);
		const auto lambda = getInSys()->lambda();
		// Number of units used for first fragments of the 3-d row
		const int x0 = *pPrevSolution->solution(usedSolIndex);

		if (lambda == 2 * x0) {
			// We do have following system:
			// 0*a0 + 1*a1 + 2*a2 + ... + n * an = lambda * (k -2)
			//   a0 +   a1 +   a2 + ... +     an = v - 2
			//
			// Substract second from the first one
			// (n - 1) * an + (n - 2) * a(n-1) + ... + a2 - a0 = lambda * (k - 2) - (v - 2)
			// since ai >= 0,
			// (n - 1) * an + (n - 2) * a(n-1) + ... + a2 >= lambda * (k - 2) - (v - 2)
			// an + a(n-1) *(n-2)/(n-1) + .... + a2 * 1/(n-1) >= (lambda * (k - 2) - (v - 2)) / (n - 1)

			if (x0 == 2) {
				const auto k = getInSys()->GetK();
				const auto v = rowNumb();
				// The difference represent minimal number of descendants of the solution used for 3-d row
				if (lambda * (k - 2) - v + 2 > 2) {
					// This should work for BIBD(22, 8,4)
					// In this case at least 2 descendants of the solution used for 3-d row
					// will be used for remaining rows
					// It means that the current rows solution should be the descendant of the solution used for 3-d row
					// Suppose that it's not.  Then at lest two rows with the numbers nRow >= 4 will form following
					// configuration:
					//    0 0 1 1
					//    0 0 1 1
					// When we change 3-d and 4-th row with these two AND columns 1 -2 with 3-4 we will increase the code
					// of the matrix. Because of that this matrix cannot be canonical.
					const bool useCanonGroup = USE_CANON_GROUP && canonChecker()->groupOrder() > 1;
					CRowSolution *pSol = pCurrSolution;
					while ((pSol = pSol->NextSolution(useCanonGroup)) != NULL) {
						const VECTOR_ELEMENT_TYPE *pSolValue = pSol->currSolution();
						if (*pSolValue + *(pSolValue + 1) >= x0)  // Why x0 here ???
							break;
					}

					if (!pSol)
						return false;

					pSol->prevSolutionIndex();
				}
			}
		}
	}

	return pCurrSolution->checkChoosenSolution(colOrbit(nRow), matrix()->rowNumb() - nRow, getInSys()->GetK());
}

bool CBIBD_Enumerator::TestFeatures(CEnumInfo *pEnumInfo)
{
	char buff[256];
	const auto paramR = getInSys()->GetR();
	const auto paramLambda = getInSys()->lambda();
	const auto iMax = rowNumb();
	for (size_t i = 0; i < iMax; i++) {
		const MATRIX_ELEMENT_TYPE *pRow = matrix()->GetRow(i);
		int r = 0;
		for (auto j = colNumb(); j--;)
			r += *(pRow + j);

		if (r != paramR) {
			matrix()->printOut(outFile());
			SPRINTF(buff, "Wrong number of units in the row # %lu: %d != %lu\n", i, r, getInSys()->GetR());
			outString(buff, outFile());
			throw;
			return false;
		}

		for (size_t j = 0; j < i; j++) {
			int lambda = 0;
			const MATRIX_ELEMENT_TYPE *pRowCurr = matrix()->GetRow(j);
			for (auto k = colNumb(); k--;)
				lambda += *(pRowCurr + k) * *(pRow + k);

			if (lambda != paramLambda) {
				SPRINTF(buff, "Wrong number of common units in the rows (%lu, %lu): %d != %lu\n", i, j, lambda, getInSys()->lambda());
				outString(buff, outFile());
				throw;
				return false;
			}
		}
	}

	const auto paramK = getInSys()->GetK();
	bool noReplicatedBlockFound = true;
	for (auto i = colNumb(); i--;) {
		int k = 0;
		for (auto j = rowNumb(); j--;)
			k += *(matrix()->GetRow(j) + i);

		if (k != paramK) {
			SPRINTF(buff, "Wrong number of units in the column # %lu: %d != %lu\n", i, k, paramK);
			outString(buff, outFile());
			throw;
			return false;
		}

		if (noReplicatedBlockFound && i) {
			const auto iPrev = i - 1;
			auto j = rowNumb();
			while (j-- && *(matrix()->GetRow(j) + i) == *(matrix()->GetRow(j) + iPrev));
			noReplicatedBlockFound = j != -1;
		}
	}

	pEnumInfo->setSimpleMatrFlag(noReplicatedBlockFound);

#if USE_THREADS < 2
    // For multithread case it's not soo easy to define that we will not construct the BIBDs with no replacated blocks
    // Since this flag is optioanl and it's used for information only, we can skip this procedure.
	if (!(noReplicatedBlockFound || pEnumInfo->constructedAllNoReplBlockMatrix())) {
		// This flag was not set yet
		const VECTOR_ELEMENT_TYPE *pSolution = rowStuff(paramK - 1)->currSolution(); // Is it correct for multi-thread case ???
		pEnumInfo->setNoReplBlockFlag(*pSolution > 1);
	}
#endif
    
	return noReplicatedBlocks() ? noReplicatedBlockFound : true;
}

bool CBIBD_Enumerator::makeFileName(char *buffer, int lenBuffer, const char *ext) const
{
	sprintf_s(buffer, lenBuffer, "%lu_%lu_%lu.%s", rowNumb(), getInSys()->GetK(), getInSys()->lambda(), ext? ext : "txt");
	return true;
}

bool CBIBD_Enumerator::makeJobTitle(char *buffer, int lenBuffer, const char *comment) const
{
	const auto v = rowNumb();
	const auto b = matrix()->colNumb();
	const auto k = getInSys()->GetK();
	sprintf_s(buffer, lenBuffer, "BIBD(%3lu, %3lu, %2lu, %2lu, %2lu)%s", v, b, b * k / v, k, getInSys()->lambda(), comment);
	return true;
}
