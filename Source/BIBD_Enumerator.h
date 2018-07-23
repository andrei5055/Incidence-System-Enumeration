#pragma once
#include "InsSysEnumerator.h"

template<class T>
class CBIBD_Enumerator : public C_InSysEnumerator<T>
{
public:
	CK CBIBD_Enumerator(const C_BIBD<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		C_InSysEnumerator(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
#if !CONSTR_ON_GPU
	virtual bool makeJobTitle(char *buffer, int lenBuffer, const char *comment = "") const;
#endif
	CK virtual bool isTDesign_enumerator(size_t t)  const { return t <= 2; }

protected:
	CK virtual bool sortSolutions(CRowSolution<T> *ptr, PERMUT_ELEMENT_TYPE idx);
	virtual int unforcedElement(const CColOrbit<T> *p, int nRow) const;
	CK virtual bool solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, size_t nRow) const;
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
#endif
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL) const;
private:
	CK bool checkChoosenSolution(CRowSolution<T> *pPrevSolution, size_t nRow, PERMUT_ELEMENT_TYPE usedSolIndex) const;
	CK virtual bool checkForcibleLambda(size_t fLambda) const { return fLambda == getInSys()->lambda(); }
};

template<class T>
bool CBIBD_Enumerator<T>::sortSolutions(CRowSolution<T> *pSolution, PERMUT_ELEMENT_TYPE idx)
{
	if (currentRowNumb() + 1 == rowNumb())
		return true;        // We will be here when lambda > 1 AND one of the colOrbits was splitted into 2 parts  
							// in previous row. (6,10,5,3,2) is one such example.

	if (!pSolution->numSolutions())
		return false;

	pSolution->sortSolutions(useCanonGroup() ? this : NULL);
	if (!pSolution->findFirstValidSolution(inSysRowEquation()->variableMaxLimitPntr(), GetData()))
		return false;

	return checkChoosenSolution(pSolution, currentRowNumb(), idx);
}

#if CONSTR_ON_GPU
#define OUT_STRING(buf, len, ...)
#else
#define OUT_STRING(buf, len, ...) { char buf[len]; SPRINTF(buf,  __VA_ARGS__); outString(buf, outFile()); }
#endif

template<class T>
bool CBIBD_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags) const
{
	const auto paramR = getInSys()->GetR();
	const auto paramLambda = getInSys()->lambda();
	const auto iMax = rowNumb();
	for (T i = 0; i < iMax; i++) {
		const auto pRow = pMatrix->GetRow(i);
		T r = 0;
		for (auto j = colNumb(); j--;)
			r += *(pRow + j);

		if (r != paramR) {
#if !CONSTR_ON_GPU
			(static_cast<const CMatrix<T> *>(pMatrix))->printOut(outFile());
#endif
			OUT_STRING(buff, 256, "Wrong number of units in the row # "ME_FRMT": "ME_FRMT" != "ME_FRMT"\n", i, r, getInSys()->GetR());
			THROW();
			return false;
		}

		for (T j = 0; j < i; j++) {
			T lambda = 0;
			const auto pRowCurr = pMatrix->GetRow(j);
			for (auto k = colNumb(); k--;)
				lambda += *(pRowCurr + k) * *(pRow + k);

			if (lambda != paramLambda) {
				OUT_STRING(buff, 256, "Wrong number of common units in the rows ("ME_FRMT", "ME_FRMT"): "ME_FRMT" != "ME_FRMT"\n", i, j, lambda, getInSys()->lambda());
				THROW();
				return false;
			}
		}
	}

	const auto paramK = getInSys()->GetK();
	bool noReplicatedBlockFound = true;
	for (auto i = colNumb(); i--;) {
		T k = 0;
		for (auto j = rowNumb(); j--;)
			k += *(pMatrix->GetRow(j) + i);

		if (k != paramK) {
			OUT_STRING(buff, 256, "Wrong number of units in the column # "ME_FRMT": "ME_FRMT" != "ME_FRMT"\n", i, k, paramK);
			THROW();
			return false;
		}

		if (noReplicatedBlockFound && i) {
			const auto iPrev = i - 1;
			auto j = rowNumb();
			while (j-- && *(pMatrix->GetRow(j) + i) == *(pMatrix->GetRow(j) + iPrev));
			noReplicatedBlockFound = j != MATRIX_ELEMENT_MAX;
		}
	}

	if (pMatrFlags)
		*pMatrFlags = noReplicatedBlockFound ? 1 : 0;

	//	pEnumInfo->setSimpleMatrFlag(noReplicatedBlockFound);

#if USE_THREADS < 2
	// For multithread case it's not so easy to define that we will not construct the BIBDs with no replacated blocks
	// Since this flag is optioanl and it's used for information only, we can skip this procedure.
	if (!(noReplicatedBlockFound || pEnumInfo->constructedAllNoReplBlockMatrix())) {
		// This flag was not set yet
		const auto *pSolution = rowStuff(paramK - 1)->currSolution(); // Is it correct for multi-thread case ???
		pEnumInfo->setNoReplBlockFlag(*pSolution > 1);
	}
#endif

	return noReplicatedBlocks() ? noReplicatedBlockFound : true;
}

template<class T>
bool CBIBD_Enumerator<T>::solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, size_t nRow) const
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

template<class T>
bool CBIBD_Enumerator<T>::checkChoosenSolution(CRowSolution<T> *pCurrSolution, size_t nRow, PERMUT_ELEMENT_TYPE usedSolIndex) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allow to skip testing of some solutions for some rows
	if (nRow == 3) {
		const CRowSolution<T> *pPrevSolution = rowStuff(nRow - 1);
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
				if (lambda * (k - 2) > v) {
					// lambda * (k - 2) - v + 2 > 2 should be use as: lambda * (k - 2) > v
					// Othervise we do have problem when v,k,l are defined as unsigned
					//
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
					const bool useCanonGroup = USE_CANON_GROUP && groupOrder() > 1;
					auto pSol = pCurrSolution;
					while ((pSol = pSol->NextSolution(useCanonGroup)) != NULL) {
						const auto *pSolValue = pSol->currSolution();
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


