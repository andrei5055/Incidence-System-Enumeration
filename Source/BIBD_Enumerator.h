#pragma once
#include "InsSysEnumerator.h"

#if CONSTR_ON_GPU
#define OUT_STRING(buf, len, ...)
#else
#define OUT_STRING(buf, len, ...) { char buf[len]; SPRINTF(buf,  __VA_ARGS__); outString(buf, this->outFile()); }
#endif

Class2Def(CBIBD_Enumerator) : public Class2(C_InSysEnumerator)
{
public:
	CK CBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		Class2(C_InSysEnumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {
		setR(getInSys()->GetR(this->lenStabilizer()));
		C_InSysEnumerator::setFirstUnforcedRow(pBIBD->rowNumb() - pBIBD->GetK());
	}
#if !CONSTR_ON_GPU
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment = "") const;
#endif
protected:
	CK virtual bool checkSolutions(RowSolutionPntr ptr, S nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting = true);
	virtual int unforcedElement(const CColOrbit<S> *p, int nRow) const;
	CK virtual bool checkSolutionsForRight(S nRow, S nPart) const	{ return nRow == firtstNonfixedRowNumber() + 1 && !nPart; }
	CK virtual bool solutionsForRightSideNeeded(const S *pRighPart, const S *pCurrSolution, const VectorPntr pLambdaSet) const;
	bool isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const;
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
	virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
#endif
	CK virtual bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags = NULL, EnumeratorPntr pEnum = NULL) const;
	CK virtual bool checkLambda(size_t val) const					{ return val == lambda(); }
	CK virtual void ReportLamdaProblem(S i, S j, size_t lambda) const {
		OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): %zu != " ME_FRMT "\n",
			i, j, lambda, this->getInSys()->lambda());
	}
	CK virtual const char *getObjName() const                       { return "BIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char *pFrmt = NULL, int *pLambdaSetSize = NULL) const
		{ return SNPRINTF(buffer, lenBuffer, pFrmt, lambda()); }
	int addLambdaInform(const Class1(CVector)* lambdaSet, char* buffer, size_t lenBuffer, int *pLambdaSetSize) const;
	CK virtual void setFirstUnforcedRow(size_t rowNum = 0)			{}
	CK virtual void resetFirstUnforcedRow()							{}
	CK inline auto getR() const										{ return m_r; }
private:
	virtual void getEnumerationObjectKey(char* pInfo, int len) const;
	CK bool checkChoosenSolution(RowSolutionPntr pPrevSolution, S nRow, S nPart, PERMUT_ELEMENT_TYPE usedSolIndex) const;
	CK virtual bool checkForcibleLambda(size_t fLambda) const		 { return checkLambda(fLambda); }
	CK inline auto lambda() const									 { return this->getInSys()->lambda(); }
	CK inline void setR(S val)										 { m_r = val; }
	virtual const char *getTopLevelDirName() const					 { return "BIBDs"; }

	S m_r;
};

FClass2(CBIBD_Enumerator, bool)::checkSolutions(RowSolutionPntr pSolution, S nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting)
{
	if (this->currentRowNumb() + 1 == this->rowNumb())
		return true;        // We will be here when lambda > 1 AND one of the colOrbits was split into 2 parts  
							// in previous row. (6,10,5,3,2) is one such example.

	if (!pSolution->numSolutions())
		return false;

	pSolution->sortSolutions(doSorting, this->useCanonGroup() ? this : NULL);
	if (!pSolution->findFirstValidSolution(this->inSysRowEquation()->variableMaxLimitPntr(), this->GetData()))
		return false;

	return checkChoosenSolution(pSolution, this->currentRowNumb(), nPart, idx);
}

FClass2(CBIBD_Enumerator, bool)::TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags, EnumeratorPntr pEnum) const
{
	const auto iMin = this->lenStabilizer();
	const auto paramR = this->getInSys()->GetR(iMin);
	const auto iMax = this->rowNumb();
	for (S i = iMin; i < iMax; i++) {
		const auto pRow = pMatrix->GetRow(i);
		S r = 0;
		for (auto j = this->colNumb(); j--;)
			r += *(pRow + j);

		if (r != paramR) {
#if !CONSTR_ON_GPU
			pMatrix->printOut(this->outFile());
#endif
			OUT_STRING(buff, 256, "Wrong number of units in the row # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", i, r, paramR);
			THROW();
			return false;
		}

		for (S j = 0; j < i; j++) {
			S lambda = 0;
			const auto pRowCurr = pMatrix->GetRow(j);
			for (auto k = this->colNumb(); k--;)
				lambda += *(pRowCurr + k) * *(pRow + k);

			if (!checkLambda(lambda)) {
				ReportLamdaProblem(i, j, lambda);
				THROW();
				return false;
			}
		}
	}

	const auto paramK = this->getInSys()->GetK();
	bool noReplicatedBlockFound = true;
	for (auto j = this->colNumb(); j--;) {
		S k = 0;
		for (auto i = iMin; i < iMax; i++)
			k += *(pMatrix->GetRow(i) + j);

		if (k != paramK) {
			OUT_STRING(buff, 256, "Wrong number of units in the column # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", j, k, paramK);
			THROW();
			return false;
		}

		if (noReplicatedBlockFound && j) {
			S i = iMin - 1;
			while (++i < iMax) {
				const auto pRow = pMatrix->GetRow(i) + j;
				if (*pRow != *(pRow - 1))
					break;
			}

			noReplicatedBlockFound = i < iMax;
		}
	}

	if (pMatrFlags && noReplicatedBlockFound)
		*pMatrFlags = t_noReplicatedBlock;

#if USE_THREADS < 2
	// For multithread case it's not so easy to define that we will not construct the BIBDs with no replacated blocks
	// Since this flag is optional and it's used for information only, we can skip this procedure.
	if (!(noReplicatedBlockFound || pEnumInfo->constructedAllNoReplBlockMatrix())) {
		// This flag was not set yet
		const auto *pSolution = rowStuff(paramK - 1)->currSolution(); // Is it correct for multi-thread case ???
		pEnumInfo->setNoReplBlockFlag(*pSolution > 1);
	}
#endif

	return this->noReplicatedBlocks() ? noReplicatedBlockFound : true;
}

FClass2(CBIBD_Enumerator, bool)::solutionsForRightSideNeeded(const S *pRighPart, const S *pCurrSolution, const VectorPntr pLambdaSet) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allows to do not consider the solutions for some right parts
	// NOTE: this method should be consistent with CBIBD_Enumerator::checkSolutionsForRight(S nRow, S nPart)

	// Condition should eliminate right part with x1 = 0 for some BIBD, and (8, 4, 3) is one of them.
	// Becase only solutions (and their descendants) with first coordinate equal to 1 could be used for next v - 2 rows.
	// Otherwise, the condition for k will not be satisfied.
	if (!*pRighPart && *pCurrSolution == 1) {
		const auto k = this->getInSys()->GetK();
		const auto lambda = getLambda(pLambdaSet);
		const auto v = this->rowNumb();
		if ((k - 2) * lambda == v - 2)
			return false;
	}

	return true;
}

FClass2(CBIBD_Enumerator, bool)::checkChoosenSolution(RowSolutionPntr pCurrSolution, S nRow, S nPart, PERMUT_ELEMENT_TYPE usedSolIndex) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allows to skip testing of some solutions for some rows
	if (nRow == 3) {
		const auto pPrevSolution = this->rowStuff(nRow - 1);
		const auto lambda = this->getInSys()->lambda();
		// Number of units used for first fragments of the 3-d row
		const auto x0 = *pPrevSolution->solution(usedSolIndex);

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
				const auto k = this->getInSys()->GetK();
				const auto v = this->rowNumb();
				// The difference represent minimal number of descendants of the solution used for 3-d row
				if (lambda * (k - 2) > v) {
					// lambda * (k - 2) - v + 2 > 2 should be use as: lambda * (k - 2) > v
					// Otherwise we do have problem when v,k,l are defined as unsigned
					//
					// This should work for BIBD(22, 8, 4)
					// In this case at least 2 descendants of the solution used for 3-d row
					// will be used for remaining rows
					// It means that the current rows solution should be the descendant of the solution used for 3-d row
					// Suppose that it's not.  Then at lest two rows with the numbers nRow >= 4 will form following
					// configuration:
					//    0 0 1 1
					//    0 0 1 1
					// When we change 3-d and 4-th row with these two AND columns 1 -2 with 3-4 we will increase the code
					// of the matrix. Therefore, this matrix cannot be canonical.
					const bool useCanonGroup = USE_CANON_GROUP && this->groupOrder() > 1;
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

	return pCurrSolution->checkChoosenSolution(this->colOrbit(nRow, nPart), this->matrix()->rowNumb() - nRow, this->getInSys()->GetK());
}


