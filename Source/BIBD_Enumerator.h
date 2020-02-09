#pragma once
#include "InsSysEnumerator.h"

#if CONSTR_ON_GPU
#define OUT_STRING(buf, len, ...)
#else
#define OUT_STRING(buf, len, ...) { char buf[len]; SPRINTF(buf,  __VA_ARGS__); outString(buf, this->outFile()); }
#endif

template<class T>
class CBIBD_Enumerator : public C_InSysEnumerator<T>
{
public:
	CK CBIBD_Enumerator(const C_InSys<T> *pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		C_InSysEnumerator<T>(pBIBD, enumFlags, treadIdx, nCanonChecker) {
		setR(getInSys()->GetR());
		m_firstUnforced = pBIBD->rowNumb() - pBIBD->GetK();
	}
#if !CONSTR_ON_GPU
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment = "") const;
#endif
protected:
	CK virtual bool sortSolutions(CRowSolution<T> *ptr, PERMUT_ELEMENT_TYPE idx);
	virtual int unforcedElement(const CColOrbit<T> *p, int nRow) const;
	CK virtual bool solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, size_t nRow) const;
	bool isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const;
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
	virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
#endif
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const;
	CK virtual bool checkLambda(size_t val) const					{ return val == lambda(); }
	CK virtual void ReportLamdaProblem(T i, T j, size_t lambda) const {
		OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): %zu != " ME_FRMT "\n",
			i, j, lambda, this->getInSys()->lambda());
	}
	CK virtual const char *getObjName() const                        { return "BIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, int *pLambdaSetSize = NULL) const 
		{ return SNPRINTF(buffer, lenBuffer, "%2" _FRMT, lambda()); }
	CK virtual size_t firstUnforcedRow() const						{ return m_firstUnforced; }
	CK virtual void setFirstUnforcedRow(size_t rowNum = 0)			{}
	CK virtual void resetFirstUnforcedRow()							{}
private:
	virtual void getEnumerationObjectKey(char* pInfo, int len) const;
	CK bool checkChoosenSolution(CRowSolution<T> *pPrevSolution, size_t nRow, PERMUT_ELEMENT_TYPE usedSolIndex) const;
	CK virtual bool checkForcibleLambda(size_t fLambda) const		 { return checkLambda(fLambda); }
	CK inline auto lambda() const									 { return this->getInSys()->lambda(); }
	CK inline void setR(T val)										 { m_r = val; }
	CK inline auto getR() const                                      { return m_r; }
	virtual const char *getTopLevelDirName() const					 { return "BIBDs"; }

	T m_r;
	size_t m_firstUnforced;
};

template<class T>
bool CBIBD_Enumerator<T>::sortSolutions(CRowSolution<T> *pSolution, PERMUT_ELEMENT_TYPE idx)
{
	if (this->currentRowNumb() + 1 == this->rowNumb())
		return true;        // We will be here when lambda > 1 AND one of the colOrbits was splitted into 2 parts  
							// in previous row. (6,10,5,3,2) is one such example.

	if (!pSolution->numSolutions())
		return false;

	pSolution->sortSolutions(this->useCanonGroup() ? this : NULL);
	if (!pSolution->findFirstValidSolution(this->inSysRowEquation()->variableMaxLimitPntr(), this->GetData()))
		return false;

	return checkChoosenSolution(pSolution, this->currentRowNumb(), idx);
}

template<class T>
bool CBIBD_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags, CEnumerator<T> *pEnum) const
{
	const auto paramR = this->getInSys()->GetR();
	const auto iMax = this->rowNumb();
	for (T i = 0; i < iMax; i++) {
		const auto pRow = pMatrix->GetRow(i);
		T r = 0;
		for (auto j = this->colNumb(); j--;)
			r += *(pRow + j);

		if (r != paramR) {
#if !CONSTR_ON_GPU
			(static_cast<const CMatrix<T> *>(pMatrix))->printOut(this->outFile());
#endif
			OUT_STRING(buff, 256, "Wrong number of units in the row # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", i, r, this->getInSys()->GetR());
			THROW();
			return false;
		}

		for (T j = 0; j < i; j++) {
			size_t lambda = 0;
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
	for (auto i = this->colNumb(); i--;) {
		T k = 0;
		for (auto j = this->rowNumb(); j--;)
			k += *(pMatrix->GetRow(j) + i);

		if (k != paramK) {
			OUT_STRING(buff, 256, "Wrong number of units in the column # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", i, k, paramK);
			THROW();
			return false;
		}

		if (noReplicatedBlockFound && i) {
			const auto iPrev = i - 1;
			auto j = this->rowNumb();
			while (j-- && *(pMatrix->GetRow(j) + i) == *(pMatrix->GetRow(j) + iPrev));
			noReplicatedBlockFound = j != MATRIX_ELEMENT_MAX;
		}
	}

	if (pMatrFlags && noReplicatedBlockFound)
		*pMatrFlags = t_noReplicatedBlock;

	//	pEnumInfo->setSimpleMatrFlag(noReplicatedBlockFound);

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

template<class T>
bool CBIBD_Enumerator<T>::solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, size_t nRow) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allow to do not consider the solutions for some right parts
	if (nRow == 3) {
		// Condition should eliminate right part with x1 = 0 for some BIBD, and (8, 4, 3) is one of them
		if (!*pRighPart && *pCurrSolution == 1) {
			const auto k = this->getInSys()->GetK();
			const auto lambda = this->getInSys()->lambda();
			const auto v = this->rowNumb();
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
	// allows to skip testing of some solutions for some rows
	if (nRow == 3) {
		const CRowSolution<T> *pPrevSolution = this->rowStuff(nRow - 1);
		const auto lambda = this->getInSys()->lambda();
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

	return pCurrSolution->checkChoosenSolution(this->colOrbit(nRow), this->matrix()->rowNumb() - nRow, this->getInSys()->GetK());
}


