#pragma once
#include "InsSysEnumerator.h"
#include "IntersectionStorage.h"
#include "designParam.h"

#if CONSTR_ON_GPU
#define OUT_STRING(buf, len, ...)
#else
#define OUT_STRING(buf, len, ...) { char buf[len]; SPRINTF(buf,  __VA_ARGS__); outString(buf, this->outFile()); }
#endif

Class2Def(CBIBD_Enumerator) : public Class2(C_InSysEnumerator)
{
public:
	CK CBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~CBIBD_Enumerator()											{ delete rowIntersection(); }

#if !CONSTR_ON_GPU
	virtual void makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment = "") const;
	CK bool outNonCombinedDesigns(designParam* pParam, const CDesignDB* designDB, std::string& outputInfo, bool addInfo = false);

#endif
	CK inline auto getR() const										{ return m_r; }

protected:
	virtual VariableMappingPntr prepareCheckSolutions(size_t nVar);
	CK virtual bool checkSolutions(RowSolutionPntr ptr, T nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting = true);
	virtual int unforcedElement(const CColOrbit<S> *p, int nRow) const;
	CK virtual bool checkSolutionsForRight(T nRow, T nPart) const	{ return nRow == firtstNonfixedRowNumber() + 1 && !nPart; }
	CK virtual bool solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, const VectorPntr pLambdaSet) const;
	bool isValidSolution(const T* pSol, T λ) const override;
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
	virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
#endif
	CK virtual bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags = NULL, const EnumeratorPntr pEnum = NULL) const;
	CK virtual bool checkLambda(T val) const						{ return val == lambda(); }
	CK virtual void ReportLamdaProblem(T i, T j, T lambda) const {
		OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): " ME_FRMT " != " ME_FRMT "\n",
			i, j, lambda, this->getInSys()->lambda());
	}
	CK virtual const char *getObjName() const                       { return "BIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char *pFrmt = NULL, size_t *pLambdaSetSize = NULL) const
		{ return SNPRINTF(buffer, lenBuffer, pFrmt, lambda()); }
	int addLambdaInform(const Class1(CVector)* lambdaSet, char* buffer, size_t lenBuffer, size_t *pLambdaSetSize) const;
	CK void CreateForcedRows() override;
	CK void CreateAuxiliaryStructures(EnumeratorPntr pMaster) override;
	CK virtual void setFirstUnforcedRow(T rowNum = 0)				{}
	CK virtual void resetFirstUnforcedRow()							{}
	CK inline void setR(T val)										{ m_r = val; }
	CK virtual bool check_X0_3(T nPart) const						{ return !nPart; }
	CK void initDesignDB(const EnumeratorPntr pMaster, size_t rowAdj = 0);
	CK virtual bool outputMaster() const							{ return false; }
	CK void ConstructedDesignProcessing() const override			{ AddMatrixToDB(this); }
	CK void AddMatrixToDB(const CMatrixCanonChecker* pCanonChecker, int rowAdj = 0) const;
	CK virtual bool checkRightParts(T nRow)                         { return nRow == 3; }
	CK virtual bool useAsRightPart(CRowSolution<TDATA_TYPES>* pRowSol, PERMUT_ELEMENT_TYPE idx);
	CK bool createNewFile(const char* fName) const override			{
		return designParams()->find_all_2_decomp != 1 || designParams()->logFile().empty();
	}
	CK bool outFileIsValid(const struct stat& info, const char* pFileName=NULL) const override;

private:
	virtual void getEnumerationObjectKey(char* pInfo, int len) const;
	CK bool checkChoosenSolution(RowSolutionPntr pPrevSolution, T nRow, T nPart, PERMUT_ELEMENT_TYPE usedSolIndex) const;
	CK virtual bool checkForcibleLambda(T fLambda, T nRows, T numPart) const { return nRows == 2? checkLambda(fLambda) : fLambda <= lambda(); }
	CK inline auto lambda() const									 { return this->getInSys()->lambda(); }
	virtual const char *getTopLevelDirName() const					 { return designParams()->find_all_2_decomp == 1? "Quasi-Combined_BIBDs" : 
																		(enumFlags() & t_kSystems? "k-systems" : "BIBDs"); }
	CK bool sharedDB() const										 { return designParams()->threadNumb > 1 && !designParams()->thread_master_DB; }
	CK inline void setRecordLen(size_t val)						     { m_recordLength = val; }
	CK inline auto recordLen() const                                 { return m_recordLength; }
	CK virtual bool lambdaCond(T currentK) const					 { return useLambdaCond() && lambda() == currentK; }
	CK void setUseLambdaCond(T k) {
		// Condition λ <= 2 uses the Hall-Connor Theorems
		m_bUseLambdaCond = lambda() <= 2 || getR() == k;
	}
	CK inline bool useLambdaCond() const							 { return m_bUseLambdaCond; }
	const auto *rowIntersection() const								 { return m_pRowIntersection; }
	inline bool useFilterFor_3d_RowSolutions() const				 { return m_bUseFilterFor_3d_RowSolutions; }

	T m_r;
	size_t m_recordLength = 0;
	bool m_bUseLambdaCond = false;
	CIntersection<T, S>* m_pRowIntersection = NULL;
	bool m_bUseFilterFor_3d_RowSolutions = false;
	T m_Num_3[2];
};

FClass2(CBIBD_Enumerator, bool)::checkSolutions(RowSolutionPntr pSolution, T nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting)
{
	if (this->currentRowNumb() + 1 == this->rowNumb())
		return true;        // We will be here when lambda > 1 AND one of the colOrbits was split into 2 parts
							// in previous row. (6,10,5,3,2) is one such example.

	if (!pSolution->numSolutions())
		return false;

	// NOTES: 
	// 1. If to enumerate combined BIBDs we use a method that allows us not always to build 
	// solutions for all parts (and to use previously constructed sets of solutions for the part 
	// that has not been changed), we can use the group of automorphisms of a partially constructed 
	// design to order sets of solutions ONLY for the first part of the design
	// 2. Besides that, we need to restore the whole sets of the solutions for the first part and 
	// reorder them according to the Aut(M'), where M' differs from M and so Aut(M) differs from Aut(M').
	// First is base for 9 and 10
//	pSolution->sortSolutions(true, this->useCanonGroup() && !nPart ? permStorage(nPart) : NULL);
	// CBIBD(9, 3, { 2, 1, 1 }) +6 designs  CBIBD(8, 4, { 6, 3 })    +1215 designs
	pSolution->sortSolutions(doSorting, numParts() == 1 && this->useCanonGroup() && !nPart ? permStorage(nPart) : NULL);
	if (!pSolution->findFirstValidSolution())
		return false;

	return checkChoosenSolution(pSolution, this->currentRowNumb(), nPart, idx);
}

FClass2(CBIBD_Enumerator, bool)::TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags, const EnumeratorPntr pEnum) const
{
	const auto partsInfo = pMatrix->partsInfo();
	const auto iMin = this->lenStabilizer();
	const auto paramK = this->getInSys()->GetK();
	const auto iMax = this->rowNumb();
	const auto pLambdaSet = this->paramSet(t_lSet);
	char buf[32];
	char* pBuf = NULL;
	for (T n = 0; n < pMatrix->numParts(); n++) {
		const auto colNumb = partsInfo ? partsInfo->colNumb(n) : this->colNumb();
		const auto paramR = partsInfo ? colNumb * paramK / (iMax - iMin) : this->getR();
		const auto lambdaCmp = this->getLambda(pLambdaSet, 0, n);
		for (T i = iMin; i < iMax; i++) {
			T len;
			const auto pRow = pMatrix->GetRow(i, n, &len);
			T r = 0;
			for (auto j = colNumb; j--;)
				r += *(pRow + j);

			if (r != paramR) {
#if !CONSTR_ON_GPU
				pMatrix->printOut(this->outFile());
#endif
				if (partsInfo)
					snprintf(pBuf = buf, sizeof(buf), "Matrix Part # " ME_FRMT ": ", n);

				OUT_STRING(buff, 256, "%sWrong number of units in the row # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", pBuf? pBuf : "", i, r, paramR);

				THROW();
				FCLOSE(outFile());
				return false;
			}

			for (T j = iMin; j < i; j++) {
				T lambda = 0;
				const auto pRowCurr = pMatrix->GetRow(j, n, &len);
				for (auto k = colNumb; k--;)
					lambda += *(pRowCurr + k) * *(pRow + k);

				if (partsInfo) {
					if (blockIdx()) {
						// Testing Kirkman Triple System
						if (lambda <= 1)
							continue;
					}
					else
					if (lambda == lambdaCmp)
						continue;

					snprintf(buf, sizeof(buf), "Matrix Part # " ME_FRMT ": ", n);
					OUT_STRING(buff, 256, "%sWrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): " ME_FRMT " != " ME_FRMT "\n",
						buf, i, j, lambda, lambdaCmp);

				} else {
					if (checkLambda(lambda))
						continue;
					ReportLamdaProblem(i, j, lambda);
				}
				THROW();
				return false;
			}
		}
	}

	bool noReplicatedBlockFound = true;
	for (auto j = this->colNumb(); j--;) {
		T k = 0;
		for (auto i = iMin; i < iMax; i++)
			k += *(pMatrix->GetRow(i) + j);

		if (k != paramK) {
			OUT_STRING(buff, 256, "Wrong number of units in the column # " ME_FRMT ": " ME_FRMT " != " ME_FRMT "\n", j, k, paramK);
			THROW();
			return false;
		}

		if (noReplicatedBlockFound && j) {
			T i = iMin - 1;
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

FClass2(CBIBD_Enumerator, bool)::solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, const VectorPntr pLambdaSet) const
{
	// Collection of conditions to be tested for specific BIBDs, which
	// allows to do not consider the solutions for some right parts
	// NOTE: this method should be consistent with CBIBD_Enumerator::checkSolutionsForRight(T nRow, T nPart)

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

FClass2(CBIBD_Enumerator, bool)::checkChoosenSolution(RowSolutionPntr pCurrSolution, T nRow, T nPart, PERMUT_ELEMENT_TYPE usedSolIndex) const
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
					auto pSol = pCurrSolution;
					while ((pSol = pSol->NextSolution(useCanonGroup())) != NULL) {
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
