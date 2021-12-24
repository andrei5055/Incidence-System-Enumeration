#include "C_tDesignEnumerator.h"
#include "IntersectionStorage.h"
#include "VariableMapping.h"

template class C_tDesignEnumerator<TDATA_TYPES> ;

FClass2(C_tDesignEnumerator)::C_tDesignEnumerator(const TDesignPntr pBIBD, uint enumFlags, int treadIdx, uint nCanonChecker) :
	Class2(CBIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker)
#if USE_EXRA_EQUATIONS
		, CEquSystem(matrix()->colNumb(), matrix()->rowNumb(), pBIBD->getT()), COrbToVar(matrix()->colNumb())
#endif
{
	m_pIntersectionStorage = NULL;

#if USE_EXRA_EQUATIONS
	const size_t numElem = COrbToVar::numElement();
	OrbToVarMapping *pntr = new OrbToVarMapping [numElem];
	for (size_t i = 0; i < numElem; i++)
		COrbToVar::addElement(pntr + i);
#endif
}

FClass2(C_tDesignEnumerator)::~C_tDesignEnumerator()
{
    delete intersectionStorage();
#if USE_EXRA_EQUATIONS
	delete [] COrbToVar::GetAt(0);
#endif
}

#if !CONSTR_ON_GPU
FClass2(C_tDesignEnumerator, bool)::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer);
	const auto t = tDesign()->getT();
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, ME_FRMT"-(" ME_FRMT"_" ME_FRMT"_" ME_FRMT")%s",
						t, this->rowNumb(), this->getInSys()->GetK(), tDesign()->lambda(), ext ? ext : FILE_NAME(""));
	return true;
}

FClass2(C_tDesignEnumerator, bool)::makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment) const
{
	const auto t = tDesign()->getT();
	const auto k = tDesign()->GetK();
	SNPRINTF(buffer, lenBuffer, ME_FRMT"-(%3"  _FRMT", %3"  _FRMT", %2"  _FRMT")%s",
						t, this->rowNumb(), k, tDesign()->lambda(), comment);
	return true;
}
#endif

FClass2(C_tDesignEnumerator, void)::outputTitle(FILE* file) const {
	fprintf(file, "%12s:        %9s:  %9s:       %9s: %9s:      %9s:\n", this->getTopLevelDirName(), "Total #", "Simple #", "Run Time", "Date", "Comments");
}

FClass2(C_tDesignEnumerator, void)::prepareToTestExtraFeatures()
{
	m_pIntersectionStorage = new CIntersectionStorage<S>(tDesign()->getT(), this->rowNumb(), tDesign()->GetNumSet(t_lSet));
}

FClass2(C_tDesignEnumerator, PERMUT_ELEMENT_TYPE *)::getIntersectionParam(const size_t **ppNumb) const
{
	const auto *pPrev = intersectionStorage()->rowsIntersection(this->currentRowNumb());
    *ppNumb = pPrev->numbIntersection();
	return pPrev->rowIntersectionPntr();
}

FClass2(C_tDesignEnumerator, CVariableMapping<T> *)::prepareCheckSolutions(size_t nVar)
{
	if (this->currentRowNumb() <= 1)
		return NULL;			// Nothing to test

	const auto lastRowIdx = this->currentRowNumb() - 1;
	S t = tDesign()->getT() - 2;
	if (t >= this->currentRowNumb())
		t = lastRowIdx;

    const auto *pCurrRow = this->matrix()->GetRow(0);
	const auto *pLastRow = this->matrix()->GetRow(lastRowIdx);
	S tuple[10];
	T *matrixRowPntr[10];
	auto *pTuple = t <= countof(tuple) ? tuple : new S[t];
	auto pMatrixRowPntr = t <= countof(matrixRowPntr) ? matrixRowPntr : new T *[t];

	// Create indices of block containing last element
	const auto r = this->getR();		// TO DO: Need to be modified to support combined t-designs construction
	auto ppBlockIdx = new S[r];
	S idx = 0;
	const auto nCol = this->colNumb();
	for (S j = 0; j < nCol; j++) {
		if (*(pLastRow + j)) {
			*(ppBlockIdx + idx++) = j;
			if (idx == r)
				break;   // No more blocks, containing last element
		}
	}

	// Create indices of block containing 0-th and last, 1-st and last etc element.
	const size_t *pNumb;
	auto *pIntersection = getIntersectionParam(&pNumb);
	const auto lambda = tDesign()->GetNumSet(t_lSet)->GetAt(0);
	for (auto k = pNumb[0]; k--; pCurrRow += nCol, pIntersection += lambda) {
		size_t i = 0;
		for (size_t j = 0; j  < r; j++) {
			if (*(pCurrRow + ppBlockIdx[j])) {
				*(pIntersection + i) = ppBlockIdx[j];
				if (++i == lambda)
					break;	// No more blocks, containing two last elements
			}
		}
	}

	for (uint i = 1; i < t; i++) {
		// Construct all (i+1)-subsets of first currentRowNumb() elements
		uint k = 0;
		pTuple[0] = -1;
		while (k != -1) {
			auto n = pTuple[k];
			for (; k <= i; k++)
				pMatrixRowPntr[k] = this->matrix()->GetRow(pTuple[k] = ++n);

			for (auto j = idx; j--;) {
				const auto blockIdx = ppBlockIdx[j];
				uint m = -1;
				while (++m <= i && *(pMatrixRowPntr[m] + blockIdx));
				if (m > i)
					*pIntersection++ = blockIdx;
			}

			// Construct next (i+1) subset
			k = i + 1;
			n = lastRowIdx;
			while (k-- && pTuple[k] == --n);
		}
	}

	if (pTuple != tuple)
		delete[] pTuple;

	if (pMatrixRowPntr != matrixRowPntr)
		delete[] pMatrixRowPntr;

	delete[] ppBlockIdx;

#if USE_EXRA_EQUATIONS
	return constructExtraEquations(t, nVar);
#else
	return NULL;
#endif
}

FClass2(C_tDesignEnumerator, bool)::isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const
{
#if USE_EXRA_EQUATIONS == 0
	// Check if solution is valid (for elimination of invalid solutions)
	if (this->currentRowNumb() <= 1)
		return true;		// Nothing to test

	auto t = tDesign()->getT();
	if (t >= this->currentRowNumb() + 2)
		t = this->currentRowNumb() + 1;

	CMatrixCanonChecker::MakeRow(pSol, false);

	const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
	const auto *pCurrRow = this->matrix()->GetRow(this->currentRowNumb());
	const size_t *pNumb;
	auto *pIntersection = getIntersectionParam(&pNumb);
	auto lambda = pLambdaSet->GetAt(0);
    for (uint i = 2; i < t; i++) {
		const auto lambdaPrev = lambda;
        lambda = pLambdaSet->GetAt(i-1);
		for (size_t k = pNumb[i - 2]; k--; pIntersection += lambdaPrev) {
            size_t val = 0;
			for (auto j = lambdaPrev; j--;) {
				if (*(pCurrRow + *(pIntersection + j)))
                    val++;
            }

            if (val != lambda)
                return false;
        }
    }
    
//	OUTPUT_MATRIX(matrix(), outFile(), currentRowNumb() + 1);
#endif
    return true;
}

#if USE_EXRA_EQUATIONS
FClass2(C_tDesignEnumerator, CVariableMapping *)::constructExtraEquations(size_t t, size_t nVar)
{
	const CColOrbit *pColOrbit = colOrbit(currentRowNumb());
	if (!pColOrbit)		// This could happend for last 1-2 rows
		return NULL;		

	const CColOrbit *pColOrbitIni = colOrbitIni(currentRowNumb());

	resetVarPtr(nVar);

	CColOrbit *pUnforcedIni = forcibleLambda(currentRowNumb()) != -1 ? unforcedOrbits(currentRowNumb())[0] : NULL;

	const auto pKSet = tDesign()->GetNumSet(t_kSet);
	const auto k = pKSet->GetAt(0);
	const auto pLambdaSet = tDesign()->GetNumSet(t_lSet);
	const size_t *pNumb;
	auto *pIntersection = getIntersectionParam(&pNumb);
	auto lambda = pLambdaSet->GetAt(0);
	VECTOR_ELEMENT_TYPE equIdx = 0;
	for (size_t i = 0; i < t; i++) {
		const auto lambdaPrev = lambda;
		lambda = pLambdaSet->GetAt(i + 1);
		// Loop over the equations for current lambda as their lambda right parts
		for (auto n = pNumb[i]; n--; equIdx++) {
			
			CVariable *pFirstVar, *pCurrVar;
			pFirstVar = pCurrVar = NULL;

			size_t j, nVarInEquation, idxMap, idx, adjLambda;
			j = nVarInEquation = idxMap = idx = adjLambda = 0;

			while (true) {
				// idx - index of the first block of the next block's orbit, containing current (i+2)-tuple of elements
				// nCol - starting column (first block number) of the next block's orbit
				const size_t nCol = *(pIntersection + idx);
				const CColOrbit *pColOrbitLast = (CColOrbit *)((char *)pColOrbitIni + nCol * colOrbitLen());

				const OrbToVarMapping *pVarMap;
				do {
					pVarMap = COrbToVar::GetAt(idxMap++);
				} while (pVarMap->pColOrb < pColOrbitLast);

				if (pVarMap->pColOrb == pColOrbitLast) {			
					const size_t numVar = pVarMap->nVar;
					pCurrVar = new CVariable(equIdx, numVar, pCurrVar);

					if (!pFirstVar)
						pFirstVar = pCurrVar;

					addVariable(pCurrVar, numVar);
					nVarInEquation++;
				}
				else {
					// Forcible orbits found
					if (pColOrbitLast->columnWeight() != k)
						adjLambda++;	// it's a all last 1's orbit

					// to use previous orb-to-var mapping one more time 
					idxMap--;    
				}

				if ((j += pColOrbitLast->length()) < lambdaPrev)
					idx += pColOrbitLast->length();
				else
					break;
			}

			if (pFirstVar) {
				// Make a loop of variables appearing in the current equation
				pFirstVar->linkVariable(pCurrVar);
				addEquation(pFirstVar, nVarInEquation, lambda - adjLambda);
			}

			pIntersection += lambdaPrev;
		}
	}

	closeVarloop();
	return solveExtraEquations();
}

FClass2(C_tDesignEnumerator, void)::addColOrbitForVariable(S nVar, CColOrbit *pColOrb)
{
	OrbToVarMapping *pOrbToVar = COrbToVar::GetAt(COrbToVar::numb());
	pOrbToVar->nVar = nVar;
	pOrbToVar->pColOrb = pColOrb;
	COrbToVar::incNumb();
}

#endif

FClass2(C_tDesignEnumerator, void)::copyInfoFromMaster(const EnumeratorPntr pMaster)
{
	prepareToTestExtraFeatures();
	S t = tDesign()->getT() - 2;
	if (t == 1)
		return;			// Nothing to copy for 3-design

	if (t >= this->currentRowNumb())
		t = this->currentRowNumb() - 1;

	t--;
	const size_t *pNumb;
	// Note: since currentRowNumb() is different for this and pMaster, we need
	//   a) call getIntersectionParam differently
	//   b) call pMaster->getIntersectionParam() last
	const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
	PERMUT_ELEMENT_TYPE *pIntersectionTo = getIntersectionParam(&pNumb);
	const auto *pIntersectionFrom = dynamic_cast<const Class2(C_tDesignEnumerator) *>(pMaster)->getIntersectionParam(&pNumb);
	for (uint i = 0; i < t; i++) {
		const size_t len = pNumb[i] * pLambdaSet->GetAt(i);
		memcpy(pIntersectionTo, pIntersectionFrom, len * sizeof(*pIntersectionTo));
		pIntersectionTo += len;
		pIntersectionFrom += len;
	}
}

FClass2(C_tDesignEnumerator, bool)::TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags, const EnumeratorPntr pEnum) const
{
	if (!Class2(CBIBD_Enumerator)::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	const auto t = tDesign()->getT();
	const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
	T *rowPntr[10];
	const auto pRow = t <= countof(rowPntr)? rowPntr : new T *[t];
	S tupleIdx[10];
	auto pTuple = t <= countof(tupleIdx)? tupleIdx : new S[t];
	bool retVal = true;
	size_t i = 2; 
	while (retVal && ++i <= t) {
		const size_t t_lambda = pLambdaSet->GetAt(i - 2);
		size_t k = 0;
		pTuple[0] = -1;
		do {
			// Create next tuple
			auto val = pTuple[k];
			for (size_t j = k; j < i; j++)
				pRow[j] = pMatrix->GetRow(pTuple[j] = ++val);

			// Check constructed tuple:
			size_t lambda = 0;
			for (auto n = this->colNumb(); n--;) {
				auto j = i;
				while (j-- && *(pRow[j] + n));
				if (j == -1)
					lambda++;
			}

			if (lambda != t_lambda) {
				retVal = false;
				break;
			}
				
			// Find first element of the tuple to be changed
			k = i;
			val = this->rowNumb();
			while (k-- && pTuple[k] == --val);
		} while (k != -1);
	}

	if (pTuple != tupleIdx)
		delete[] pTuple;

	if (pRow != rowPntr)
		delete[] pRow;

	return retVal;
}

/*
1. System of equations for 2-d (starting from 0) row of 3-design is:
x0 + x1 = labda_2
x2 + x3 = r - lambda_2
x0 + x2 = lambda_2
x0      = lambda_3

First three equations are regular equations for 2-design, last one is the conditions for intersection of first three row.
The solution of (1) is:
x0 = lamda_3
x1 = lamba_2 - lambda_3
x2 = lamba_2 - lambda_3
x3 = r - 2 * lambda_2 + lambda_3

2. System of equations for 3-d (starting from 0) row of 3-design is:
x0 + x1 = labda_3
x2 + x3 = lambda_2 - lambda_3
x4 + x5 = lambda_2 - lambda_3
x6 + x7 = r - 2 * lambda_2 + lambda_3
x0 + x2 + x4 + x6 = lambda2
x0 + x1 = lambda_3
x0 + x2 = lambda_3
x0 + x4 = lambda_3

First five equations are regular equations for 2-design, last three is the conditions for intersection of 3-d row with (0,1), (0,2) and (1,2)
a) We can remove equation 6, since 1st is exactly the same.
x1 = x2 = x4 ==> from last 3 equation

Now we can re-write 2-d and 5-th equations as (2a)
x1 + x3 = lambda_2 - lambda_3
x1 + x6 = lambda_2 - lambda_3

Therefore, x3 = x5 = x6 ==> from 2-d and 3-d of (2) and (2a)
Now (2) has the same solutions as (2b)
x0 + x1 = lambda_3
x1 + x3 = lambda_2 - lambda_3
x3 + x7 = r - 2 * lambda_2 + lambda_3

Obviously, 0 <= x7 <= b - (r  + r - lambda_2 + r - 2 * lambda_2 + lambda_3) = b - 3 * (r - lambda_2) - lambda_3
We know that 
b = v*(v-1)*lambda_2/(k*(k-1)) 
r = (v-1)*lambda_2/(k-1)

And we can see that x3 = x5 = x6
x0 = a
x1 = x2 = x4 = lambda_3 - a
x3 = x5 = x6 = lambda_2 - 2 * lambda_3 + a
x7 = r - 3 * (lambda_2 - lambda_3) - a

*/