#include "C_tDesignEnumerator.h"
#include "IntersectionStorage.h"
#include "VariableMapping.h"

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

C_tDesignEnumerator::C_tDesignEnumerator(const C_tDesign *pBIBD, bool matrOwner, bool noReplicatedBlocks) : CBIBD_Enumerator(pBIBD, matrOwner, noReplicatedBlocks)
{
	m_pIntersectionStorage = NULL;
	m_pVarPntr = new CVariable *[matrix()->colNumb() << 1];

	// For 3 design it will be enough and we won't realloc that memory
	resetEquNumb();
	allocateMemoryForEquations(matrix()->rowNumb());
}

C_tDesignEnumerator::~C_tDesignEnumerator()
{
    delete intersectionStorage();
	releaseVariables();
	delete [] varPntr();
	for (auto i = equNumbMax(); i--;)
		delete equation(i);

	delete [] equPntr();
	delete [] equIndx();
}

bool C_tDesignEnumerator::makeFileName(char *buffer, int lenBuffer, const char *ext) const
{
	const auto t = tDesign()->getT();
	sprintf_s(buffer, lenBuffer, "%lu-(%lu_%lu_%lu).%s", t, rowNumb(), getInSys()->GetK(), tDesign()->lambda(), ext? ext : "txt");
	return true;
}

bool C_tDesignEnumerator::makeJobTitle(char *buffer, int lenBuffer, const char *comment) const
{
	const auto t = tDesign()->getT();
	const auto k = tDesign()->GetK();
	sprintf_s(buffer, lenBuffer, "%lu-(%3lu, %3lu, %2lu)%s", t, rowNumb(), k, tDesign()->lambda(), comment);
	return true;
}

void C_tDesignEnumerator::prepareToTestExtraFeatures()
{
	m_pIntersectionStorage = new CIntersectionStorage(tDesign()->getT(), rowNumb(), tDesign()->GetNumSet(t_lSet));
}

PERMUT_ELEMENT_TYPE *C_tDesignEnumerator::getIntersectionParam(const size_t **ppNumb) const
{
    const auto *pIntStorage = intersectionStorage();
    const auto *pPrev = pIntStorage->rowsIntersection(currentRowNumb());
    *ppNumb = pPrev->numbIntersection();
	return pPrev->rowIntersectionPntr();
}

void C_tDesignEnumerator::prepareCheckSolutions(size_t nVar)
{
	if (currentRowNumb() <= 1)
		return;			// Nothing to test

	const auto lastRowIdx = currentRowNumb() - 1;
	auto t = tDesign()->getT() - 2;
	if (t >= currentRowNumb())
		t = lastRowIdx;

    const auto *pCurrRow = matrix()->GetRow(0);
	const auto *pLastRow = matrix()->GetRow(lastRowIdx);
	const size_t *pNumb;
	PERMUT_ELEMENT_TYPE tuple[10];
	auto *pTuple = t <= countof(tuple) ? tuple : new PERMUT_ELEMENT_TYPE[t];
	MATRIX_ELEMENT_TYPE *matrixRowPntr[10];
	auto *pMatrixRowPntr = t <= countof(matrixRowPntr) ? matrixRowPntr : new MATRIX_ELEMENT_TYPE *[t];
	auto *pIntersection = getIntersectionParam(&pNumb);

	// Create indices of block containing last element
	PERMUT_ELEMENT_TYPE *ppBlockIdx = new PERMUT_ELEMENT_TYPE[tDesign()->GetR()];
	size_t idx = 0;
	for (auto j = colNumb(); j--;)
		if (*(pLastRow + j))
			*(ppBlockIdx + idx++) = j;

	const auto nCol = colNumb();
    for (size_t i = 0; i < t; i++) {
		if (i == 0) {
			for (auto k = pNumb[i]; k--; pCurrRow += nCol) {
				for (auto j = idx; j--;) {
					if (*(pCurrRow + ppBlockIdx[j]))
						*pIntersection++ = ppBlockIdx[j];
				}
			}
		} else {
			// Construct all (i+1)-subsets of first currentRowNumb() elements
			size_t k = 0;
			pTuple[0] = -1;
			while (k != -1) {
				size_t n = pTuple[k];
				for (; k <= i; k++) {
					pTuple[k] = ++n;
					pMatrixRowPntr[k] = matrix()->GetRow(n);
				}

				for (auto j = idx; j--;) {
					const auto blockIdx = ppBlockIdx[j];
					size_t k = -1;
					while (++k <= i && *(pMatrixRowPntr[k] + blockIdx));
					if (k > i)
						*pIntersection++ = blockIdx;
				}

				// Construct next (i+1) subset
				k = i + 1;
				n = lastRowIdx;
				while (k-- && pTuple[k] == --n);
			}
		}
    }

	if (pTuple != tuple)
		delete[] pTuple;

	if (pMatrixRowPntr != matrixRowPntr)
		delete[] pMatrixRowPntr;

	delete[] ppBlockIdx;

	constructAdditionalEquations(t, nVar);
}

bool C_tDesignEnumerator::isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const 
{ 
	if (currentRowNumb() <= 1)
		return true;		// Nothing to test

	auto t = tDesign()->getT() - 2;
	if (t >= currentRowNumb())
		t = currentRowNumb() - 1;

    MakeRow(pSol);

    const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
    const auto *pCurrRow = matrix()->GetRow(currentRowNumb());
	const size_t *pNumb;
	auto *pIntersection = getIntersectionParam(&pNumb);
	auto lambda = pLambdaSet->GetAt(0);
    for (size_t i = 0; i < t; i++) {
		const auto lambdaPrev = lambda;
        lambda = pLambdaSet->GetAt(i+1);
		for (size_t k = pNumb[i]; k--; pIntersection += lambdaPrev) {
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
    return true;
}

void C_tDesignEnumerator::releaseVariables()
{
	for (auto i = equNumb(); i--;)
		equation(i)->releaseEquation();

	resetEquNumb();
}

void C_tDesignEnumerator::setEquation(CEquation *pEqu, size_t idx)
{ 
	*(equPntr() + idx) = pEqu;

	// Change equIndx to be able find this equation by it's initial index  
	const auto pVar = pEqu->firstVariable();
	*(equIndx() + pVar->equIdx()) = idx;
}

void C_tDesignEnumerator::setEquation(CEquation *pEquFrom, const CEquation *pEquTo)
{
	const auto pVar = pEquTo->firstVariable();
	setEquation(pEquFrom, *(equIndx() + pVar->equIdx()));
}

void C_tDesignEnumerator::constructAdditionalEquations(size_t t, size_t nVar)
{
	const CColOrbit *pColOrbit = colOrbit(currentRowNumb());
	if (!pColOrbit)		// This could happend for last 1-2 rows
		return;		

	const CColOrbit *pColOrbitIni = colOrbitIni(currentRowNumb());

	releaseVariables();
	memset(varPntr(), 0, nVar * sizeof(*varPntr()));

	const auto pLambdaSet = tDesign()->GetNumSet(t_lSet);
	const size_t *pNumb;
	auto *pIntersection = getIntersectionParam(&pNumb);
	auto lambda = pLambdaSet->GetAt(0);
	VECTOR_ELEMENT_TYPE equIdx = 0;
	for (size_t i = 0; i < t; i++) {
		const auto lambdaPrev = lambda;
		lambda = pLambdaSet->GetAt(i + 1);
		// Loop over the equations for current lambda as their right parts
		for (auto k = pNumb[i]; k--; equIdx++) {
			CVariable *pFirstVar, *pCurrVar, *pPrevVar;
			pFirstVar = pCurrVar = pPrevVar = NULL;
			const CColOrbit *pColOrbitTmp = pColOrbit;
			size_t numVar, nVarInEquation;
			nVarInEquation = numVar = 0;
			size_t j = 0;
			while (j < lambdaPrev) {
				const size_t nCol = *(pIntersection + j);
				const void *pColOrbitLast = (char *)pColOrbitIni + nCol * colOrbitLen();
				while (pColOrbitTmp < pColOrbitLast) {
					numVar++;
					pColOrbitTmp = pColOrbitTmp->next();
				}

				nVarInEquation++;
				pCurrVar = new CVariable(equIdx, numVar, pCurrVar);

				if (!pFirstVar)
					pFirstVar = pCurrVar;

				if (*(varPntr() + numVar))
					pCurrVar->setSameVarNextEqu(*(varPntr() + numVar + nVar));
				else
					*(varPntr() + numVar) = pCurrVar;

				*(varPntr() + numVar + nVar) = pCurrVar;
				j += pColOrbitTmp->lenght();
			}

			// Make a loop of variables appearing in current equation
			pFirstVar->linkVariable(pCurrVar);
			addEquation(pFirstVar, nVarInEquation, lambda);
		}
	}

	// Make the loop of all appearance of given variable in different equations
	for (auto i = nVar; i--;) {
		auto pVar = varPntr()[i];
		if (pVar)
			pVar->setSameVarNextEqu(varPntr()[i + nVar]);
	}

	solveAdditionalEquations();
}

int C_tDesignEnumerator::solveAdditionalEquations()
{
	CVariableMapping VarValue(10), *pVarValue = &VarValue;
	pVarValue->resetMapping();
	// Loop over all constructed equations

	auto iMax = equNumb();
	size_t i = 0;
	while (i < iMax) {
		CEquation *pEquation = equation(i);
		if (pEquation->numbVar() == 1) {
			// We found trivial equation
			const auto pVariable = pEquation->firstVariable();
			const auto varValue = pEquation->rightPart();
			
			pVarValue->addMapping(pVariable->varIndex(), varValue);

			// Get next equation with the current variable
			// (we needd to do it BEFORE releasing this variable)
			auto pVariableNext = pVariable->sameVarNextEqu();
			// Release current equation and copy last equation of the system to the current slot
			pEquation->releaseEquation();
			if (i < --iMax)
				setEquation(equation(iMax), i);

			CVariable *pVariableTmp, *pVarTmp, *pTmp;
			while (pVariableNext != pVariable) {
				pVariableNext = (pVariableTmp = pVariableNext)->sameVarNextEqu();

				CEquation *pCurrEquation = equation(pVariableTmp);
				if (pCurrEquation->rightPart() < varValue)
					return -1;    // no solution
				
				auto *pVarNext = pVariableTmp->nextVarSameEqu();
				if (pCurrEquation->rightPart() == varValue) {
					// The values of all remaining variables of current equation are 0's
					while (pVarNext != pVariableTmp) {
						pVarNext = (pVarTmp = pVarNext)->nextVarSameEqu();
						// Map variable with it's value 0
						pVarValue->addMapping(pVarTmp->varIndex(), 0);
						// ... and remove this variable from remaining equations
						auto *pNxt = pVarTmp->sameVarNextEqu();
						while (pNxt != pVarTmp) {
							pNxt = (pTmp = pNxt)->sameVarNextEqu();
							pTmp->prevVarSameEqu()->linkVariable(pTmp->nextVarSameEqu());
							delete pTmp;
						}
						delete pVarTmp;
					}

					// Current equation already solved and we need to delete it and move the last equation to its place
					setEquation(equation(--iMax), pCurrEquation);
				}
				else {
					if (pVarNext == pVariableTmp)  // Only one variable in this equation
						return -1;     

					// Adjust right part
					pCurrEquation->adjustRightPart(varValue);
					// Remove variable
					pTmp->prevVarSameEqu()->linkVariable(pVarNext);
					delete pVariableTmp;
				}
			}
		} else 
			i++;
	}

	setEquNumb(iMax);
	return 0;
}

void C_tDesignEnumerator::allocateMemoryForEquations(size_t nEqu)
{
	setEquNumbMax(nEqu);
	m_pEquPntr = new CEquation *[equNumbMax()];
	m_pEquIndx = new VECTOR_ELEMENT_TYPE[equNumbMax()];
	for (auto i = equNumb(); i < nEqu; i++)
		m_pEquPntr[i] = new CEquation();
}

void C_tDesignEnumerator::addEquation(CVariable *pFirstVar, size_t nVar, size_t lambda)
{
	if (equNumb() == equNumbMax()) {
		auto *pTmp = equPntr();
		auto *pIdx = equIndx();
		allocateMemoryForEquations(equNumbMax() << 1);
		memcpy(equPntr(), pTmp, equNumb() * sizeof(*pTmp));
		memcpy(equIndx(), pIdx, equNumb() * sizeof(*pIdx));
		delete[] pTmp;
		delete[] pIdx;
	}

	equPntr()[m_nEquNumb]->initEquation(pFirstVar, nVar, lambda);
	equIndx()[m_nEquNumb] = m_nEquNumb;
	m_nEquNumb++;
}

void C_tDesignEnumerator::copyInfoFromMaster(const CEnumerator *pMaster)
{
	prepareToTestExtraFeatures();
	auto t = tDesign()->getT() - 2;
	if (t == 1)
		return;			// Nothing to copy for 3-design

	if (t >= currentRowNumb())
		t = currentRowNumb() - 1;

	t--;
	const size_t *pNumb;
	// Note: since currentRowNumb() is different for this and pMaster, we need
	//   a) call getIntersectionParam differently
	//   b) call pMaster->getIntersectionParam() last
	const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
	PERMUT_ELEMENT_TYPE *pIntersectionTo = getIntersectionParam(&pNumb);
	const auto *pIntersectionFrom = dynamic_cast<const C_tDesignEnumerator *>(pMaster)->getIntersectionParam(&pNumb);
	for (size_t i = 0; i < t; i++) {
		const size_t len = pNumb[i] * pLambdaSet->GetAt(i);;
		memcpy(pIntersectionTo, pIntersectionFrom, len * sizeof(*pIntersectionTo));
		pIntersectionTo += len;
		pIntersectionFrom += len;
	}
}

bool C_tDesignEnumerator::TestFeatures(CEnumInfo *pEnumInfo)
{
	if (!CBIBD_Enumerator::TestFeatures(pEnumInfo))
		return false;

	const auto t = tDesign()->getT();
	const auto *pLambdaSet = tDesign()->GetNumSet(t_lSet);
	const auto pRow = new MATRIX_ELEMENT_TYPE *[t];
	auto *pTuple = new PERMUT_ELEMENT_TYPE[t];
	bool retVal = true;
	size_t i = 2; 
	while (retVal && ++i <= t) {
		const size_t t_lambda = pLambdaSet->GetAt(i - 2);
		size_t k = 0;
		pTuple[0] = -1;
		do {
			// Create next tuple
			size_t val = pTuple[k];
			for (size_t j = k; j < i; j++)
				pRow[j] = matrix()->GetRow(pTuple[j] = ++val);

			// Check constructed tuple:
			size_t lambda = 0;
			for (auto n = colNumb(); n--;) {
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
			val = rowNumb();
			while (k-- && pTuple[k] == --val);
		} while (k != -1);
	}

	delete[] pTuple;
	delete[] pRow;
	return retVal;
}