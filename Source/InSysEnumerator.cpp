#include "InsSysEnumerator.h"
#include "RightPartFilter.h"
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
C_InSysEnumerator::C_InSysEnumerator(const C_InSys *pInSys, bool matrOwner, bool noReplicatedBlocks) : 
					m_bNoReplBlock(noReplicatedBlocks),
					CEnumerator(pInSys, matrOwner), CInSysSolver(pInSys->colNumb() >> 1), CVector(pInSys->colNumb())
{
	const auto nCol = pInSys->colNumb();
	setRowEquation(new CInSysRowEquation(nCol));
	setX0_3(nCol);
	m_pRightPartFilter = new CRightPartFilter(nCol);
	const auto nRow = pInSys->rowNumb();
	setForcibleLambdaPntr(new size_t[nRow]);
	memset(forcibleLambdaPntr(), 0, nRow * sizeof(*forcibleLambdaPntr()));
	setFirstUnforcedRow();
	setForcibleLambda(nRow - 1, getInSys()->lambda());
}

C_InSysEnumerator::~C_InSysEnumerator()
{
	delete rightPartFilter();
	delete[] forcibleLambdaPntr();
}

CRowSolution *C_InSysEnumerator::setFirstRowSolutions()
{
	CRowSolution *pSolutions = rowStuff(0);
	const C_InSys *pISys = getInSys();
	const CVector *pR_set = pISys->GetNumSet(t_rSet);
	pSolutions->InitSolutions(pR_set->GetSize(), 1);
	for (size_t i = 0; i < pR_set->GetSize(); i++)
		pSolutions->AddElement(pR_set->GetAt(i));

	return pSolutions;
}

size_t C_InSysEnumerator::MakeSystem()
{
	size_t nVar = 0;
	// Total number of equations (some of them corresponds to the forcibly constructed columns)
	VECTOR_ELEMENT_TYPE equationIdx = -1;
	// Number of equations corresponding only to the colums which ARE NOT forcibly constructed
	VECTOR_ELEMENT_TYPE eqIdx = 0;

	CInSysRowEquation *pRowEquation = inSysRowEquation();

	pRowEquation->resetMappings();
	const auto nRow = currentRowNumb();

	// When we are using the group of canonical matrix, we need to adjust previously constructed
	// generators of the automorphism group on columns (because some of them could be forcibly constructed
	setUseCanonGroup(USE_CANON_GROUP && canonChecker()->groupOrder() > 1);

	int buffer[256], *pColGroupIdx = NULL;
    size_t lenPerm;
	if (useCanonGroup()) {
		lenPerm = canonChecker()->lenPerm();
		pColGroupIdx = lenPerm <= countof(buffer) ? buffer : new int[lenPerm];
	} else
		lenPerm = 0;

	const C_InSys *pIS = getInSys();
	const auto nRowToBuild = pIS->rowNumb() - nRow;
	const auto k = pIS->GetK();
	rightPartFilter()->reset();

	int colGroupIdx = 0;
	const CColOrbit *pColOrbPrev = colOrbit(nRow - 1);
	CColOrbit *pColOrbit, *pColOrbitNext = colOrbit(nRow);
	CColOrbit *pPrev = NULL;
	while ((pColOrbit = pColOrbitNext) != NULL) {
		pColOrbitNext = pColOrbit->next();
		equationIdx++;
		const auto diffWeight = k - pColOrbit->colomnWeight();
		if (!diffWeight || diffWeight == nRowToBuild) {
			colGroupIdx++;			// We need to skip this column's orbit
			// All remaining elements of all columns
			// of the current orbit should be equal to 0 or 1, respectively
			const auto currLen = pColOrbit->lenght();
			if (noReplicatedBlocks() && currLen > 1)
				break;

			addForciblyConstructedColOrbit(pColOrbit, pPrev, diffWeight ? 1 : 0);
			if (currLen == pColOrbPrev->lenght()) {
				rightPartFilter()->addFilter(equationIdx, diffWeight ? currLen : 0);
				pColOrbPrev = pColOrbPrev->next();
				continue;
			}

			pColOrbitNext = (pColOrbit = pColOrbitNext)->next();

			// When we are here, only first column(s) of pColOrbPrev orbits should be forcible and diffWeight = 0
			// Limit for right part should be established here 
			rightPartFilter()->addFilter(equationIdx, 0, nVar);
		}

		if (pColGroupIdx)
			pColGroupIdx[nVar] = colGroupIdx++;

		int mapSetIdx;
		const auto currLen = pColOrbit->lenght();
		setVariableLimit(nVar, currLen, nRowToBuild, k, pColOrbit->colomnWeight());
		if (pColOrbit->colomnWeight() != pColOrbPrev->colomnWeight()) {
			auto len = pColOrbPrev->lenght() - currLen;
			// Column weight was changed. It means that x[i] > 0
			if (len) {
				// The length of the current orbit is NOT the same.
				// It means that current x[i] < xMax[i] (orbit was splited in two sub-orbits)
				// Check if second sub-orbit should contain all unit
				if (diffWeight + 1 == nRowToBuild) {
					if (noReplicatedBlocks() && len > 1)
						break;

					CColOrbit *pTmp;
					pColOrbitNext = (pTmp = pColOrbitNext)->next();
					pTmp->Init(len, pColOrbitNext);

					addForciblyConstructedColOrbit(pTmp, pColOrbit, 1);
					rightPartFilter()->addFilter(equationIdx, len, nVar);
					colGroupIdx++;    // This group ol columns will be skipped
					mapSetIdx = 1;
				} else {
					setVariableLimit(++nVar, len, nRowToBuild, k, pColOrbPrev->colomnWeight());
					pColOrbitNext = (pColOrbit = pColOrbitNext)->next();
					if (pColGroupIdx)
						pColGroupIdx[nVar] = colGroupIdx++;

					mapSetIdx = 2;
				}
			}
			else {
				// The length of the current orbit is the same.
				// It means that current x[i] = xMax[i] (orbit was not splited in two sub-orbits)
				mapSetIdx = 1;
			}
		} else
			mapSetIdx = 0;

		pRowEquation->addVarMapping(mapSetIdx, nVar++, eqIdx++);
		pColOrbPrev = pColOrbPrev->next();
		pPrev = pColOrbit;
	}

	if (!pColOrbit) {
		if (nRowToBuild > 1) {
			auto fLambda = forcibleLambda(nRow - 1);
			CColOrbit *pTmp = *(currUnforcedOrbPtr() + 1);
			while (pTmp) {
				fLambda += pTmp->lenght();
				pTmp = pTmp->next();
			}

			if (nRowToBuild != 2 || checkForcibleLambda(fLambda))
				setForcibleLambda(nRow, fLambda);
			else
				return nVar = -1;
		}
	} else
		nVar = SIZE_MAX;

	if (useCanonGroup()) {
		if (nVar < lenPerm && nRowToBuild > 1) {
			// We realy need to adjust the generators of the automorphism group on columns
			canonChecker()->adjustGenerators(pColGroupIdx, nVar);
		}

		if (pColGroupIdx != buffer)
			delete[] pColGroupIdx;
	}

	return nVar;
}

CRowSolution *C_InSysEnumerator::FindSolution(size_t nVar, PERMUT_ELEMENT_TYPE lastRightPartIndex)
{
	const auto nRow = currentRowNumb();
	const auto nRowPrev = nRow - 1;
	const CRowSolution *pPrevRowSolution = rowStuff(nRowPrev);
	CRowSolution *pCurrRowSolution = rowStuff(nRow);
	pCurrRowSolution->setSolutionSize(nVar);
	pCurrRowSolution->resetSolution();
	const VECTOR_ELEMENT_TYPE *pFirstSol = pPrevRowSolution->firstSolution();
	const CSolutionPerm *pSolPerm = pPrevRowSolution->solutionPerm();
	const PERMUT_ELEMENT_TYPE *pPerm = pSolPerm ? pSolPerm->GetData() : NULL;
	const auto len = pPrevRowSolution->solutionSize();
	CInSysRowEquation *pRowEquation = inSysRowEquation();

	// We will keep the variation intervals of "lambda" variables from the 2-variable equations
	// Because there are no more than nVar/2 of them, we will use
	initSolver(pCurrRowSolution, pRowEquation->variableMinValPntr());

	const auto pLambdaSet = getInSys()->GetNumSet(t_lSet);
	// For all possible right parts of the systems
	// (which are lexicographically not greater, than the vector, used for the current row)
	VECTOR_ELEMENT_TYPE *pResult = pCurrRowSolution->newSolution();
	VECTOR_ELEMENT_TYPE *buffer = new VECTOR_ELEMENT_TYPE[nVar];
	// Solution used for last constracted matrix's row:
	const VECTOR_ELEMENT_TYPE *pCurrSolution = pPrevRowSolution->solution(lastRightPartIndex);
	const auto forcibleLambdaValue = forcibleLambda(nRowPrev);
    const size_t nLambdas = constructing_t_Design()? 1 : pLambdaSet->GetSize();
	bool readyToCheckSolution = false;
	for (PERMUT_ELEMENT_TYPE j = 0; j <= lastRightPartIndex; j++) {
		if (!pPrevRowSolution->isValidSolution(j))
			continue;

		const PERMUT_ELEMENT_TYPE idx = pPerm ? *(pPerm + j) : j;
		const VECTOR_ELEMENT_TYPE *pRightSide = rightPartFilter()->getRightPart(pFirstSol + idx * len, inSysRowEquation()->variableMaxValPntr(), len, buffer);
		if (!pRightSide)
			continue;

		if (!solutionsForRightSideNeeded(pRightSide, pCurrSolution, nRow))
			continue;

		// Resolve trivial equations and define lambdaMin (part of any current lambda
		// which is already splited among the lambda-variables
		resetMapping();
		int lambdaMin = pRowEquation->resolveTrivialEquations(pRightSide, pResult, nVar, this);
		if (lambdaMin < 0)
			continue;

		lambdaMin += forcibleLambdaValue;

		if (!readyToCheckSolution) {
			readyToCheckSolution = true;
			prepareCheckSolutions(nVar);
		}

		// For all possible values of right part of "lambda" equation
		for (size_t i = 0; i < nLambdas; i++) {
			const int lambdaToSplit = pLambdaSet->GetAt(i) - lambdaMin;
			if (lambdaToSplit < 0)
				continue;

			VECTOR_ELEMENT_TYPE *pResultTmp = findAllSolutionsForLambda(pResult, lambdaToSplit);
			if (pResultTmp) {
				pResult = pResultTmp;
				continue;
			}

			// There are no solution for current lambda from pLambdaSet
			// Since lambda set is ordered, we will not find solutions for
			// remaining lambdas as well
			break;
		}
	}

	delete[] buffer;
	return  pCurrRowSolution->getSolution();
}

bool C_InSysEnumerator::sortSolutions(CRowSolution *pSolution, size_t i)
{
	if (!pSolution->numSolutions())
		return false;

	return pSolution->findFirstValidSolution(inSysRowEquation()->variableMaxValPntr());
}

void C_InSysEnumerator::setVariableLimit(size_t nVar, size_t len, size_t nRowToBuild, size_t k, int colWeight)
{
	// Minimal value for the nVar-th element which could be used as a first valid candidate for next row
	SetAt(nVar, ((k - colWeight) * len + nRowToBuild - 1) / nRowToBuild);
	inSysRowEquation()->setVariableMaxLimit(nVar, len);

	// Maximal value of any Xi cannot be bigger then X0_3, when i-th group of column
	// already have at least 2 inits in each column. Otherwise, these two rows with the row
	// which would use the soution we eliminating by using X0_3 as a maximum,
	// would be bigger than current 3 rows
	if (len > getX0_3() && colWeight >= 2)
		len = getX0_3();

	if (noReplicatedBlocks() && getInSys()->GetK() - colWeight == 1)
		len = 1;

	inSysRowEquation()->setVariableMaxVal(nVar, len);
}

void C_InSysEnumerator::addForciblyConstructedColOrbit(CColOrbit *pColOrbit, CColOrbit *pPrev, int idx)
{
	if (pPrev)
		pPrev->setNext(pColOrbit->next());
	else
		setColOrbitCurr(pColOrbit->next());

	// All remaining elements of all columns
	// of the current orbit should be equal to 0 or 1, respectively
	CColOrbitManager::addForciblyConstructedColOrbit(pColOrbit, idx);

	// if it's possible, define first unforced row number
	if (idx > 0 && !firstUnforcedRow())
		setFirstUnforcedRow(currentRowNumb());
}

void C_InSysEnumerator::resetFirstUnforcedRow()
{	
	if (firstUnforcedRow() == currentRowNumb())
		setFirstUnforcedRow(0);
}
