#pragma once
#include "Enumerator.h"
#include "RightPartFilter.h"

IClass2Def(_InSysEnumerator) : public IClass2(Enumerator), public IClass2(InSysSolver), public CVector<S>
{
public:
	CK C_InSysEnumerator(const InSysPntr pInSysm, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_InSysEnumerator();
	CK virtual S getX0_3() const								{ return m_x0_3; }
	CK virtual size_t firstUnforcedRow() const                  { return m_firstUnforcedRow; }
	CK virtual void setFirstUnforcedRow(size_t rowNum = 0)      { m_firstUnforcedRow = rowNum; }
	CK virtual S *forcibleLambdaPntr() const					{ return m_pForsibleLambda; }
	CK virtual bool noReplicatedBlocks() const					{ return m_bNoReplBlock; }
	CK virtual bool isPBIB_enumerator() const					{ return false;  }
	CK int define_MT_level(int v) const							{ return v / 2; }
	CK int define_MT_level(const designParam *pParam) const		{ return pParam->lambda()[0] == 1?
																		 pParam->v / pParam->k : define_MT_level(pParam->v); }
protected:
	CK virtual void setX0_3(S value)							{ m_x0_3 = value; }
	CK virtual bool sortSolutions(RowSolutionPntr ptr, size_t idx);
	CK inline CInSysRowEquation<S> * inSysRowEquation() const	{ return (CInSysRowEquation<S> *)this->rowEquation(); }
	CK virtual bool solutionsForRightSideNeeded(const S *pRighPart, const S *pCurrSolution, size_t nRow) const
																{ return true; }
	CK virtual CEquSystem *equSystem()							{ return NULL;  }
	CK CColOrbit<S> **unforcedOrbits(size_t nRow) const			{ return getUnforcedColOrbPntr() + this->rank() * nRow; }
	CK virtual CColOrbit<S> **getUnforcedColOrbPntr() const		{ return forcibleLambda(this->currentRowNumb()) != -1 ? this->unforcedColOrbPntr() : NULL; }
	CK virtual S forcibleLambda(size_t i) const					{ return *(forcibleLambdaPntr() + i); }
	CK virtual void setColOrbitForCurrentRow(CColOrbit<S> *pColOrb){}
	CK virtual void addColOrbitForVariable(size_t nVar, CColOrbit<S> *pColOrb)	{}
	CK virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix);
	virtual void CanonizeByColumns(MatrixDataPntr pMatrix, S *pColIdxStorage, CanonicityCheckerPntr pCanonChecker = NULL, bool permCol = false) const;

private:
	CK void addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, int idx);
	CK virtual RowSolutionPntr setFirstRowSolutions();
	CK virtual size_t MakeSystem();
	CK virtual RowSolutionPntr FindSolution(size_t nVar, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX);
	CK void setVariableLimit(size_t nVar, S len, size_t nRowToBuild, size_t k, int colWeight);
	CK virtual bool checkForcibleLambda(size_t fLambda) const   { return true; }
	CK virtual void resetFirstUnforcedRow()						{ if (firstUnforcedRow() == this->currentRowNumb())
																	setFirstUnforcedRow(0);
																}
	virtual CVariableMapping<T> *prepareCheckSolutions(size_t n){ return NULL; }
    CK virtual bool constructing_t_Design()                     { return false; }
	CK inline CRightPartFilter<S> *rightPartFilter()			{ return m_pRightPartFilter; }
	CK inline void setForcibleLambdaPntr(S *p)					{ m_pForsibleLambda = p; }
	CK inline void setForcibleLambda(S i, S val)				{ *(forcibleLambdaPntr() + i) = val; } // i < nCol, val <= nCol

	S m_x0_3;
	S *m_pForsibleLambda;
	size_t m_firstUnforcedRow;
	CRightPartFilter<S> *m_pRightPartFilter;
	const bool m_bNoReplBlock;
};

TClass2(_InSysEnumerator)::C_InSysEnumerator(const InSysPntr pInSys, uint enumFlags, int treadIdx, uint nCanonChecker) :
	m_bNoReplBlock(enumFlags & t_noReplicatedBlocks), IClass2(Enumerator)(pInSys, enumFlags | t_IS_enumerator, treadIdx, nCanonChecker),
	IClass2(InSysSolver)(pInSys->colNumb() >> 1, pInSys->GetT()), CVector<S>(pInSys->colNumb())
{
	const auto nCol = pInSys->colNumb();
	this->setRowEquation(new CInSysRowEquation<S>(nCol, pInSys->GetT() > 2));
	setX0_3(nCol);
	m_pRightPartFilter = new CRightPartFilter<S>(nCol);
	const auto nRow = pInSys->rowNumb();
	setForcibleLambdaPntr(new S[nRow]);
	memset(forcibleLambdaPntr(), 0, nRow * sizeof(*forcibleLambdaPntr()));
	setFirstUnforcedRow();
	setForcibleLambda(nRow - 1, this->getInSys()->lambda());
}

TClass2(_InSysEnumerator)::~C_InSysEnumerator() {
	delete rightPartFilter();
	delete[] forcibleLambdaPntr();
}

TClass2(_InSysEnumerator, bool)::sortSolutions(RowSolutionPntr pSolution, size_t i) {
	if (!pSolution->numSolutions())
		return false;

	return pSolution->findFirstValidSolution(inSysRowEquation()->variableMaxValPntr());
}

TClass2(_InSysEnumerator, RowSolutionPntr)::setFirstRowSolutions() {
	auto *pSolutions = this->rowStuff(0);
	const auto *pISys = this->getInSys();
	const auto *pR_set = pISys->GetNumSet(t_rSet);
	pSolutions->InitSolutions(pR_set->GetSize(), 1);
	for (size_t i = 0; i < pR_set->GetSize(); i++)
		pSolutions->AddElement(pR_set->GetAt(i));

	return pSolutions;
}

TClass2(_InSysEnumerator, size_t)::MakeSystem()
{
	VECTOR_ELEMENT_TYPE nVar = 0;
	// Total number of equations (some of them corresponds to the forcibly constructed columns)
	VECTOR_ELEMENT_TYPE equationIdx = -1;// (1 << sizeof(S) - 1;
	// Number of equations corresponding only to the colums which ARE NOT forcibly constructed
	VECTOR_ELEMENT_TYPE eqIdx = 0;

	auto pRowEquation = inSysRowEquation();

	pRowEquation->resetMappings();
	const auto nRow = this->currentRowNumb();

	// When we are using the group of canonical matrix, we need to adjust previously constructed
	// generators of the automorphism group on columns (because some of them could be forcibly constructed
	this->setUseCanonGroup(USE_CANON_GROUP && this->groupOrder() > 1);

	int buffer[256], *pColGroupIdx = NULL;
	size_t lenPermut;
	if (this->useCanonGroup()) {
		lenPermut = this->lenPerm();
		pColGroupIdx = lenPermut <= countof(buffer) ? buffer : new int[lenPermut];
	}
	else
		lenPermut = 0;

	const auto *pIS = this->getInSys();
	const auto nRowToBuild = pIS->rowNumb() - nRow;
	const auto k = pIS->GetK();
	rightPartFilter()->reset();

	int colGroupIdx = 0;
	const auto *pColOrbPrev = this->colOrbit(nRow - 1);
	CColOrbit<S> *pColOrbit, *pColOrbitNext = this->colOrbit(nRow);
#if USE_EXRA_EQUATIONS
	setColOrbitForCurrentRow(pColOrbitNext);
	const size_t nextRowOrbShift = pIS->colNumb() * colOrbitLen();
#endif

	CColOrbit<S> *pPrev = NULL;
	while ((pColOrbit = pColOrbitNext) != NULL) {
		pColOrbitNext = pColOrbit->next();
		equationIdx++;
		const auto diffWeight = k - pColOrbit->columnWeight();
		if (!diffWeight || diffWeight == nRowToBuild) {
			colGroupIdx++;			// We need to skip this column's orbit
									// All remaining elements of all columns
									// of the current orbit should be equal to 0 or 1, respectively
			const auto currLen = pColOrbit->length();
			if (noReplicatedBlocks() && currLen > 1)
				break;

			addForciblyConstructedColOrbit(pColOrbit, pPrev, diffWeight ? 1 : 0);
			if (currLen == pColOrbPrev->length()) {
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
		const auto currLen = pColOrbit->length();
		setVariableLimit(nVar, currLen, nRowToBuild, k, pColOrbit->columnWeight());
		if (pColOrbit->columnWeight() != pColOrbPrev->columnWeight()) {
			addColOrbitForVariable(nVar, pColOrbit);
			auto len = pColOrbPrev->length() - currLen;
			// Column weight was changed. It means that x[i] > 0
			if (len) {
				// The length of the current orbit is NOT the same.
				// It means that current x[i] < xMax[i] (the orbit was splitted in two sub-orbits)
				// Check if second sub-orbit should contain all units
				if (diffWeight + 1 == nRowToBuild) {
					if (noReplicatedBlocks() && len > 1)
						break;

					auto pTmp = pColOrbitNext;
					pColOrbitNext = pTmp->next();
					pTmp->Init(len, pColOrbitNext);

					addForciblyConstructedColOrbit(pTmp, pColOrbit, 1);
					rightPartFilter()->addFilter(equationIdx, len, nVar);
					colGroupIdx++;    // This group of columns will be skipped
					mapSetIdx = 1;
#if USE_EXRA_EQUATIONS
					// Assign the length of current orbits to all its descendant
					const size_t lenOrb = pTmp->length();
					size_t i = 0;
					while (++i < nRowToBuild) {
						pTmp = (CColOrbit *)((char *)pTmp + nextRowOrbShift);
						pTmp->setLength(lenOrb);
					}
#endif
				}
				else {
					setVariableLimit(++nVar, len, nRowToBuild, k, pColOrbPrev->columnWeight());
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
		}
		else
			mapSetIdx = 0;

		pRowEquation->addVarMapping(mapSetIdx, nVar++, eqIdx++);
		pColOrbPrev = pColOrbPrev->next();
		pPrev = pColOrbit;
	}

	if (!pColOrbit) {
		if (nRowToBuild > 1) {
			auto fLambda = forcibleLambda(nRow - 1);
			auto *pTmp = *(this->currUnforcedOrbPtr() + 1);
			while (pTmp) {
				fLambda += pTmp->length();
				pTmp = pTmp->next();
			}

			if (nRowToBuild != 2 || checkForcibleLambda(fLambda))
				setForcibleLambda(nRow, fLambda);
			else
				return (size_t)-1;
		}
	}
	else
		nVar = ELEMENT_MAX;

	if (this->useCanonGroup()) {
		if (nVar < lenPermut && nRowToBuild > 1) {
			// We realy need to adjust the generators of the automorphism group on columns
			this->adjustGenerators(pColGroupIdx, nVar);
		}

		if (pColGroupIdx != buffer)
			delete[] pColGroupIdx;
	}

	return nVar;
}

TClass2(_InSysEnumerator, RowSolutionPntr)::FindSolution(size_t nVar, PERMUT_ELEMENT_TYPE lastRightPartIndex)
{
	const auto nRow = this->currentRowNumb();
	const auto nRowPrev = nRow - 1;
	const auto pPrevRowSolution = this->rowStuff(nRowPrev);
	auto pCurrRowSolution = this->rowStuff(nRow);
	pCurrRowSolution->setSolutionSize(nVar);
	pCurrRowSolution->resetSolution();
	const auto pFirstSol = pPrevRowSolution->firstSolution();
	const CSolutionPerm *pSolPerm = pPrevRowSolution->solutionPerm();
	const PERMUT_ELEMENT_TYPE *pPerm = pSolPerm ? pSolPerm->GetData() : NULL;
	const auto len = pPrevRowSolution->solutionSize();
	auto pRowEquation = inSysRowEquation();
	pRowEquation->setEquSystem(equSystem());

	// We will keep the variation intervals of "lambda" variables from the 2-variable equations
	// Because there are no more than nVar/2 of them, we will use
	this->initSolver(pCurrRowSolution, pRowEquation->variableMinValPntr());

	const auto pLambdaSet = this->getInSys()->GetNumSet(t_lSet);
	// For all possible right parts of the systems
	// (which are lexicographically not greater, than the vector, used for the current row)
	auto pResult = pCurrRowSolution->newSolution();
	VECTOR_ELEMENT_TYPE buffer[256];
	auto pBuffer = nVar <= countof(buffer)? buffer : new VECTOR_ELEMENT_TYPE[nVar];
	// Solution used for last constructed matrix's row:
	const auto pCurrSolution = pPrevRowSolution->solution(lastRightPartIndex);
	const auto forcibleLambdaValue = forcibleLambda(nRowPrev);
	const size_t nLambdas = constructing_t_Design() ? 1 : pLambdaSet->GetSize();
	const auto *pMaxVal = inSysRowEquation()->variableMaxValPntr();
	bool readyToCheckSolution = false;
	for (PERMUT_ELEMENT_TYPE j = 0; j <= lastRightPartIndex; j++) {
		if (!pPrevRowSolution->isValidSolution(j))
			continue;

		const auto idx = pPerm ? *(pPerm + j) : j;
		const auto pRightSide = rightPartFilter()->getRightPart(pFirstSol + idx * len, pMaxVal, len, pBuffer);
		if (!pRightSide)
			continue;

		if (!solutionsForRightSideNeeded(pRightSide, pCurrSolution, nRow))
			continue;

		// Resolve trivial equations and define lambdaMin (part of any current lambda
		// which is already splited among the lambda-variables)
		this->resetMapping();
		int lambdaMin = pRowEquation->resolveTrivialEquations(pRightSide, pResult, nVar, this);
		if (lambdaMin < 0)
			continue;

		lambdaMin += (int)forcibleLambdaValue;

#if USE_EXRA_EQUATIONS
		resetExtra();
#endif
		CVariableMapping<T> *pVarMapping = NULL;
		if (!readyToCheckSolution) {
			readyToCheckSolution = true;
			pVarMapping = prepareCheckSolutions(nVar);
#if USE_EXRA_EQUATIONS
			// we need to revert the list of the forcibly constructed orbits
			//	addForciblyConstructedColOrbit(pColOrbit, pPrev, diffWeight ? 1 : 0);
			//			CColOrbitManager::addForciblyConstructedColOrbit(pColOrbit, idx);
#endif
		}

#if USE_EXRA_EQUATIONS
		const size_t shiftToBaseVar = getLastMapping() - getMappingPntr() + 2;
		const size_t nTotal = 0;
		size_t nElem = pVarMapping ? pVarMapping->nElement() : 0;
		int diffLambda = 0;
		CVariableMapping *pVarValue = pVarMapping;

		CVariableMapping *pRightPartAdjValues = NULL;
		bool flag = true;
		size_t shift = 0;
		while (flag || nElem) {
			const size_t n = getVarValues()->nElement();
			if (nElem) {
				diffLambda = pRowEquation->resolveExtraMapping(pVarValue, pMaxVal, pResult, this, getVarValues(), shift);

				if (diffLambda < 0)
					break;

				lambdaMin += diffLambda; // Solution constructed
			}

			if (flag) {
				flag = false;
				// Adjust the right parts of the extra equations
				pRightPartAdjValues = pRowEquation->AdjustExtraEquationRightParts();
				pRowEquation->addDefinedVariables(getVarValues(), pResult, this);
			}

			if (n == getVarValues()->nElement())
				break;

			shift = nElem = pRowEquation->excludeVariables(pVarValue = getVarValues());
		}

		if (diffLambda >= 0) {
			/*
			const size_t nAddedVar = pVarMapping->nElement() - nTotal;
			if (false) //nAddedVar)
			addVarIdxToStack(nAddedVar, shiftToBaseVar);
			*/

			setEquSystem(equSystem());
#endif

			// For all possible values of right part of "lambda" equation
			for (size_t i = 0; i < nLambdas; i++) {
				const int lambdaToSplit = pLambdaSet->GetAt(i) - lambdaMin;
				if (lambdaToSplit < 0)
					continue;

#if USE_EXRA_EQUATIONS
				equSystem()->resetVariables(nVar);
#endif
				VECTOR_ELEMENT_TYPE *pResultTmp = this->findAllSolutionsForLambda(pResult, lambdaToSplit);
				if (!pResultTmp) {
					// There are no solution for current lambda from pLambdaSet
					// Since lambda set is ordered, we will not find solutions for
					// remaining lambdas as well
					break;
				}

				pResult = pResultTmp;
			}

#if USE_EXRA_EQUATIONS
		}

		// We need to restore right parts of extra equations
		if (pRightPartAdjValues) {
			// First restore them for those leading variables from t_dual equations, 
			// which have minimal values > 0
			printResults(NULL, lambdaMin, -1, -1);
			pRowEquation->AdjustExtraEquationRightParts(pRightPartAdjValues);
			delete pRightPartAdjValues;
			printFinalResultForRightPart(getVarValues());
		}

		// Finaly restore them from all other variables
		pRowEquation->AdjustExtraEquationRightParts(getVarValues(), true);
		printResults(NULL, lambdaMin, -1, -1);
#endif
	}

	if (pBuffer != buffer)
		delete[] pBuffer;

	return  pCurrRowSolution->getSolution();
}

TClass2(_InSysEnumerator, void)::addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, int idx)
{
	if (pPrev)
		pPrev->setNext(pColOrbit->next());
	else
		this->setColOrbitCurr(pColOrbit->next());

	// All remaining elements of all columns
	// of the current orbit should be equal to 0 or 1, respectively
	CColOrbitManager<S>::addForciblyConstructedColOrbit(pColOrbit, idx);
#if USE_EXRA_EQUATIONS
	equSystem()->addForcibleOrb(pColOrbit);
#endif

	// if it's possible, define first unforced row number
	if (idx > 0 && !firstUnforcedRow())
		setFirstUnforcedRow(this->currentRowNumb());
}

TClass2(_InSysEnumerator, void)::setVariableLimit(size_t nVar, S len, size_t nRowToBuild, size_t k, int colWeight)
{
	// Minimal value for the nVar-th element which could be used as a first valid candidate for next row
	this->SetAt(nVar, static_cast<S>(((k - colWeight) * len + nRowToBuild - 1) / nRowToBuild));
	inSysRowEquation()->setVariableMaxLimit(nVar, len);

	// Maximal value of any Xi cannot be bigger than X0_3, when i-th group of columns
	// already have at least 2 inits in each column. Otherwise, these two rows with the row
	// which would use the solution, would be bigger than current 3 rows
	if (len > getX0_3() && colWeight >= 2)
		len = getX0_3();

	if (noReplicatedBlocks() && this->getInSys()->GetK() - colWeight == 1)
		len = 1;

	inSysRowEquation()->setVariableMaxVal(nVar, len);
}

