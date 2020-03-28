#pragma once
#include "Enumerator.h"
#include "RightPartFilter.h"

Class2Def(C_InSysEnumerator) : public Class2(CEnumerator), public Class2(CInSysSolver), public CVector<S>
{
public:
	CK C_InSysEnumerator(const InSysPntr pInSysm, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_InSysEnumerator();
	CK virtual S getX0_3() const								{ return m_x0_3; }
	CK virtual S firstUnforcedRow() const						{ return m_firstUnforcedRow; }
	CK virtual void setFirstUnforcedRow(S rowNum = 0)			{ m_firstUnforcedRow = rowNum; }
	CK virtual S *forcibleLambdaPntr(S nRow = 0) const			{ return m_pForsibleLambda + nRow * numParts(); }
	CK virtual bool noReplicatedBlocks() const					{ return m_bNoReplBlock; }
	CK virtual bool isPBIB_enumerator() const					{ return false;  }
	CK int define_MT_level(int v) const							{ return v / 2; }
	CK int define_MT_level(const designParam *pParam) const		{ return pParam->lambda()[0] == 1?
																		 pParam->v / pParam->k : define_MT_level(pParam->v); }
protected:
	CK virtual void setX0_3(S value)							{ m_x0_3 = value; }
	CK virtual bool checkSolutions(RowSolutionPntr ptr, PERMUT_ELEMENT_TYPE idx, bool doSorting = true);
	CK inline auto inSysRowEquation() const						{ return (CInSysRowEquation<S> *)this->rowEquation(); }
	CK virtual S getLambda(const VectorPntr pLambdaSet, S idx = 0, S numPart = 0) const { return pLambdaSet->GetAt(idx); }
	CK virtual bool checkSolutionsForRight(S row, S part) const	{ return false; }
	CK virtual bool solutionsForRightSideNeeded(const S *pRighPart, const S *pCurrSolution, const VectorPntr pLambdaSet) const
																{ return true; }
	CK virtual CEquSystem *equSystem()							{ return NULL;  }
	CK CColOrbit<S> **unforcedOrbits(size_t nRow, S iPart = 0) const	{ return getUnforcedColOrbPntr(iPart) + this->rank() * nRow; }
	CK virtual CColOrbit<S> **getUnforcedColOrbPntr(S iPart = 0) const {
			return forcibleLambda(this->currentRowNumb(), iPart) != ELEMENT_MAX ? this->unforcedColOrbPntr(iPart) : NULL;
	}
	CK virtual S forcibleLambda(S nRow, S nPart) const			{ return *(forcibleLambdaPntr(nRow) + nPart); }
	CK virtual void setColOrbitForCurrentRow(CColOrbit<S> *pColOrb){}
	CK virtual VectorPntr paramSet(t_numbSetType idx) const		{ return this->getInSys()->GetNumSet(idx); }
#if USE_EXRA_EQUATIONS
	CK virtual void addColOrbitForVariable(S nVar, CColOrbit<S> *pColOrb)	{}
#else
#define addColOrbitForVariable(...)
#endif
	CK virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix);
	virtual void CanonizeByColumns(MatrixDataPntr pMatrix, S *pColIdxStorage, CanonicityCheckerPntr pCanonChecker = NULL, bool permCol = false) const;
private:
	CK void addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, int idx);
	CK virtual RowSolutionPntr setFirstRowSolutions();
	CK virtual S MakeSystem(S numPart);
	CK virtual RowSolutionPntr FindSolution(S nVar, S nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX);
	CK void setVariableLimit(S nVar, S len, S nRowToBuild, S k, S colWeight);
	CK virtual bool checkForcibleLambda(size_t fLambda) const   { return true; }
	CK virtual void resetFirstUnforcedRow()						{ if (firstUnforcedRow() == this->currentRowNumb())
																	setFirstUnforcedRow(0);
																}
	virtual CVariableMapping<T> *prepareCheckSolutions(size_t n){ return NULL; }
	CK virtual size_t numLambdas()								{ return this->paramSet(t_lSet)->GetSize(); }
	CK inline auto rightPartFilter()							{ return m_pRightPartFilter; }
	CK inline void setForcibleLambdaPntr(S *p)					{ m_pForsibleLambda = p; }
	CK inline void setForcibleLambda(S nRow, S val, S nPart)	{ *(forcibleLambdaPntr(nRow) + nPart) = val; }

	S m_x0_3;
	S *m_pForsibleLambda;
	S m_firstUnforcedRow;
	CRightPartFilter<S> *m_pRightPartFilter;
	const bool m_bNoReplBlock;
};

FClass2(C_InSysEnumerator)::C_InSysEnumerator(const InSysPntr pInSys, uint enumFlags, int treadIdx, uint nCanonChecker) :
	m_bNoReplBlock(enumFlags & t_noReplicatedBlocks), Class2(CEnumerator)(pInSys, enumFlags | t_IS_enumerator, treadIdx, nCanonChecker),
	Class2(CInSysSolver)(pInSys->colNumb() >> 1, pInSys->GetT()), CVector<S>(pInSys->colNumb())
{
	const auto nCol = pInSys->colNumb();
	const auto tDesign = pInSys->GetT() > 2;
	const auto numParts = this->numParts();
	if (false && numParts > 1) { // Will not use for now  (search for this comment)
		auto pRowEquation = new CInSysRowEquation<S>[numParts];
		this->setRowEquation(pRowEquation);
		for (S i = 0; i < numParts; i++)
			pRowEquation[i].InitRowEquation(nCol, tDesign);
	} else
		this->setRowEquation(new CInSysRowEquation<S>(nCol, tDesign));

	setX0_3(nCol);
	m_pRightPartFilter = new CRightPartFilter<S>(nCol);
	const auto nRow = pInSys->rowNumb();
	setForcibleLambdaPntr(new S[nRow * numParts]);
	memset(forcibleLambdaPntr(), 0, nRow * numParts * sizeof(*forcibleLambdaPntr()));
	setFirstUnforcedRow();
	//setForcibleLambda(nRow - 1, this->getInSys()->lambda(), 0); // It looks like we dont need this
}

FClass2(C_InSysEnumerator)::~C_InSysEnumerator() {
	delete rightPartFilter();
	delete[] forcibleLambdaPntr();
}

FClass2(C_InSysEnumerator, bool)::checkSolutions(RowSolutionPntr pSolution, PERMUT_ELEMENT_TYPE i, bool doSorting) {
	if (!pSolution->numSolutions())
		return false;

	return pSolution->findFirstValidSolution(inSysRowEquation()->variableMaxValPntr());
}

FClass2(C_InSysEnumerator, RowSolutionPntr)::setFirstRowSolutions() {
	auto pSolutions = this->rowStuff(0);
	const auto* pR_set = this->getInSys()->GetNumSet(t_rSet);
	auto i = pR_set->GetSize();
	pSolutions->InitSolutions(1, i);
	while (i--)
		pSolutions->AddElement(pR_set->GetAt(i));

	return pSolutions;
}

FClass2(C_InSysEnumerator, S)::MakeSystem(S numPart)
{
	S nVar = 0;
	// Total number of equations (some of them corresponds to the forcibly constructed columns)
	S equationIdx = ELEMENT_MAX;
	// Number of equations corresponding only to the colums which ARE NOT forcibly constructed
	S eqIdx = 0;

	auto pRowEquation = inSysRowEquation();

	pRowEquation->resetMappings();
	const auto nRow = this->currentRowNumb();

	// When we are using the group of canonical matrix, we need to adjust previously constructed
	// generators of the automorphism group on columns (because some of them could be forcibly constructed
	this->setUseCanonGroup(USE_CANON_GROUP && this->groupOrder() > 1);

	int buffer[256], *pColGroupIdx = NULL;
	S lenPermut = 0;
	if (this->useCanonGroup()) {
		lenPermut = this->lenPerm();
		pColGroupIdx = lenPermut <= countof(buffer) ? buffer : new int[lenPermut];
	}

	const auto *pIS = this->getInSys();
	const auto nRowToBuild = pIS->rowNumb() - nRow;
	const auto k = pIS->GetK();
	rightPartFilter()->reset();

	int colGroupIdx = 0;
	setCurrentNumPart(numPart);
	const auto *pColOrbPrev = this->colOrbit(nRow - 1, numPart);
	CColOrbit<S> *pColOrbit, *pColOrbitNext = this->colOrbit(nRow, numPart);
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
			const auto len = pColOrbPrev->length() - currLen;
			// Column weight was changed. It means that x[i] > 0
			if (len) {
				// The length of the current orbit is NOT the same.
				// It means that current x[i] < xMax[i] (the orbit was split in two sub-orbits)
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

	if (pColOrbit)
		return ELEMENT_MAX;

	if (nRowToBuild > 1) {
		auto fLambda = forcibleLambda(nRow - 1, numPart);
		auto *pTmp = *(this->currUnforcedOrbPtr() + 1);
		while (pTmp) {
			fLambda += pTmp->length();
			pTmp = pTmp->next();
		}

		if (nRowToBuild != 2 || checkForcibleLambda(fLambda))
			setForcibleLambda(nRow, fLambda, numPart);
		else
			return ELEMENT_MAX;
	}

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

FClass2(C_InSysEnumerator, RowSolutionPntr)::FindSolution(S nVar, S nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex)
{
	const auto nRow = this->currentRowNumb();
	const auto nRowPrev = nRow - 1;
	const auto pPrevRowSolution = this->rowStuff(nRowPrev, nPart);
	auto pCurrRowSolution = this->rowStuff(nRow, nPart);
	pCurrRowSolution->setSolutionLength(nVar);
	pCurrRowSolution->resetSolution();
	const auto pFirstSol = pPrevRowSolution->firstSolution();
	const CSolutionPerm *pSolPerm = pPrevRowSolution->solutionPerm();
	const PERMUT_ELEMENT_TYPE *pPerm = pSolPerm ? pSolPerm->GetData() : NULL;
	const auto len = pPrevRowSolution->solutionLength();
	auto pRowEquation = inSysRowEquation();
	pRowEquation->setEquSystem(equSystem());

	// We will keep the variation intervals of "lambda" variables from the 2-variable equations
	// Because there are no more than nVar/2 of them, we will use
	this->initSolver(pCurrRowSolution, pRowEquation->variableMinValPntr());

	const auto pLambdaSet = this->paramSet(t_lSet);
	// For all possible right parts of the systems
	// (which are lexicographically not greater, than the vector, used for the current row)
	auto pResult = pCurrRowSolution->newSolution();
	S buffer[256];
	auto pBuffer = nVar <= countof(buffer)? buffer : new S[nVar];
	// Solution used for last constructed matrix's row:
	const auto pCurrSolution = pPrevRowSolution->solution(lastRightPartIndex);
	const auto forcibleLambdaValue = forcibleLambda(nRowPrev, nPart);
	const auto nLambdas = numLambdas();
	const auto *pMaxVal = inSysRowEquation()->variableMaxValPntr();
	bool readyToCheckSolution = false;
	for (PERMUT_ELEMENT_TYPE j = 0; j <= lastRightPartIndex; j++) {
		if (!pPrevRowSolution->isValidSolution(j))
			continue;

		const auto idx = pPerm ? *(pPerm + j) : j;
		const auto pRightSide = rightPartFilter()->getRightPart(pFirstSol + idx * len, pMaxVal, len, pBuffer);
		if (!pRightSide)
			continue;

		if (checkSolutionsForRight(nRow, nPart)) {
			if (!solutionsForRightSideNeeded(pRightSide, pCurrSolution, pLambdaSet))
				continue;
		}

		// Resolve trivial equations and define lambdaMin (part of any current lambda
		// which is already splitted among the lambda-variables)
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
			for (S i = 0; i < nLambdas; i++) {
				const int lambdaToSplit = this->getLambda(pLambdaSet, i, nPart) - lambdaMin;
				if (lambdaToSplit < 0)
					continue;

#if USE_EXRA_EQUATIONS
				equSystem()->resetVariables(nVar);
#endif
				auto pResultTmp = this->findAllSolutionsForLambda(pResult, lambdaToSplit);
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

	return pCurrRowSolution->getSolution();
}

FClass2(C_InSysEnumerator, void)::addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, int idx)
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

FClass2(C_InSysEnumerator, void)::setVariableLimit(S nVar, S len, S nRowToBuild, S k, S colWeight)
{
	// Minimal value for the nVar-th element which could be used as a first valid candidate for next row
	this->SetAt(nVar, static_cast<S>(((k - colWeight) * len + nRowToBuild - 1) / nRowToBuild));
	inSysRowEquation()->setVariableMaxLimit(nVar, len);

	// Maximal value of any Xi cannot be bigger than X0_3, when i-th group of columns
	// already have at least 2 inits in each column. Otherwise, these two rows with the row
	// which would use the solution, would be bigger than current 3 rows
	if (colWeight >= 2 && len > getX0_3())
		len = getX0_3();

	if (noReplicatedBlocks() && k - colWeight == 1)
		len = 1;

	inSysRowEquation()->setVariableMaxVal(nVar, len);
}

