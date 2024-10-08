﻿#pragma once
#include "Enumerator.h"
#include "RightPartFilter.h"
#include "designParam.h"

Class2Def(C_InSysEnumerator) : public Class2(CEnumerator), public Class2(CInSysSolver), public CVector<S>
{
public:
	CK C_InSysEnumerator(const InSysPntr pInSysm, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_InSysEnumerator();
	CK virtual T getX0_3() const								{ return m_x0_3; }
	CK virtual T firstUnforcedRow() const						{ return m_firstUnforcedRow; }
	CK virtual void setFirstUnforcedRow(T rowNum = 0)			{ m_firstUnforcedRow = rowNum; }
	CK virtual T *forcibleLambdaPntr(T nRow = 0) const			{ return m_pForsibleLambda + nRow * numParts(); }
	CK virtual bool noReplicatedBlocks() const					{ return m_bNoReplBlock; }
	CK virtual bool isPBIB_enumerator() const					{ return false;  }
	CK int define_MT_level(int v) const							{ return v / 2; }
	CK int define_MT_level(const designParam *pParam) const		{ return pParam->lambda()[0] == 1?
																		 pParam->v / pParam->k : define_MT_level(pParam->v); }
	CK void ConstructCanonicalMatrix(int k = 0) {
		CanonicityCheckerPntr pClassGroup = NULL;
		CanonizeMatrix(k, &pClassGroup);
		matrix()->printOut(outFile(), matrix()->rowNumb(), 0, this);
		if (pClassGroup) {
			pClassGroup->outputAutomorphismInfo(outFile());
			delete pClassGroup;
		}
	}
protected:
	CK virtual void setX0_3(T value)							{ m_x0_3 = value; }
	CK virtual bool checkSolutions(RowSolutionPntr ptr, PERMUT_ELEMENT_TYPE idx, bool doSorting = true);
	CK inline auto inSysRowEquation() const						{ return (CInSysRowEquation<S> *)this->rowEquation(); }
	CK virtual T getLambda(const VectorPntr pLambdaSet, T idx = 0, T numPart = 0) const { return pLambdaSet->GetAt(idx); }
	CK virtual bool checkSolutionsForRight(T row, T part) const	{ return false; }
	CK virtual bool solutionsForRightSideNeeded(const T *pRighPart, const T *pCurrSolution, const VectorPntr pLambdaSet) const
																{ return true; }
	CK virtual CEquSystem *equSystem()							{ return NULL;  }
	CK CColOrbit<S>** unforcedOrbits(T nRow, T iPart) const		{ return getUnforcedColOrbPntr(iPart) + shiftToUnforcedOrbit(nRow); }
	CK virtual CColOrbit<S> **getUnforcedColOrbPntr(T iPart) const {
		return forcibleLambda(this->currentRowNumb(), iPart) != ELEMENT_MAX ? this->unforcedColOrbPntr(iPart) : NULL;
	}
	CK virtual T forcibleLambda(T nRow, T nPart) const			{ return *(forcibleLambdaPntr(nRow) + nPart); }
	CK virtual void setColOrbitForCurrentRow(CColOrbit<S> *pColOrb){}
	CK virtual VectorPntr paramSet(t_numbSetType idx) const		{ return this->getInSys()->GetNumSet(idx); }
	CK virtual bool check_X0_3(T nPart) const					{ return false; }
	CK virtual bool checkRightParts(T nRow)						{ return false; }
	CK virtual bool useAsRightPart(CRowSolution<TDATA_TYPES> *pRowSol, PERMUT_ELEMENT_TYPE idx)		{ return true;  }
	CK void copyLimits(RowSolutionPntr pRowSolution, bool saveValues) const override;
#if !CONSTR_ON_GPU
	bool makeFileName(char* buffer, size_t len, const char* ext = NULL) const override;
#endif
#if USE_EXRA_EQUATIONS
	CK virtual void addColOrbitForVariable(S nVar, CColOrbit<S> *pColOrb)	{}
#else
#define addColOrbitForVariable(...)
#endif
	CK virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix);
	virtual void CanonizeByColumns(InSysPntr pMatrix, T *pColIdxStorage, CanonicityCheckerPntr pCanonChecker = NULL, bool permCol = false) const;
private:
	CK void addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, S nPart, S idx);
	CK virtual RowSolutionPntr setFirstRowSolutions();
	CK virtual T MakeSystem(T numPart);
	CK virtual RowSolutionPntr FindSolution(T nVar, T nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX);
	CK void setVariableLimit(T nVar, T len, T nRowToBuild, T currentK, T weightDeficit);
	CK virtual bool checkForcibleLambda(T fLambda, T nRows, T numPart) const		{ return true; }
	CK virtual void resetFirstUnforcedRow()						{ if (firstUnforcedRow() == this->currentRowNumb())
																	setFirstUnforcedRow(0);
																}
	virtual CVariableMapping<T> *prepareCheckSolutions(size_t n){ return NULL; }
	CK virtual T numLambdas() const								{ return static_cast<T>(paramSet(t_lSet)->GetSize()); }
	CK inline auto rightPartFilter()							{ return m_pRightPartFilter; }
	CK inline void setForcibleLambdaPntr(T *p)					{ m_pForsibleLambda = p; }
	CK virtual void setForcibleLambda(T nRow, T val, T nPart)	{ *(forcibleLambdaPntr(nRow) + nPart) = val; }
	CK virtual bool lambdaCond(T currentK) const				{ return false; }

	T m_x0_3;
	T *m_pForsibleLambda;
	T m_firstUnforcedRow;
	CRightPartFilter<T> *m_pRightPartFilter;
	const bool m_bNoReplBlock;
};

Class2Def(C_InSysCanonizator) : public Class2(C_InSysEnumerator) {
public:
	C_InSysCanonizator(Class2(C_InSys) * pInsSys, uint enumFlags) : Class2(C_InSysEnumerator)(pInsSys, enumFlags) {
		m_pGroupOrder = new CGroupOrder<T>;
		setGroupOrder(1);
	}
	~C_InSysCanonizator()					{ delete m_pGroupOrder; }
	inline void setGroupOrder(size_t val)   { m_pGroupOrder->setGroupOrder(val); }
	virtual void makeJobTitle(const designParam* pParam, char* buffer, int len, const char* comment = "") const {
		SNPRINTF(buffer, len, "Canonization of %" _FRMT "x%" _FRMT " matrix", this->rowNumb(),this->colNumb());
	}
	CK CGroupOrder<T>* extraGroupOrder() const override { return m_pGroupOrder; }
	CGroupOrder<T>* m_pGroupOrder;
};

FClass2(C_InSysEnumerator)::C_InSysEnumerator(const InSysPntr pInSys, uint enumFlags, int treadIdx, uint nCanonChecker) :
	m_bNoReplBlock(enumFlags & t_noReplicatedBlocks), Class2(CEnumerator)(pInSys, enumFlags | t_IS_enumerator, treadIdx, nCanonChecker),
	Class2(CInSysSolver)(pInSys->colNumb() >> 1, pInSys->GetT()), CVector<S>(pInSys->colNumb())
{
	const auto nCol = pInSys->colNumb();
	const auto tDesign = pInSys->GetT() > 2;
	const auto numParts = this->numParts();
	this->setRowEquation(new CInSysRowEquation<S>(nCol, tDesign));

	setX0_3(nCol);

	m_pRightPartFilter = new CRightPartFilter<T>(nCol);
	const auto nRow = pInSys->rowNumb();
	setForcibleLambdaPntr(new T[nRow * numParts]);
	memset(forcibleLambdaPntr(), 0, nRow * numParts * sizeof(*forcibleLambdaPntr()));

	setFirstUnforcedRow();
}

FClass2(C_InSysEnumerator)::~C_InSysEnumerator() {
	delete rightPartFilter();
	delete[] forcibleLambdaPntr();
}

FClass2(C_InSysEnumerator, bool)::checkSolutions(RowSolutionPntr pSolution, PERMUT_ELEMENT_TYPE i, bool doSorting) {
	if (!pSolution->numSolutions())
		return false;

	return pSolution->findFirstValidSolution(false);
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

FClass2(C_InSysEnumerator, T)::MakeSystem(T numPart)
{
	T nVar = 0;
	// Total number of equations (some of them corresponds to the forcibly constructed columns)
	T equationIdx = ELEMENT_MAX;
	// Number of equations corresponding only to the columns which ARE NOT forcibly constructed
	T eqIdx = 0;

	auto pRowEquation = inSysRowEquation();

	pRowEquation->resetMappings();
	const auto nRow = this->currentRowNumb();
	const auto pPermStorage = permStorage(numPart);

	// When we are using the group of canonical matrix, we need to adjust previously constructed
	// generators of the automorphism group on columns (because some of them could be forcibly constructed
	int buffer[256], *pColGroupIdx;
	T lenPermut = 0;
	if (this->useCanonGroup()) {
		lenPermut = pPermStorage->lenPerm();
		pColGroupIdx = lenPermut <= countof(buffer) ? buffer : new int[lenPermut];
	}
	else
		pColGroupIdx = NULL;

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
	const char* const pRowFirstOrbit = (char *)(this->colOrbitIni(nRow, numPart));
	CColOrbit<S> *pPrev = NULL;
	while ((pColOrbit = pColOrbitNext) != NULL) {
		pColOrbitNext = pColOrbit->next();
		equationIdx++;
		const auto columnWeight = pColOrbit->columnWeight();
		const auto nColll = ((char *)pColOrbit - pRowFirstOrbit)/ colOrbitLen();
		auto weightDeficit = k - columnWeight;
		if (!weightDeficit || weightDeficit == nRowToBuild) {
			colGroupIdx++;			// We need to skip this column's orbit
									// All remaining elements of all columns
									// of the current orbit should be equal to 0 or 1, respectively
			const auto currLen = pColOrbit->length();
			if (noReplicatedBlocks() && currLen > 1)
				break;

			addForciblyConstructedColOrbit(pColOrbit, pPrev, numPart, weightDeficit ? 1 : 0);
			if (currLen == pColOrbPrev->length()) {
				rightPartFilter()->addFilter(equationIdx, weightDeficit ? currLen : 0);
				pColOrbPrev = pColOrbPrev->next();
				continue;
			}

			// We are here only if orbit of blocks for previous row is split AND
			// the blocks of the first suborbits just got the kth unit.
			assert(!weightDeficit);
			weightDeficit = 1;
			pColOrbitNext = (pColOrbit = pColOrbitNext)->next();

			// When we are here, only first column(s) of pColOrbPrev orbits should be forcible and diffWeight = 0
			// Limit for right part should be established here 
			rightPartFilter()->addFilter(equationIdx, 0, nVar);
		}

		if (pColGroupIdx)
			pColGroupIdx[nVar] = colGroupIdx++;

		int mapSetIdx = 0;
		const auto currLen = pColOrbit->length();
		setVariableLimit(nVar, currLen, nRowToBuild, columnWeight, weightDeficit);
		const auto columnWeightPrev = pColOrbPrev->columnWeight();
		if (pColOrbit->columnWeight() != columnWeightPrev) {
			// Column weight was changed. It means that x[i] > 0
			addColOrbitForVariable(nVar, pColOrbit);
			const auto len = pColOrbPrev->length() - currLen;

			// The length of the current orbit is the same.
			// We could have two options here:
			//    (a) orbit was not splited in two sub-orbits (current x[i] = xMax[i]), in that lase len is equal to 0
			//    (b) orbit was splited in two sub-orbits, but second sub-orbit is forcibly constructed by 1's
			mapSetIdx = 1;
			if (len) {
				// The length of the current orbit is NOT the same.
				// Check if second sub-orbit should contain all units
				if (++weightDeficit == nRowToBuild) {
					if (noReplicatedBlocks() && len > 1)
						break;

					auto pTmp = pColOrbitNext;
					pColOrbitNext = pTmp->next();
					pTmp->InitOrb(len, pColOrbitNext);

					addForciblyConstructedColOrbit(pTmp, pColOrbit, numPart, 1);
					rightPartFilter()->addFilter(equationIdx, len, nVar);
					colGroupIdx++;    // This group of columns will be skipped
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
					setVariableLimit(++nVar, len, nRowToBuild, columnWeightPrev, weightDeficit);
					pColOrbitNext = (pColOrbit = pColOrbitNext)->next();
					if (pColGroupIdx)
						pColGroupIdx[nVar] = colGroupIdx++;

					mapSetIdx = 2;
				}
			}
		}

		pRowEquation->addVarMapping(mapSetIdx, nVar++, eqIdx++);
		pColOrbPrev = pColOrbPrev->next();
		pPrev = pColOrbit;
	}

	if (pColOrbit)
		return ELEMENT_MAX;

	if (nRowToBuild > 1) {
		auto fLambda = forcibleLambda(nRow - 1, numPart);
		// Take a pointer to the first orbit, forcibly constructed with 1's 
		const auto *pTmp = *(this->currUnforcedOrbPtr(numPart) + 1);
		const auto X0_3 = getX0_3();
		const auto check_X0_3_flg = X0_3 && check_X0_3(numPart);
		while (pTmp) {
			const auto len = pTmp->length();
			if (check_X0_3_flg && len > X0_3 && (pTmp->columnWeight() || nRowToBuild > 2))
				 return ELEMENT_MAX;

			fLambda += len;
			pTmp = pTmp->next();
		}

		if (!checkForcibleLambda(fLambda, nRowToBuild, numPart))
			return ELEMENT_MAX;

		setForcibleLambda(nRow, fLambda, numPart);
	}

	if (this->useCanonGroup()) {
		if (nVar < lenPermut && nRowToBuild > 1) {
			// We realy need to adjust the generators of the automorphism group on columns
			pPermStorage->adjustGenerators(pColGroupIdx, nVar);
		}

		if (pColGroupIdx != buffer)
			delete[] pColGroupIdx;
	}

	// For BIBDs when we construct system for last row (nRowToBuild=1), nVar > 0 
	// iff at least one column orbits was splitted by currently used solution.
	return nVar;
}

FClass2(C_InSysEnumerator, RowSolutionPntr)::FindSolution(T nVar, T nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex)
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

	const auto* pLambdaSet = getInSys()->objectType() != t_objectType::t_Kirkman_Triple?
		this->paramSet(t_lSet) : NULL;
	// For all possible right parts of the systems
	// (which are lexicographically not greater, than the vector, used for the current row)
	auto pResult = pCurrRowSolution->newSolution();
	T buffer[256];
	auto pBuffer = nVar <= countof(buffer)? buffer : new T[nVar];
	// Solution used for last constructed matrix's row:
	const auto pCurrSolution = pPrevRowSolution->solution(lastRightPartIndex);
	const auto forcibleLambdaValue = (int)forcibleLambda(nRowPrev, nPart);
	const auto nLambdas = numLambdas();
	const auto *pMaxVal = inSysRowEquation()->variableMaxValPntr();
	bool readyToCheckSolution = false;
	const bool flg = checkRightParts(nRow);
	for (PERMUT_ELEMENT_TYPE j = 0; j <= lastRightPartIndex; j++) {
		if (flg && !useAsRightPart(pPrevRowSolution, j))
			continue;

#if !HARD_REMOVE
		if (!pPrevRowSolution->validSolution(j))
			continue;
#endif
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

		lambdaMin += forcibleLambdaValue;

#if USE_EXRA_EQUATIONS
		resetExtra();
#endif
		CVariableMapping<T> *pVarMapping = NULL;
		if (!readyToCheckSolution) {
			readyToCheckSolution = true;
			pVarMapping = prepareCheckSolutions(nVar);
#if USE_EXRA_EQUATIONS
			// we need to revert the list of the forcibly constructed orbits
			//	addForciblyConstructedColOrbit(pColOrbit, pPrev, nPart, diffWeight ? 1 : 0);
			//			CColOrbitManager::addForciblyConstructedColOrbit(pColOrbit, nPart, idx);
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
			for (T i = 0; i < nLambdas; i++) {
				const auto λ = pLambdaSet ? static_cast<int>(this->getLambda(pLambdaSet, i, nPart)) : i;
				const int lambdaToSplit = λ - lambdaMin;
				if (lambdaToSplit < 0)
					continue;

#if USE_EXRA_EQUATIONS
				equSystem()->resetVariables(nVar);
#endif
				auto *pResultTmp = this->findAllSolutionsForLambda(pResult, lambdaToSplit, λ);
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

FClass2(C_InSysEnumerator, void)::addForciblyConstructedColOrbit(CColOrbit<S> *pColOrbit, CColOrbit<S> *pPrev, S nPart, S idx)
{
	if (pPrev)
		pPrev->setNext(pColOrbit->next());
	else
		this->setColOrbitCurr(pColOrbit->next(), nPart);
	// All remaining elements of all columns
	// of the current orbit should be equal to 0 or 1, respectively
	CColOrbitManager<S>::addForciblyConstructedColOrbit(pColOrbit, nPart, idx);
#if USE_EXRA_EQUATIONS
	equSystem()->addForcibleOrb(pColOrbit);
#endif

	// if it's possible, define first unforced row number
	if (!idx && !firstUnforcedRow())
		setFirstUnforcedRow(this->currentRowNumb());
}

FClass2(C_InSysEnumerator, void)::setVariableLimit(T nVar, T len, T nRowToBuild, T currentK, T weightDeficit)
{
	// Minimal value for the nVar-th element which could be used as a first valid candidate for next row
	this->SetAt(nVar, static_cast<T>((weightDeficit * len + nRowToBuild - 1) / nRowToBuild));
	inSysRowEquation()->setVariableMaxLimit(nVar, len);

	if (noReplicatedBlocks() && (lambdaCond(currentK) || weightDeficit == 1))
		len = 1;
	else {
		// Maximal value of any Xi cannot be bigger than X0_3, when i-th group of columns
		// already have at least 2 inits in each column. Otherwise, these two rows with the row
		// which would use the solution, would be bigger than current 3 rows

		if (!currentNumPart() && currentK > 2 && len > getX0_3())
			len = getX0_3();
	}

	inSysRowEquation()->setVariableMaxVal(nVar, len);
}

#if !CONSTR_ON_GPU
FClass2(C_InSysEnumerator, bool)::makeFileName(char* buffer, size_t lenBuffer, const char* ext) const
{
	const auto inSys = this->getInSys();
	std::string inputFile(designParams()->logFile());
	auto pos = inputFile.find_last_of("/");
	if (pos != std::string::npos)
		inputFile = inputFile.substr(0, pos + 1);
	else
		inputFile.clear();

	const auto len = SNPRINTF(buffer, lenBuffer, "%s", inputFile.c_str());
	SNPRINTF(buffer + len, lenBuffer - len, ME_FRMT"_" ME_FRMT, rowNumb(), colNumb());
	return true;
}
#endif