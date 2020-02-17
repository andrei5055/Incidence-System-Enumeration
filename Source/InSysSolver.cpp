//
//  InSysSolver.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/2/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "InSysSolver.h"
#include "RowSolution.h"
#if USE_EXRA_EQUATIONS
#include "EquSystem.h"
#endif

template class CInSysSolver<MATRIX_ELEMENT_TYPE, SIZE_TYPE>;


#if PRINT_RESULTS
	CEquation *pMyEquA[PR_EQU_NUMB];
	int ppp, ggg;
	size_t nSolutions;
	size_t nSolutionsPrev;
	size_t printResRowNumb;
	size_t printResNumVar;

	void setPrintResultRowNumber(size_t nRow)	{ printResRowNumb = nRow; }
	void setPrintResultNumVar(size_t numVar)	{ printResNumVar = numVar; }

	void assignPrintResultInfo(int &myID, CEquation *pntr)
	{
		static int id;
		if (PRINT_RES_ROW_NUMB == printResRowNumb && ggg == PRINT_RES_SOL_NUMB)
			pMyEquA[myID = id++] = pntr;
		else
			myID = -1;

	}

	void assignOrbID(int &myOrbID)
	{
		static int orbID;
		myOrbID = orbID++;
	}

	void printResults(FILE *file, VECTOR_ELEMENT_TYPE *pResult, size_t len, int varIdx)
	{
		const size_t lenBuffer = len * 3 + 2;
		char *buffer = new char[lenBuffer];

		char *pntr = buffer;
		for (size_t i = 0; i < len; i++)
			pntr += sprintf_s(pntr, lenBuffer - (pntr - buffer), "%3d", i);

		fprintf(file, "%s\n", buffer);

		pntr = buffer;
		for (size_t i = 0; i < len; i++)
			pntr += sprintf_s(pntr, lenBuffer - (pntr - buffer), "%3d", pResult[i]);

		fprintf(file, "%s\n", buffer);

		if (0 <= varIdx && varIdx < (int)len) {
			memset(buffer, ' ', strlen(buffer));
			buffer[3 * varIdx + 2] = '^';
			fprintf(file, "%s\n", buffer);
		}

		delete[] buffer;
	}

	FILE *openOutFile()
	{
		if (!pMyEquA[0])
			return NULL;

		FILE *file;
		if (fopen_s(&file, "C:/ttt/TTT.txt", ppp++ ? "a" : "w"))
			return NULL;

		return file;
	}

	void printResults(VECTOR_ELEMENT_TYPE *pResult, int lambda, int idx, int varIdx)
	{
		FILE *file = openOutFile();
		if (!file)
			return;

		const char *pTag = idx >= 0 ? "+++++++++++++++++++" : (pResult ? "-------------------" : "===================");
		fprintf(file, "Total: %2d: %s %s = %d  ppp = %2d\n", nSolutions, pTag, pResult ? "lambdaToSplit" : "lambdaMin", lambda, ppp);
		for (int i = 0; i < PRINT_RES_EQU_NUMB; i++) {
			CEquation *pMyEqu = pMyEquA[i];
			if (!pMyEqu)
				continue;

			fprintf(file, "%d:  m_nVar = %d  right part = %d  solved = %3s\n", i, pMyEqu->numbVar(),
				pMyEqu->rightPart(), pMyEqu->solved() ? "Yes" : "No");
		}

		if (idx >= 0)
			nSolutions++;

		if (pResult)
			printResults(file, pResult, printResNumVar, varIdx);

		fclose(file);
	}

	void printMapInfo(CVariableMapping *pVarValue, VECTOR_ELEMENT_TYPE *pResult)
	{
		if (!ppp)
			return;

		FILE *file = openOutFile();
		if (!file)
			return;

		const auto pMapping = pVarValue->getMappingPntr();
		size_t nElem = pVarValue->nElement();
		int varIdx = *(pMapping + (nElem << 1) - 2);
		fprintf(file, "nElem = %2d  varIdx = %2d  varVal = %2d\n", nElem, varIdx, pResult[varIdx]);
		fclose(file);
	}

	void printAddEquationVariable(const CEquation *pEqu, VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE varVal, bool addFlg)
	{
		if (pEqu->myID < 0 || EQUATION_ID_TO_PRINT >= 0 && pEqu->myID != EQUATION_ID_TO_PRINT)
			return;

		FILE *file = openOutFile();
		if (!file)
			return;

		fprintf(file, "ppp = %2d:  [%d] %s  m_nVar = %d  right part = %d  solved = %3s  varIdx = %2d  valVar = %d\n", ppp, pEqu->myID, addFlg ? "+++" : "---", pEqu->numbVar(), pEqu->rightPart(), pEqu->solved() ? "Yes" : "No", varIdx, varVal);
		fclose(file);
	}

	void printDefinedVariable(VECTOR_ELEMENT_TYPE varID, VECTOR_ELEMENT_TYPE varValue)
	{
		FILE *file = openOutFile();
		if (!file)
			return;

		fprintf(file, "varID = %2d = %d  (%3d)\n", varID, varValue , ppp);
		fclose(file);
	}

	void printFinalResultForRightPart(CVariableMapping *pVarValue)
	{
		FILE *file = openOutFile();
		if (!file)
			return;

		fprintf(file, "===>>> Constructed %d solutions for current right part, Map Pos: %s \n", nSolutions - nSolutionsPrev, pVarValue->nElement() == 0 ? "OK" : "WRONG");
		nSolutionsPrev = nSolutions;
		fclose(file);
	}

#endif

#if USE_EXRA_EQUATIONS
template<class T>
int CInSysSolver<T>::unassignVariable(int nVar, VECTOR_ELEMENT_TYPE *pResult, size_t adj, bool adjustSecond, int varIdx, bool excludeVar, bool undoAssignNeeded)
{
	if (excludeVar) {
		// Unassign nVar variables which are in getVarValues() stack
		// varIdx -  is the index of leading variable. It will be unassigned last.
		printMapInfo(getVarValues(), pResult);
		equSystem()->excludeVariables(getVarValues(), 1, false, nVar, NULL, undoAssignNeeded);
	}

	if (varIdx >= 0) {
		// We need to undo assigning for all, but first variable of the group, 
		// if we plan to adjust it here
		nVar--;
	}

	int adjLambda = 0;
	if (nVar) {
		if (undoAssignNeeded)
			adjLambda = undoAssignVariableValues(pResult, nVar - adj, adj);

		getVarValues()->removeLastMapping(nVar);  // Remove just added variables

		if (!excludeVar && varIdx >= 0) //???? varIdx == 11 /*
			equSystem()->excludeVariables(getVarValues(), 1, false, 1);
	}

	if (varIdx < 0)
		return adjLambda;

	const size_t idx = (getVarValues()->nElement() << 1) - 1;
	getVarValues()->adjustElement(idx, 1);
	if (adjustSecond)  // Second variable of the pair is in one of extra equations 
		getVarValues()->adjustElement(idx + 2, -1); // we need to adjust it's value

	--*(pResult + varIdx);
	++*(pResult + varIdx + 1);
	return adjLambda + 1;
}

template<class T>
int CInSysSolver<T>::removeVariableFromExtraEquations(int varIdx, VECTOR_ELEMENT_TYPE *pResult, VECTOR_ELEMENT_TYPE variation, int &lambdaToSplit) CONST
{
	int diffLambda = 0;
	const size_t nTotal = getVarValues()->nElement();
	do {
		int nAddedVar = equSystem()->excludeVariables(pResult, varIdx, variation, varMinValPntr(), getVarValues());
		bool exclFlag = false;
		const bool adjustSecond = !tDesignEnum() && getVarValues()->nElement() - nAddedVar == 2;  // Define if second variable of the pair is in one of extra equations
		// NOTE: It is easy to prove that it's never TRUE for t-design
		while (nAddedVar > 0) {
			exclFlag = true;
			diffLambda = assignVariableValues(pResult, nAddedVar, getLastMapping() - 2, lambdaToSplit);
			if (diffLambda < 0)
				break;

			lambdaToSplit -= diffLambda;
			nAddedVar = equSystem()->excludeVariables(getVarValues(), 1, true, nAddedVar);
		}

		if (!nAddedVar)
			return 1;

		if (*(pResult + varIdx) == *(varMinValPntr() + varIdx))
			varIdx = -1;

		variation = 0;

#if VAR_1
		const bool exclVar = nAddedVar > 0;
#else
		const bool exclVar = exclFlag; 
#endif
		lambdaToSplit += unassignVariable(getVarValues()->nElement() - nTotal/* + varIdx < 0? 1 : 0 */, pResult, exclVar ? 1 : 0, adjustSecond, varIdx, exclVar, diffLambda != -2);
	} while (varIdx >= 0);

	return -1;
}

template<class T>
int CInSysSolver<T>::assignVariableValues(VECTOR_ELEMENT_TYPE *pResult, size_t nVar, const VECTOR_ELEMENT_TYPE *pVariation, int lambdaToSplit) const
{
	int diffLambda = 0;
	auto *pVarDefined = equSystem()->varDefinedPtr();
	const auto pMapping = getVarValues()->getMappingPntr();
	int j = getVarValues()->nElement() << 1;
	size_t n = nVar;
	while (n--) {
		j -= 2;
		const int variation = *(pMapping + j + 1);
		if (!variation)
			continue;   // Both variables already have correct values 

		// Check validity of current variation 
		if (variation < 0)
			break;		//  New value of first variable would be less than its minimum

		if (variation > lambdaToSplit - diffLambda)
			break;		// Negative amount of units to be splitted for remaining variables

		const auto idxVar = *(pMapping + j);
		if (variation > *(pResult + idxVar + 1))
			break;      // New value of second variable would be < 0

		// Searching for the maximum value of variation for current variable.
		// Notes: 
		//     a) variations  a stored in the reverse order of the variable indices
		//     b) binary search could be implemented here.   
		const VECTOR_ELEMENT_TYPE *pVariationCurr = pVariation;
		while (*pVariationCurr < idxVar)
			pVariationCurr -= 2;

		if (variation > *(pVariationCurr + 1))
			break;		// New value of first variable would be bigger than its maximum

		// We need to modify the variable's value
		*(pResult + idxVar) += variation;

		// NOTE: This is correct for t-designs: current variable is always the first variable in some pair
		//       Let's modify the second variable of the pair
		*(pResult + idxVar + 1) -= variation;
		// This should be done ONLY in removeVariableFroExtraEquations
		// *(pVarDefined + idxVar) = *(pVarDefined + idxVar+ 1) = 1;
		diffLambda += variation;
	}

	if (n == -1)
		return diffLambda;

	if (n + 1 == nVar)  // No need to call undoAssignVariableValues, the method will do nothing
		return -1;		// 

	// Assignment was failed. We need to undo all assignments just made.
	undoAssignVariableValues(pResult, nVar, 0, n, j, false);

	return -2;
}

template<class T>
size_t CInSysSolver<T>::undoAssignVariableValues(VECTOR_ELEMENT_TYPE *pResult, size_t nVar, size_t adj, size_t n, int j, bool backStep) const
{
	auto *pVarDefined = equSystem()->varDefinedPtr();
	const auto *pMapping = getVarValues()->getMappingPntr();
	if (n == -1)
		j = (getVarValues()->nElement() - adj) << 1;

	int step = backStep ? -2 : 2;
	size_t diffLambda = 0;
	while (++n < nVar)
		diffLambda += undoAssignVariableValue(*(pMapping + (j += step)), pVarDefined, pResult);

	return diffLambda;
}

template<class T>
size_t CInSysSolver<T>::undoAssignVariableValue(VECTOR_ELEMENT_TYPE idxVar, char *pVarDefined, VECTOR_ELEMENT_TYPE *pResult) const
{
	*(pVarDefined + idxVar) = *(pVarDefined + idxVar + 1) = 0;
	const size_t diff = *(pResult + idxVar) - *(varMinValPntr() + idxVar);
	if (diff) {
		*(pResult + idxVar) -= diff;
		*(pResult + idxVar + 1) += diff;
	}

	return diff;
}

template<class T>
void CInSysSolver<T>::addVarIdxToStack(size_t nAddedVarIdx, size_t varShift)
{ 
	// nAddedVarIdx - number of variables added to the stack at that stage
	// varShift     - shift to the mapping of the leading variable of that group
	addedVarStack()->addMapping(nAddedVarIdx, varShift); 
}

template<class T>
size_t CInSysSolver<T>::getVarIdxFromStack()
{ 
	addedVarStack()->removeLastMapping();
	return *addedVarStack()->getLastMapping();
}

#endif

