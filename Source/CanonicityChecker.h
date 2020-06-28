//
//  CanonicityChecker.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__CanonicityChecker__
#define __BIBD_Mac__CanonicityChecker__
#pragma once

#include "DataTypes.h"
#include "PermutStorage.h"
#include "ColOrbits.h"

#define CPermut				CSimpleArray<S>
#define CColNumbStorage		CContainer<S>

typedef enum {
	t_saveNothing			= 0,
	t_saveRowPermutations	= (1 << 0),
	t_saveColPermutations	= (1 << 1),
	t_saveRowToChange		= (1 << 2)
} t_canonOutInfo;

Class2Def(CMatrix);
Class2Def(CMatrixCol);
Class2Def(CRowSolution);

Class2Def(CCanonicityChecker) {
public:
	CC CCanonicityChecker(S nRow, S nCol, int rank = 2, uint enumFlags = t_enumDefault, S numParts = 1);
	CC ~CCanonicityChecker();
	void InitCanonicityChecker(S nRow, S nCol, int rank, S *pMem);
	CC bool TestCanonicity(S nRowMax, const MatrixColPntr pEnum, uint outInfo, S *pPartNumb = NULL, S *pRowOut = NULL, RowSolutionPntr pRowSolution = NULL);
	void outputAutomorphismInfo(FILE *file, const MatrixDataPntr pMatrix = NULL) const;
	CC auto enumFlags() const						{ return m_enumFlags; }
	CC inline auto groupOrder() const				{ return m_nGroupOrder; }
	CC inline void setGroupOrder(uint val)			{ m_nGroupOrder = val; }
	CC inline auto permRow() const					{ return m_pPermutRow->elementPntr(); }
	CC inline auto permCol() const					{ return m_pPermutCol->elementPntr(); }
	CC inline auto permStorage(S nPart) const		{ return permStorage() + nPart; }
	CC inline auto permStorage() const				{ return m_pPermutStorage[0]; }
	CC inline auto permColStorage() const			{ return m_pPermutStorage[1]; }
	CC inline auto permRowStorage() const			{ return m_pPermutStorage[2]; }
	CC inline auto improvedSolution() const			{ return m_pImprovedSol; }
	CC bool groupIsTransitive() const;
	bool printMatrix(const designParam *pParam) const;
	CC S stabiliserLengthExt() const				{ return m_nStabExtern; }
protected:
	CC inline int rank() const						{ return m_rank; }
	CC virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix)		{}
	virtual void CanonizeByColumns(MatrixDataPntr pMatrix, S *pColIdxStorage = NULL, CanonicityCheckerPntr pCanonChecker = NULL) const	{}
	CC inline auto getRowOrbits(int idx) const		{ return m_pObits[0][idx]; }
	CC inline auto getColOrbits(int idx) const		{ return m_pObits[1][idx]; }
	inline bool checkProperty(uint flag) const		{ return enumFlags() & flag; }
	CC inline auto numRow() const					{ return m_nNumRow; }
	CC void setStabiliserLengthExt(S len)			{ m_nStabExtern = len; }
	CC inline S numParts() const					{ return m_numParts; }
	CC virtual S lenStabilizer() const				{ return 0; }
private:
	CC void init(S nRow, bool savePerm);
	CC S next_permutation(S idx = ELEMENT_MAX, S lenStab = 0);
	CC void addAutomorphism(bool rowPermut = true);
	CC int checkColOrbit(S orbLen, S nColCurr, const T *pRow, const T *pRowPerm) const;
	CC inline void setStabilizerLength(S len)		{ m_nStabLength = len; }
	CC inline auto stabilizerLength() const			{ return m_nStabLength; }
	CC inline void setStabilizerLengthAut(S l)		{ m_nStabLengthAut = l; }
	CC inline auto stabilizerLengthAut() const		{ return m_nStabLengthAut; }
	CC inline void setNumRow(S nRow)				{ m_nNumRow = nRow; }
	CC inline auto numCol() const					{ return static_cast<S>(m_pPermutCol->numElement()); }
#define orbits()	m_pObits[0][0] 
	CC void updateGroupOrder();
	CC inline auto colNumbStorage() const			{ return m_nColNumbStorage; }
	CC inline auto counter() const					{ return m_pCounter; }
	CC inline void setColIndex(S *p)				{ m_pColIndex = p; }
	CC inline auto colIndex() const					{ return m_pColIndex; }
	CC void revert(S i);
	CC S constructColIndex(const ColOrbPntr pColOrbit, const ColOrbPntr pColOrbitIni, size_t colOrbLen, S shift = 0) const;
	CC inline void setPermStorage(PermutStoragePntr p, int idx = 0)	{ m_pPermutStorage[idx] = p; }
	CC S rowToChange(S nRow) const;
	void reconstructSolution(const ColOrbPntr pColOrbitStart, const ColOrbPntr pColOrbit,
		size_t colOrbLen, const ColOrbPntr pColOrbitIni, const T *pRowPerm, const S *pRowSolution, size_t solutionSize);
	CC void UpdateOrbits(const S *permut, S lenPerm, S *pOrbits, bool rowPermut, bool updateGroupOrder = false);
#if USE_STRONG_CANONICITY
	inline void setSolutionStorage(CSolutionStorage *p) { m_pSolutionStorage = p; }
	inline CSolutionStorage *solutionStorage() const { return m_pSolutionStorage; }

	CSolutionStorage *m_pSolutionStorage;
#else
#define setSolutionStorage(x)
#define solutionStorage()		(CSolutionStorage *)NULL
#endif

	S m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	S m_nStabLength;
	S m_nStabLengthAut;
    int m_rank;
    CPermut *m_pPermutRow;
    CPermut *m_pPermutCol;
    CColNumbStorage **m_nColNumbStorage;
	PermutStoragePntr m_pPermutStorage[3];
	S *m_pObits[2][2];
	CCounter<int> *m_pCounter;
	S *m_pColIndex;
	S *m_pImprovedSol;
	const uint m_enumFlags;
	uint m_nGroupOrder;
	S m_nNumRow;
	const S m_numParts;
};

CanonicityChecker()::CCanonicityChecker(S nRow, S nCol, int rank, uint enumFlags, S numParts) : m_rank(rank), m_enumFlags(enumFlags), m_numParts(numParts)
{
	m_pPermutRow = new CPermut(nRow);
	m_pPermutCol = new CPermut(nCol);
	setColIndex(new S[nCol << 1]);
	m_pCounter = new CCounter<int>(rank);
	m_nColNumbStorage = new CColNumbStorage *[rank];
	for (int i = rank; i--;)
		m_nColNumbStorage[i] = new CColNumbStorage(nCol);

	m_pImprovedSol = new S[(rank - 1) * nCol];
	setSolutionStorage(new CSolutionStorage());
	const auto mult = enumFlags & t_outStabilizerOrbit ? 2 : 1;

	const auto outColumnOrbits = enumFlags & t_outColumnOrbits;
	const auto keepRowPermute = enumFlags & t_alwaysKeepRowPermute;
	const auto nPermutStorages = numParts + (outColumnOrbits ? 1 : 0) + (keepRowPermute ? 1 : 0);
	auto pPermStorage = new CPermutStorage<T, S>[nPermutStorages]();
	setPermStorage(pPermStorage);

	memset(m_pObits, 0, sizeof(m_pObits));
	const auto len = nRow + (outColumnOrbits ? nCol : 0);
	auto pntr = m_pObits[0][0] = new S[mult * len];	// memory to keep orbits of different types
	int permStorageIdx = numParts;
	if (outColumnOrbits) {
		setPermStorage(pPermStorage + permStorageIdx++, 1);
		m_pObits[1][0] = pntr + mult * nRow;
		if (mult > 1)
			m_pObits[1][1] = m_pObits[1][0] + nCol;
	}
	else
		setPermStorage(NULL, 1);

	setPermStorage(keepRowPermute ? pPermStorage + permStorageIdx : NULL, 2);

	if (mult > 1)
		m_pObits[0][1] = pntr + nRow;
}

CanonicityChecker()::~CCanonicityChecker()
{
	delete m_pPermutRow;
	delete m_pPermutCol;
	delete counter();
	delete [] permStorage();
	for (int i = rank(); i--;)
		delete m_nColNumbStorage[i];

	delete[] colNumbStorage();
	delete[] colIndex();
	delete[] improvedSolution();
	delete[] m_pObits[0][0];
#if USE_STRONG_CANONICITY
	delete solutionStorage();
#endif
}

CanonicityChecker(bool)::TestCanonicity(S nRowMax, const MatrixColPntr pEnum, uint outInfo, S *pPartNumb, S *pRowOut, RowSolutionPntr pRowSolution)
{
	// Construct trivial permutations for rows and columns
	const bool rowPermut = outInfo & t_saveRowPermutations;
	init(nRowMax--, rowPermut);

	const auto *pMatr = pEnum->matrix();
	const auto colOrbLen = pEnum->colOrbitLen();
	const auto colNumb = pEnum->colNumb();
	auto **colOrbit = pEnum->colOrbits();
	auto **colOrbitIni = pEnum->colOrbitsIni();
#ifndef USE_CUDA
	const S *pCurrSolution;
	size_t solutionSize;
	if (pRowSolution) {
		// Because the position of current solution (pRowSolution->solutionIndex()) 
		// could be changed, let's take pointer here
		pCurrSolution = pRowSolution->currSolution();
		solutionSize = pRowSolution->solutionLength();
#if USE_STRONG_CANONICITY
		solutionStorage()->clear();
		pRowSolution->setLenOrbitOfSolution(0);
#endif
	}

	size_t startIndex = 0;
#endif

	const auto pPartInfo = pMatr->partsInfo();
	PREPARE_PERM_OUT(permColStorage());
	const auto lenStab = stabiliserLengthExt();

	bool retVal = true;
	S nRow = ELEMENT_MAX;
	const auto* permColumn = permCol();
	while (true) {
	next_permut:
		nRow = next_permutation(nRow, lenStab);
		if (nRow == ELEMENT_MAX || nRow < lenStabilizer())
			break;

		OUT_PERM(permRowStorage(), permRow(), nRowMax + 1);
		OUTPUT_PERMUTATION(permColStorage(), pEnum->outFile(), orbits(), numRow());

		// Loop for all remaining matrix's rows
		for (; nRow <= nRowMax; nRow++) {
			const auto* pRow = pMatr->GetRow(nRow);
			const auto* pRowPerm = pMatr->GetRow(*(permRow() + nRow));
			for (S nPart = 0; nPart < numParts(); nPart++) {
				const auto* pColOrbitIni = pEnum->colOrbitIni(nRow, nPart);
				const auto* pColOrbit = pEnum->colOrbit(nRow, nPart);
				const auto shift = nPart? pPartInfo->getShift(nPart) : 0;
				while (pColOrbit) {
					// Define the number of column to start with
					const auto nColCurr = shift + static_cast<S>(((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen);
					const auto orbLen = pColOrbit->length();
					int diff;
					if (orbLen > 1)
						diff = checkColOrbit(orbLen, nColCurr, pRow, pRowPerm);
					else
						diff = (int) * (pRow + nColCurr) - *(pRowPerm + *(permColumn + nColCurr));

					if (diff > 0)
						goto next_permut;

					if (!diff) {
						pColOrbit = pColOrbit->next();
						continue;
					}

					if (!nPart && pRowOut) {
						if (outInfo & t_saveRowToChange)
							*pRowOut = rowToChange(nRow);
#ifndef USE_CUDA
						else {
							if (pRowSolution && !pRowSolution->isLastSolution()) {
								// NOTE: No need to remove last solution, we will leave this level anyway
								if (nRow < nRowMax) {
									// We can remove all solutions which are in the same group with the solution just tested.
									// We don't need them even as the right parts of our equations, because the usage of everyone 
									// will make the matrix non canonical
									pRowSolution->removeNoncanonicalSolutions(startIndex);
								}
								else {
									reconstructSolution(pEnum->colOrbit(nRow, nPart), pColOrbit, colOrbLen, pColOrbitIni, pRowPerm, pCurrSolution, solutionSize);
#if USE_STRONG_CANONICITY
									if (solutionStorage()) {
										for (auto* pSolution : *solutionStorage()) {
											if (!MEMCMP(pSolution, improvedSolution(), solutionSize * sizeof(*pSolution)))
												goto next_permut; // We already saw this solution
										}
									}

									startIndex = pRowSolution->moveNoncanonicalSolutions(improvedSolution(), startIndex, solutionStorage(), pRowOut);
									if (startIndex != SIZE_MAX) {
										retVal = false;
										// we will try to find better alternative
										goto next_permut;
									}
#else
									pRowSolution->moveNoncanonicalSolutions(improvedSolution(), 0, solutionStorage());
#endif
								}
							}
						}
#endif				
					}

					if (pPartNumb)
						*pPartNumb = nPart;

					return false;
				}
			}
		}

		// Automorphism found:
#if (PRINT_CURRENT_MATRIX && PRINT_PERMUTATION)
		outString("-----\n", pEnum->outFile());
#endif
		if (!rowPermut) {
			// We are here to define the canonicity of partially constructed matrix AND we just found the non-trivial automorphism.
			// Let's construct the permutation of the column's orbits which corresponds to just found automorphism
			for (S nPart = 0; nPart < numParts(); nPart++) {
				auto pPermStorage = permStorage(nPart);
				const auto* pColOrbitIni = pEnum->colOrbitIni(nRow, nPart);
				const auto* pColOrbit = pEnum->colOrbit(nRow, nPart);
				const auto shift = nPart ? pPartInfo->getShift(nPart) : 0;
				// Saving permutations, acting on the orbits of columns
				S lenPermut;
				if (pPermStorage->isEmpty()) {		// Index for columns was not yet constructed
					lenPermut = constructColIndex(pColOrbit, pColOrbitIni, colOrbLen, shift);
					pPermStorage->allocateMemoryForPermut(lenPermut);	// Memory for identity permutation just in case we will need it
						                                                // NOTE: as of 04/05/2020, we do not use it
				} else
					lenPermut = pPermStorage->lenPerm();

				auto pVarPerm = pPermStorage->allocateMemoryForPermut(lenPermut);
				while (pColOrbit) {
					const auto nColCurr = shift + ((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen;
					*pVarPerm++ = *(colIndex() + *(permCol() + nColCurr));
					pColOrbit = pColOrbit->next();
				}
			}
		}

//#ifndef USE_CUDA
		// We need the permutations on columns AND
		//  (a) matrix is completely constructed OR
		//  (b) we will need to analyse the group on partially constructed matrix
		// (As of Oct.11, 2018 this is used only for IGraphs (semi-symmetric graphs)
		if (permColStorage() && (rowPermut || permRowStorage()))
			ConstructColumnPermutation(pEnum->matrix());
//#endif

		addAutomorphism(rowPermut);
		nRow = ELEMENT_MAX - 1;
	}

	if (rowPermut)
		updateGroupOrder();

	return retVal;
}

#endif /* defined(__BIBD_Mac__CanonicityChecker__) */
