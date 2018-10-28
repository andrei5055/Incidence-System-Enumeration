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

#define CPermut				CSimpleArray<T>
#define CColNumbStorage		CContainer<T>

typedef enum {
	t_saveNothing			= 0,
	t_saveRowPermutations	= (1 << 0),
	t_saveColPermutations	= (1 << 1),
	t_saveRowToChange		= (1 << 2)
} t_canonOutInfo;

template<class T> class CMatrix;
template<class T> class CMatrixData;
template<class T> class CMatrixCol;
template<class T> class CRowSolution;

template<class T>
class CCanonicityChecker
{
public:
    CC CCanonicityChecker(T nRow, T nCol, int rank = 2, uint enumFlags = t_enumDefault);
    CC ~CCanonicityChecker();
	void InitCanonicityChecker(T nRow, T nCol, int rank, T *pMem);
	CC bool TestCanonicity(T nRowMax, const CMatrixCol<T> *pEnum, int outInfo = 0, T *pRowOut = NULL, CRowSolution<T> *pRowSolution = NULL);
    void outputAutomorphismInfo(FILE *file, const CMatrixData<T> *pMatrix = NULL) const;
	CC uint enumFlags() const						{ return m_enumFlags; }
	CC inline uint groupOrder() const				{ return m_nGroupOrder; }
	CC inline T numColOrb() const					{ return m_nNumColOrb; }
	CC inline void setGroupOrder(uint val)			{ m_nGroupOrder = val; }
	CK inline T lenPerm() const						{ return permStorage()->lenPerm(); }
	CK inline void adjustGenerators(int *idx, T len){ permStorage()->adjustGenerators(idx, len); }
	CK inline size_t constructGroup()				{ return permStorage()->constructGroup(); }
	CK inline size_t findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon)
													{ return permStorage()->findSolutionIndex(pFirst, idx, pMem, pCanonIdx, nCanon); }
    CC inline T * permRow() const					{ return m_pPermutRow->elementPntr(); }
	CC inline T * permCol() const					{ return m_pPermutCol->elementPntr(); }
	CC inline CPermutStorage<T> *permStorage() const		{ return m_pPermutStorage[0]; }
	CC inline CPermutStorage<T> *permColStorage() const		{ return m_pPermutStorage[1]; }
	CC inline CPermutStorage<T> *permRowStorage() const		{ return m_pPermutStorage[2]; }
	CC inline VECTOR_ELEMENT_TYPE *improvedSolution() const	{ return m_pImprovedSol; }
	CC bool groupIsTransitive() const;
	bool printMatrix(const designParam *pParam) const;

protected:
	CC inline int rank() const						{ return m_rank; }
	virtual void ConstructColumnPermutation(const CMatrixData<T> *pMatrix)		{}
	virtual void CanonizeByColumns(CMatrixData<T> *pMatrix, T *pColIdxStorage = NULL, CCanonicityChecker *pCanonChecker = NULL) const	{}
	CC inline T *getRowOrbits(int idx) const		{ return m_pObits[0][idx]; }
	CC inline T *getColOrbits(int idx) const		{ return m_pObits[1][idx]; }
	inline bool checkProperty(uint flag) const		{ return enumFlags() & flag; }
	CC inline T numRow() const						{ return m_nNumRow; }
private:
    CC void init(T nRow, bool savePerm);
	CC T next_permutation(T idx = MATRIX_ELEMENT_MAX);
	CC void addAutomorphism(bool rowPermut = true);
    CC int checkColOrbit(size_t orbLen, size_t nColCurr, const T *pRow, const T *pRowPerm) const;
    CC inline void setStabilizerLength(T len)		{ m_nStabLength = len; }
    CC inline T stabilizerLength() const			{ return m_nStabLength; }
    CC inline void setStabilizerLengthAut(T l)		{ m_nStabLengthAut = l; }
    CC inline T stabilizerLengthAut() const			{ return m_nStabLengthAut; }
	CC inline void setNumRow(T nRow)				{ m_nNumRow = nRow; }
    CC inline T numCol() const						{ return static_cast<T>(m_pPermutCol->numElement()); }
#define orbits()	m_pObits[0][0] 
	CC void updateGroupOrder();
    CC inline CColNumbStorage **colNumbStorage() const	{ return m_nColNumbStorage; }
    CC inline CCounter<int> *counter() const			{ return m_pCounter; }
	CC inline void setColIndex(T *p)				{ m_pColIndex = p; }
	CC inline T *colIndex() const					{ return m_pColIndex; }
	CC inline void setNumColOrb(T v)				{ m_nNumColOrb = v; }
    CC void revert(T i);
	CC T getLenPermutCol(T **permCol) const;
	CC T * constructColIndex(const CColOrbit<T> *pColOrbit, const CColOrbit<T> *pColOrbitIni, size_t colOrbLen);
	CC inline void setPermStorage(CPermutStorage<T> *p, int idx = 0)	{ m_pPermutStorage[idx] = p; }
	CC T rowToChange(T nRow) const;
	void reconstructSolution(const CColOrbit<T> *pColOrbitStart, const CColOrbit<T> *pColOrbit, 
		size_t colOrbLen, const CColOrbit<T> *pColOrbitIni, const T *pRowPerm, const VECTOR_ELEMENT_TYPE *pRowSolution, size_t solutionSize);
	void UpdateOrbits(const T *permut, T lenPerm, T *pOrbits, bool rowPermut, bool updateGroupOrder = false);
#if USE_STRONG_CANONICITY
	inline void setSolutionStorage(CSolutionStorage *p) { m_pSolutionStorage = p; }
	inline CSolutionStorage *solutionStorage() const { return m_pSolutionStorage; }

	CSolutionStorage *m_pSolutionStorage;
#else
#define setSolutionStorage(x)
#define solutionStorage()		(CSolutionStorage *)NULL
#endif

	T m_nStabLength;
	T m_nStabLengthAut;
    int m_rank;
    CPermut *m_pPermutRow;
    CPermut *m_pPermutCol;
    CColNumbStorage **m_nColNumbStorage;
	CPermutStorage<T> *m_pPermutStorage[3];
	T *m_pObits[2][2];
    CCounter<int> *m_pCounter;
	T *m_pColIndex;
	VECTOR_ELEMENT_TYPE *m_pImprovedSol;
	const uint m_enumFlags;
	uint m_nGroupOrder;
	T m_nNumRow;
	T m_nNumColOrb;
};

template<class T>
CCanonicityChecker<T>::CCanonicityChecker(T nRow, T nCol, int rank, uint enumFlags) : 
	        m_rank(rank), m_enumFlags(enumFlags)
{
	m_pPermutRow = new CPermut(nRow);
	m_pPermutCol = new CPermut(nCol);
	setColIndex(new T[nCol << 1]);
	m_pCounter = new CCounter<int>(rank);
	m_nColNumbStorage = new CColNumbStorage *[rank];
	for (int i = rank; i--;)
		m_nColNumbStorage[i] = new CColNumbStorage(nCol);

	m_pImprovedSol = new VECTOR_ELEMENT_TYPE[(rank - 1) * nCol];
	setSolutionStorage(new CSolutionStorage());
	const auto mult = enumFlags & t_outStabilizerOrbit ? 2 : 1;
	setPermStorage(new CPermutStorage<T>());

	memset(m_pObits, 0, sizeof(m_pObits));
	const auto outColumnOrbits = enumFlags & t_outColumnOrbits;
	const auto len = nRow + (outColumnOrbits ? nCol : 0);
	auto pntr = m_pObits[0][0] = new T[mult * len];	// memory to keep orbits of different types
	if (outColumnOrbits) {
		setPermStorage(new CPermutStorage<T>(), 1);
		m_pObits[1][0] = pntr + mult * nRow;
		if (mult > 1)
			m_pObits[1][1] = m_pObits[1][0] + nCol;
	}
	else
		setPermStorage(NULL, 1);

	setPermStorage(enumFlags & t_alwaisKeepRowPermute ? new CPermutStorage<T>() : NULL, 2);

	if (mult > 1)
		m_pObits[0][1] = pntr + nRow;
}

template<class T>
CCanonicityChecker<T>::~CCanonicityChecker()
{
	delete m_pPermutRow;
	delete m_pPermutCol;
	delete counter();
	delete permStorage();
	delete permColStorage();
	delete permRowStorage();
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

template<class T>
bool CCanonicityChecker<T>::TestCanonicity(T nRowMax, const CMatrixCol<T> *pEnum, int outInfo, T *pRowOut, CRowSolution<T> *pRowSolution)
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
	const VECTOR_ELEMENT_TYPE *pCurrSolution;
	size_t solutionSize;
	if (pRowSolution) {
		// Because the position of current solution (pRowSolution->solutionIndex()) 
		// could be changed, let's take pointer here
		// Don't worry about that and similar compilation warnings: For now we don't plan to call this method with pRowSolution != NULL
		pCurrSolution = pRowSolution->currSolution();
		solutionSize = pRowSolution->solutionSize();
#if USE_STRONG_CANONICITY
		solutionStorage()->clear();
		pRowSolution->setLenOrbitOfSolution(0);
#endif
	}

	size_t startIndex = 0;
#endif

	T *pColIndex = NULL;
	T *pVarPerm = NULL;
	bool retVal = true;
	T nRow = MATRIX_ELEMENT_MAX;
	while (true) {
	next_permut:
		nRow = next_permutation(nRow);
		if (nRow == MATRIX_ELEMENT_MAX)
			break;

		OUTPUT_PERMUTATION(this, pEnum->outFile(), orbits(), numRow());

		// Loop for all remaining matrix's rows
		for (; nRow <= nRowMax; nRow++) {
			const auto *pRow = pMatr->GetRow(nRow);
			const auto *pRowPerm = pMatr->GetRow(*(permRow() + nRow));
			const auto *pColOrbitIni = colOrbitIni[nRow];
			const auto *pColOrbit = colOrbit[nRow];
			while (pColOrbit) {
				// Define the number of column to start with
				const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
				const auto orbLen = pColOrbit->length();
				int diff;
				if (orbLen > 1)
					diff = checkColOrbit(orbLen, nColCurr, pRow, pRowPerm);
				else
					diff = (int)*(pRow + nColCurr) - *(pRowPerm + *(permCol() + nColCurr));

				if (diff > 0)
					goto next_permut;

				if (!diff) {
					pColOrbit = pColOrbit->next();
					continue;
				}

#if RECURSIVE_CANON
				if (orbLen > 1)
					return false;

				// Construct matrix, for just found permRow() and permCol()
				CMatrix *pMatr = matrix();
				CMatrix matr(pMatr, permRow(), permCol());
				CCanonicityChecker canon(rowNumb(), colNumb(), rank());
				setMatrix(&matr);
				std::cout << "########";
				matr.printOut();
				setCanonChecker(&canon);

				// Check, if this matrix is canonical
				int level;
				bool retVal = TestCanonicity(level);
				setMatrix(pMatr);
#else
				if (pRowOut) {
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
								reconstructSolution(colOrbit[nRow], pColOrbit, colOrbLen, pColOrbitIni, pRowPerm, pCurrSolution, solutionSize);
#if USE_STRONG_CANONICITY
								if (solutionStorage()) {
									for (auto *pSolution : *solutionStorage()) {
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

				return false;
#endif
			}
		}

		// Automorphism found:
#if (PRINT_CURRENT_MATRIX && PRINT_PERMUTATION)
		outString("-----\n", pEnum->outFile());
#endif
		if (!rowPermut) {
			// We are here to define the canonicity of partially constructed 
			// matrix AND we just found the non-trivial automorphism.
			// Let's construct the permutation of the column's orbits 
			// which corresponds to just found automorphism
			const auto *pColOrbitIni = colOrbitIni[nRow - (rowPermut? 2 : 0)];
			const auto *pColOrbit = colOrbit[nRow - (rowPermut ? 2 : 0)];
			if (!pColIndex) // Index for columns was not yet constructed
				pVarPerm = (pColIndex = constructColIndex(pColOrbit, pColOrbitIni, colOrbLen)) + colNumb;

			size_t varIdx = 0;
			while (pColOrbit) {
				const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
				*(pVarPerm + varIdx++) = *(pColIndex + *(permCol() + nColCurr));
				pColOrbit = pColOrbit->next();
			}
		}
		
		// We need the permutations on columns AND
		//  (a) matrix is completely constructed OR
		//  (b) we will need to analyse the group on partially constructed matrix
		// (As of Oct.11, 2018 this is used only for IGraphs (inconsistent graphs)
		if (permColStorage() && (rowPermut || permRowStorage()))
			ConstructColumnPermutation(pEnum->matrix());

		addAutomorphism(rowPermut);
		nRow = MATRIX_ELEMENT_MAX - 1;
	}

	if (rowPermut)
		updateGroupOrder();

	return retVal;
}

#endif /* defined(__BIBD_Mac__CanonicityChecker__) */
