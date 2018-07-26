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
	t_saveRowPermutations	= 1,
	t_saveRowToChange		= (1 << 1)
} t_canonOutInfo;

template<class T> class CMatrixCol;
template<class T> class CRowSolution;

template<class T>
class CCanonicityChecker : public CPermut
{
public:
    CC CCanonicityChecker(T nRow, size_t nCol, int rank);
    CC ~CCanonicityChecker();
	void InitCanonicityChecker(T nRow, size_t nCol, int rank, T *pMem);
	CC bool TestCanonicity(T nRowMax, const CMatrixCol<T> *pEnum, int outInfo = 0, T *pRowOut = NULL, CRowSolution<T> *pRowSolution = NULL);
    void outputAutomorphismInfo(FILE *file) const;
	CC inline uint groupOrder() const				{ return m_nGroupOrder; }
	CC inline T numColOrb() const					{ return m_nNumColOrb; }
	CC inline void setGroupOrder(uint val)			{ m_nGroupOrder = val; }
	CK inline T lenPerm() const						{ return permStorage()->lenPerm(); }
	CK inline void adjustGenerators(int *idx, T len){ permStorage()->adjustGenerators(idx, len); }
	CK inline size_t constructGroup()				{ return permStorage()->constructGroup(); }
	CK inline size_t findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon)
													{ return permStorage()->findSolutionIndex(pFirst, idx, pMem, pCanonIdx, nCanon); }
	CC inline CPermutStorage<T> *permStorage() const				{ return m_pPermutStorage; }
	CC inline VECTOR_ELEMENT_TYPE *improvedSolution() const			{ return m_pImprovedSol; }
	CC bool groupIsTransitive() const;
	bool printMatrix(const designRaram *pParam) const;

protected:
	CC inline int rank() const						{ return m_rank; }
private:
    CC void init(T nRow, bool savePerm);
	CC T next_permutation(T idx = MATRIX_ELEMENT_MAX);
    CC void addAutomorphism(bool rowPermut = true);
    CC int checkColOrbit(size_t orbLen, size_t nColCurr, const T *pRow, const T *pRowPerm) const;
    CC inline T * permRow() const					{ return m_pPermutRow->elementPntr(); }
    CC inline T * permCol() const					{ return m_pPermutCol->elementPntr(); }
    CC inline void setStabilizerLength(T len)		{ m_nStabLength = len; }
    CC inline T stabilizerLength() const			{ return m_nStabLength; }
    CC inline void setStabilizerLengthAut(T l)		{ m_nStabLengthAut = l; }
    CC inline T stabilizerLengthAut() const			{ return m_nStabLengthAut; }
	CC inline void setNumRow(T nRow)				{ m_nNumRow = nRow; }
	CC inline T numRow() const						{ return m_nNumRow; }
    CC inline T numCol() const						{ return static_cast<T>(m_pPermutCol->numElement()); }
#define orbits()	this->elementPntr()
	CC void updateGroupOrder();
    CC inline CColNumbStorage **colNumbStorage() const			{ return m_nColNumbStorage; }
    CC inline CCounter<int> *counter() const					{ return m_pCounter; }
	CC inline void setColIndex(T *p)				{ m_pColIndex = p; }
	CC inline T *colIndex() const					{ return m_pColIndex; }
	CC inline void setNumColOrb(T v)				{ m_nNumColOrb = v; }
    CC void revert(T i);
	CC T getLenPermutCol(T **permCol) const;
	CC T * constructColIndex(const CColOrbit<T> *pColOrbit, const CColOrbit<T> *pColOrbitIni, size_t colOrbLen);
	CC inline void setPermStorage(CPermutStorage<T> *p){ m_pPermutStorage = p; }
	CC T rowToChange(T nRow) const;
	void reconstructSolution(const CColOrbit<T> *pColOrbitStart, const CColOrbit<T> *pColOrbit, 
		size_t colOrbLen, const CColOrbit<T> *pColOrbitIni, const T *pRowPerm, const VECTOR_ELEMENT_TYPE *pRowSolution, size_t solutionSize);
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
	CPermutStorage<T> *m_pPermutStorage;
    CCounter<int> *m_pCounter;
	T *m_pColIndex;
	VECTOR_ELEMENT_TYPE *m_pImprovedSol;

	uint m_nGroupOrder;
	T m_nNumRow;
	T m_nNumColOrb;
};

template<class T>
CCanonicityChecker<T>::CCanonicityChecker(T nRow, size_t nCol, int rank) : CPermut(nRow), m_rank(rank)
{
	m_pPermutRow = new CPermut(nRow);
	m_pPermutCol = new CPermut(nCol);
	setColIndex(new T[nCol << 1]);
	m_pCounter = new CCounter<int>(rank);
	setPermStorage(new CPermutStorage<T>());
	m_nColNumbStorage = new CColNumbStorage *[rank];
	for (int i = rank; i--;)
		m_nColNumbStorage[i] = new CColNumbStorage(nCol);

	m_pImprovedSol = new VECTOR_ELEMENT_TYPE[(rank - 1) * nCol];
	setSolutionStorage(new CSolutionStorage());
}

template<class T>
CCanonicityChecker<T>::~CCanonicityChecker()
{
	delete m_pPermutRow;
	delete m_pPermutCol;
	delete counter();
	delete permStorage();
	for (int i = rank(); i--;)
		delete m_nColNumbStorage[i];

	delete[] colNumbStorage();
	delete[] colIndex();
	delete[] improvedSolution();
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
			// We are here to define the canonicity of partially constructed matrix AND
			// we just found the nontrivial automorphism.
			const auto *pColOrbitIni = colOrbitIni[nRow];
			const auto *pColOrbit = colOrbit[nRow];
			if (!pColIndex) // Index for columns was not yet constructed
				pVarPerm = (pColIndex = constructColIndex(pColOrbit, pColOrbitIni, colOrbLen)) + colNumb;

			size_t varIdx = 0;
			while (pColOrbit) {
				const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
				*(pVarPerm + varIdx++) = *(pColIndex + *(permCol() + nColCurr));
				pColOrbit = pColOrbit->next();
			}
		}

		addAutomorphism(rowPermut);
		nRow = MATRIX_ELEMENT_MAX - 1;
	}

	if (rowPermut)
		updateGroupOrder();

	return retVal;
}

#endif /* defined(__BIBD_Mac__CanonicityChecker__) */
