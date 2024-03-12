//
//  CanonicityChecker.h
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#ifndef __BIBD_Mac__CanonicityChecker__
#define __BIBD_Mac__CanonicityChecker__
#endif
#pragma once

#include "DataTypes.h"
#include "PermutStorage.h"
#include "ColOrbits.h"
#include "GroupOnParts.h"
#include "CheckCanon.h"
#include "GroupOrder.h"


#define CHECK_PERMUTS		0 // Check permutations used for combined BIBD enumeration

#define CPermut				CSimpleArray<S>
#define CColNumbStorage		CContainer<S>
#define IDX_MAX				(ELEMENT_MAX - 1)

typedef enum {
	t_saveNothing			= 0,
	t_saveRowPermutations	= (1 << 0),
	t_saveColPermutations	= (1 << 1),
	t_saveRowToChange		= (1 << 2)
} t_canonOutInfo;

Class2Def(CMatrix);
Class2Def(CMatrixCol);
Class2Def(CRowSolution);

template <typename T, typename S>
struct TestCanonParams {
	const MatrixColPntr pEnum;
	const MatrixDataPntr pMatrix;
	T numParts;
	T* pPartNumb;
	T* pRowOut;
	CGroupOnParts<T>* pGroupOnParts;
	MatrixDataPntr pSpareMatrix;         // used when pGroupOnParts != NULL
	T* pPermCol;
	T startingRowNumb;                   // starting row for the loop in TestCanonicity (used when pPermCol != NULL)
};


Class2Def(CCanonicityChecker) : public CGroupOrder<T>, public CRank {
public:
	CC CCanonicityChecker(T nRow, T nCol, T rank = 2, uint enumFlags = t_enumDefault, T numParts = 1);
	CC ~CCanonicityChecker();
	void InitCanonicityChecker(T nRow, T nCol, int rank, char *pMem);
	CC bool TestCanonicity(T nRowMax, const TestCanonParams<T, S> *pCanonParam, uint outInfo = t_saveNothing, RowSolutionPntr pRowSolution = NULL);
	void outputAutomorphismInfo(FILE *file, const MatrixDataPntr pMatrix = NULL) const;
	CC auto enumFlags() const						{ return m_enumFlags; }
	CC inline T *permRow() const					{ return m_pPermutRow->elementPntr(); }
	CC inline T *permCol() const					{ return m_pPermutCol->elementPntr(); }
	CC inline auto permStorage(T nPart) const		{ return permStorage() + nPart; }
	CC inline auto permStorage() const				{ return m_pPermutStorage[0]; }
	CC inline auto permColStorage() const			{ return m_pPermutStorage[1]; }
	CC inline auto permRowStorage() const			{ return m_pPermutStorage[2]; }
	CC inline auto improvedSolution() const			{ return m_pImprovedSol; }
	CC bool groupIsTransitive() const;
	bool printMatrix(const designParam *pParam) const;
	CC auto stabiliserLengthExt() const				{ return m_nStabExtern; }
	CK virtual CGroupOrder<T>* extraGroupOrder() const { return NULL; }
	inline bool CheckCanonicity(const T* result, int nLines, int *pGrpNumb, T* bResult = NULL) {
		return m_pCheckerKSystemCanon->CheckCanonicity(result, nLines, pGrpNumb, bResult);
	}
	inline bool CheckPermutations(const T* inputMatrix, const T* bestResult, int nRows) {
		return m_pCheckerKSystemCanon->CheckPermutations(inputMatrix, bestResult, nRows);
	}
	inline bool improvedResultIsReady(t_bResultFlags flag = t_bResultFlags::t_readyCompletely) const {
		return m_pCheckerKSystemCanon->improvedResultIsReady(flag);
	}
	inline char* comment() const {
		return m_pCheckerKSystemCanon->comment();
	}
	inline void setPreordered(bool flag) {
		m_pCheckerKSystemCanon->setPreordered(flag);
	}
protected:
	void updateCanonicityChecker(T rowNumb, T colNumb);
	CC virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix)		{}
	virtual void CanonizeByColumns(MatrixDataPntr pMatrix, T *pColIdxStorage = NULL, CanonicityCheckerPntr pCanonChecker = NULL) const	{}
	CC inline auto getRowOrbits(int idx) const		{ return m_pObits[0][idx]; }
	CC inline auto getColOrbits(int idx) const		{ return m_pObits[1][idx]; }
	inline bool checkProperty(uint flag) const		{ return enumFlags() & flag; }
	CC inline auto numRow() const					{ return m_nNumRow; }
	CC void setStabiliserLengthExt(T len)			{ m_nStabExtern = len; }
	CC inline auto numParts() const					{ return m_numParts; }
	CC virtual T lenStabilizer() const				{ return 0; }
	CK inline auto shiftToUnforcedOrbit(T nRow) const { return m_pShift[nRow]; }
	CC virtual void resetGroupOrder()				{}
	CC virtual void incGroupOrder()					{}
	CK inline void setGroupOnParts(CGroupOnParts<T>* pntr) { m_pGroupOnParts = pntr; }
	CK inline auto getGroupOnParts() const			{ return m_pGroupOnParts; }
	CK virtual CGroupOnParts<T>* makeGroupOnParts(const CCanonicityChecker *owner) { return NULL; }
private:
	CC T *init(T nRow, T numParts, bool savePerm, T *pOrbits, T **pPermRows, bool groupOnParts, T* pPermCol = NULL);
	CC void addAutomorphism(const T nRow, const T *pRowPerm, T *pOrbits, bool rowPermut = true, bool savePermut = false, bool calcGroupOrder = true);
	CC int checkColOrbit(T orbLen, T nColCurr, const S *pRow, const T *pRowPerm, T *pColPerm) const;
	CC inline void setNumRow(T nRow)				{ m_nNumRow = nRow; }
	CC inline auto numCol() const					{ return static_cast<T>(m_pPermutCol->numElement()); }
#define orbits()	m_pObits[0][0] 
	CC inline auto colNumbStorage() const			{ return m_nColNumbStorage; }
	CC inline auto counter() const					{ return m_pCounter; }
	CC inline void setColIndex(T *p)				{ m_pColIndex = p; }
	CC inline auto colIndex() const					{ return m_pColIndex; }

	CC T constructColIndex(const ColOrbPntr pColOrbit, const ColOrbPntr pColOrbitIni, size_t colOrbLen, T shift = 0) const;
	CC inline void setPermStorage(PermutStoragePntr p, int idx = 0)	{ m_pPermutStorage[idx] = p; }
	CC T rowToChange(T nRow) const;
	void reconstructSolution(const ColOrbPntr pColOrbitStart, const ColOrbPntr pColOrbit,
		size_t colOrbLen, const ColOrbPntr pColOrbitIni, const T *pRowPerm, const T *pRowSolution, size_t solutionSize);
//	CC void UpdateOrbits(const T *permut, const T lenPerm, T *pOrbits, bool rowPermut, bool updateGroupOrder = false);
#if CHECK_PERMUTS
	CC void check_permut(const T* permut, T len_perm) const;
#else
#define check_permut(permut, len_perm)
#endif

#if USE_STRONG_CANONICITY
	inline void setSolutionStorage(CSolutionStorage *p) { m_pSolutionStorage = p; }
	inline CSolutionStorage *solutionStorage() const { return m_pSolutionStorage; }

	CSolutionStorage *m_pSolutionStorage;
#else
#define setSolutionStorage(x)
#define solutionStorage()		(CSolutionStorage *)NULL
#endif

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	CPermut *m_pPermutRow;
	CPermut *m_pPermutCol;
	CPermut* m_pPermutSparse = NULL;

	CColNumbStorage **m_nColNumbStorage;
	PermutStoragePntr m_pPermutStorage[3];
	T *m_pObits[2][2];
	CCounter<int> *m_pCounter;
	T *m_pColIndex;
	T *m_pImprovedSol;
	const uint m_enumFlags;
	T m_nNumRow;
	const T m_numParts;
	CGroupOnParts<T>* m_pGroupOnParts = NULL;
	T* m_pTrivialPermutCol = NULL;     // Trivial permutation on columns
	const T m_numElem;				   // The number of elements that will be the same for all partially constructed objects
	                                   // (it is equal nCol for combinatorial designs or number of players for k-system) 
	CCheckerCanon<SIZE_TYPE> *m_pCheckerKSystemCanon = NULL;
};

CanonicityChecker()::CCanonicityChecker(T nRow, T nCol, T rank, uint enumFlags, T numParts) : 
	CRank(nRow, rank), m_numElem(nCol), m_enumFlags(enumFlags), m_numParts(numParts)
{
	m_pPermutRow = new CPermut(enumFlags & t_EnumeratorFlags::t_kSystems? nCol : nRow);
	m_pPermutCol = new CPermut(nCol);
	setColIndex(new T[2 * nCol]);
	m_pCounter = new CCounter<int>(rank);
	m_nColNumbStorage = new CColNumbStorage *[rank];
	for (auto i = rank; i--;)
		m_nColNumbStorage[i] = new CColNumbStorage(nCol);

	m_pImprovedSol = new T[(rank - 1) * nCol];
	setSolutionStorage(new CSolutionStorage());
	const auto mult = enumFlags & t_outStabilizerOrbit ? 2 : 1;

	const auto outColumnOrbits = enumFlags & t_outColumnOrbits;
	const auto keepRowPermute = enumFlags & (t_alwaysKeepRowPermute | t_outRowPermute);
	const auto nPermutStorages = numParts + (outColumnOrbits ? 1 : 0) + (keepRowPermute ? 1 : 0);
	auto pPermStorage = new CPermutStorage<T, S>[nPermutStorages]();
	setPermStorage(pPermStorage);

	memset(m_pObits, 0, sizeof(m_pObits));

	// length of memory to keep orbits of spare matrix, if it is used
	const auto extraOrbLen = enumFlags & t_useGroupOnParts ? nRow : 0;
	// memory to keep orbits of different types
	const auto len = nRow + (outColumnOrbits ? nCol : 0);
	auto pntr = m_pObits[0][0] = new T[mult * len + extraOrbLen];
	pntr += extraOrbLen;
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

	m_pTrivialPermutCol = new T[nCol];
	for (auto i = nCol; i--;)
		m_pTrivialPermutCol[i] = i;

	if (enumFlags & t_EnumeratorFlags::t_kSystems &&
		(!m_pCheckerKSystemCanon || m_pCheckerKSystemCanon->numDays() != nRow)) {
		delete m_pCheckerKSystemCanon;
		m_pCheckerKSystemCanon = new CCheckerCanon<unsigned char>(nRow, m_numElem, 3);
	}
}

CanonicityChecker()::~CCanonicityChecker()
{
	delete m_pPermutRow;
	delete m_pPermutCol;
	delete[] m_pPermutSparse;

	delete counter();
	delete [] permStorage();
	for (int i = rank(); i--;)
		delete m_nColNumbStorage[i];

	delete[] colNumbStorage();
	delete[] colIndex();
	delete[] improvedSolution();
	delete[] m_pObits[0][0];
	delete[] m_pTrivialPermutCol;
	if (getGroupOnParts() && getGroupOnParts()->owner() == this)
		delete getGroupOnParts();

	delete m_pCheckerKSystemCanon;
#if USE_STRONG_CANONICITY
	delete solutionStorage();
#endif
}

CanonicityChecker(void)::updateCanonicityChecker(T rowNumb, T colNumb)
{
	// When we use automorphism on parts, we need two copies each of the following
	m_pPermutSparse = new CPermut[2];
	m_pPermutSparse[0].Init(rowNumb, new T[rowNumb]);
	m_pPermutSparse[1].Init(rowNumb, new T[colNumb]);
}

#if CHECK_PERMUTS
CanonicityChecker(void)::check_permut(const T* permut, T len_perm) const {
	for (auto i = len_perm; i--;) {
		for (auto j = len_perm; j-- && permut[j] != i;);
		assert(j != ELEMENT_MAX);
	}
}
#endif

CanonicityChecker(bool)::TestCanonicity(T nRowMax, const TestCanonParams<T, S>* pCanonParam, uint outInfo, RowSolutionPntr pRowSolution)
{
	const auto lenStab = stabiliserLengthExt();
	if (nRowMax == lenStab - 1)
		return true;

	const auto* pEnum = pCanonParam->pEnum;
	auto *pPartNumb = pCanonParam->pPartNumb;
	auto *pRowOut = pCanonParam->pRowOut;
	const auto* pGroupOnParts = pCanonParam->pGroupOnParts;
	// Construct trivial permutations for rows and columns
	const auto rowPermut = outInfo & t_saveRowPermutations;
	const auto* pMatr = pCanonParam->pMatrix;
	const auto* pMatrPerm = pMatr;
	const auto numParts = pCanonParam->numParts;
	bool check_trivial_row_perm = pCanonParam->pPermCol != NULL;
	auto savePermut = rowPermut && (enumFlags() & t_outRowPermute);

	const auto colOrbLen = pEnum->colOrbitLen();
	const auto* pPartInfo = pMatr->partsInfo();
	const auto nonCombinedDesign = numParts == 1;
	
	// Reset permutation counters, if used
	PREPARE_PERM_OUT(permColStorage());
	PREPARE_PERM_OUT(permRowStorage());

	bool calcGroupOrder = true;
	bool retVal = true;
	bool usingGroupOnBlocks = false;
#if !USE_COL_PERMUT
	size_t lenMatr;
	S* pMatrTo, * pMatrFrom;
	T colNumb;
#endif
	size_t numGroups = 0;
	T startingRowNumb = pCanonParam->startingRowNumb;
	size_t idxPerm[16] = {};	// to keep the indices of currently used permutation
	size_t* pIndxPerms = NULL;	// of i-th symmetrical group acting on the parts

	T idxPartSrc[16] = {};		// to keep the initial (source) indices of the parts
	T* pPartSrc = numParts <= countof(idxPartSrc) ? idxPartSrc : new T[numParts];
	for (auto i = numParts; i--;)
		pPartSrc[i] = i;

	const T* pCurrSolution;
	size_t solutionSize;
	if (pRowSolution) {
		// Because the position of current solution (pRowSolution->solutionIndex()) 
		// can be changed, take pointer here
		pCurrSolution = pRowSolution->currSolution();
		solutionSize = pRowSolution->solutionLength();
#if USE_STRONG_CANONICITY
		solutionStorage()->clear();
		pRowSolution->setLenOrbitOfSolution(0);
#endif
	}

	size_t startIndex = 0;
	T* permColumn = pCanonParam->pPermCol;
	const auto len_stab = lenStabilizer();
	auto pOrbits = orbits();
	while (true) {
#ifndef USE_CUDA
		T* permRows = NULL;
		permColumn = init(nRowMax, numParts, rowPermut, pOrbits, &permRows, usingGroupOnBlocks, permColumn);

		T nRow = ELEMENT_MAX;
		if (check_trivial_row_perm) {
			nRow = startingRowNumb;
			goto try_permut;
		}
		while (true) {

		next_permut:
			nRow = next_permutation(permRows, pOrbits, numRow(), nRow, lenStab);
			if (nRow == ELEMENT_MAX || nRow < len_stab)
				break;

		try_permut:
#if 0  // Activate to keep/print tested permutations  
			OUT_PERM(permRowStorage(), permRows, nRowMax);
#endif
			OUTPUT_PERMUTATION(permColStorage(), pEnum->outFile(), pOrbits, numRow());

			// Loop for all remaining matrix's rows
			for (; nRow < nRowMax; nRow++) {
				const auto* pRow = pMatr->GetRow(nRow);

				const auto* pRowPerm = pMatrPerm->GetRow(*(permRows + nRow));
				for (T nPart = 0; nPart < numParts; nPart++) {
					const auto nPartSrc = pPartSrc[nPart];
					const auto* pColOrbitIni = pEnum->colOrbitIni(nRow, nPartSrc);
					const auto* pColOrbit = pEnum->colOrbit(nRow, nPartSrc);
					const auto shift = nPart ? pPartInfo->getShift(nPart) : 0;
					while (pColOrbit) {
						// Define the number of column to start with
						const auto nColCurr = shift + static_cast<S>(((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen);
						const auto orbLen = pColOrbit->length();
						int diff;
						if (orbLen > 1)
							diff = checkColOrbit(orbLen, nColCurr, pRow, pRowPerm, permColumn);
						else
							diff = static_cast<int>(*(pRow + nColCurr)) - *(pRowPerm + *(permColumn + nColCurr));

						if (diff > 0)
							goto next_permut;

						if (!diff) {
							pColOrbit = pColOrbit->next();
							continue;
						}

						if (!usingGroupOnBlocks && nonCombinedDesign && pRowOut) {
							if (outInfo & t_saveRowToChange)
								*pRowOut = rowToChange(nRow);
#ifndef USE_CUDA
							else {
								if (pRowSolution && !pRowSolution->isLastSolution()) {
									// NOTE: No need to remove last solution, we will leave this level anyway
									if (nRow + 1 < nRowMax) {
										// We can remove all solutions which are in the same group with the solution just tested.
										// We don't need them even as the right parts of our equations, because the usage of everyone 
										// will make the matrix non canonical
										pRowSolution->removeNoncanonicalSolutions(startIndex);
									}
									else {
										// The current solution is not canonical.
										// Let's define an equivalent solution that is lexicographically greater than the current one.
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
										// Move the current solution and its orbit closer to just found canonical solution
										pRowSolution->moveNoncanonicalSolutions(improvedSolution(), 0, solutionStorage());
#endif
									}
								}
							}
#endif				
						}

						if (pPartNumb)
							*pPartNumb = -1; // usingGroupOnBlocks ? -1 : nPart;

						if (pIndxPerms != idxPerm)
							delete[] pIndxPerms;

						if (pPartSrc != idxPartSrc)
							delete[] pPartSrc;

						return false;
					}
				}
			}

			// Automorphism found:
			if (!usingGroupOnBlocks) {
				// We are not using group on blocks yet. If pGroupOnParts != NULL, we will use it on next path.
#if (PRINT_CURRENT_MATRIX && PRINT_PERMUTATION)
				outString("-----\n", pEnum->outFile());
#endif
				if (!rowPermut) {
					// We are here to define the canonicity of partially constructed matrix AND we just found the non-trivial automorphism.
					// Let's construct the permutation of the column's orbits which corresponds to just found automorphism
					for (T nPart = 0; nPart < numParts; nPart++) {
						auto pPermStorage = permStorage(nPart);
						const auto* pColOrbitIni = pEnum->colOrbitIni(nRow, nPart);
						const auto* pColOrbit = pEnum->colOrbit(nRow, nPart);
						const auto shift = nPart ? pPartInfo->getShift(nPart) : 0;
						// Saving permutations, acting on the orbits of columns
						T lenPermut;
						if (pPermStorage->isEmpty()) {		// Index for columns was not yet constructed
							lenPermut = constructColIndex(pColOrbit, pColOrbitIni, colOrbLen, shift);
							pPermStorage->allocateMemoryForPermut(lenPermut);	// Memory for identity permutation just in case we will need it
						}
						else
							lenPermut = pPermStorage->lenPerm();

						auto pVarPerm = pPermStorage->allocateMemoryForPermut(lenPermut);
						while (pColOrbit) {
							const auto nColCurr = shift + ((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen;
							*pVarPerm++ = *(colIndex() + *(permColumn + nColCurr));
							pColOrbit = pColOrbit->next();
						}
					}
				}

				//#ifndef USE_CUDA
						// We need the permutations on columns AND
						//  (a) matrix is completely constructed OR
						//  (b) we will need to analyse the group on partially constructed matrix
						// (As of Oct.11, 2018 this is used only for  Semi-Symmetric Graphs (IGraphs)
				if (permColStorage() && (rowPermut || permRowStorage()))
					ConstructColumnPermutation(pEnum->matrix());
				//#endif
			}

			if (!calcGroupOrder) {
				// If we are here then it is possible to get the SAME matrix by some permutation of its parts and rows
				// So we don't need to continue try different permutations of rows, because if for some permutation of
				// the rows we would find non-canonicity, we would find it earlier without doing permutation of parts.
				if (rowPermut) {
					// We are checking completely constructed matrix with the nontrivial group action on its parts.
					incGroupOrder();
				}

				break;
			}

			addAutomorphism(nRowMax, permRows, pOrbits, rowPermut, savePermut, calcGroupOrder);
			nRow = IDX_MAX;
		}

		if (rowPermut && calcGroupOrder)
			updateGroupOrder(nRowMax, pOrbits);

		if (!retVal || !pGroupOnParts || !pGroupOnParts->useGroupOnParts(nRowMax))
			break;    // We don't have to test on groups of blocks

		if (!usingGroupOnBlocks) {
			// Corresponding data were not initialized yet
			usingGroupOnBlocks = true;
			numGroups = pGroupOnParts->numGroups();
			// Initialization of variables used with the group acting on the parts of the matrix
			pIndxPerms = numGroups < countof(idxPerm) ? idxPerm : new size_t[numGroups];
			memset(pIndxPerms, 0, numGroups * sizeof(*pIndxPerms));
#if USE_COL_PERMUT
			startingRowNumb = pGroupOnParts ? 1 : 0;
#else
			colNumb = (pMatrPerm = pCanonParam->pSpareMatrix)->colNumb();
			startingRowNumb = pGroupOnParts->getStartingRowNumb();
			// copy first startingRowNumb rows
			lenMatr = colNumb * startingRowNumb;
			pMatrTo = pMatrPerm->GetDataPntr();
			pMatrFrom = pMatr->GetDataPntr();
			memcpy(pMatrTo, pMatrFrom, lenMatr * sizeof(*pMatrTo));

			// Calculate pointers and length for copying remaining part of the matrix
			pMatrTo += lenMatr;
			pMatrFrom += lenMatr;
			lenMatr = (nRowMax - startingRowNumb) * colNumb * sizeof(*pMatrTo);
#endif
			calcGroupOrder = savePermut = false;	// Permutations and group order should not be saved for Spare Matrix
			pOrbits += numRow();
			resetGroupOrder();
		}

		// Get next set of indices of permutations used for each group
		size_t idx = 0;
		while (idx < numGroups) {
			if (++*(pIndxPerms + idx) >= pGroupOnParts->groupHandle(idx)->groupOrder())
				*(pIndxPerms + idx++) = 0;
			else break;
		}
		if (idx >= numGroups)
			break;

		memcpy(permColumn, m_pTrivialPermutCol, sizeof(permColumn[0]) * pMatr->colNumb());
#if !USE_COL_PERMUT
		// Permuting the parts of the matrix in accordance with the current set of permutations
		memcpy(pMatrTo, pMatrFrom, lenMatr);
#endif
		for (auto i = numParts; i--;)
			pPartSrc[i] = i;

		const auto* pPartsInfo = pMatr->partsInfo();
		for (size_t idx = 0; idx < numGroups; idx++) {
			const auto idxPerm = *(pIndxPerms + idx);
			// For the current set of parts, check if all parts remain in place.
			if (!idxPerm)
				continue;

			const auto pntr = pGroupOnParts->groupHandle(idx);
			// Address of next permutation to try for current group
			const auto* pPerm = pntr->getPermutation(idxPerm);
			const auto idxPart = pntr->partIdx();  // absolute index of the first part moved by current symmetrical group
			const auto nCol = pPartsInfo->colNumb(idxPart);
#if !USE_COL_PERMUT
			const size_t lenPart = nCol * sizeof(*pMatrTo);
#endif
			for (auto jFrom = 0; jFrom < pntr->permLength(); ++jFrom) {
				// For the current set of parts, check if the j-th part remains in place.
				auto jTo = *(pPerm + jFrom);
				if (jFrom == jTo)
					continue;

				// Calculate absolute indices and save the source index
				const auto idxTo = pPartsInfo->getShift(jTo += idxPart);
				const auto idxFrom = pPartsInfo->getShift(pPartSrc[jTo] = jFrom + idxPart);
#if USE_COL_PERMUT
				// Copying a set of column indices corresponding to the moved part of the Combined BIBD
				memcpy(permColumn + idxTo, m_pTrivialPermutCol + idxFrom, sizeof(permColumn[0]) * nCol);
#else
				auto* pPartTo = pMatrTo + idxTo;
				const auto* pPartFrom = pMatrFrom + idxFrom;
				for (auto j = startingRowNumb; j < nRowMax; j++) {
					memcpy(pPartTo, pPartFrom, lenPart);
					pPartTo += colNumb;
					pPartFrom += colNumb;
				}
#endif
			}
		}

		check_permut(pPartSrc, numParts);
		check_permut(permColumn, pMatr->colNumb());

		// After permitation of parts we also need 
		// to try trivial permutation on rows
		check_trivial_row_perm = true; 
		continue;
	}

	if (pIndxPerms != idxPerm)
		delete[] pIndxPerms;

	if (pPartSrc != idxPartSrc)
		delete[] pPartSrc;

	return retVal;
}

#endif /* defined(__BIBD_Mac__CanonicityChecker__) */
