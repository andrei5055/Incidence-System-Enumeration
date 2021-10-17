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
Class1Def(CGroupOnParts);

template <typename T, typename S>
struct TestCanonParams {
	const MatrixColPntr pEnum;
	S* pPartNumb;
	S* pRowOut;
	CGroupOnParts<T>* pGroupOnParts;
	MatrixDataPntr pSpareMatrix;         // used when pGroupOnParts != NULL
};

Class2Def(CCanonicityChecker) : public CRank {
public:
	CC CCanonicityChecker(T nRow, T nCol, int rank = 2, uint enumFlags = t_enumDefault, S numParts = 1);
	CC ~CCanonicityChecker();
	void InitCanonicityChecker(T nRow, T nCol, int rank, char *pMem);
	CC bool TestCanonicity(T nRowMax, const TestCanonParams<T, S> *pCanonParam, uint outInfo, RowSolutionPntr pRowSolution = NULL);
	void outputAutomorphismInfo(FILE *file, const MatrixDataPntr pMatrix = NULL) const;
	CC auto enumFlags() const						{ return m_enumFlags; }
	CC inline auto groupOrder() const				{ return m_nGroupOrder; }
	CC inline void setGroupOrder(uint val)			{ m_nGroupOrder = val; }
	CC inline T *permRow() const					{ return m_pPermutRow->elementPntr(); }
	CC inline T *permCol() const					{ return m_pPermutCol->elementPntr(); }
	CC inline auto permStorage(S nPart) const		{ return permStorage() + nPart; }
	CC inline auto permStorage() const				{ return m_pPermutStorage[0]; }
	CC inline auto permColStorage() const			{ return m_pPermutStorage[1]; }
	CC inline auto permRowStorage() const			{ return m_pPermutStorage[2]; }
	CC inline auto improvedSolution() const			{ return m_pImprovedSol; }
	CC bool groupIsTransitive() const;
	bool printMatrix(const designParam *pParam) const;
	CC auto stabiliserLengthExt() const				{ return m_nStabExtern; }
protected:
	void updateCanonicityChecker(T rowNumb, T colNumb);
	CC virtual void ConstructColumnPermutation(const MatrixDataPntr pMatrix)		{}
	virtual void CanonizeByColumns(MatrixDataPntr pMatrix, T *pColIdxStorage = NULL, CanonicityCheckerPntr pCanonChecker = NULL) const	{}
	CC inline auto getRowOrbits(int idx) const		{ return m_pObits[0][idx]; }
	CC inline auto getColOrbits(int idx) const		{ return m_pObits[1][idx]; }
	inline bool checkProperty(uint flag) const		{ return enumFlags() & flag; }
	CC inline auto numRow() const					{ return m_nNumRow; }
	CC void setStabiliserLengthExt(T len)			{ m_nStabExtern = len; }
	CC inline T numParts() const					{ return m_numParts; }
	CC virtual T lenStabilizer() const				{ return 0; }
	CK inline auto shiftToUnforcedOrbit(T nRow) const { return m_pShift[nRow]; }
private:
	CC T *init(T nRow, bool savePerm, T *pOrbits, T **pPermRows, bool groupOnParts);
	CC T next_permutation(T *perm, const T *pOrbits, T idx = ELEMENT_MAX, T lenStab = 0);
	CC void addAutomorphism(const T *pRowPerm, T* pOrbits, bool rowPermut = true, bool savePermut = false, bool calcGroupOrder = true);
	CC int checkColOrbit(T orbLen, T nColCurr, const S *pRow, const T *pRowPerm, T *pColPerm) const;
	CC inline void setStabilizerLength(T len)		{ m_nStabLength = len; }
	CC inline auto stabilizerLength() const			{ return m_nStabLength; }
	CC inline void setStabilizerLengthAut(T l)		{ m_nStabLengthAut = l; }
	CC inline auto stabilizerLengthAut() const		{ return m_nStabLengthAut; }
	CC inline void setNumRow(T nRow)				{ m_nNumRow = nRow; }
	CC inline auto numCol() const					{ return static_cast<T>(m_pPermutCol->numElement()); }
#define orbits()	m_pObits[0][0] 
	CC void updateGroupOrder();
	CC inline auto colNumbStorage() const			{ return m_nColNumbStorage; }
	CC inline auto counter() const					{ return m_pCounter; }
	CC inline void setColIndex(T *p)				{ m_pColIndex = p; }
	CC inline auto colIndex() const					{ return m_pColIndex; }
	CC inline void revert(T *perm, T j, T i) const {
		while (++i < --j) perm[i] ^= (perm[j] ^= (perm[i] ^= perm[j]));	
	}
	CC T constructColIndex(const ColOrbPntr pColOrbit, const ColOrbPntr pColOrbitIni, size_t colOrbLen, T shift = 0) const;
	CC inline void setPermStorage(PermutStoragePntr p, int idx = 0)	{ m_pPermutStorage[idx] = p; }
	CC T rowToChange(T nRow) const;
	void reconstructSolution(const ColOrbPntr pColOrbitStart, const ColOrbPntr pColOrbit,
		size_t colOrbLen, const ColOrbPntr pColOrbitIni, const T *pRowPerm, const T *pRowSolution, size_t solutionSize);
	CC void UpdateOrbits(const T *permut, T lenPerm, T *pOrbits, bool rowPermut, bool updateGroupOrder = false);
#if USE_STRONG_CANONICITY
	inline void setSolutionStorage(CSolutionStorage *p) { m_pSolutionStorage = p; }
	inline CSolutionStorage *solutionStorage() const { return m_pSolutionStorage; }

	CSolutionStorage *m_pSolutionStorage;
#else
#define setSolutionStorage(x)
#define solutionStorage()		(CSolutionStorage *)NULL
#endif

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T m_nStabLength;
	T m_nStabLengthAut;
	CPermut *m_pPermutRow;
	CPermut *m_pPermutCol;
	CPermut* m_pPermutSparse;
	CColNumbStorage **m_nColNumbStorage;
	PermutStoragePntr m_pPermutStorage[3];
	T *m_pObits[2][2];
	CCounter<int> *m_pCounter;
	T *m_pColIndex;
	T *m_pImprovedSol;
	const uint m_enumFlags;
	uint m_nGroupOrder;
	T m_nNumRow;
	const S m_numParts;
};

CanonicityChecker()::CCanonicityChecker(T nRow, T nCol, int rank, uint enumFlags, S numParts) : CRank(nRow, rank), m_enumFlags(enumFlags), m_numParts(numParts)
{
	m_pPermutRow = new CPermut(nRow);
	m_pPermutCol = new CPermut(nCol);
	m_pPermutSparse = NULL;

	setColIndex(new T[nCol << 1]);
	m_pCounter = new CCounter<int>(rank);
	m_nColNumbStorage = new CColNumbStorage *[rank];
	for (int i = rank; i--;)
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
#if USE_STRONG_CANONICITY
	delete solutionStorage();
#endif
}

CanonicityChecker(void)::updateCanonicityChecker(T rowNumb, T colNumb)
{
	// When we use automorphism on parts, we need two copies aech opf the following
	m_pPermutSparse = new CPermut[2];
	m_pPermutSparse[0].Init(rowNumb, new T[rowNumb]);
	m_pPermutSparse[1].Init(rowNumb, new T[colNumb]);
}

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
	const auto* pMatr = pEnum->matrix();
	const auto* pMatrPerm = pMatr;
	auto savePermut = rowPermut && (enumFlags() & t_outRowPermute);

	const auto colOrbLen = pEnum->colOrbitLen();
	const auto* pPartInfo = pMatr->partsInfo();
	PREPARE_PERM_OUT(permColStorage());

	bool calcGroupOrder = true;
	bool retVal = true;
	bool usingGroupOnBlocks = false;
	size_t lenMatr, numGroups = 0;

	T colNumb, startingRowNumb = 0;
	S *pMatrTo, *pMatrFrom;
	size_t idxPerm[16];			// to keep the indices of currently used permutation
	size_t* pIndxPerms = NULL;	// of i-th symmetrical group acting on the parts
	bool check_trivial_row_perm = false;

	T idxPartSrc[16];			// to keep the initial (source) indices of the parts
	T* pPartSrc = numParts() <= countof(idxPartSrc) ? idxPartSrc : new T[numParts()];
	for (auto i = numParts(); i--;)
		pPartSrc[i] = i;

	const S* pCurrSolution;
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
#endif
	auto pOrbits = orbits();
	while (true) {
#ifndef USE_CUDA
		T* permRows;
		T *permColumn = init(nRowMax, rowPermut, pOrbits, &permRows, usingGroupOnBlocks);

		T nRow = ELEMENT_MAX;
		while (true) {
			if (check_trivial_row_perm) {
				nRow = startingRowNumb;
				goto try_permut;
			}
		next_permut:
			nRow = next_permutation(permRows, pOrbits, nRow, lenStab);
			if (nRow == ELEMENT_MAX || nRow < lenStabilizer())
				break;

			check_trivial_row_perm = false;

		try_permut:
			OUT_PERM(permRowStorage(), permRows, nRowMax);
			OUTPUT_PERMUTATION(permColStorage(), pEnum->outFile(), pOrbits, numRow());

			// Loop for all remaining matrix's rows
			for (; nRow < nRowMax; nRow++) {
				const auto* pRow = pMatr->GetRow(nRow);
				const auto nRW = *(permRows + nRow);
				const auto* pRowPerm = pMatrPerm->GetRow(*(permRows + nRow));
				for (T nPart = 0; nPart < numParts(); nPart++) {
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

						if (!usingGroupOnBlocks && !nPart && pRowOut) {
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
							*pPartNumb = usingGroupOnBlocks? -1 : nPart;

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
					for (T nPart = 0; nPart < numParts(); nPart++) {
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

			addAutomorphism(permRows, pOrbits, rowPermut, savePermut, calcGroupOrder);
			if (check_trivial_row_perm) {
				check_trivial_row_perm = false;
				nRow = ELEMENT_MAX;
			}
			else
				nRow = ELEMENT_MAX - 1;
		}

		if (rowPermut && calcGroupOrder)
			updateGroupOrder();

		if (!retVal || !pGroupOnParts || !pGroupOnParts->useGroupOnParts(nRowMax))
			break;    // We don't have to test on groups of blocks

		if (!usingGroupOnBlocks) {
			// Corresponding data were not initialized yet
			usingGroupOnBlocks = true;
			numGroups = pGroupOnParts->numGroups();
			// Initialization of variables to use the group acting on the parts of the matrix
			pIndxPerms = numGroups < countof(idxPerm) ? idxPerm : new size_t[numGroups];
			memset(pIndxPerms, 0, numGroups * sizeof(*pIndxPerms));
			colNumb = (pMatrPerm = pCanonParam->pSpareMatrix)->colNumb();
			startingRowNumb = pGroupOnParts->getStartingRowNumb();
			// copy first startingRowNumb rows
			lenMatr = colNumb * startingRowNumb;
			pMatrTo = pMatrPerm->GetDataPntr();
			pMatrFrom = pMatr->GetDataPntr();
			memcpy(pMatrTo, pMatrFrom, lenMatr * sizeof(S));

			// Calculate pointers and length for copying remaining part of the matrix
			pMatrTo += lenMatr;
			pMatrFrom += lenMatr;
			lenMatr = (nRowMax - startingRowNumb) * colNumb * sizeof(S);
 
			calcGroupOrder = savePermut = false;	// Permutations and group order should not be saved for Spare Matrix
			pOrbits += numRow();
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

		// Permuting the parts of the matrix in accordance with the current set of permutations
		const auto* pPartsInfo = pMatr->partsInfo();
		memcpy(pMatrTo, pMatrFrom, lenMatr);

		for (auto i = numParts(); i--;)
			pPartSrc[i] = i;

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
			const size_t lenPart = nCol * sizeof(S);
			for (auto jFrom = 0; jFrom < pntr->permLength(); ++jFrom) {
				// For the current set of parts, check if the j-th part remains in place.
				auto jTo = *(pPerm + jFrom);
				if (jFrom == jTo)
					continue;

				// Calculate absolute indices and save source index 
				auto *pPartTo = pMatrTo + pPartsInfo->getShift(jTo += idxPart);
				const auto *pPartFrom = pMatrFrom + pPartsInfo->getShift(jFrom += idxPart);
				pPartSrc[jTo] = jFrom;

				for (auto j = startingRowNumb; j < nRowMax; j++) {
					memcpy(pPartTo, pPartFrom, lenPart);
					pPartTo += colNumb;
					pPartFrom += colNumb;
				}
			}
		}

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
