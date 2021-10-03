#pragma once
#include "DataTypes.h"
#include <stdio.h>

Class2Def(CMatrixData);
#define PermutStoragePntr Class2(CPermutStorage) *

Class2Def(CPermutStorage)
{
public:
	CC CPermutStorage();
	CC ~CPermutStorage()							{ delete[] permutMemory(); }
	void outputAutomorphismInfo(FILE *file, const T *pRowOrbits, const PermutStoragePntr pPermColumn = NULL,
		const T *pColOrbits = NULL, const MatrixDataPntr pMatrix = NULL, bool outPermut = true) const;
	void outputPerm(FILE *file, const S *perm, S lenPerm, S lenPerm2 = 0, const char *pColPerm = NULL, char **pBuffer = NULL, size_t *pLenBuffer = NULL,
			char **ppFormat = NULL) const;
	CK void adjustGenerators(int *pIdx, S lenIdx);
	CK size_t constructGroup();
	CK size_t findSolutionIndex(const S *pFirst, size_t idx, S *pMem, size_t *pCanonIdx, int &nCanon);
	CC inline auto lenPerm() const					{ return m_nLenPerm; }
	CC inline void initPermutStorage()				{ setLenMemUsed(0); }
	CC void savePermut(const S lenPermut, const S *perm = NULL);
	CK inline auto getPermutByIndex(size_t i) const	{ return permutMemory() + i * lenPerm(); }
	CK inline auto numPerm() const					{ return lenMemUsed() / lenPerm(); }
	CC S *allocateMemoryForPermut(S lenPermut);
	CC inline void setLenPerm(S val)				{ m_nLenPermByte = (m_nLenPerm = val) * sizeof(m_pPermutMem[0]); }
	CC void UpdateOrbits(const S *permut, S lenPerm, S *pOrbits, S idx = 0) const;
	S *CreateOrbits(const PermutStoragePntr pPermColumn, const  MatrixDataPntr pMatrix, S *pRowOrbits = NULL, S *pColOrbits = NULL, int firstpermIdx = 1) const;
	CC inline bool isEmpty() const					{ return !lenMemUsed(); }
protected:
private:
	CK inline auto lenPermByte() const				{ return m_nLenPermByte;  }
	CC inline void setPermutMemory(S *pntr)			{ m_pPermutMem = pntr; }
	CC inline auto permutMemory() const				{ return m_pPermutMem; }
	CC inline auto lenMemMax() const				{ return m_nLenMax; }
	CC inline auto lenMemUsed() const				{ return m_nLenUsed; }
	CC inline void setLenMemMax(size_t val)         { m_nLenMax = val; }
	CC inline void setLenMemUsed(size_t val)        { m_nLenUsed = val; }
	CC inline void deallocateMemoryForPermut()		{ m_nLenUsed -= lenPerm(); }
	void orderPermutations(size_t *pPermPerm);
	CK S *multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, S *pMultRes = NULL, size_t *pToIdx = NULL);
	CK void multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm);
	void outputOrbits(FILE *file, const T *pOrbits, T len, const PermutStoragePntr pPermColumn = NULL) const;
	void outputOrbits(FILE *file, const PermutStoragePntr pPermColumn,
		const MatrixDataPntr pMatrix, const S *pRowOrbits, const S *pColOrbits) const;
	void outputPermutations(FILE *file, T len, const PermutStoragePntr pPermColumn = NULL,
		const T *permutMemoryCol = NULL, const T *permutMemoryRow = NULL, int nOrbs = 0) const;

	S *m_pPermutMem;
	size_t m_nLenMax;
	size_t m_nLenUsed;
	S m_nLenPerm;
	size_t m_nLenPermByte;

#if OUT_PERMUTATION
	mutable int m_cntr = 0;
public:
	void resetGroupCntr()								{ if (this) m_cntr = 0; }
	void printPerm(const S* pPerm, bool savePerm = false, int add = 1, S permLen = 0) const;

#define PREPARE_PERM_OUT(x)		x->resetGroupCntr()
#define OUT_PERM(x, y, len)		x->printPerm(y, false, 1, len)
#else
#define resetGroupCntr()
#define printPerm(x,...)
#define PREPARE_PERM_OUT(x)
#define OUT_PERM(x,...)
#endif
};

PermutStorage()::CPermutStorage()
{
	setPermutMemory(NULL);
	setLenMemUsed(0);
	setLenMemMax(0);
}

PermutStorage(S *)::allocateMemoryForPermut(S lenPermut)
{
	const auto lenUsed = lenMemUsed();
	if (!lenUsed)
		setLenPerm(lenPermut);

	const auto newLength = lenUsed + lenPermut;
	if (lenMemMax() < newLength) {
		const auto len = 2 * newLength;
		setLenMemMax(len);
		auto *permutMem = new S[len];
		if (permutMemory()) {
			memcpy(permutMem, permutMemory(), lenUsed * sizeof(S));
			delete[] permutMemory();
		}

		setPermutMemory(permutMem);
	}

	setLenMemUsed(newLength);
	return permutMemory() + lenUsed;
}

PermutStorage(void)::savePermut(const S lenPermut, const S *perm)
{
	auto *pMem = allocateMemoryForPermut(lenPermut);
	printPerm(perm, true);
	if (perm)
		memcpy(pMem, perm, lenPermut * sizeof(S));
}

PermutStorage(void)::outputPerm(FILE *file, const S *perm, S lenPerm, S lenRowPerm, const char *pColPerm,
				char **ppBuffer, size_t *pLenBuffer, char **ppFormat) const {
	// Output of the last element in permutation storage. It is not necessary a permutation. 
	// Perhaps it's a result of some unsuccessful attempt to construct it.
	size_t len;
	char *pBuffer, *pFormat;
	if (!ppBuffer || *ppBuffer == NULL) {
		len = lenPerm + lenRowPerm;
		const size_t lenFrmt = len < 10 ? 2 : len < 100 ? 3 : len < 1000 ? 4 : 5;
		len *= lenFrmt;
		len += 4;
		pBuffer = new char[len + 8];
		SPRINTF(pFormat = pBuffer + len, "%%%zud", lenFrmt);
		if (ppBuffer) {
			*ppBuffer = pBuffer;
			*pLenBuffer = len;
			*ppFormat = pFormat;
		}
	}
	else {
		len = *pLenBuffer;
		pBuffer = *ppBuffer;
		pFormat = *ppFormat;
	}

	char *pBuf = pBuffer;
	for (size_t i = 0; i < lenPerm; i++)
		pBuf += SNPRINTF(pBuf, len - (pBuf - pBuffer), pFormat, *(perm + i));

	if (file) {
		outString(pBuffer, file);
		if (pColPerm)
			outString(pColPerm, file);

		outString("\n", file);
	}
	else
		printf("%s\n", pBuffer);

	if (!ppBuffer)
		delete[] pBuffer;
}

PermutStorage(void)::adjustGenerators(int *pIdx, S lenIdx)
{
	// Adjustment of the generators of the automorphism group on columns according to 
	// non-unforcible group of columns defined by indices in pIdx 
	size_t i, iMax = numPerm();
	if (!iMax)
		return;		// trivial group

	int buff[256];
	int *pFrom = iMax <= countof(buff) ? buff : new int[iMax];
	int j = 0;
	for (i = 0; i < lenIdx; i++)
		pFrom[pIdx[i]] = j++;

	auto *pTo = permutMemory();
	for (i = 0; i < iMax; i++, pTo += lenIdx) {
		auto *pColPerm = getPermutByIndex(i);
		for (size_t j = 0; j < lenIdx; j++)
			pTo[j] = pFrom[pColPerm[pIdx[j]]];
	}

	if (pFrom != buff)
		delete pFrom;

	setLenPerm(lenIdx);
	setLenMemUsed(lenIdx * iMax);
}

PermutStorage(size_t)::constructGroup()
{
	// Algorithm from http://www.sipria.ru/pdf/dt24114.pdf is implemented here.
	size_t buffer[16];
	const auto nPerm = numPerm();
	if (!nPerm)
		return 0;

	auto nextStart = nPerm < countof(buffer) ? buffer : new size_t[nPerm + 1];

	// Loop for all permutations except the trivial one
	resetGroupCntr();
	auto lastIdx = nPerm;
	S *pTmpPerm = NULL;
	for (size_t i = 1; i < nPerm; i++) {
		// Construct all degrees of current permutations
		auto tmpPermIdx = i;
		auto lastIdxTmp = lastIdx;
		printPerm(getPermutByIndex(i));
		while (!(pTmpPerm = multiplyPermutations(tmpPermIdx, i, pTmpPerm, &lastIdxTmp)))
			tmpPermIdx = lastIdxTmp - 1;

		// Save index for next round of multiplication
		const auto jMax = lastIdx;
		nextStart[i] = lastIdx = lastIdxTmp;

		// Find left and right multiplication of current permutation and  all permultations from previous round
		for (auto j = i + 1; j < jMax; j++) {
			pTmpPerm = multiplyPermutations(j, i, pTmpPerm, &lastIdx);
			pTmpPerm = multiplyPermutations(i, j, pTmpPerm, &lastIdx);
		}
	}

	size_t lastPrev = SIZE_MAX;
	while (lastPrev != lastIdx) {
		lastPrev = lastIdx;
		for (size_t i = 1; i < nPerm; i++) {
			for (auto j = nextStart[i]; j < lastIdx; j++) {
				pTmpPerm = multiplyPermutations(j, i, pTmpPerm, &lastIdx);
				pTmpPerm = multiplyPermutations(i, j, pTmpPerm, &lastIdx);
			}
		}
	}

	if (pTmpPerm)
		deallocateMemoryForPermut();

	if (nextStart != buffer)
		delete[] nextStart;

	return numPerm();
}

PermutStorage(void)::multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm)
{
	for (size_t j = fromIdx; j < toIdx; j++) {
		const auto permIdx = pPermPerm ? pPermPerm[j] : j;
		multiplyPermutations(permIdx, currPermIdx);

		for (size_t p = 0; p < permOrder; p++)
			multiplyPermutations(permIdx, lastIdx + p);
	}
}
#if CONSTR_ON_GPU
CK int cudaMemCmp(const unsigned char* left, const unsigned char* right, size_t length);
#define MEMCMP(left, right, length)		cudaMemCmp(left, right, length)
#else
#define MEMCMP(left, right, length)		memcmp(left, right, length)
#endif

PermutStorage(S *)::multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, S *pMultRes, size_t *pToIdx)
{
	if (!pMultRes)
		pMultRes = allocateMemoryForPermut(lenPerm());

	// Since the permutations could be moved in allocateMemoryForPermut,
	// we need to access the permuts by their indices here
	const auto *pFirst = getPermutByIndex(firstPermIdx);
	const auto *pSecond = getPermutByIndex(secondPermIdx);
	printPerm(pFirst, false, 0);
	// Do multiplication of two permutations here:
	for (auto k = lenPerm(); k--;)
		*(pMultRes + k) = *(pFirst + *(pSecond + k));

	// Compare with the previously constructed permutations
	for (auto i = *pToIdx; i--;) {
		if (!MEMCMP(getPermutByIndex(i), pMultRes, lenPermByte()))
			return pMultRes;   // this memory will be reused
	}

	printPerm(getPermutByIndex(*pToIdx));
	++*pToIdx;
	return NULL;
}

PermutStorage(size_t)::findSolutionIndex(const S *pFirst, size_t idx, S *pMem, size_t *pCanonIdx, int &nCanon)
{
	const S *pCurrSolution = pFirst + lenPerm() * idx;
	S *pCanonical = (S *)pCurrSolution;
	const auto len = lenPerm() * sizeof(*pMem);

	// Loop for all permutations, except the trivial one
	const auto *pPermut = permutMemory();
	const auto *pPermutEnd = pPermut + lenMemUsed();
	while ((pPermut += lenPerm()) < pPermutEnd) {
		for (auto i = lenPerm(); i--;)
			pMem[i] = pCurrSolution[*(pPermut + i)];

		if (MEMCMP(pMem, pCanonical, len) <= 0)
			continue;

		// Copy current solution as canonical
		if (pCanonical != pCurrSolution) {
			auto pTmp = pMem;
			pMem = pCanonical;
			pCanonical = pTmp;
		}
		else {
			pCanonical = pMem;
			pMem += lenPerm();
		}
	}

	if (pCanonical == pCurrSolution)
		return pCanonIdx[nCanon++] = idx; // Current solution is the canonical one

										  // Find corresponding canonical solution
	int i = nCanon;
	while (i-- && MEMCMP(pFirst + lenPerm() * pCanonIdx[i], pCanonical, len));

	if (i < 0) {
		// This could happen when we use reordering of solutions by the group of automorphisms
		// BIBD(12, 44, 11, 3, 2) provides many such examples
		return SIZE_MAX;
	}

	return pCanonIdx[i];
}

PermutStorage(void)::UpdateOrbits(const S* permut, S lenPerm, S* pOrb, S idx) const {
	// Update orbits of elements
	do {
		auto i = *(pOrb + idx);
		auto j = *(pOrb + permut[idx]);
		if (j == i)
			continue;

		if (j < i) {
			i ^= j;
			i ^= j ^= i;
		}

		for (auto k = lenPerm; k--;) {
			if (*(pOrb + k) == j)
				*(pOrb + k) = i;
		}
	} while (++idx < lenPerm);
}

PermutStorage(void)::outputAutomorphismInfo(FILE* file, const T* pRowOrbits,
	const  PermutStoragePntr pPermColumn, const T* pColOrbits, const MatrixDataPntr pMatrix, bool outPermut) const {
	const auto nRows = lenPerm();
	if (pPermColumn)
		outputOrbits(file, pPermColumn, pMatrix, pRowOrbits, pColOrbits);
	else
		outputOrbits(file, pRowOrbits, nRows);

	if (outPermut)
	    outputPermutations(file, nRows, pPermColumn);
}

PermutStorage(void)::outputOrbits(FILE* file, const T* pOrbits, T lenPerm, const PermutStoragePntr pPermColumn) const {
	outputPerm(file, pOrbits, lenPerm, pPermColumn ? pPermColumn->lenPerm() : 0);
}

PermutStorage(void)::outputOrbits(FILE* file, const PermutStoragePntr pPermColumn,
	const MatrixDataPntr pMatrix, const S* pRowOrbits, const S* pColOrbits) const {
	const auto nRows = lenPerm();
	const auto constructOrbits = !pColOrbits || !pRowOrbits;
	if (constructOrbits) {
		if (!pMatrix)
			return;

		pRowOrbits = CreateOrbits(pPermColumn, pMatrix);
		pColOrbits = pRowOrbits + 2 * nRows;
	}

	outputPermutations(file, nRows, pPermColumn, pColOrbits, pRowOrbits, 2);
	if (constructOrbits)
		delete[] pRowOrbits;
}

PermutStorage(void)::outputPermutations(FILE* file, T lenPerm, const PermutStoragePntr pPermColumn,
	const T* permutMemoryCol, const T* permutMemoryRow, int nOrbs) const
{
	char* pFormat;
	char* pBuffer = NULL;
	char* pBufferRows = NULL;
	size_t lenBuffer, lenBufferCol;
	T lenPermCol = 0;
	if (pPermColumn) {
		lenPermCol = pPermColumn->lenPerm();
		if (!permutMemoryCol)
			permutMemoryCol = pPermColumn->permutMemory();
	}

	size_t lenMax;
	if (!permutMemoryRow) {
		lenMax = lenMemUsed();
		permutMemoryRow = permutMemory();
	}
	else
		lenMax = nOrbs * lenPerm;


	for (size_t len = 0; len < lenMax; len += lenPerm) {
		if (permutMemoryCol) {
			outputPerm(NULL, permutMemoryCol, lenPermCol, lenPerm, NULL, &pBuffer, &lenBufferCol, &pFormat);
			permutMemoryCol += lenPermCol;
			if (!len) {
				lenBuffer = lenBufferCol;
				lenBufferCol = (lenBufferCol - 4) / (lenPermCol + lenPerm) * lenPermCol + 2;
				pBufferRows = pBuffer + lenBufferCol;
				lenBuffer -= lenBufferCol;
			}
		}

		outputPerm(file, permutMemoryRow + len, lenPerm, 0, pBuffer, &pBufferRows, &lenBuffer, &pFormat);
	}

	if (pBuffer)
		delete[] pBuffer;
	else
		delete[] pBufferRows;
}

PermutStorage(S*)::CreateOrbits(const  Class2(CPermutStorage)* pPermColumn,
	const MatrixDataPntr pMatrix, S* pOrbits, S* pColOrbits, int firstpermIdx) const {
	const auto nRows = lenPerm();
	const auto nCols = pPermColumn->lenPerm();
	auto* pRowOrbits = pOrbits ? pOrbits : new S[2 * (nRows + (pColOrbits ? 0 : nCols))];
	if (!pColOrbits)
		pColOrbits = pRowOrbits + 2 * nRows;

	// Orbits will be printed first and stabilizer will be second 
	auto* pColOrbitsTmp = pColOrbits;
	auto* pRowOrbitsTmp = pRowOrbits;
	pColOrbitsTmp[0] = 0;
	auto jPrev = 0;
	for (S j = 1; j < nCols; ++j) {
		// Compare i-th column with the previous one
		T* pElem = pMatrix->GetDataPntr() + jPrev;
		S i = 0;
		while (i < nRows && *pElem == *(pElem + 1)) {
			pElem += nCols;
			++i;
		}

		pColOrbitsTmp[j] = i == nRows ? pColOrbitsTmp[jPrev] : j;
		jPrev = j;
	}

	for (S j = 0; j < nRows; ++j)
		pRowOrbitsTmp[j] = j;

	// Identical permutation will be skipped
	const auto iMax = pPermColumn->numPerm();
	for (size_t i = firstpermIdx; i < iMax; ++i) {
		const auto pRowPermut = getPermutByIndex(i);
		if (pRowPermut[0] && pColOrbitsTmp == pColOrbits) {
			// Copying the orbits of stabilizer of first elements
			memcpy(pRowOrbitsTmp += nRows, pRowOrbits, nRows * sizeof(S));
			memcpy(pColOrbitsTmp += nCols, pColOrbits, nCols * sizeof(S));
		}

		UpdateOrbits(pRowPermut, nRows, pRowOrbitsTmp);
		UpdateOrbits(pPermColumn->getPermutByIndex(i), nCols, pColOrbitsTmp);
	}

	return pRowOrbits;
}