#include "DataTypes.h"
#include <stdio.h>

template<class T>
class CPermutStorage
{
public:
	CC CPermutStorage();
	CC ~CPermutStorage()							{ delete[] permutMemory(); }
	void outputPermutations(FILE *file, T len, const CPermutStorage<T> *pPermColumn = NULL) const;
	void outputPerm(FILE *file, const T *perm, size_t lenPerm, size_t lenPerm2 = 0, char **pBuffer = NULL, size_t *pLenBuffer = NULL,
			char **ppFormat = NULL, const char *pColPerm = NULL) const;
	CK void adjustGenerators(int *pIdx, T lenIdx);
	CK size_t constructGroup();
	CK size_t findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon);
	CC inline T lenPerm() const						{ return m_nLenPerm; }
	CC inline void initPermutStorage()				{ setLenMemUsed(0); setLenPerm(MATRIX_ELEMENT_MAX); }
	CC void savePermut(const T lenPermut, const T *perm = NULL);
	CK inline T *getPermutByIndex(size_t i) const	{ return permutMemory() + i * lenPerm(); }
	CK inline size_t numPerm() const				{ return lenMemUsed() / lenPerm(); }
	CC T *allocateMemoryForPermut(T lenPermut);
	CC inline void setLenPerm(T val)				{ m_nLenPermByte = (m_nLenPerm = val) * sizeof(m_pPermutMem[0]); }
protected:
private:
	CK inline size_t lenPermByte() const			{ return m_nLenPermByte;  }
	CC inline void setPermutMemory(T *pntr)			{ m_pPermutMem = pntr; }
	CC inline T *permutMemory() const				{ return m_pPermutMem; }
	CC inline size_t lenMemMax() const              { return m_nLenMax; }
	CC inline size_t lenMemUsed() const             { return m_nLenUsed; }
	CC inline void setLenMemMax(size_t val)         { m_nLenMax = val; }
	CC inline void setLenMemUsed(size_t val)        { m_nLenUsed = val; }
	CC inline void deallocateMemoryForPermut()		{ m_nLenUsed -= lenPerm(); }
	void orderPermutations(size_t *pPermPerm);
	CK T *multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, T *pMultRes = NULL, size_t *pToIdx = NULL);
	CK void multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm);

	T *m_pPermutMem;
	size_t m_nLenMax;
	size_t m_nLenUsed;
	T m_nLenPerm;
	size_t m_nLenPermByte;
};


template<class T>
CPermutStorage<T>::CPermutStorage()
{
	setPermutMemory(NULL);
	setLenMemUsed(0);
	setLenMemMax(0);
}

template<class T>
T *CPermutStorage<T>::allocateMemoryForPermut(T lenPermut)
{
	const size_t lenUsed = lenMemUsed();
	const size_t newLength = lenUsed + lenPermut;

	if (lenMemMax() < newLength) {
		const size_t len = 2 * newLength;
		setLenMemMax(len);
		auto *permutMem = new T[len];
		if (permutMemory()) {
			memcpy(permutMem, permutMemory(), lenUsed * sizeof(*permutMemory()));
			delete[] permutMemory();
		}

		setPermutMemory(permutMem);
	}

	setLenMemUsed(newLength);
	return permutMemory() + lenUsed;
}

template<class T>
void CPermutStorage<T>::savePermut(const T lenPermut, const T *perm)
{
	auto *pMem = allocateMemoryForPermut(lenPermut);

	if (lenPerm() == MATRIX_ELEMENT_MAX)
		setLenPerm(lenPermut);

	if (!perm) {
		for (T i = lenPermut; i--;)
			*(pMem + i) = i;
	}
	else
		memcpy(pMem, perm, lenPermut * sizeof(*permutMemory()));
}

template<class T>
void CPermutStorage<T>::outputPerm(FILE *file, const T *perm, size_t lenPerm, size_t lenRowPerm, char **ppBuffer, size_t *pLenBuffer, 
									char **ppFormat, const char *pColPerm) const
{
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

	if (!ppBuffer)
		delete[] pBuffer;
}

template<class T>
void CPermutStorage<T>::outputPermutations(FILE *file, T lenPerm, const CPermutStorage<T> *pPermColumn) const
{
	char *pFormat;
	char *pBuffer = NULL;
	char *pBufferRows = NULL;
	size_t lenBuffer, lenBufferCol;
	for (size_t len = 0; len < lenMemUsed(); len += lenPerm) {
		if (pPermColumn) {
			auto lenPermCol = pPermColumn->lenPerm();
			outputPerm(NULL, pPermColumn->permutMemory() + lenPermCol, lenPermCol, lenPerm, &pBuffer, &lenBufferCol, &pFormat);
			if (!len) {
				lenBuffer = lenBufferCol;
				lenBufferCol = (lenBufferCol - 4) / (lenPermCol + lenPerm) * lenPermCol + 2;
				pBufferRows = pBuffer + lenBufferCol;
				lenBuffer -= lenBufferCol;
			}
		}

		outputPerm(file, permutMemory() + len, lenPerm, 0, &pBufferRows, &lenBuffer, &pFormat, pBuffer);
	}

	if (pBuffer)
		delete[] pBuffer;
	else
		delete[] pBufferRows;
}

template<class T>
void CPermutStorage<T>::adjustGenerators(int *pIdx, T lenIdx)
{
	// Adjustement of the generators of the automorphism group on columns according to 
	// non-unforcible group of columns defined by indeces in pIdx 
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

template<class T>
size_t CPermutStorage<T>::constructGroup()
{
	// Algorithm from http://www.sipria.ru/pdf/dt24114.pdf is implemented here.
	size_t buffer[16];
	const auto nPerm = numPerm();
	if (!nPerm)
		return 0;

	size_t *nextStart = nPerm < countof(buffer) ? buffer : new size_t[nPerm + 1];

	// Loop for all permutations except the trivial one
	auto lastIdx = nPerm;
	T *pTmpPerm = NULL;
	for (size_t i = 1; i < nPerm; i++) {
		// Construct all degrees of current permutations
		auto tmpPermIdx = i;
		auto lastIdxTmp = lastIdx;
		while (!(pTmpPerm = multiplyPermutations(tmpPermIdx, i, pTmpPerm, &lastIdxTmp)))
			tmpPermIdx = lastIdxTmp - 1;

		// Save index for next round of multiplication
		const auto jMax = lastIdx;
		nextStart[i] = lastIdx = lastIdxTmp;

		// Find left and right multiplication of current permutation and  all permultations from previous round
		for (size_t j = i + 1; j < jMax; j++) {
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

template<class T>
void CPermutStorage<T>::multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm)
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

template<class T>
T *CPermutStorage<T>::multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, T *pMultRes, size_t *pToIdx)
{
	if (!pMultRes)
		pMultRes = allocateMemoryForPermut(lenPerm());

	// Since the permutations could be moved in allocateMemoryForPermut,
	// we need to access the permuts by their indices here
	const auto *pFirst = getPermutByIndex(firstPermIdx);
	const auto *pSecond = getPermutByIndex(secondPermIdx);

	// Do multiplication of two permutations here:
	for (auto k = lenPerm(); k--;)
		*(pMultRes + k) = *(pFirst + *(pSecond + k));

	// Compare with the previously constructed permutations
	for (auto i = *pToIdx; i--;) {
		if (!MEMCMP(getPermutByIndex(i), pMultRes, lenPermByte()))
			return pMultRes;   // this memory will be reused
	}

	++*pToIdx;
	return NULL;
}

template<class T>
size_t CPermutStorage<T>::findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon)
{
	const VECTOR_ELEMENT_TYPE *pCurrSolution = pFirst + lenPerm() * idx;
	VECTOR_ELEMENT_TYPE *pCanonical = (VECTOR_ELEMENT_TYPE *)pCurrSolution;
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
			VECTOR_ELEMENT_TYPE *pTmp = pMem;
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
		// It could happen when we are using reordering of solutions by the automorphism group
		// BIBD(12, 44, 11, 3, 2) provides a lot's of such examples
		return SIZE_MAX;
	}

	return pCanonIdx[i];
}
