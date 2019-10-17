#include "PermutStorage.h"
#include "matrix.h"

template class CPermutStorage<MATRIX_ELEMENT_TYPE>;

template<class T>
void CPermutStorage<T>::orderPermutations(size_t *pPermPerm)
{
	const auto nPerm = numPerm();
	assert(nPerm > 0);
	pPermPerm[0] = 0;
	pPermPerm[1] = 1;
	if (nPerm <= 2)
		return;

	for (size_t i = nPerm; i-- > 2;)
		pPermPerm[i] = i;

	for (size_t i = 0; i < nPerm; i++) {
		auto jBest = UINT64_MAX;
        size_t idx = pPermPerm[i];
		const auto *pFirst = permutMemory() + lenPerm() * idx;
		for (auto j = i + 1; j < nPerm; j++) {
			auto jdx = pPermPerm[j];
			const auto *pSecnd = permutMemory() + lenPerm() * jdx;
			for (int k = 0; true; k++) {
				const int diff = (int)*(pFirst + k) - (int)*(pSecnd + k);
				if (!diff)
					continue;

				if (diff > 0) {
					jBest = j;
					idx = jdx;
					pFirst = pSecnd;
				}

				break;
			}
		}

		if (idx == pPermPerm[i])
			continue;

		pPermPerm[jBest] = pPermPerm[i];
		pPermPerm[i] = idx;
	}
}

template<class T>
void CPermutStorage<T>::UpdateOrbits(const T *permut, T lenPerm, T *pOrb, T idx) const
{
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


template<class T>
void CPermutStorage<T>::outputAutomorphismInfo(FILE *file, const T *pRowOrbits, 
	const CPermutStorage<T> *pPermColumn, const T *pColOrbits, const CMatrixData<T> *pMatrix) const
{
	const auto nRows = lenPerm();
	if (pPermColumn) {
		outputOrbits(file, pPermColumn, pMatrix, pRowOrbits, pColOrbits);
	} else	
		outputOrbits(file, pRowOrbits, nRows);

	outputPermutations(file, nRows, pPermColumn);
}

template<class T>
T *CPermutStorage<T>::CreateOrbits(const CPermutStorage<T> *pPermColumn,
									const CMatrixData<T> *pMatrix, T *pOrbits, T *pColOrbits, int firstpermIdx) const
{
	const auto nRows = lenPerm();
	const auto nCols = pPermColumn->lenPerm();
	T *pRowOrbits = pOrbits ? pOrbits : new T[2 * (nRows + (pColOrbits? 0 : nCols))];
	if (!pColOrbits)
		pColOrbits = pRowOrbits + 2 * nRows;

	// Orbits will be printed first and stabilizer will be second 
	T *pColOrbitsTmp = pColOrbits;
	T *pRowOrbitsTmp = pRowOrbits;
	pColOrbitsTmp[0] = 0;
	auto jPrev = 0;
	for (int j = 1; j < nCols; ++j) {
		// Compare i-th column with the previous one
		T *pElem = pMatrix->GetDataPntr() + jPrev;
		int i = 0;
		while (i < nRows && *pElem == *(pElem + 1)) {
			pElem += nCols;
			++i;
		}

		pColOrbitsTmp[j] = i == nRows ? pColOrbitsTmp[jPrev] : j;
		jPrev = j;
	}

	for (int j = 0; j < nRows; ++j)
		pRowOrbitsTmp[j] = j;

	// Identical permutation will be skipped
	const auto iMax = pPermColumn->numPerm();
	for (size_t i = firstpermIdx; i < iMax; ++i) {
		const auto pRowPermut = getPermutByIndex(i);
		if (pRowPermut[0] && pColOrbitsTmp == pColOrbits) {
			// Copying the orbits of stabilizer of first elements
			memcpy(pRowOrbitsTmp += nRows, pRowOrbits, nRows * sizeof(T));
			memcpy(pColOrbitsTmp += nCols, pColOrbits, nCols * sizeof(T));
		}

		UpdateOrbits(pRowPermut, nRows, pRowOrbitsTmp);
		UpdateOrbits(pPermColumn->getPermutByIndex(i), nCols, pColOrbitsTmp);
	}

	return pRowOrbits;
}

template<class T>
void CPermutStorage<T>::outputOrbits(FILE *file, const T *pOrbits, T lenPerm, const CPermutStorage<T> *pPermColumn) const
{
	outputPerm(file, pOrbits, lenPerm, pPermColumn ? pPermColumn->lenPerm() : 0);
}

template<class T>
void CPermutStorage<T>::outputOrbits(FILE *file, const CPermutStorage<T> *pPermColumn, 
	const CMatrixData<T> *pMatrix, const T *pRowOrbits, const T *pColOrbits) const
{
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

template<class T>
void CPermutStorage<T>::outputPermutations(FILE *file, T lenPerm, const CPermutStorage<T> *pPermColumn,
	const T *permutMemoryCol, const T *permutMemoryRow, int nOrbs) const
{
	char *pFormat;
	char *pBuffer = NULL;
	char *pBufferRows = NULL;
	size_t lenBuffer, lenBufferCol;
	size_t lenPermCol = 0;
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
