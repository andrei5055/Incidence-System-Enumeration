#include "PermutStorage.h"
#include "matrix.h"

template class CPermutStorage<MATRIX_ELEMENT_TYPE, SIZE_TYPE>;

PermutStorage(void)::orderPermutations(size_t *pPermPerm)
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

PermutStorage(void)::UpdateOrbits(const S *permut, S lenPerm, S *pOrb, S idx) const {
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


PermutStorage(void)::outputAutomorphismInfo(FILE *file, const S *pRowOrbits,
	const  IClass2(PermutStorage) *pPermColumn, const S *pColOrbits, const IClass2(MatrixData) *pMatrix) const {
	const auto nRows = lenPerm();
	if (pPermColumn) {
		outputOrbits(file, pPermColumn, pMatrix, pRowOrbits, pColOrbits);
	} else	
		outputOrbits(file, pRowOrbits, nRows);

	outputPermutations(file, nRows, pPermColumn);
}

PermutStorage(S *)::CreateOrbits(const  IClass2(PermutStorage)*pPermColumn,
								const IClass2(MatrixData) *pMatrix, S *pOrbits, S *pColOrbits, int firstpermIdx) const {
	const auto nRows = lenPerm();
	const auto nCols = pPermColumn->lenPerm();
	auto *pRowOrbits = pOrbits ? pOrbits : new S[2 * (nRows + (pColOrbits? 0 : nCols))];
	if (!pColOrbits)
		pColOrbits = pRowOrbits + 2 * nRows;

	// Orbits will be printed first and stabilizer will be second 
	auto *pColOrbitsTmp = pColOrbits;
	auto *pRowOrbitsTmp = pRowOrbits;
	pColOrbitsTmp[0] = 0;
	auto jPrev = 0;
	for (S j = 1; j < nCols; ++j) {
		// Compare i-th column with the previous one
		T *pElem = pMatrix->GetDataPntr() + jPrev;
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
			memcpy(pRowOrbitsTmp += nRows, pRowOrbits, nRows * sizeof(T));
			memcpy(pColOrbitsTmp += nCols, pColOrbits, nCols * sizeof(T));
		}

		UpdateOrbits(pRowPermut, nRows, pRowOrbitsTmp);
		UpdateOrbits(pPermColumn->getPermutByIndex(i), nCols, pColOrbitsTmp);
	}

	return pRowOrbits;
}

PermutStorage(void)::outputOrbits(FILE *file, const S *pOrbits, S lenPerm, const IClass2(PermutStorage) *pPermColumn) const {
	outputPerm(file, pOrbits, lenPerm, pPermColumn ? pPermColumn->lenPerm() : 0);
}

PermutStorage(void)::outputOrbits(FILE *file, const IClass2(PermutStorage) *pPermColumn,
	const IClass2(MatrixData) *pMatrix, const S *pRowOrbits, const S *pColOrbits) const {
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

PermutStorage(void)::outputPermutations(FILE *file, S lenPerm, const IClass2(PermutStorage) *pPermColumn,
	const S *permutMemoryCol, const S *permutMemoryRow, int nOrbs) const
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
