#include "IntersectionStorage.h"

template class CIntersection<TDATA_TYPES>;

FClass2(CIntersection, VariableMappingPntr)::prepareRowIntersections(const InSysPntr pMatrix, T currRowNumb, T lambda, T t) const
{
	if (currRowNumb <= 1)
		return NULL;			// Nothing to test

	const auto lastRowIdx = currRowNumb - 1;

	const auto* pCurrRow = pMatrix->GetRow(0);
	const auto* pLastRow = pMatrix->GetRow(lastRowIdx);

	// Create the set of the indices of block containing last element
	const auto r = pMatrix->GetR(0);		// TO DO: Need to be modified to support combined t-designs construction

	size_t blockIdx[256];
	auto ppBlockIdx = r <= countof(blockIdx) ? blockIdx : new size_t[r];
	size_t idx = 0;
	const auto nCol = pMatrix->colNumb();
	for (size_t j = 0; j < nCol; j++) {
		if (*(pLastRow + j)) {
			*(ppBlockIdx + idx++) = j;
			if (idx == r)
				break;   // No more blocks, containing last element
		}
	}

	// Create indices of block containing 0-th and last, 1-st and last etc element.
	const size_t* pNumb;
	auto* pIntersection = intersectionParam(&pNumb, currRowNumb);
	for (auto k = pNumb[0]; k--; pCurrRow += nCol, pIntersection += lambda) {
		size_t i = 0;
		for (size_t j = 0; j < r; j++) {
			if (*(pCurrRow + ppBlockIdx[j])) {
				*(pIntersection + i) = ppBlockIdx[j];
				if (++i == lambda)
					break;	// No more blocks, containing two last elements
			}
		}
	}

	t -= 2;
	if (t >= currRowNumb)
		t = lastRowIdx;

	S tuple[10];
	T* matrixRowPntr[10];
	auto* pTuple = t <= countof(tuple) ? tuple : new S[t];
	auto pMatrixRowPntr = t <= countof(matrixRowPntr) ? matrixRowPntr : new T * [t];
	for (T i = 1; i < t; i++) {
		// Construct all (i+1)-subsets of first currentRowNumb() elements
		uint k = 0;
		pTuple[0] = -1;
		while (k != -1) {
			auto n = pTuple[k];
			for (; k <= i; k++)
				pMatrixRowPntr[k] = pMatrix->GetRow(pTuple[k] = ++n);

			for (auto j = idx; j--;) {
				const auto blockIdx = ppBlockIdx[j];
				uint m = -1;
				while (++m <= i && *(pMatrixRowPntr[m] + blockIdx));
				if (m > i)
					*pIntersection++ = blockIdx;
			}

			// Construct next (i+1) subset
			k = i + 1;
			n = lastRowIdx;
			while (k-- && pTuple[k] == --n);
		}
	}

	if (pTuple != tuple)
		delete[] pTuple;

	if (pMatrixRowPntr != matrixRowPntr)
		delete[] pMatrixRowPntr;

	if (ppBlockIdx != blockIdx)
		delete[] ppBlockIdx;

#if USE_EXRA_EQUATIONS
	return constructExtraEquations(t, nVar);
#else
	return NULL;
#endif
}
